// PhantomMidSide.cpp
// A simple stereo mid-side processor using JACK.
// It takes a stereo input (in_left and in_right), converts it to mid and side,
// applies independent gain controls (midGain and sideGain), then reconstructs the stereo signal.
// The pan processing follows an equal-power logic by using the formulas:
//   mid = (L + R) / 2
//   side = (L - R) / 2
//   L_out = midGain * mid + sideGain * side
//   R_out = midGain * mid - sideGain * side
//
// Real-time control via a console allows updating midGain and sideGain.
// Compile with:
//   g++ -std=c++11 PhantomMidSide.cpp -ljack -lpthread -o PhantomMidSide

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <jack/jack.h>
#include <iostream>
#include <atomic>
#include <thread>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <chrono>
#include <string>

using namespace std;

class PhantomMidSide {
private:
    jack_client_t* client;
    jack_port_t* in_left;
    jack_port_t* in_right;
    jack_port_t* out_left;
    jack_port_t* out_right;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // User-controllable parameters: mid and side gain.
    // When both are set to 1.0, the output is identical to the input.
    atomic<float> midGain;   // Default 1.0.
    atomic<float> sideGain;  // Default 1.0.

    // JACK process callback: processes each block of audio.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomMidSide* pms = static_cast<PhantomMidSide*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(pms->in_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(pms->in_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(pms->out_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(pms->out_right, nframes));

        float currentMidGain = pms->midGain.load();
        float currentSideGain = pms->sideGain.load();

        // Process each sample.
        for (jack_nframes_t i = 0; i < nframes; i++) {
            float L = inL[i];
            float R = inR[i];
            // Convert stereo to mid-side.
            float mid = (L + R) * 0.5f;
            float side = (L - R) * 0.5f;
            // Apply gain adjustments.
            float L_out = currentMidGain * mid + currentSideGain * side;
            float R_out = currentMidGain * mid - currentSideGain * side;
            // Write outputs.
            outL[i] = L_out;
            outR[i] = R_out;
        }
        return 0;
    }

    // Control thread: allows real-time updates of mid and side gain.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomMidSide] Enter new parameters (midGain sideGain), or type 'q' to quit:" << endl;
                cout << "For example: \"1.0 0.5\" (mid at unity, side reduced to 0.5): ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newMidGain, newSideGain;
            if (!(iss >> newMidGain >> newSideGain)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomMidSide] Invalid input. Please try again." << endl;
                continue;
            }
            midGain.store(newMidGain);
            sideGain.store(newSideGain);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomMidSide] Updated parameters: midGain = " << newMidGain
                    << ", sideGain = " << newSideGain << endl;
            }
        }
    }

public:
    PhantomMidSide(const char* client_name = "PhantomMidSide")
        : client(nullptr), running(true), midGain(1.0f), sideGain(1.0f)
    {
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomMidSide: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_left = jack_port_register(client, "in_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        in_right = jack_port_register(client, "in_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_left || !in_right || !out_left || !out_right) {
            jack_client_close(client);
            throw runtime_error("PhantomMidSide: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomMidSide: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomMidSide: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomMidSide::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomMidSide] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomMidSide] Default midGain = " << midGain.load() << ", sideGain = " << sideGain.load() << " (center)" << endl;
    }

    ~PhantomMidSide() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomMidSide] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomMidSide] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomMidSide midSide;
        midSide.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomMidSide] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

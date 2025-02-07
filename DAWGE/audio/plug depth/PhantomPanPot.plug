// PhantomPanner.cpp
// A simple stereo panner plugin using JACK.
// It takes a mono input and pans it to stereo using an equal-power panning law.
// The pan parameter ranges from -1.0 (full left) to +1.0 (full right).
//
// Compile with:
//   g++ -std=c++11 PhantomPanner.cpp -ljack -lpthread -o PhantomPanner

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

class PhantomPanner {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_left;
    jack_port_t* out_right;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Pan parameter: -1.0 (full left) to +1.0 (full right); 0.0 is center.
    atomic<float> pan;

    // JACK process callback: reads mono input and outputs stereo.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomPanner* panner = static_cast<PhantomPanner*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(panner->in_port, nframes));
        float* left_out = static_cast<float*>(jack_port_get_buffer(panner->out_left, nframes));
        float* right_out = static_cast<float*>(jack_port_get_buffer(panner->out_right, nframes));

        float currentPan = panner->pan.load();
        // Map pan from [-1, 1] to angle between 0 and π/2.
        float angle = (currentPan + 1.0f) * (M_PI / 4.0f);
        float leftGain = cos(angle);
        float rightGain = sin(angle);

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float sample = in[i];
            left_out[i] = sample * leftGain;
            right_out[i] = sample * rightGain;
        }
        return 0;
    }

    // Control thread: allows updating the pan parameter.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomPanner] Enter pan value (-1.0 for full left to 1.0 for full right), or 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newPan;
            if (!(iss >> newPan)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomPanner] Invalid input. Please try again." << endl;
                continue;
            }
            // Clamp newPan to [-1.0, 1.0].
            if (newPan < -1.0f)
                newPan = -1.0f;
            if (newPan > 1.0f)
                newPan = 1.0f;
            pan.store(newPan);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomPanner] Updated pan value: " << newPan << endl;
            }
        }
    }

public:
    PhantomPanner(const char* client_name = "PhantomPanner")
        : client(nullptr), running(true), pan(0.0f)  // Default pan is center.
    {
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomPanner: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_left || !out_right) {
            jack_client_close(client);
            throw runtime_error("PhantomPanner: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomPanner: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomPanner: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomPanner::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomPanner] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomPanner] Default pan: " << pan.load() << " (center)" << endl;
    }

    ~PhantomPanner() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomPanner] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomPanner] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomPanner panner;
        panner.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomPanner] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

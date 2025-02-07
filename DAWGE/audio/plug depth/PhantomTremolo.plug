// PhantomTremolo.cpp
// A simple mono tremolo effect using JACK.
// The effect modulates the amplitude of the input signal with an LFO.
// Real-time adjustable parameters:
//   - LFO Frequency (Hz)
//   - Depth (0.0 to 1.0)
//   - Mix (0.0 = dry, 1.0 = fully modulated)
// Compile with:
//   g++ -std=c++11 PhantomTremolo.cpp -ljack -lpthread -o PhantomTremolo

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

using namespace std;

class PhantomTremolo {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Tremolo parameters.
    atomic<float> lfoFreq;  // LFO frequency in Hz.
    atomic<float> depth;    // Depth (0.0 to 1.0); 0 means no modulation, 1 means full modulation.
    atomic<float> mix;      // Mix between dry and processed (0.0 = dry, 1.0 = fully modulated).

    // LFO phase accumulator.
    float lfoPhase;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomTremolo* trem = static_cast<PhantomTremolo*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(trem->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(trem->out_port, nframes));

        float currentFreq = trem->lfoFreq.load();
        float currentDepth = trem->depth.load();
        float currentMix = trem->mix.load();

        // Calculate LFO increment (radians per sample).
        float dt = 1.0f / trem->sample_rate;
        float lfoInc = 2.0f * M_PI * currentFreq * dt;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Compute the modulation factor:
            // When sin(LFO_phase) is 1, modFactor = 1 (no attenuation).
            // When sin(LFO_phase) is -1, modFactor = 1 - depth.
            float modFactor = 1.0f - currentDepth + currentDepth * (0.5f * (1.0f + sinf(trem->lfoPhase)));
            // Apply amplitude modulation.
            float modulated = dry * modFactor;
            // Blend dry and modulated signals.
            out[i] = (1.0f - currentMix) * dry + currentMix * modulated;

            // Advance LFO phase.
            trem->lfoPhase += lfoInc;
            if (trem->lfoPhase >= 2.0f * M_PI)
                trem->lfoPhase -= 2.0f * M_PI;
        }
        return 0;
    }

    // Control thread: update parameters in real time.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomTremolo] Enter parameters: LFO Frequency (Hz), Depth (0.0-1.0), Mix (0.0-1.0)" << endl;
                cout << "e.g., \"5 0.8 0.7\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newFreq, newDepth, newMix;
            if (!(iss >> newFreq >> newDepth >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomTremolo] Invalid input. Please try again." << endl;
                continue;
            }
            if (newFreq < 0.0f) newFreq = 0.0f;
            if (newDepth < 0.0f) newDepth = 0.0f;
            if (newDepth > 1.0f) newDepth = 1.0f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            lfoFreq.store(newFreq);
            depth.store(newDepth);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomTremolo] Updated parameters:" << endl;
                cout << "  LFO Frequency = " << newFreq << " Hz" << endl;
                cout << "  Depth = " << newDepth << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomTremolo(const char* client_name = "PhantomTremolo")
        : client(nullptr), running(true), lfoPhase(0.0f)
    {
        // Set default parameters.
        lfoFreq.store(5.0f);   // 5 Hz LFO.
        depth.store(0.8f);     // 80% modulation depth.
        mix.store(0.7f);       // 70% mix.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomTremolo: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomTremolo: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomTremolo: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomTremolo: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomTremolo::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomTremolo] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomTremolo] Default parameters:" << endl;
        cout << "  LFO Frequency = " << lfoFreq.load() << " Hz" << endl;
        cout << "  Depth = " << depth.load() << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomTremolo() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomTremolo] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomTremolo] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomTremolo tremolo;
        tremolo.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomTremolo] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

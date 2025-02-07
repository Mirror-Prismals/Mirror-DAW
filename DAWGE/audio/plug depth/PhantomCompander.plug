// PhantomCompander.cpp
// A simple mono compander plugin using JACK and standard C++.
// It applies compression when the input exceeds a threshold and expansion when it falls below,
// with separate mix controls for the compression and expansion effects.
//
// Real-time adjustable parameters (via console):
//   - Threshold (dB): The reference level (e.g., -20 dB)
//   - Compression Ratio: For signals above the threshold (e.g., 4.0)
//   - Expansion Ratio: For signals below the threshold (e.g., 2.0)
//   - compMix: Mix for the compression branch (0.0 = dry, 1.0 = fully processed)
//   - expMix: Mix for the expansion branch (0.0 = dry, 1.0 = fully processed)
//
// Compile with:
//   g++ -std=c++11 PhantomCompander.cpp -ljack -lpthread -o PhantomCompander

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

class PhantomCompander {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Parameters (set via control thread)
    // Threshold in dB (e.g., -20 dB); will be converted to linear inside process.
    atomic<float> threshold_dB;
    // Compression ratio (for signals above threshold).
    atomic<float> compRatio;
    // Expansion ratio (for signals below threshold).
    atomic<float> expRatio;
    // Mix for the compression branch (0.0 = dry, 1.0 = fully processed).
    atomic<float> compMix;
    // Mix for the expansion branch.
    atomic<float> expMix;

    // Utility: sign function declared as static.
    static inline float sign(float x) {
        return (x >= 0.0f) ? 1.0f : -1.0f;
    }

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomCompander* compander = static_cast<PhantomCompander*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(compander->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(compander->out_port, nframes));

        // Convert threshold from dB to linear.
        float thresh_lin = powf(10.0f, compander->threshold_dB.load() / 20.0f);
        float cRatio = compander->compRatio.load();
        float eRatio = compander->expRatio.load();
        float mixComp = compander->compMix.load();
        float mixExp = compander->expMix.load();

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float x = in[i];
            float absX = fabs(x);
            float processed = x;  // Default to dry signal.
            // Avoid division by zero.
            if (absX < 1e-9f) {
                out[i] = x;
                continue;
            }
            if (absX > thresh_lin) {
                // Compression branch.
                float excess = absX - thresh_lin;
                float desired = thresh_lin + excess / cRatio;
                float gain = desired / absX;
                processed = sign(x) * absX * gain;
                out[i] = (1.0f - mixComp) * x + mixComp * processed;
            }
            else if (absX < thresh_lin) {
                // Expansion branch.
                float deficit = thresh_lin - absX;
                float desired = thresh_lin - deficit * eRatio;
                if (desired < 0.0f) desired = 0.0f;
                float gain = (absX > 0) ? (desired / absX) : 1.0f;
                processed = sign(x) * absX * gain;
                out[i] = (1.0f - mixExp) * x + mixExp * processed;
            }
            else {
                out[i] = x;
            }
        }
        return 0;
    }

    // Control thread: update parameters via console.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomCompander] Enter parameters:" << endl;
                cout << "Format: <Threshold_dB> <CompRatio> <ExpRatio> <compMix> <expMix>" << endl;
                cout << "e.g., \"-20 4.0 2.0 1.0 1.0\" for -20 dB threshold, 4:1 compression, 2:1 expansion, full effect," << endl;
                cout << "or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newThresh, newCompRatio, newExpRatio, newCompMix, newExpMix;
            if (!(iss >> newThresh >> newCompRatio >> newExpRatio >> newCompMix >> newExpMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomCompander] Invalid input. Please try again." << endl;
                continue;
            }
            // Clamp mix values between 0 and 1.
            if (newCompMix < 0.0f) newCompMix = 0.0f;
            if (newCompMix > 1.0f) newCompMix = 1.0f;
            if (newExpMix < 0.0f) newExpMix = 0.0f;
            if (newExpMix > 1.0f) newExpMix = 1.0f;
            threshold_dB.store(newThresh);
            compRatio.store(newCompRatio);
            expRatio.store(newExpRatio);
            compMix.store(newCompMix);
            expMix.store(newExpMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomCompander] Updated parameters:" << endl;
                cout << "  Threshold = " << newThresh << " dB" << endl;
                cout << "  Compression Ratio = " << newCompRatio << endl;
                cout << "  Expansion Ratio = " << newExpRatio << endl;
                cout << "  Compression Mix = " << newCompMix << endl;
                cout << "  Expansion Mix = " << newExpMix << endl;
            }
        }
    }

public:
    PhantomCompander(const char* client_name = "PhantomCompander")
        : client(nullptr), running(true)
    {
        // Set default parameters.
        threshold_dB.store(-20.0f);  // -20 dB threshold.
        compRatio.store(4.0f);       // 4:1 compression.
        expRatio.store(2.0f);        // 2:1 expansion.
        compMix.store(1.0f);         // Fully apply compression effect for signals above threshold.
        expMix.store(1.0f);          // Fully apply expansion effect for signals below threshold.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomCompander: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomCompander: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomCompander: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomCompander: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomCompander::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomCompander] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomCompander] Default parameters:" << endl;
        cout << "  Threshold = " << threshold_dB.load() << " dB" << endl;
        cout << "  Compression Ratio = " << compRatio.load() << endl;
        cout << "  Expansion Ratio = " << expRatio.load() << endl;
        cout << "  Compression Mix = " << compMix.load() << endl;
        cout << "  Expansion Mix = " << expMix.load() << endl;
    }

    ~PhantomCompander() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomCompander] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomCompander] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomCompander compander;
        compander.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomCompander] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

// PhantomDuck.cpp
// A simple sidechain duck compressor using JACK and standard C++.
// It takes two mono inputs: "main" (the signal to be ducked) and "side" (the key signal).
// When the sidechain envelope exceeds a threshold, the main signal is attenuated.
// The output is a blend between the original main signal and the compressed (ducked) main signal.
//
// Real-time adjustable parameters (via console):
//   - Threshold (dB) [e.g., -30 dB]
//   - Ratio (e.g., 4.0)
//   - Attack Time (ms) [e.g., 10 ms]
//   - Release Time (ms) [e.g., 50 ms]
//   - Mix (0.0 = dry, 1.0 = fully ducked)
//
// Compile with:
//   g++ -std=c++11 PhantomDuck.cpp -ljack -lpthread -o PhantomDuck

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

class PhantomDuck {
private:
    jack_client_t* client;
    // Two input ports: one for the main signal, one for the sidechain.
    jack_port_t* main_in_port;
    jack_port_t* side_in_port;
    // One output port.
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Compressor parameters.
    // Threshold (dB) at which to start ducking.
    atomic<float> threshold_dB;  // e.g., -30 dB.
    atomic<float> ratio;         // e.g., 4.0.
    atomic<float> attackTime;    // in ms, e.g., 10.
    atomic<float> releaseTime;   // in ms, e.g., 50.
    atomic<float> mix;           // 0.0 = dry main, 1.0 = fully ducked main.

    // Envelope for the sidechain signal.
    float envelope;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomDuck* duck = static_cast<PhantomDuck*>(arg);
        float* mainIn = static_cast<float*>(jack_port_get_buffer(duck->main_in_port, nframes));
        float* sideIn = static_cast<float*>(jack_port_get_buffer(duck->side_in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(duck->out_port, nframes));

        float dt = 1.0f / duck->sample_rate;
        float dt_ms = dt * 1000.0f;
        float att = duck->attackTime.load();
        float rel = duck->releaseTime.load();
        float currentThreshold_dB = duck->threshold_dB.load();
        float currentThreshold = powf(10.0f, currentThreshold_dB / 20.0f); // convert threshold dB to linear.
        float currentRatio = duck->ratio.load();
        float currentMix = duck->mix.load();

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float mainSample = mainIn[i];
            float sideSample = sideIn[i];
            float absSide = fabs(sideSample);
            // Update envelope using attack/release exponential smoothing.
            if (absSide > duck->envelope)
                duck->envelope = expf(-dt_ms / att) * duck->envelope + (1.0f - expf(-dt_ms / att)) * absSide;
            else
                duck->envelope = expf(-dt_ms / rel) * duck->envelope + (1.0f - expf(-dt_ms / rel)) * absSide;
            // Avoid log of zero by ensuring envelope is at least a tiny value.
            float effectiveEnv = (duck->envelope > 1e-6f) ? duck->envelope : 1e-6f;
            // Convert envelope to dB.
            float env_dB = 20.0f * log10f(effectiveEnv);
            float gain;
            if (env_dB > currentThreshold_dB) {
                // Compute gain reduction in dB:
                float reduction_dB = (env_dB - currentThreshold_dB) * (1.0f - 1.0f / currentRatio);
                // Convert reduction dB to linear gain.
                gain = powf(10.0f, -reduction_dB / 20.0f);
            }
            else {
                gain = 1.0f;
            }
            // Apply gain reduction to the main signal.
            float duckedSample = mainSample * gain;
            // Blend the dry main signal with the ducked version.
            out[i] = (1.0f - currentMix) * mainSample + currentMix * duckedSample;
        }
        return 0;
    }

    // Control thread: allows updating parameters in real time.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomDuck] Enter parameters:" << endl;
                cout << "Threshold (dB), Ratio, Attack (ms), Release (ms), Mix (0.0-1.0)" << endl;
                cout << "e.g., \"-30 4.0 10 50 1.0\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newThreshold, newRatio, newAttack, newRelease, newMix;
            if (!(iss >> newThreshold >> newRatio >> newAttack >> newRelease >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomDuck] Invalid input. Please try again." << endl;
                continue;
            }
            // Clamp mix to [0, 1].
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            threshold_dB.store(newThreshold);
            ratio.store(newRatio);
            attackTime.store(newAttack);
            releaseTime.store(newRelease);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomDuck] Updated parameters:" << endl;
                cout << "  Threshold = " << newThreshold << " dB" << endl;
                cout << "  Ratio = " << newRatio << endl;
                cout << "  Attack = " << newAttack << " ms" << endl;
                cout << "  Release = " << newRelease << " ms" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomDuck(const char* client_name = "PhantomDuck")
        : client(nullptr), running(true), envelope(0.0f)
    {
        // Set default parameters.
        threshold_dB.store(-30.0f);
        ratio.store(4.0f);
        attackTime.store(10.0f);
        releaseTime.store(50.0f);
        mix.store(1.0f); // Fully ducked by default.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomDuck: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Register input ports.
        main_in_port = jack_port_register(client, "main", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        side_in_port = jack_port_register(client, "side", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        // Register mono output port.
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!main_in_port || !side_in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomDuck: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomDuck: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomDuck: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomDuck::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomDuck] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomDuck] Default parameters:" << endl;
        cout << "  Threshold = " << threshold_dB.load() << " dB" << endl;
        cout << "  Ratio = " << ratio.load() << endl;
        cout << "  Attack = " << attackTime.load() << " ms" << endl;
        cout << "  Release = " << releaseTime.load() << " ms" << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomDuck() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomDuck] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomDuck] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomDuck duck;
        duck.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomDuck] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

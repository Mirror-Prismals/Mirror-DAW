// PhantomGlitch.cpp
// A simple mono glitch (stutter) plugin using JACK and standard C++.
// It randomly triggers a stutter effect that captures a short segment of audio
// (of specified duration) and repeatedly replays it until the segment is finished.
// The output is a blend between the dry signal and the stuttered (glitched) signal.
// Real-time adjustable parameters:
//   - Stutter Duration (ms): Length of the stutter segment.
//   - Stutter Probability (per second): Chance of triggering a stutter.
//   - Mix: Blend between dry and stuttered signals (0.0 = dry, 1.0 = fully stuttered).
//
// Compile with:
//   g++ -std=c++11 PhantomGlitch.cpp -ljack -lpthread -o PhantomGlitch

#include <jack/jack.h>
#include <iostream>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <string>

using namespace std;

class PhantomGlitch {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Plugin parameters:
    atomic<float> stutterDuration_ms;   // Duration of the stutter segment in milliseconds.
    atomic<float> stutterProbability;   // Chance per second to trigger a stutter (0.0 to 1.0).
    atomic<float> mix;                  // Mix between dry and stuttered signal (0.0 = dry, 1.0 = full stutter).

    // Internal stutter state:
    bool stutterActive;                 // True if currently in stutter mode.
    vector<float> stutterBuffer;        // Holds the captured stutter segment.
    int stutterSamples;                 // Number of samples in the stutter segment.
    int stutterIndex;                   // Current read position in the stutter buffer.
    int stutterRemaining;               // Samples remaining to output from stutterBuffer.

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomGlitch* pg = static_cast<PhantomGlitch*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(pg->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(pg->out_port, nframes));

        // Compute per-sample stutter trigger probability.
        // If stutterProbability is, say, 0.5 per second, then per sample probability is (0.5 / sample_rate).
        float perSampleProb = pg->stutterProbability.load() / pg->sample_rate;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            float processed = dry; // Default: pass through.

            if (!pg->stutterActive) {
                // In normal mode, check if a stutter should trigger.
                float randVal = static_cast<float>(rand()) / RAND_MAX;  // Random in [0,1].
                if (randVal < perSampleProb) {
                    // Trigger stutter: capture stutterBuffer.
                    // Compute stutter duration in samples.
                    pg->stutterSamples = static_cast<int>(pg->stutterDuration_ms.load() * pg->sample_rate / 1000.0f);
                    if (pg->stutterSamples < 1) pg->stutterSamples = 1;
                    pg->stutterBuffer.resize(pg->stutterSamples);
                    // For simplicity, capture the current sample and the following (stutterSamples-1) samples from the input.
                    // If not enough samples remain in this process callback, fill the rest with dry.
                    for (int j = 0; j < pg->stutterSamples; j++) {
                        if (i + j < nframes)
                            pg->stutterBuffer[j] = in[i + j];
                        else
                            pg->stutterBuffer[j] = dry;
                    }
                    pg->stutterActive = true;
                    pg->stutterIndex = 0;
                    pg->stutterRemaining = pg->stutterSamples;
                    // Immediately use the captured sample.
                    processed = pg->stutterBuffer[pg->stutterIndex];
                    pg->stutterIndex = (pg->stutterIndex + 1) % pg->stutterSamples;
                    pg->stutterRemaining--;
                    // Skip processing additional samples from input for stutter if desired.
                    // (Here, we simply process the current sample as stutter and let subsequent samples continue to be processed normally.)
                }
            }
            else {
                // If a stutter is active, output the captured stutterBuffer.
                processed = pg->stutterBuffer[pg->stutterIndex];
                pg->stutterIndex = (pg->stutterIndex + 1) % pg->stutterSamples;
                pg->stutterRemaining--;
                if (pg->stutterRemaining <= 0) {
                    // End stutter mode.
                    pg->stutterActive = false;
                }
            }
            // Blend dry and processed signals.
            out[i] = (1.0f - pg->mix.load()) * dry + pg->mix.load() * processed;
        }
        return 0;
    }

    // Control thread: allow real-time updates via console.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomGlitch] Enter parameters: stutterDuration (ms), stutterProbability (per second, 0.0-1.0), mix (0.0-1.0)" << endl;
                cout << "e.g., \"100 0.3 1.0\" for 100 ms stutter, 30% chance per second, and full stutter effect, or 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newDuration, newProb, newMix;
            if (!(iss >> newDuration >> newProb >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomGlitch] Invalid input. Please try again." << endl;
                continue;
            }
            if (newDuration < 1.0f) newDuration = 1.0f;
            if (newProb < 0.0f) newProb = 0.0f;
            if (newProb > 1.0f) newProb = 1.0f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            stutterDuration_ms.store(newDuration);
            stutterProbability.store(newProb);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomGlitch] Updated parameters:" << endl;
                cout << "  Stutter Duration = " << newDuration << " ms" << endl;
                cout << "  Stutter Probability = " << newProb << " per second" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomGlitch(const char* client_name = "PhantomGlitch")
        : client(nullptr), running(true), stutterActive(false), stutterSamples(0),
        stutterIndex(0), stutterRemaining(0), mix(0.0f)
    {
        // Set default parameters.
        stutterDuration_ms.store(100.0f);  // 100 ms stutter duration.
        stutterProbability.store(0.3f);    // 30% chance per second.
        mix.store(1.0f);                   // Full stutter effect (100% processed).

        // Open JACK client.
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomGlitch: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);
        // Allocate a delay buffer if needed (here we don't need a long buffer because we capture the stutter from the current block).
        // However, we can allocate a minimal buffer to hold the current input block.
        // For simplicity, we won't use a global delay buffer—our stutter capture is done on-the-fly.

        // Register mono input and output ports.
        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomGlitch: Failed to register JACK ports");
        }
        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomGlitch: Failed to set process callback");
        }
        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomGlitch: Failed to activate JACK client");
        }

        // Initialize the stutter buffer (will be resized upon trigger).
        stutterBuffer.resize(0);

        // Launch control thread.
        control_thread = thread(&PhantomGlitch::control_loop, this);

        // Print initialization info.
        {
            lock_guard<mutex> lock(print_mutex);
            cout << "[PhantomGlitch] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
            cout << "[PhantomGlitch] Default parameters:" << endl;
            cout << "  Stutter Duration = " << stutterDuration_ms.load() << " ms" << endl;
            cout << "  Stutter Probability = " << stutterProbability.load() << " per second" << endl;
            cout << "  Mix = " << mix.load() << endl;
        }
    }

    ~PhantomGlitch() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomGlitch] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomGlitch] Shutting down." << endl;
    }
};

int main() {
    try {
        // Seed the random number generator.
        srand(static_cast<unsigned int>(time(nullptr)));
        PhantomGlitch glitch;
        glitch.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomGlitch] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

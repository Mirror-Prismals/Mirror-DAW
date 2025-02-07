// PhantomDeNoiser.cpp
// A simple mono noise reduction (de-noiser) plugin using JACK and standard C++.
// It adaptively estimates a noise floor from the input (when the signal is very low)
// and subtracts a scaled version of that noise estimate from the incoming signal.
// The output is blended with the dry signal using a mix parameter.
//
// Real-time adjustable parameters:
//   - Threshold (dB): The level below which the signal is considered noise (e.g., -60 dB).
//   - Reduction: A factor (0.0 to 1.0) specifying how much of the noise estimate to subtract.
//   - Learning Time (ms): Time constant for updating the noise estimate (e.g., 100 ms).
//   - Mix: Blend between dry and processed signals (0.0 = dry, 1.0 = fully noise-reduced).
//
// Compile with:
//   g++ -std=c++11 PhantomDeNoiser.cpp -ljack -lpthread -o PhantomDeNoiser

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

class PhantomDeNoiser {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Parameters:
    atomic<float> threshold_dB;  // Noise threshold in dB (e.g., -60 dB).
    atomic<float> reduction;     // Reduction factor (0.0 to 1.0).
    atomic<float> learningTime_ms; // Learning time for noise estimate (ms).
    atomic<float> mix;           // Mix between dry and processed (0.0 to 1.0).

    // Internal state:
    float noiseEstimate;  // Running noise floor estimate (in linear amplitude).

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomDeNoiser* nd = static_cast<PhantomDeNoiser*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(nd->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(nd->out_port, nframes));

        // Convert threshold from dB to linear.
        float currentThreshold_lin = powf(10.0f, nd->threshold_dB.load() / 20.0f);
        float currentReduction = nd->reduction.load();
        float currentLearningTime = nd->learningTime_ms.load();
        float currentMix = nd->mix.load();

        // Calculate dt in ms per sample.
        float dt_ms = 1000.0f / nd->sample_rate;
        // Compute exponential smoothing coefficient for noise estimation.
        // We use: alpha = exp(-dt / T)
        float alpha = expf(-dt_ms / currentLearningTime);

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float x = in[i];
            float absX = fabs(x);

            // Update noise estimate only if the current absolute value is below the noise threshold.
            if (absX < currentThreshold_lin) {
                nd->noiseEstimate = alpha * nd->noiseEstimate + (1.0f - alpha) * absX;
            }
            // Compute processed sample: subtract a scaled noise estimate.
            // Preserve the sign of x.
            float processed = x;
            if (x > 0)
                processed = x - currentReduction * nd->noiseEstimate;
            else if (x < 0)
                processed = x + currentReduction * nd->noiseEstimate;
            // Optionally, you might want to clamp the result to avoid inversion.
            // For this simple implementation, we leave it as is.
            // Blend processed with dry signal.
            out[i] = (1.0f - currentMix) * x + currentMix * processed;
        }
        return 0;
    }

    // Control thread: allows updating parameters in real time.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomDeNoiser] Enter parameters: threshold (dB), reduction (0.0-1.0), learning time (ms), mix (0.0-1.0)" << endl;
                cout << "e.g., \"-60 1.0 100 1.0\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newThreshold_dB, newReduction, newLearningTime, newMix;
            if (!(iss >> newThreshold_dB >> newReduction >> newLearningTime >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomDeNoiser] Invalid input. Please try again." << endl;
                continue;
            }
            // Optionally clamp mix between 0 and 1.
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            threshold_dB.store(newThreshold_dB);
            reduction.store(newReduction);
            learningTime_ms.store(newLearningTime);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomDeNoiser] Updated parameters:" << endl;
                cout << "  Threshold = " << newThreshold_dB << " dB" << endl;
                cout << "  Reduction = " << newReduction << endl;
                cout << "  Learning Time = " << newLearningTime << " ms" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomDeNoiser(const char* client_name = "PhantomDeNoiser")
        : client(nullptr), running(true), noiseEstimate(0.0f)
    {
        // Set default parameters.
        threshold_dB.store(-60.0f);    // Default threshold: -60 dB.
        reduction.store(1.0f);         // Default reduction: full subtraction.
        learningTime_ms.store(100.0f); // Default learning time: 100 ms.
        mix.store(1.0f);               // Fully processed by default (100% noise-reduced).

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomDeNoiser: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomDeNoiser: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomDeNoiser: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomDeNoiser: Failed to activate JACK client");
        }

        // Initialize the noise estimate.
        noiseEstimate = 0.0f;

        control_thread = thread(&PhantomDeNoiser::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomDeNoiser] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomDeNoiser] Default parameters:" << endl;
        cout << "  Threshold = " << threshold_dB.load() << " dB" << endl;
        cout << "  Reduction = " << reduction.load() << endl;
        cout << "  Learning Time = " << learningTime_ms.load() << " ms" << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomDeNoiser() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomDeNoiser] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomDeNoiser] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomDeNoiser deNoiser;
        deNoiser.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomDeNoiser] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

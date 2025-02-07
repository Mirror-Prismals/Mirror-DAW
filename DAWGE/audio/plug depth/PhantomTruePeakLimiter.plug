// PhantomTruePeakLimiter.cpp
// A simple mono true-peak limiter plugin using JACK and standard C++.
// It approximates inter-sample (true) peaks by oversampling each input sample (via linear interpolation),
// computes a desired gain to ensure that the peak does not exceed a specified ceiling (in dB),
// and then smooths this gain using attack and release time constants.
// The processed signal is optionally mixed with the dry input.
//
// Real-time adjustable parameters:
//   - Ceiling (dB): Desired output ceiling (e.g., 0 dB for 0 dBFS)
//   - Attack Time (ms): How quickly the limiter applies gain reduction.
//   - Release Time (ms): How quickly the limiter releases gain reduction.
//   - Mix: Blend between dry and limited signal (0.0 = dry, 1.0 = fully limited)
//
// Compile with:
//   g++ -std=c++11 PhantomTruePeakLimiter.cpp -ljack -lpthread -o PhantomTruePeakLimiter

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

class PhantomTruePeakLimiter {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Limiter parameters (all atomic for real-time control):
    // Ceiling in dB (e.g., 0.0 dB is unity; user can set, say, -0.5, 0.0, +1.0, etc.)
    atomic<float> ceiling_dB;
    // Attack and Release times in milliseconds.
    atomic<float> attackTime;  // e.g., 10 ms
    atomic<float> releaseTime; // e.g., 50 ms
    // Mix: 0.0 = dry, 1.0 = fully limited.
    atomic<float> mix;

    // Current gain (smoothed gain factor, linear). Not atomic; updated in real-time callback.
    float currentGain;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomTruePeakLimiter* limiter = static_cast<PhantomTruePeakLimiter*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(limiter->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(limiter->out_port, nframes));

        // Retrieve parameters.
        float ceiling_dB_val = limiter->ceiling_dB.load();
        // Convert ceiling from dB to linear. (0 dB -> 1.0; negative values yield lower ceilings.)
        float linearCeiling = powf(10.0f, ceiling_dB_val / 20.0f);
        float attTime = limiter->attackTime.load();    // in ms
        float relTime = limiter->releaseTime.load();     // in ms
        float mixVal = limiter->mix.load();

        // Compute dt in ms per sample.
        float dt_ms = 1000.0f / limiter->sample_rate;
        // Compute smoothing coefficients.
        // For attack: if desired gain < currentGain, use attack coefficient.
        // For release: if desired gain > currentGain, use release coefficient.
        // alpha = exp(-dt_ms / T)
        float attackCoeff = expf(-dt_ms / attTime);
        float releaseCoeff = expf(-dt_ms / relTime);

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];

            // Approximate true peak by oversampling between current sample and next sample.
            // For the last sample, we assume no interpolation.
            float maxVal = fabs(dry);
            if (i < nframes - 1) {
                float next = in[i + 1];
                // Compute three interpolated values.
                float v1 = dry + 0.25f * (next - dry);
                float v2 = dry + 0.50f * (next - dry);
                float v3 = dry + 0.75f * (next - dry);
                maxVal = fmax(maxVal, fabs(v1));
                maxVal = fmax(maxVal, fabs(v2));
                maxVal = fmax(maxVal, fabs(v3));
            }
            // Determine desired gain for this sample.
            float desiredGain = 1.0f;
            if (maxVal > linearCeiling && linearCeiling > 0.0f) {
                desiredGain = linearCeiling / maxVal;
            }
            // Smooth currentGain toward desiredGain.
            if (desiredGain < limiter->currentGain) {
                // Attack: reduce gain quickly.
                limiter->currentGain = attackCoeff * limiter->currentGain + (1.0f - attackCoeff) * desiredGain;
            }
            else {
                // Release: increase gain more slowly.
                limiter->currentGain = releaseCoeff * limiter->currentGain + (1.0f - releaseCoeff) * desiredGain;
            }
            // Apply currentGain to dry signal.
            float limitedSample = limiter->currentGain * dry;
            // Blend with dry signal using mix parameter.
            out[i] = (1.0f - mixVal) * dry + mixVal * limitedSample;
        }
        return 0;
    }

    // Control thread: read parameter updates from the console.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomTruePeakLimiter] Enter parameters:" << endl;
                cout << "Ceiling (dB), Attack Time (ms), Release Time (ms), Mix (0.0-1.0)" << endl;
                cout << "e.g., \"0.0 10 50 1.0\" (0.0 dB ceiling, 10 ms attack, 50 ms release, full limiting) or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newCeiling_dB, newAttack, newRelease, newMix;
            if (!(iss >> newCeiling_dB >> newAttack >> newRelease >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomTruePeakLimiter] Invalid input. Please try again." << endl;
                continue;
            }
            // Clamp mix.
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            ceiling_dB.store(newCeiling_dB);
            attackTime.store(newAttack);
            releaseTime.store(newRelease);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomTruePeakLimiter] Updated parameters:" << endl;
                cout << "  Ceiling = " << newCeiling_dB << " dB" << endl;
                cout << "  Attack Time = " << newAttack << " ms" << endl;
                cout << "  Release Time = " << newRelease << " ms" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomTruePeakLimiter(const char* client_name = "PhantomTruePeakLimiter")
        : client(nullptr), running(true), currentGain(1.0f)
    {
        // Set default parameters.
        ceiling_dB.store(0.0f);      // 0 dB ceiling (unity).
        attackTime.store(10.0f);     // 10 ms attack.
        releaseTime.store(50.0f);    // 50 ms release.
        mix.store(1.0f);             // Fully limited (100% wet).

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomTruePeakLimiter: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomTruePeakLimiter: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomTruePeakLimiter: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomTruePeakLimiter: Failed to activate JACK client");
        }

        // Initialize currentGain.
        currentGain = 1.0f;

        control_thread = thread(&PhantomTruePeakLimiter::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomTruePeakLimiter] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomTruePeakLimiter] Default parameters:" << endl;
        cout << "  Ceiling = " << ceiling_dB.load() << " dB" << endl;
        cout << "  Attack Time = " << attackTime.load() << " ms" << endl;
        cout << "  Release Time = " << releaseTime.load() << " ms" << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomTruePeakLimiter() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomTruePeakLimiter] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomTruePeakLimiter] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomTruePeakLimiter limiter;
        limiter.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomTruePeakLimiter] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

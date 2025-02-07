// PhantomPlateReverb.cpp
// A simple mono plate reverb plugin using JACK.
// It uses 4 parallel comb filters followed by 2 cascaded all-pass filters.
// The feedback of each comb filter is computed from its delay time and a user‑specified RT60 value.
// The output is a mix between the dry signal and the reverb (wet) signal.
// Real-time adjustable parameters: RT60 (seconds) and mix (0.0 = dry, 1.0 = fully wet).

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

using namespace std;

// ----------------------------
// Comb Filter Structure
struct CombFilter {
    vector<float> buffer;
    int bufferSize;
    int writeIndex;
    float delaySamples; // delay time in samples
    float feedback;     // feedback coefficient

    CombFilter(int delay_ms, int fs, float rt60) : writeIndex(0) {
        // Convert delay from ms to samples.
        delaySamples = (delay_ms * fs) / 1000.0f;
        // Use a buffer size slightly larger than delaySamples.
        bufferSize = static_cast<int>(delaySamples) + 1;
        buffer.resize(bufferSize, 0.0f);
        // Compute feedback coefficient from RT60.
        // feedback = 10^(-3 * delay / RT60)  where delay is in seconds.
        float delay_s = delaySamples / static_cast<float>(fs);
        feedback = pow(10.0, (-3.0 * delay_s) / rt60);
    }

    // Process one sample through the comb filter.
    float process(float input) {
        float output = buffer[writeIndex];
        buffer[writeIndex] = input + output * feedback;
        writeIndex = (writeIndex + 1) % bufferSize;
        return output;
    }

    // Update feedback based on a new RT60 (in seconds)
    void updateFeedback(float rt60, int fs) {
        float delay_s = delaySamples / static_cast<float>(fs);
        feedback = pow(10.0, (-3.0 * delay_s) / rt60);
    }
};

// ----------------------------
// All-Pass Filter Structure
struct AllPassFilter {
    vector<float> buffer;
    int bufferSize;
    int writeIndex;
    float feedback;   // feedback coefficient

    AllPassFilter(int delay_ms, int fs, float fb) : writeIndex(0), feedback(fb) {
        bufferSize = (delay_ms * fs) / 1000 + 1;
        buffer.resize(bufferSize, 0.0f);
    }

    // Process one sample.
    float process(float input) {
        float buffered = buffer[writeIndex];
        float output = -feedback * input + buffered;
        buffer[writeIndex] = input + output * feedback;
        writeIndex = (writeIndex + 1) % bufferSize;
        return output;
    }
};

// ----------------------------
// PhantomPlateReverb Class
class PhantomPlateReverb {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Reverb parameters.
    atomic<float> rt60;  // RT60 in seconds (e.g., 3.0 seconds)
    atomic<float> mix;   // Dry/Wet mix (0.0 = dry, 1.0 = fully wet)

    // Comb filters (4 in parallel).
    vector<CombFilter> combs;
    // Preset comb delays in ms (typical for plate reverb)
    const int combDelays_ms[4] = { 50, 56, 61, 68 };

    // All-pass filters (2 in series).
    vector<AllPassFilter> allpasses;
    // Preset all-pass delay in ms and feedback.
    const int allpassDelay_ms = 12;
    const float allpassFeedback = 0.7f;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomPlateReverb* reverb = static_cast<PhantomPlateReverb*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(reverb->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(reverb->out_port, nframes));

        float currentRT60 = reverb->rt60.load();
        float currentMix = reverb->mix.load();

        // Update comb filter feedbacks based on RT60.
        for (auto& comb : reverb->combs) {
            comb.updateFeedback(currentRT60, reverb->sample_rate);
        }

        // Process each sample.
        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            float combSum = 0.0f;
            // Process through all comb filters in parallel.
            for (auto& comb : reverb->combs) {
                combSum += comb.process(dry);
            }
            // Average the comb outputs.
            float combOutput = combSum / reverb->combs.size();
            // Pass through all all-pass filters in series.
            float apOutput = combOutput;
            for (auto& ap : reverb->allpasses) {
                apOutput = ap.process(apOutput);
            }
            // Mix with dry signal.
            out[i] = (1.0f - currentMix) * dry + currentMix * apOutput;
        }
        return 0;
    }

    // Control thread: allows real-time parameter updates.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomPlateReverb] Enter parameters: RT60 (sec) and mix (0.0-1.0)" << endl;
                cout << "e.g., \"3.0 0.7\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newRT60, newMix;
            if (!(iss >> newRT60 >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomPlateReverb] Invalid input. Please try again." << endl;
                continue;
            }
            if (newRT60 <= 0.0f) newRT60 = 0.1f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            rt60.store(newRT60);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomPlateReverb] Updated parameters:" << endl;
                cout << "  RT60 = " << newRT60 << " sec" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomPlateReverb(const char* client_name = "PhantomPlateReverb")
        : client(nullptr), running(true)
    {
        // Set default parameters.
        rt60.store(3.0f);  // 3 seconds decay.
        mix.store(0.7f);   // 70% wet.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomPlateReverb: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Initialize comb filters.
        for (int i = 0; i < 4; i++) {
            // Create a comb filter with delay from combDelays_ms and initial RT60.
            CombFilter comb(combDelays_ms[i], sample_rate, rt60.load());
            combs.push_back(comb);
        }

        // Initialize all-pass filters.
        // We'll use 2 all-pass filters in series.
        for (int i = 0; i < 2; i++) {
            AllPassFilter ap(allpassDelay_ms, sample_rate, allpassFeedback);
            allpasses.push_back(ap);
        }

        // Register mono JACK ports.
        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomPlateReverb: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomPlateReverb: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomPlateReverb: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomPlateReverb::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomPlateReverb] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomPlateReverb] Default parameters:" << endl;
        cout << "  RT60 = " << rt60.load() << " sec" << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomPlateReverb() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomPlateReverb] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomPlateReverb] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomPlateReverb reverb;
        reverb.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomPlateReverb] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

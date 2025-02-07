// PhantomTransientShaper.cpp
// A simple mono transient shaper plugin using JACK and standard C++.
// It detects transients by computing an envelope and its derivative, then applies a boost during attack
// and a reduction during sustain. The processed signal is blended with the dry signal via a mix parameter.
//
// Real-time adjustable parameters:
//   - Attack Time (ms): how quickly the envelope rises (e.g., 10 ms)
//   - Release Time (ms): how quickly the envelope falls (e.g., 50 ms)
//   - Attack Boost: gain factor applied during the attack (e.g., 2.0)
//   - Sustain Factor: gain factor applied during sustain (e.g., 0.8)
//   - Mix: blend between dry and processed signal (0.0 = dry, 1.0 = fully processed)
//
// Compile with:
//   g++ -std=c++11 PhantomTransientShaper.cpp -ljack -lpthread -o PhantomTransientShaper

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

class PhantomTransientShaper {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Transient shaper parameters (all real-time adjustable).
    atomic<float> attackTime;    // in ms (e.g., 10)
    atomic<float> releaseTime;   // in ms (e.g., 50)
    atomic<float> attackBoost;   // multiplier during attack (e.g., 2.0)
    atomic<float> sustainFactor; // multiplier during sustain (e.g., 0.8)
    atomic<float> mix;           // mix between dry and processed (0.0 to 1.0)

    // Envelope detector state.
    float envelope;      // current envelope
    float prevEnvelope;  // previous envelope value

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomTransientShaper* pts = static_cast<PhantomTransientShaper*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(pts->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(pts->out_port, nframes));

        // Retrieve parameters.
        float attTime = pts->attackTime.load();
        float relTime = pts->releaseTime.load();
        float atkBoost = pts->attackBoost.load();
        float susFactor = pts->sustainFactor.load();
        float currentMix = pts->mix.load();

        // Calculate dt in ms per sample.
        float dt = 1000.0f / pts->sample_rate;  // ms per sample

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Rectify input to get instantaneous amplitude.
            float rect = fabs(dry);
            // Update envelope with separate attack and release coefficients.
            if (rect > pts->envelope)
                pts->envelope = expf(-dt / attTime) * pts->envelope + (1.0f - expf(-dt / attTime)) * rect;
            else
                pts->envelope = expf(-dt / relTime) * pts->envelope + (1.0f - expf(-dt / relTime)) * rect;

            // Compute derivative of the envelope.
            float delta = pts->envelope - pts->prevEnvelope;
            pts->prevEnvelope = pts->envelope;

            // Decide gain multiplier based on the derivative.
            float gainMult = 1.0f;
            if (delta > 0) {
                // Rising envelope: transient detected; boost.
                gainMult = atkBoost;
            }
            else {
                // Falling or steady: reduce sustain.
                gainMult = susFactor;
            }
            // Apply gain multiplier to the input.
            float processed = dry * gainMult;
            // Mix processed with dry signal.
            out[i] = (1.0f - currentMix) * dry + currentMix * processed;
        }
        return 0;
    }

    // Control thread: listens for parameter updates via console.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomTransientShaper] Enter parameters:" << endl;
                cout << "Attack Time (ms), Release Time (ms), Attack Boost, Sustain Factor, Mix (0.0-1.0)" << endl;
                cout << "e.g., \"10 50 2.0 0.8 0.7\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newAttack, newRelease, newAtkBoost, newSusFactor, newMix;
            if (!(iss >> newAttack >> newRelease >> newAtkBoost >> newSusFactor >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomTransientShaper] Invalid input. Please try again." << endl;
                continue;
            }
            if (newAttack < 1.0f) newAttack = 1.0f;
            if (newRelease < 1.0f) newRelease = 1.0f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            attackTime.store(newAttack);
            releaseTime.store(newRelease);
            attackBoost.store(newAtkBoost);
            sustainFactor.store(newSusFactor);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomTransientShaper] Updated parameters:" << endl;
                cout << "  Attack Time = " << newAttack << " ms" << endl;
                cout << "  Release Time = " << newRelease << " ms" << endl;
                cout << "  Attack Boost = " << newAtkBoost << endl;
                cout << "  Sustain Factor = " << newSusFactor << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomTransientShaper(const char* client_name = "PhantomTransientShaper")
        : client(nullptr), running(true), envelope(0.0f), prevEnvelope(0.0f)
    {
        // Set default parameters.
        attackTime.store(10.0f);    // 10 ms attack.
        releaseTime.store(50.0f);   // 50 ms release.
        attackBoost.store(2.0f);    // Double the amplitude during attack.
        sustainFactor.store(0.8f);  // Reduce sustain to 80% of original.
        mix.store(0.7f);            // 70% processed signal.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomTransientShaper: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomTransientShaper: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomTransientShaper: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomTransientShaper: Failed to activate JACK client");
        }

        // Initialize envelope variables.
        envelope = 0.0f;
        prevEnvelope = 0.0f;

        control_thread = thread(&PhantomTransientShaper::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomTransientShaper] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomTransientShaper] Default parameters:" << endl;
        cout << "  Attack Time = " << attackTime.load() << " ms" << endl;
        cout << "  Release Time = " << releaseTime.load() << " ms" << endl;
        cout << "  Attack Boost = " << attackBoost.load() << endl;
        cout << "  Sustain Factor = " << sustainFactor.load() << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomTransientShaper() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomTransientShaper] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomTransientShaper] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomTransientShaper ts;
        ts.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomTransientShaper] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

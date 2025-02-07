// PhantomDynamicEQ.cpp
// A simple mono dynamic EQ plugin using JACK and standard C++.
// It splits the input into three bands (low, mid, high) using simple first-order filters,
// computes an envelope for each band, and if the envelope exceeds a specified threshold,
// reduces that band's gain according to a compression ratio. The processed bands are then
// summed and blended with the dry input using a mix parameter.
//
// Real-time adjustable parameters (via console):
//   Low Band:    Threshold (dB) and Ratio
//   Mid Band:    Threshold (dB) and Ratio
//   High Band:   Threshold (dB) and Ratio
//   Global:      Attack Time (ms), Release Time (ms), Mix (0.0 = dry, 1.0 = fully processed)
//
// Compile with:
//   g++ -std=c++11 PhantomDynamicEQ.cpp -ljack -lpthread -o PhantomDynamicEQ

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

// Simple first-order low-pass filter.
class LPF {
public:
    float cutoff;  // cutoff frequency in Hz.
    float a;       // smoothing coefficient.
    float prevY;
    LPF(float cutoffHz, int fs) : cutoff(cutoffHz), prevY(0.0f) {
        float dt = 1.0f / fs;
        float RC = 1.0f / (2.0f * M_PI * cutoff);
        a = dt / (RC + dt);
    }
    float process(float x) {
        float y = a * x + (1.0f - a) * prevY;
        prevY = y;
        return y;
    }
    void update(float cutoffHz, int fs) {
        cutoff = cutoffHz;
        float dt = 1.0f / fs;
        float RC = 1.0f / (2.0f * M_PI * cutoff);
        a = dt / (RC + dt);
    }
    void reset() { prevY = 0.0f; }
};

// Simple first-order high-pass filter.
class HPF {
public:
    float cutoff;  // cutoff frequency in Hz.
    float beta;    // smoothing coefficient.
    float prevX, prevY;
    HPF(float cutoffHz, int fs) : cutoff(cutoffHz), prevX(0.0f), prevY(0.0f) {
        float dt = 1.0f / fs;
        float RC = 1.0f / (2.0f * M_PI * cutoff);
        beta = RC / (RC + dt);
    }
    float process(float x) {
        float y = beta * (prevY + x - prevX);
        prevX = x;
        prevY = y;
        return y;
    }
    void update(float cutoffHz, int fs) {
        cutoff = cutoffHz;
        float dt = 1.0f / fs;
        float RC = 1.0f / (2.0f * M_PI * cutoff);
        beta = RC / (RC + dt);
    }
    void reset() { prevX = 0.0f; prevY = 0.0f; }
};

// PhantomDynamicEQ plugin class.
class PhantomDynamicEQ {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Per-band parameters (for Low, Mid, High bands):
    atomic<float> thresholdLow_dB;   // e.g., -30 dB.
    atomic<float> ratioLow;          // e.g., 2.0.
    atomic<float> thresholdMid_dB;   // e.g., -25 dB.
    atomic<float> ratioMid;          // e.g., 2.5.
    atomic<float> thresholdHigh_dB;  // e.g., -20 dB.
    atomic<float> ratioHigh;         // e.g., 3.0.

    // Global parameters.
    atomic<float> attackTime;    // in ms.
    atomic<float> releaseTime;   // in ms.
    atomic<float> mix;           // 0.0 (dry) to 1.0 (fully processed).

    // Filters for band splitting.
    LPF lowFilter;    // Low-pass filter for low band (default cutoff ~300 Hz).
    HPF highFilter;   // High-pass filter for high band (default cutoff ~3000 Hz).
    // For mid band, we'll compute: mid = input - (low + high).

    // Envelope variables for each band.
    float envLow;
    float envMid;
    float envHigh;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomDynamicEQ* deq = static_cast<PhantomDynamicEQ*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(deq->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(deq->out_port, nframes));

        // Get global parameters.
        float att = deq->attackTime.load();
        float rel = deq->releaseTime.load();
        float mixVal = deq->mix.load();

        // Precompute dt (ms per sample).
        float dt_ms = 1000.0f / deq->sample_rate;
        float attCoeff = expf(-dt_ms / att);
        float relCoeff = expf(-dt_ms / rel);

        // Convert band thresholds from dB to linear.
        auto dBToLinear = [](float dB) -> float {
            return powf(10.0f, dB / 20.0f);
        };
        float threshLow_lin = dBToLinear(deq->thresholdLow_dB.load());
        float threshMid_lin = dBToLinear(deq->thresholdMid_dB.load());
        float threshHigh_lin = dBToLinear(deq->thresholdHigh_dB.load());

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float x = in[i];
            // Split signal into bands.
            float low = deq->lowFilter.process(x);
            float high = deq->highFilter.process(x);
            float mid = x - (low + high);

            // Update envelopes.
            float absLow = fabs(low);
            if (absLow > deq->envLow)
                deq->envLow = attCoeff * deq->envLow + (1.0f - attCoeff) * absLow;
            else
                deq->envLow = relCoeff * deq->envLow + (1.0f - relCoeff) * absLow;

            float absMid = fabs(mid);
            if (absMid > deq->envMid)
                deq->envMid = attCoeff * deq->envMid + (1.0f - attCoeff) * absMid;
            else
                deq->envMid = relCoeff * deq->envMid + (1.0f - relCoeff) * absMid;

            float absHigh = fabs(high);
            if (absHigh > deq->envHigh)
                deq->envHigh = attCoeff * deq->envHigh + (1.0f - attCoeff) * absHigh;
            else
                deq->envHigh = relCoeff * deq->envHigh + (1.0f - relCoeff) * absHigh;

            // Compute gain reduction for each band.
            auto computeGain = [&](float env, float thresh, float ratio) -> float {
                if (env > thresh && thresh > 0)
                {
                    float desired = thresh + (env - thresh) / ratio;
                    return desired / env;
                }
                return 1.0f;
            };

            float gainLow = computeGain(deq->envLow, threshLow_lin, deq->ratioLow.load());
            float gainMid = computeGain(deq->envMid, threshMid_lin, deq->ratioMid.load());
            float gainHigh = computeGain(deq->envHigh, threshHigh_lin, deq->ratioHigh.load());

            // Apply gain reduction.
            float procLow = low * gainLow;
            float procMid = mid * gainMid;
            float procHigh = high * gainHigh;
            // Recombine processed bands.
            float procSignal = procLow + procMid + procHigh;
            // Mix with dry signal.
            out[i] = (1.0f - mixVal) * x + mixVal * procSignal;
        }
        return 0;
    }

    // Control thread: allows real-time updating of parameters.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomDynamicEQ] Enter parameters:" << endl;
                cout << "Low Threshold (dB), Low Ratio, Mid Threshold (dB), Mid Ratio, High Threshold (dB), High Ratio, Attack (ms), Release (ms), Mix (0.0-1.0)" << endl;
                cout << "e.g., \"-30 2.0 -25 2.5 -20 3.0 10 50 1.0\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float lt, lr, mt, mr, ht, hr, att, rel, m;
            if (!(iss >> lt >> lr >> mt >> mr >> ht >> hr >> att >> rel >> m)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomDynamicEQ] Invalid input. Please try again." << endl;
                continue;
            }
            if (m < 0.0f) m = 0.0f;
            if (m > 1.0f) m = 1.0f;
            thresholdLow_dB.store(lt);
            ratioLow.store(lr);
            thresholdMid_dB.store(mt);
            ratioMid.store(mr);
            thresholdHigh_dB.store(ht);
            ratioHigh.store(hr);
            attackTime.store(att);
            releaseTime.store(rel);
            mix.store(m);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomDynamicEQ] Updated parameters:" << endl;
                cout << "  Low:    Threshold = " << lt << " dB, Ratio = " << lr << endl;
                cout << "  Mid:    Threshold = " << mt << " dB, Ratio = " << mr << endl;
                cout << "  High:   Threshold = " << ht << " dB, Ratio = " << hr << endl;
                cout << "  Attack = " << att << " ms, Release = " << rel << " ms" << endl;
                cout << "  Mix    = " << m << endl;
            }
        }
    }

public:
    PhantomDynamicEQ(const char* client_name = "PhantomDynamicEQ")
        : client(nullptr), running(true),
          lowFilter(300.0f, 44100),  // default low cutoff at 300 Hz.
          highFilter(3000.0f, 44100), // default high cutoff at 3000 Hz.
          envLow(0.0f), envMid(0.0f), envHigh(0.0f)
    {
        // Set default per-band parameters.
        thresholdLow_dB.store(-30.0f);
        ratioLow.store(2.0f);
        thresholdMid_dB.store(-25.0f);
        ratioMid.store(2.5f);
        thresholdHigh_dB.store(-20.0f);
        ratioHigh.store(3.0f);
        // Global parameters.
        attackTime.store(10.0f);
        releaseTime.store(50.0f);
        mix.store(1.0f); // Fully processed by default.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomDynamicEQ: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Update filters with actual sample rate.
        lowFilter.update(300.0f, sample_rate);
        highFilter.update(3000.0f, sample_rate);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomDynamicEQ: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomDynamicEQ: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomDynamicEQ: Failed to activate JACK client");
        }

        // Initialize envelope values.
        envLow = envMid = envHigh = 0.0f;

        control_thread = thread(&PhantomDynamicEQ::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomDynamicEQ] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomDynamicEQ] Default parameters:" << endl;
        cout << "  Low:    Threshold = " << thresholdLow_dB.load() << " dB, Ratio = " << ratioLow.load() << endl;
        cout << "  Mid:    Threshold = " << thresholdMid_dB.load() << " dB, Ratio = " << ratioMid.load() << endl;
        cout << "  High:   Threshold = " << thresholdHigh_dB.load() << " dB, Ratio = " << ratioHigh.load() << endl;
        cout << "  Attack = " << attackTime.load() << " ms, Release = " << releaseTime.load() << " ms" << endl;
        cout << "  Mix    = " << mix.load() << endl;
    }

    ~PhantomDynamicEQ() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomDynamicEQ] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomDynamicEQ] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomDynamicEQ dynEQ;
        dynEQ.run();
    } catch (const exception& e) {
        cerr << "[PhantomDynamicEQ] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

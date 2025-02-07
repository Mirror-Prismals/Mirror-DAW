// PhantomDeEsser.cpp
// A simple mono de-esser plugin using JACK and standard C++ only.
// This plugin uses a first-order high-pass filter to isolate high frequencies,
// computes an envelope for the high band, and if that envelope exceeds a threshold,
// applies gain reduction (with a given ratio) to the high band before recombining it
// with the low frequencies. Now, a mix knob is added to blend the dry and processed signals.
//
// Real-time adjustable parameters:
//   - Cutoff Frequency (Hz): for the high-pass filter (default ~5000 Hz)
//   - Threshold (dB): level above which de-essing occurs (e.g., -30 dB)
//   - Ratio: compression ratio applied to the high-band (e.g., 2.0)
//   - Attack (ms): how fast the envelope rises
//   - Release (ms): how fast the envelope falls
//   - Mix: Dry/Wet mix (0.0 = completely dry, 1.0 = fully processed)
//
// Compile with:
//   g++ -std=c++11 PhantomDeEsser.cpp -ljack -lpthread -o PhantomDeEsser

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

// Helper: Convert dB to linear amplitude.
inline float dBToLinear(float dB) {
    return powf(10.0f, dB / 20.0f);
}

// A simple first-order high-pass filter.
// Implements: y[n] = β * (y[n-1] + x[n] - x[n-1])
struct HighPass {
    float x_prev;
    float y_prev;
    HighPass() : x_prev(0.0f), y_prev(0.0f) {}
    // Process one sample with cutoff frequency fc (Hz) and sample rate fs.
    float process(float x, float fc, int fs) {
        float dt = 1.0f / fs;
        float RC = 1.0f / (2.0f * M_PI * fc);
        float beta = RC / (RC + dt);
        float y = beta * (y_prev + x - x_prev);
        x_prev = x;
        y_prev = y;
        return y;
    }
    void reset() { x_prev = 0.0f; y_prev = 0.0f; }
};

class PhantomDeEsser {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // De-esser parameters:
    std::atomic<float> cutoffHz;      // High-pass filter cutoff frequency (Hz), e.g., 5000 Hz.
    std::atomic<float> threshold_dB;  // Threshold in dB (e.g., -30 dB).
    std::atomic<float> ratio;         // Compression ratio (e.g., 2.0).
    std::atomic<float> attackTime;    // Attack time in ms.
    std::atomic<float> releaseTime;   // Release time in ms.
    std::atomic<float> mix;           // Dry/Wet mix (0.0 = completely dry, 1.0 = fully processed).

    // High-pass filter state for envelope detection.
    HighPass hpFilter;
    // Envelope value.
    float envelope;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomDeEsser* deesser = static_cast<PhantomDeEsser*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(deesser->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(deesser->out_port, nframes));

        // Retrieve current parameters.
        float currentCutoff = deesser->cutoffHz.load();
        float currentThreshold_dB = deesser->threshold_dB.load();
        float currentThreshold = dBToLinear(currentThreshold_dB); // Linear threshold.
        float currentRatio = deesser->ratio.load();
        float currentAttack = deesser->attackTime.load();   // in ms
        float currentRelease = deesser->releaseTime.load(); // in ms
        float currentMix = deesser->mix.load();             // 0.0 to 1.0

        float dt = 1.0f / deesser->sample_rate;  // seconds per sample
        float dt_ms = dt * 1000.0f;              // milliseconds per sample

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float sample = in[i];
            // Apply high-pass filter to extract high frequencies.
            float highBand = deesser->hpFilter.process(sample, currentCutoff, deesser->sample_rate);
            // Compute absolute value of high band.
            float absHigh = fabs(highBand);
            // Update envelope with attack/release smoothing.
            if (absHigh > deesser->envelope)
                deesser->envelope = expf(-dt_ms / currentAttack) * deesser->envelope + (1.0f - expf(-dt_ms / currentAttack)) * absHigh;
            else
                deesser->envelope = expf(-dt_ms / currentRelease) * deesser->envelope + (1.0f - expf(-dt_ms / currentRelease)) * absHigh;

            // Determine gain reduction factor.
            float gain = 1.0f;
            if (deesser->envelope > currentThreshold && currentThreshold > 0) {
                // Compress the excess above threshold.
                float desired = currentThreshold + (deesser->envelope - currentThreshold) / currentRatio;
                gain = desired / deesser->envelope;
            }

            // Compute low band as the original minus the high band.
            float lowBand = sample - highBand;
            // Apply gain reduction to the high band.
            float processedHigh = gain * highBand;
            // Recombine to form the processed signal.
            float processed = lowBand + processedHigh;
            // Use the mix knob to blend dry and processed signals.
            out[i] = (1.0f - currentMix) * sample + currentMix * processed;
        }
        return 0;
    }

    // Control thread: update parameters in real time.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomDeEsser] Enter parameters: cutoff (Hz), threshold (dB), ratio, attack (ms), release (ms), mix (0-1)\n"
                    << "e.g., \"5000 -30 2.0 10 50 0.8\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newCutoff, newThreshold, newRatio, newAttack, newRelease, newMix;
            if (!(iss >> newCutoff >> newThreshold >> newRatio >> newAttack >> newRelease >> newMix)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomDeEsser] Invalid input. Please try again." << std::endl;
                continue;
            }
            if (newCutoff < 20.0f) newCutoff = 20.0f;
            if (newRatio < 1.0f) newRatio = 1.0f;
            if (newAttack < 1.0f) newAttack = 1.0f;
            if (newRelease < 1.0f) newRelease = 1.0f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            cutoffHz.store(newCutoff);
            threshold_dB.store(newThreshold);
            ratio.store(newRatio);
            attackTime.store(newAttack);
            releaseTime.store(newRelease);
            mix.store(newMix);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomDeEsser] Updated parameters: cutoff = " << newCutoff
                    << " Hz, threshold = " << newThreshold << " dB, ratio = " << newRatio
                    << ", attack = " << newAttack << " ms, release = " << newRelease
                    << " ms, mix = " << newMix << std::endl;
            }
        }
    }

public:
    PhantomDeEsser(const char* client_name = "PhantomDeEsser")
        : client(nullptr), running(true), envelope(0.0f)
    {
        // Set default parameters.
        cutoffHz.store(5000.0f);      // Default cutoff frequency at 5000 Hz.
        threshold_dB.store(-30.0f);   // Threshold at -30 dB.
        ratio.store(2.0f);            // Compression ratio 2:1.
        attackTime.store(10.0f);      // 10 ms attack.
        releaseTime.store(50.0f);     // 50 ms release.
        mix.store(0.8f);              // 80% processed signal.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomDeEsser: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomDeEsser: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomDeEsser: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomDeEsser: Failed to activate JACK client");
        }

        control_thread = std::thread(&PhantomDeEsser::control_loop, this);

        hpFilter.reset();

        std::lock_guard<std::mutex> lock(print_mutex);
        std::cout << "[PhantomDeEsser] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
        std::cout << "[PhantomDeEsser] Default parameters: cutoff = " << cutoffHz.load()
            << " Hz, threshold = " << threshold_dB.load() << " dB, ratio = " << ratio.load()
            << ", attack = " << attackTime.load() << " ms, release = " << releaseTime.load()
            << " ms, mix = " << mix.load() << std::endl;
    }

    ~PhantomDeEsser() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomDeEsser] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomDeEsser] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomDeEsser deesser;
        deesser.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomDeEsser] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

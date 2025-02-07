// PhantomAutoWah.cpp
// A simple mono auto-wah plugin using JACK.
// It computes an envelope from the input, maps that envelope to a cutoff frequency,
// updates a resonant band-pass biquad filter accordingly, and blends the filtered signal with the dry signal.
//
// Real-time adjustable parameters:
//   - Attack (ms)
//   - Release (ms)
//   - Min Cutoff (Hz)
//   - Max Cutoff (Hz)
//   - Q Factor
//   - Mix (0.0 = dry, 1.0 = fully processed)
//
// Compile with:
//   g++ -std=c++11 PhantomAutoWah.cpp -ljack -lpthread -o PhantomAutoWah

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

// Structure for a biquad filter.
struct Biquad {
    // Coefficients
    float b0, b1, b2, a1, a2;
    // State variables
    float x1, x2, y1, y2;
    Biquad() : b0(1), b1(0), b2(0), a1(0), a2(0),
        x1(0), x2(0), y1(0), y2(0) {
    }
    // Process a sample.
    inline float process(float x) {
        float y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1;
        x1 = x;
        y2 = y1;
        y1 = y;
        return y;
    }
    // Reset filter state.
    inline void reset() { x1 = x2 = y1 = y2 = 0; }
};

// Function to update biquad coefficients for a bandpass filter (constant skirt gain)
// using Robert Bristow-Johnson’s formulas.
// f0: center frequency (Hz), Q: quality factor, fs: sample rate.
void updateBandpass(Biquad& bq, float f0, float Q, int fs) {
    float w0 = 2.0f * M_PI * f0 / fs;
    float cosw0 = cosf(w0);
    float sinw0 = sinf(w0);
    float alpha = sinw0 / (2.0f * Q);
    float b0 = sinw0 / 2.0f;
    float b1 = 0.0f;
    float b2 = -sinw0 / 2.0f;
    float a0 = 1.0f + alpha;
    float a1 = -2.0f * cosw0;
    float a2 = 1.0f - alpha;
    bq.b0 = b0 / a0;
    bq.b1 = b1 / a0;
    bq.b2 = b2 / a0;
    bq.a1 = a1 / a0;
    bq.a2 = a2 / a0;
}

// Auto-wah plugin class.
class PhantomAutoWah {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;
    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Parameters (all in atomic variables for real-time updates).
    std::atomic<float> attackTime;    // ms
    std::atomic<float> releaseTime;   // ms
    std::atomic<float> minCutoff;     // Hz (e.g., 500 Hz)
    std::atomic<float> maxCutoff;     // Hz (e.g., 3000 Hz)
    std::atomic<float> QFactor;       // e.g., 2.0
    std::atomic<float> mix;           // 0.0 (dry) to 1.0 (fully processed)

    // Envelope for the input signal.
    float envelope;
    // Biquad filter for processing.
    Biquad bpFilter;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomAutoWah* wah = static_cast<PhantomAutoWah*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(wah->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(wah->out_port, nframes));

        float dt = 1.0f / wah->sample_rate;      // seconds per sample
        float dt_ms = dt * 1000.0f;              // milliseconds per sample

        // Retrieve current parameters.
        float att = wah->attackTime.load();
        float rel = wah->releaseTime.load();
        float fmin = wah->minCutoff.load();
        float fmax = wah->maxCutoff.load();
        float Q = wah->QFactor.load();
        float currentMix = wah->mix.load();

        // For each sample, update the envelope and filter coefficients.
        for (jack_nframes_t i = 0; i < nframes; i++) {
            float sample = in[i];
            float absSample = fabs(sample);
            // Update envelope with attack/release smoothing.
            if (absSample > wah->envelope)
                wah->envelope = expf(-dt_ms / att) * wah->envelope + (1.0f - expf(-dt_ms / att)) * absSample;
            else
                wah->envelope = expf(-dt_ms / rel) * wah->envelope + (1.0f - expf(-dt_ms / rel)) * absSample;
            // Clamp envelope to [0,1] (assuming input is normalized).
            if (wah->envelope > 1.0f)
                wah->envelope = 1.0f;

            // Map envelope to cutoff frequency.
            float fc = fmin + (fmax - fmin) * wah->envelope;

            // Update biquad coefficients for bandpass filter at fc.
            updateBandpass(wah->bpFilter, fc, Q, wah->sample_rate);

            // Process the sample through the bandpass filter.
            float filtered = wah->bpFilter.process(sample);

            // Mix dry and processed signals.
            out[i] = (1.0f - currentMix) * sample + currentMix * filtered;
        }
        return 0;
    }

    // Control thread: allows real-time adjustment of parameters.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomAutoWah] Enter parameters: attack (ms), release (ms), minCutoff (Hz), maxCutoff (Hz), Q, mix (0-1)\n"
                    << "e.g., \"10 50 500 3000 2.0 0.8\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newAttack, newRelease, newMinCutoff, newMaxCutoff, newQ, newMix;
            if (!(iss >> newAttack >> newRelease >> newMinCutoff >> newMaxCutoff >> newQ >> newMix)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomAutoWah] Invalid input. Please try again." << std::endl;
                continue;
            }
            if (newAttack < 1.0f) newAttack = 1.0f;
            if (newRelease < 1.0f) newRelease = 1.0f;
            if (newMinCutoff < 20.0f) newMinCutoff = 20.0f;
            if (newMaxCutoff < newMinCutoff) newMaxCutoff = newMinCutoff;
            if (newQ < 0.1f) newQ = 0.1f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            attackTime.store(newAttack);
            releaseTime.store(newRelease);
            minCutoff.store(newMinCutoff);
            maxCutoff.store(newMaxCutoff);
            QFactor.store(newQ);
            mix.store(newMix);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomAutoWah] Updated parameters:" << std::endl;
                std::cout << "  Attack = " << newAttack << " ms" << std::endl;
                std::cout << "  Release = " << newRelease << " ms" << std::endl;
                std::cout << "  Min Cutoff = " << newMinCutoff << " Hz" << std::endl;
                std::cout << "  Max Cutoff = " << newMaxCutoff << " Hz" << std::endl;
                std::cout << "  Q Factor = " << newQ << std::endl;
                std::cout << "  Mix = " << newMix << std::endl;
            }
        }
    }

public:
    PhantomAutoWah(const char* client_name = "PhantomAutoWah")
        : client(nullptr), running(true), envelope(0.0f)
    {
        // Set default parameters.
        attackTime.store(10.0f);      // 10 ms attack.
        releaseTime.store(50.0f);     // 50 ms release.
        minCutoff.store(500.0f);      // 500 Hz minimum cutoff.
        maxCutoff.store(3000.0f);     // 3000 Hz maximum cutoff.
        QFactor.store(2.0f);          // Q factor of 2.
        mix.store(0.8f);              // 80% processed signal.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomAutoWah: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomAutoWah: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomAutoWah: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomAutoWah: Failed to activate JACK client");
        }

        // Reset filter state.
        bpFilter.reset();

        control_thread = std::thread(&PhantomAutoWah::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomAutoWah] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomAutoWah] Default parameters:" << std::endl;
            std::cout << "  Attack = " << attackTime.load() << " ms" << std::endl;
            std::cout << "  Release = " << releaseTime.load() << " ms" << std::endl;
            std::cout << "  Min Cutoff = " << minCutoff.load() << " Hz" << std::endl;
            std::cout << "  Max Cutoff = " << maxCutoff.load() << " Hz" << std::endl;
            std::cout << "  Q Factor = " << QFactor.load() << std::endl;
            std::cout << "  Mix = " << mix.load() << std::endl;
        }
    }

    ~PhantomAutoWah() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomAutoWah] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomAutoWah] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomAutoWah autoWah;
        autoWah.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomAutoWah] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

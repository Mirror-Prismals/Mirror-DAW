// PhantomExciter.cpp
// A simple mono exciter (harmonic enhancer) plugin using JACK.
// It boosts the high frequencies via a high-shelf filter, then applies
// nonlinear saturation to add harmonics, and finally blends the processed
// signal with the dry input.
// Real-time adjustable parameters:
//   - Drive: Controls how hard the signal is pushed into the saturation stage.
//   - High-Shelf Gain (dB): Boost (or cut) applied in the high frequency range.
//   - Mix: Blend between dry and excited signals (0.0 = dry, 1.0 = fully processed).
//   - Output Gain (dB): Overall level of the processed signal (expressed in dB, e.g. -inf to +10 dB).
//
// Compile with:
//   g++ -std=c++11 PhantomExciter.cpp -ljack -lpthread -o PhantomExciter
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

class Biquad {
public:
    float b0, b1, b2, a1, a2;
    float x1, x2, y1, y2;
    Biquad() : b0(1), b1(0), b2(0), a1(0), a2(0),
        x1(0), x2(0), y1(0), y2(0) {
    }
    inline float process(float x) {
        float y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1;
        x1 = x;
        y2 = y1;
        y1 = y;
        return y;
    }
    inline void reset() { x1 = x2 = y1 = y2 = 0; }
};

// High-shelf filter coefficient update using RBJ’s formulas.
void updateHighShelf(Biquad& bq, float f0, float gain_dB, int fs) {
    float A = powf(10.0f, gain_dB / 40.0f);
    float w0 = 2.0f * M_PI * f0 / fs;
    float cosw0 = cosf(w0);
    float sinw0 = sinf(w0);
    float S = 1.0f;  // Shelf slope
    float alpha = sinw0 / 2.0f * sqrtf((A + 1.0f / A) * (1.0f / S - 1.0f) + 2.0f);

    float b0 = A * ((A + 1.0f) - (A - 1.0f) * cosw0 + 2.0f * sqrtf(A) * alpha);
    float b1 = 2.0f * A * ((A - 1.0f) - (A + 1.0f) * cosw0);
    float b2 = A * ((A + 1.0f) - (A - 1.0f) * cosw0 - 2.0f * sqrtf(A) * alpha);
    float a0 = (A + 1.0f) + (A - 1.0f) * cosw0 + 2.0f * sqrtf(A) * alpha;
    float a1 = -2.0f * ((A - 1.0f) + (A + 1.0f) * cosw0);
    float a2 = (A + 1.0f) + (A - 1.0f) * cosw0 - 2.0f * sqrtf(A) * alpha;

    bq.b0 = b0 / a0;
    bq.b1 = b1 / a0;
    bq.b2 = b2 / a0;
    bq.a1 = a1 / a0;
    bq.a2 = a2 / a0;
}

class PhantomExciter {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;
    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Parameters (atomic for real-time updates)
    std::atomic<float> drive;         // Drive factor for saturation (≥ 1.0), e.g., 2.0.
    std::atomic<float> hsGain_dB;       // High-shelf gain in dB (boost high frequencies).
    std::atomic<float> mix;           // Mix between dry and excited signal (0.0-1.0).
    // Now output gain is expressed in dB.
    std::atomic<float> outGain_dB;    // Output gain in dB (e.g., -10 to +10).

    // High-shelf filter for extracting high frequencies.
    Biquad hsFilter;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomExciter* exciter = static_cast<PhantomExciter*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(exciter->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(exciter->out_port, nframes));

        // Retrieve parameters.
        float currentDrive = exciter->drive.load();
        float currentHsGain_dB = exciter->hsGain_dB.load();
        float currentMix = exciter->mix.load();
        float currentOutGain_dB = exciter->outGain_dB.load();
        // Convert output gain from dB to linear.
        float currentOutGain = powf(10.0f, currentOutGain_dB / 20.0f);

        // Update high-shelf filter coefficients.
        // We'll set a cutoff frequency for the high-shelf filter, e.g., 3000 Hz.
        float cutoff = 3000.0f;
        updateHighShelf(exciter->hsFilter, cutoff, currentHsGain_dB, exciter->sample_rate);

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Apply drive and then soft clipping using tanh.
            float driven = currentDrive * dry;
            float saturated = tanhf(driven);
            // Process through high-shelf filter to emphasize high frequencies.
            float excited = exciter->hsFilter.process(saturated);
            // Apply output gain (converted from dB).
            excited *= currentOutGain;
            // Mix with dry signal.
            out[i] = (1.0f - currentMix) * dry + currentMix * excited;
        }
        return 0;
    }

    // Control thread: allows real-time adjustment of parameters.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomExciter] Enter parameters: drive, high-shelf gain (dB), mix (0-1), output gain (dB)\n"
                    << "e.g., \"2.0 6.0 0.7 0.0\" (0.0 dB is unity) or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newDrive, newHsGain, newMix, newOutGain_dB;
            if (!(iss >> newDrive >> newHsGain >> newMix >> newOutGain_dB)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomExciter] Invalid input. Please try again." << std::endl;
                continue;
            }
            if (newDrive < 1.0f) newDrive = 1.0f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            drive.store(newDrive);
            hsGain_dB.store(newHsGain);
            mix.store(newMix);
            outGain_dB.store(newOutGain_dB);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomExciter] Updated parameters:" << std::endl;
                std::cout << "  Drive = " << newDrive << std::endl;
                std::cout << "  High-Shelf Gain = " << newHsGain << " dB" << std::endl;
                std::cout << "  Mix = " << newMix << std::endl;
                std::cout << "  Output Gain = " << newOutGain_dB << " dB" << std::endl;
            }
        }
    }

public:
    PhantomExciter(const char* client_name = "PhantomExciter")
        : client(nullptr), running(true)
    {
        // Set default parameters.
        drive.store(2.0f);         // Default drive: 2.0
        hsGain_dB.store(6.0f);       // Default high-shelf gain: +6 dB boost
        mix.store(0.7f);           // 70% processed signal
        outGain_dB.store(0.0f);    // 0.0 dB = unity gain (can be set from, say, -10 dB to +10 dB)

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomExciter: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomExciter: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomExciter: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomExciter: Failed to activate JACK client");
        }

        hsFilter.reset();

        control_thread = std::thread(&PhantomExciter::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomExciter] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomExciter] Default parameters:" << std::endl;
            std::cout << "  Drive = " << drive.load() << std::endl;
            std::cout << "  High-Shelf Gain = " << hsGain_dB.load() << " dB" << std::endl;
            std::cout << "  Mix = " << mix.load() << std::endl;
            std::cout << "  Output Gain = " << outGain_dB.load() << " dB" << std::endl;
        }
    }

    ~PhantomExciter() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomExciter] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomExciter] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomExciter exciter;
        exciter.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomExciter] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

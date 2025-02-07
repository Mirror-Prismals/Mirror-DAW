// PhantomResonator.cpp
// A simple real-time stereo resonator using JACK and a biquad bandpass filter.
// The resonator emphasizes a narrow band of frequencies (set by the resonant frequency and Q factor)
// and blends the resonated (wet) signal with the original dry signal.
// Parameters (adjustable in real time):
//   - Resonant Frequency (Hz)
//   - Q Factor
//   - Mix (0.0 = dry, 1.0 = fully resonated)
//   - Resonator Gain (linear multiplier)
//
// Compile with:
//   g++ -std=c++11 PhantomResonator.cpp -ljack -lpthread -o PhantomResonator

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

// --- Biquad Filter Class ---
class Biquad {
public:
    float b0, b1, b2, a1, a2;
    float x1, x2, y1, y2;
    Biquad() : b0(1), b1(0), b2(0), a1(0), a2(0),
        x1(0), x2(0), y1(0), y2(0) {
    }
    // Process one sample.
    float process(float x) {
        float y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1; x1 = x;
        y2 = y1; y1 = y;
        return y;
    }
    // Update coefficients for a bandpass filter (constant skirt gain) using RBJ's formulas.
    // f0: center frequency (Hz), Q: quality factor, fs: sample rate.
    void updateBandpass(float f0, float Q, int fs) {
        float w0 = 2.0f * M_PI * f0 / fs;
        float alpha = sinf(w0) / (2.0f * Q);
        float cosw0 = cosf(w0);
        float b0_un = sinf(w0) / 2.0f;
        float b1_un = 0.0f;
        float b2_un = -sinf(w0) / 2.0f;
        float a0 = 1.0f + alpha;
        float a1_un = -2.0f * cosw0;
        float a2_un = 1.0f - alpha;
        b0 = b0_un / a0;
        b1 = b1_un / a0;
        b2 = b2_un / a0;
        a1 = a1_un / a0;
        a2 = a2_un / a0;
    }
    // Reset filter state.
    void reset() {
        x1 = x2 = y1 = y2 = 0;
    }
};

// --- PhantomResonator Class ---
class PhantomResonator {
private:
    jack_client_t* client;
    jack_port_t* in_left, * in_right, * out_left, * out_right;
    int sample_rate;
    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Resonator parameters.
    std::atomic<float> resFreq;   // Resonant frequency in Hz.
    std::atomic<float> Q;         // Quality factor.
    std::atomic<float> mix;       // 0.0 = dry, 1.0 = fully resonated.
    std::atomic<float> resGain;   // Gain applied to the resonated signal.

    // Biquad filters (one per channel).
    Biquad leftBiquad;
    Biquad rightBiquad;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomResonator* res = static_cast<PhantomResonator*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(res->in_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(res->in_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(res->out_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(res->out_right, nframes));

        // Read current parameters.
        float currentFreq = res->resFreq.load();
        float currentQ = res->Q.load();
        float currentMix = res->mix.load();
        float currentResGain = res->resGain.load();

        // Update filter coefficients for both channels.
        res->leftBiquad.updateBandpass(currentFreq, currentQ, res->sample_rate);
        res->rightBiquad.updateBandpass(currentFreq, currentQ, res->sample_rate);

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dryL = inL[i];
            float dryR = inR[i];
            // Process through the bandpass filter.
            float resL = res->leftBiquad.process(dryL);
            float resR = res->rightBiquad.process(dryR);
            // Apply resonator gain.
            resL *= currentResGain;
            resR *= currentResGain;
            // Mix dry and resonated signals.
            outL[i] = (1.0f - currentMix) * dryL + currentMix * resL;
            outR[i] = (1.0f - currentMix) * dryR + currentMix * resR;
        }
        return 0;
    }

    // Control thread: update parameters in real time.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomResonator] Enter parameters: resonant frequency (Hz), Q factor, mix (0-1), resonator gain (linear)\n"
                    << "e.g., \"500 10 0.5 1.0\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newFreq, newQ, newMix, newResGain;
            if (!(iss >> newFreq >> newQ >> newMix >> newResGain)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomResonator] Invalid input. Please try again." << std::endl;
                continue;
            }
            if (newFreq <= 0) newFreq = 1.0f;
            if (newQ <= 0) newQ = 0.1f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            resFreq.store(newFreq);
            Q.store(newQ);
            mix.store(newMix);
            resGain.store(newResGain);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomResonator] Updated parameters: resonant frequency = " << newFreq
                    << " Hz, Q factor = " << newQ
                    << ", mix = " << newMix
                    << ", resonator gain = " << newResGain << std::endl;
            }
        }
    }

public:
    PhantomResonator(const char* client_name = "PhantomResonator")
        : client(nullptr), running(true)
    {
        // Set default parameters.
        resFreq.store(500.0f);   // 500 Hz resonant frequency.
        Q.store(10.0f);          // High Q for a narrow, pronounced resonance.
        mix.store(0.5f);         // 50% wet/dry mix.
        resGain.store(1.0f);     // Unity gain.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomResonator: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_left = jack_port_register(client, "in_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        in_right = jack_port_register(client, "in_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_left || !in_right || !out_left || !out_right) {
            jack_client_close(client);
            throw std::runtime_error("PhantomResonator: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomResonator: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomResonator: Failed to activate JACK client");
        }

        control_thread = std::thread(&PhantomResonator::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomResonator] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomResonator] Default parameters: resonant frequency = " << resFreq.load()
                << " Hz, Q factor = " << Q.load() << ", mix = " << mix.load()
                << ", resonator gain = " << resGain.load() << std::endl;
        }
    }

    ~PhantomResonator() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomResonator] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomResonator] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomResonator resonator;
        resonator.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomResonator] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

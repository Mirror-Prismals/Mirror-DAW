// PhantomRingMod.cpp
// A simple mono ring modulator effect using JACK.
// It multiplies the input signal by a sine wave oscillator at a given modulation frequency,
// and blends the modulated signal with the original dry signal according to a mix parameter.
//
// Real-time adjustable parameters:
//   - Modulation Frequency (Hz): the frequency of the sine wave oscillator.
//   - Mix: dry/wet mix (0.0 = completely dry, 1.0 = fully modulated).
//
// Compile with:
//   g++ -std=c++11 PhantomRingMod.cpp -ljack -lpthread -o PhantomRingMod

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

class PhantomRingMod {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Ring modulator parameters.
    std::atomic<float> modFrequency; // in Hz (carrier frequency)
    std::atomic<float> mix;          // dry/wet mix (0.0 = dry, 1.0 = fully modulated)

    // LFO phase accumulator.
    float lfoPhase;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomRingMod* ringMod = static_cast<PhantomRingMod*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(ringMod->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(ringMod->out_port, nframes));

        // Retrieve current parameters.
        float currentFreq = ringMod->modFrequency.load();
        float currentMix = ringMod->mix.load();

        // Calculate LFO phase increment per sample.
        float dt = 1.0f / ringMod->sample_rate;
        float lfoInc = 2.0f * M_PI * currentFreq * dt;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Compute modulated (carrier) value.
            float modSignal = sinf(ringMod->lfoPhase);
            // Ring modulate: multiply dry signal by modulator.
            float modulated = dry * modSignal;
            // Mix dry and modulated signals.
            out[i] = (1.0f - currentMix) * dry + currentMix * modulated;

            // Advance LFO phase.
            ringMod->lfoPhase += lfoInc;
            if (ringMod->lfoPhase >= 2.0f * M_PI)
                ringMod->lfoPhase -= 2.0f * M_PI;
        }
        return 0;
    }

    // Control thread: allows real-time adjustment of modulator frequency and mix.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomRingMod] Enter parameters: modFrequency (Hz) and mix (0.0-1.0)\n"
                    << "e.g., \"100 0.7\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newFreq, newMix;
            if (!(iss >> newFreq >> newMix)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomRingMod] Invalid input. Please try again." << std::endl;
                continue;
            }
            if (newFreq < 0.0f) newFreq = 0.0f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            modFrequency.store(newFreq);
            mix.store(newMix);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomRingMod] Updated parameters: modFrequency = " << newFreq
                    << " Hz, mix = " << newMix << std::endl;
            }
        }
    }

public:
    PhantomRingMod(const char* client_name = "PhantomRingMod")
        : client(nullptr), running(true), lfoPhase(0.0f)
    {
        // Set default parameters.
        modFrequency.store(100.0f); // Default modulation frequency 100 Hz.
        mix.store(0.7f);            // Default mix: 70% modulated signal.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomRingMod: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomRingMod: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomRingMod: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomRingMod: Failed to activate JACK client");
        }

        control_thread = std::thread(&PhantomRingMod::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomRingMod] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomRingMod] Default parameters: modFrequency = " << modFrequency.load()
                << " Hz, mix = " << mix.load() << std::endl;
        }
    }

    ~PhantomRingMod() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomRingMod] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomRingMod] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomRingMod ringMod;
        ringMod.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomRingMod] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

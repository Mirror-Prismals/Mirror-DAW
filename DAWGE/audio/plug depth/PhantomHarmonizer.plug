// PhantomHarmonizer.cpp
// A toy mono harmonizer pedal using JACK and standard C++.
// It pitch shifts the incoming signal by a fixed semitone interval and mixes it with the dry signal.
// The pitch shift is achieved by reading from a circular delay buffer at a rate scaled by the pitch shift ratio.
//
// Compile with:
//   g++ -std=c++11 PhantomHarmonizer.cpp -ljack -lpthread -o PhantomHarmonizer

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

// Utility: Convert semitones to pitch shift ratio.
inline float semitonesToRatio(float semitones) {
    return powf(2.0f, semitones / 12.0f);
}

class PhantomHarmonizer {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Harmonizer parameters.
    std::atomic<float> semitoneShift;  // in semitones (e.g., 4.0 for major third up)
    std::atomic<float> mix;            // dry/wet mix (0.0 = dry, 1.0 = full harmony)
    std::atomic<float> baseDelay_ms;   // base delay in ms (controls latency), e.g., 20 ms

    // Derived parameters.
    float pitchRatio;          // 2^(semitoneShift/12)
    float baseDelaySamples;    // baseDelay_ms * sample_rate / 1000

    // Delay buffer for pitch shifting.
    std::vector<float> delayBuffer;
    size_t bufferSize;         // size of delayBuffer (in samples)
    size_t writeIndex;         // integer write pointer (modulo bufferSize)
    float virtualReadIndex;    // floating-point read pointer

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomHarmonizer* ph = static_cast<PhantomHarmonizer*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(ph->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(ph->out_port, nframes));

        // Update derived parameters.
        float currentSemitone = ph->semitoneShift.load();
        ph->pitchRatio = semitonesToRatio(currentSemitone);
        float currentBaseDelay_ms = ph->baseDelay_ms.load();
        ph->baseDelaySamples = currentBaseDelay_ms * ph->sample_rate / 1000.0f;

        float currentMix = ph->mix.load();

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];

            // Write current sample to delay buffer.
            ph->delayBuffer[ph->writeIndex] = dry;

            // Read the pitch-shifted sample from the delay buffer using the virtual read pointer.
            // Use linear interpolation.
            size_t index0 = static_cast<size_t>(ph->virtualReadIndex) % ph->bufferSize;
            size_t index1 = (index0 + 1) % ph->bufferSize;
            float frac = ph->virtualReadIndex - floorf(ph->virtualReadIndex);
            float shiftedSample = (1.0f - frac) * ph->delayBuffer[index0] + frac * ph->delayBuffer[index1];

            // Mix the dry and pitch-shifted signals.
            out[i] = (1.0f - currentMix) * dry + currentMix * shiftedSample;

            // Increment write pointer.
            ph->writeIndex = (ph->writeIndex + 1) % ph->bufferSize;

            // Increment the virtual read pointer by pitchRatio.
            ph->virtualReadIndex += ph->pitchRatio;

            // Check the delay (the difference between writeIndex and virtualReadIndex, modulo buffer size).
            // We want the delay to remain roughly equal to baseDelaySamples.
            float delay = static_cast<float>(ph->writeIndex) - ph->virtualReadIndex;
            if (delay < 0)
                delay += ph->bufferSize;
            float tolerance = 5.0f; // samples tolerance.
            if (delay > ph->baseDelaySamples + tolerance) {
                // Reinitialize virtualReadIndex to maintain the desired delay.
                ph->virtualReadIndex = static_cast<float>(ph->writeIndex) - ph->baseDelaySamples;
                if (ph->virtualReadIndex < 0)
                    ph->virtualReadIndex += ph->bufferSize;
            }
        }
        return 0;
    }

    // Control thread: update harmonizer parameters in real time.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomHarmonizer] Enter parameters: semitone shift (e.g., 4.0), mix (0-1), base delay (ms) (e.g., 20)\n"
                    << "e.g., \"4.0 0.5 20\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newSemitones, newMix, newBaseDelay;
            if (!(iss >> newSemitones >> newMix >> newBaseDelay)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomHarmonizer] Invalid input. Please try again." << std::endl;
                continue;
            }
            semitoneShift.store(newSemitones);
            mix.store(newMix);
            baseDelay_ms.store(newBaseDelay);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomHarmonizer] Updated parameters: semitone shift = " << newSemitones
                    << " semitones, mix = " << newMix << ", base delay = " << newBaseDelay << " ms" << std::endl;
            }
        }
    }

public:
    PhantomHarmonizer(const char* client_name = "PhantomHarmonizer")
        : client(nullptr), running(true), writeIndex(0), virtualReadIndex(0.0f)
    {
        // Set default parameters.
        semitoneShift.store(4.0f);  // Major third up.
        mix.store(0.5f);            // 50% mix.
        baseDelay_ms.store(20.0f);   // 20 ms base delay.

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomHarmonizer: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Choose a delay buffer size that provides headroom. For a 20 ms delay at 44100 Hz, that's about 882 samples.
        // We choose a fixed size (e.g., 2048 samples) for simplicity.
        bufferSize = 2048;
        delayBuffer.resize(bufferSize, 0.0f);

        // Initialize write pointer and virtual read pointer.
        writeIndex = 0;
        float initDelaySamples = baseDelay_ms.load() * sample_rate / 1000.0f;
        virtualReadIndex = static_cast<float>(writeIndex) - initDelaySamples;
        if (virtualReadIndex < 0)
            virtualReadIndex += bufferSize;

        // Register JACK ports.
        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomHarmonizer: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomHarmonizer: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomHarmonizer: Failed to activate JACK client");
        }

        control_thread = std::thread(&PhantomHarmonizer::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomHarmonizer] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomHarmonizer] Default parameters: semitone shift = " << semitoneShift.load()
                << " semitones, mix = " << mix.load() << ", base delay = " << baseDelay_ms.load() << " ms" << std::endl;
        }
    }

    ~PhantomHarmonizer() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomHarmonizer] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomHarmonizer] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomHarmonizer harmonizer;
        harmonizer.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomHarmonizer] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

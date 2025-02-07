// PhantomCrusher.cpp
// A simple real-time bitcrusher plugin using JACK.
//
// This plugin reduces both the bit depth and the effective sample rate
// of the input signal to achieve a crunchy lo-fi effect.
// Parameters:
//   - Bit Depth: integer (1 to 16) controlling quantization resolution.
//   - Reduction Factor: integer >= 1 that holds each processed sample for N frames.
//   - Mix: dry/wet blend (0.0 = dry, 1.0 = fully bitcrushed).
//
// Compile with:
//   g++ -std=c++11 PhantomCrusher.cpp -ljack -lpthread -o PhantomCrusher

#include <jack/jack.h>
#include <iostream>
#include <atomic>
#include <thread>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <chrono>

// Helper function to quantize a sample (assumed to be in [-1,1])
// given a bit depth (1 to 16).
float quantizeSample(float x, int bitDepth) {
    // Clamp bitDepth to [1,16]
    if (bitDepth < 1) bitDepth = 1;
    if (bitDepth > 16) bitDepth = 16;
    int levels = (1 << bitDepth) - 1; // 2^bitDepth - 1
    float step = 2.0f / levels;        // Range is 2 (-1 to 1)
    float quantized = std::round((x + 1.0f) / step) * step - 1.0f;
    // Clamp to [-1, 1]
    if (quantized < -1.0f) quantized = -1.0f;
    if (quantized > 1.0f) quantized = 1.0f;
    return quantized;
}

class PhantomCrusher {
private:
    jack_client_t* client;
    jack_port_t* in_left, * in_right, * out_left, * out_right;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Bitcrusher parameters:
    std::atomic<int> bitDepth;            // e.g., default 16 (no reduction) down to lower values.
    std::atomic<int> reductionFactor;     // e.g., 1 = no sample rate reduction, 2,3,...
    std::atomic<float> mix;               // 0.0 (dry) to 1.0 (fully processed)

    // For sample rate reduction: for each channel, hold the processed sample
    // and count how many frames have passed.
    int leftCounter;
    int rightCounter;
    float leftHeldSample;
    float rightHeldSample;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomCrusher* crusher = static_cast<PhantomCrusher*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(crusher->in_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(crusher->in_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(crusher->out_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(crusher->out_right, nframes));

        // Read parameters atomically.
        int currentBitDepth = crusher->bitDepth.load();
        int currentReduction = crusher->reductionFactor.load();
        float currentMix = crusher->mix.load();

        // Process left channel.
        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = inL[i];
            float processed;
            // If counter is zero, compute a new quantized value.
            if (crusher->leftCounter == 0) {
                processed = quantizeSample(dry, currentBitDepth);
                crusher->leftHeldSample = processed;
            }
            else {
                processed = crusher->leftHeldSample;
            }
            crusher->leftCounter++;
            if (crusher->leftCounter >= currentReduction)
                crusher->leftCounter = 0;

            // Mix dry and processed signals.
            outL[i] = currentMix * processed + (1.0f - currentMix) * dry;
        }

        // Process right channel.
        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = inR[i];
            float processed;
            if (crusher->rightCounter == 0) {
                processed = quantizeSample(dry, currentBitDepth);
                crusher->rightHeldSample = processed;
            }
            else {
                processed = crusher->rightHeldSample;
            }
            crusher->rightCounter++;
            if (crusher->rightCounter >= currentReduction)
                crusher->rightCounter = 0;

            outR[i] = currentMix * processed + (1.0f - currentMix) * dry;
        }
        return 0;
    }

    // Control thread: allows real-time parameter adjustment.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomCrusher] Enter parameters: bitDepth (1-16), reductionFactor (>=1), mix (0.0-1.0)\n"
                    << "e.g., \"8 4 0.7\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            int newBitDepth, newReduction;
            float newMix;
            if (!(iss >> newBitDepth >> newReduction >> newMix)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomCrusher] Invalid input. Please try again." << std::endl;
                continue;
            }
            // Sanity checks.
            if (newBitDepth < 1) newBitDepth = 1;
            if (newBitDepth > 16) newBitDepth = 16;
            if (newReduction < 1) newReduction = 1;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            bitDepth.store(newBitDepth);
            reductionFactor.store(newReduction);
            mix.store(newMix);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomCrusher] Updated parameters: bitDepth = " << newBitDepth
                    << ", reductionFactor = " << newReduction
                    << ", mix = " << newMix << std::endl;
            }
        }
    }

public:
    PhantomCrusher(const char* client_name = "PhantomCrusher")
        : client(nullptr), leftCounter(0), rightCounter(0),
        leftHeldSample(0.0f), rightHeldSample(0.0f), running(true)
    {
        // Set default parameters.
        bitDepth.store(16);          // Default: no bit reduction.
        reductionFactor.store(1);    // Default: no sample rate reduction.
        mix.store(1.0f);             // Fully processed (bitcrushed).

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomCrusher: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_left = jack_port_register(client, "in_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        in_right = jack_port_register(client, "in_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_left || !in_right || !out_left || !out_right) {
            jack_client_close(client);
            throw std::runtime_error("PhantomCrusher: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomCrusher: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomCrusher: Failed to activate JACK client");
        }

        // Start the control thread.
        control_thread = std::thread(&PhantomCrusher::control_loop, this);

        // Initialize counters.
        leftCounter = 0;
        rightCounter = 0;

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomCrusher] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomCrusher] Default parameters: bitDepth = " << bitDepth.load()
                << ", reductionFactor = " << reductionFactor.load()
                << ", mix = " << mix.load() << std::endl;
        }
    }

    ~PhantomCrusher() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomCrusher] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomCrusher] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomCrusher crusher;
        crusher.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomCrusher] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

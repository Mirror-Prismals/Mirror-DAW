// PhantomVibrato.cpp
// A simple real-time stereo vibrato effect using JACK.
//
// This plugin implements vibrato by modulating a delay line for each channel.
// The LFO modulates the delay time (in ms) around a base delay, and linear
// interpolation is used to read the delayed sample. The wet (modulated) signal
// is mixed with the dry signal according to a mix parameter.
//
// Compile with:
//   g++ -std=c++11 PhantomVibrato.cpp -ljack -lpthread -o PhantomVibrato

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

// We'll set a maximum delay (in ms) for our delay buffers.
// For vibrato the modulation is subtle, so 50 ms is plenty.
const float MAX_DELAY_MS = 50.0f;

class PhantomVibrato {
private:
    jack_client_t* client;
    // Stereo ports
    jack_port_t* in_left, * in_right, * out_left, * out_right;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Vibrato parameters:
    // rate: LFO frequency in Hz (default: 5.0 Hz)
    // depth: modulation depth in ms (default: 2.0 ms)
    // baseDelay: base delay in ms (default: 5.0 ms) around which modulation occurs
    // mix: dry/wet mix (0.0 = fully dry, 1.0 = fully vibrato-processed) (default: 1.0)
    std::atomic<float> rate;       // in Hz
    std::atomic<float> depth;      // in ms
    std::atomic<float> baseDelay;  // in ms
    std::atomic<float> mix;        // 0.0 to 1.0

    // Delay buffers for each channel.
    std::vector<float> leftDelayBuffer;
    std::vector<float> rightDelayBuffer;
    size_t delayBufferSize; // in samples

    // Write indices for circular buffers.
    size_t leftWriteIndex;
    size_t rightWriteIndex;

    // LFO phase for each channel.
    float leftLFOPhase;
    float rightLFOPhase;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomVibrato* vib = static_cast<PhantomVibrato*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(vib->in_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(vib->in_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(vib->out_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(vib->out_right, nframes));

        // Read current parameters.
        float currentRate = vib->rate.load();
        float currentDepth = vib->depth.load();      // in ms
        float currentBaseDelay = vib->baseDelay.load(); // in ms
        float currentMix = vib->mix.load();

        // Pre-calculate LFO increment per sample.
        float dt = 1.0f / vib->sample_rate;
        float lfoInc = 2.0f * M_PI * currentRate * dt;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            // --- Process Left Channel ---
            float dryL = inL[i];
            // Write the incoming sample into the delay buffer.
            vib->leftDelayBuffer[vib->leftWriteIndex] = dryL;
            // Compute modulated delay time (in ms) for left channel.
            float modDelayL = currentBaseDelay + currentDepth * sinf(vib->leftLFOPhase);
            // Convert delay time to samples.
            float delaySamplesL = (modDelayL * vib->sample_rate) / 1000.0f;
            // Compute read index (with wrap-around).
            float readIndexL = static_cast<float>(vib->leftWriteIndex) - delaySamplesL;
            if (readIndexL < 0)
                readIndexL += vib->delayBufferSize;
            size_t indexL0 = static_cast<size_t>(readIndexL);
            size_t indexL1 = (indexL0 + 1) % vib->delayBufferSize;
            float fracL = readIndexL - static_cast<float>(indexL0);
            float wetL = (1.0f - fracL) * vib->leftDelayBuffer[indexL0] + fracL * vib->leftDelayBuffer[indexL1];
            // Mix dry and wet signals.
            float outSampleL = currentMix * wetL + (1.0f - currentMix) * dryL;
            outL[i] = outSampleL;

            // --- Process Right Channel ---
            float dryR = inR[i];
            vib->rightDelayBuffer[vib->rightWriteIndex] = dryR;
            // For right channel, add a phase offset (e.g., π/2) for a subtle stereo difference.
            float modDelayR = currentBaseDelay + currentDepth * sinf(vib->rightLFOPhase + M_PI / 2);
            float delaySamplesR = (modDelayR * vib->sample_rate) / 1000.0f;
            float readIndexR = static_cast<float>(vib->rightWriteIndex) - delaySamplesR;
            if (readIndexR < 0)
                readIndexR += vib->delayBufferSize;
            size_t indexR0 = static_cast<size_t>(readIndexR);
            size_t indexR1 = (indexR0 + 1) % vib->delayBufferSize;
            float fracR = readIndexR - static_cast<float>(indexR0);
            float wetR = (1.0f - fracR) * vib->rightDelayBuffer[indexR0] + fracR * vib->rightDelayBuffer[indexR1];
            float outSampleR = currentMix * wetR + (1.0f - currentMix) * dryR;
            outR[i] = outSampleR;

            // Increment write indices with wrap-around.
            vib->leftWriteIndex = (vib->leftWriteIndex + 1) % vib->delayBufferSize;
            vib->rightWriteIndex = (vib->rightWriteIndex + 1) % vib->delayBufferSize;

            // Advance LFO phases.
            vib->leftLFOPhase += lfoInc;
            if (vib->leftLFOPhase >= 2.0f * M_PI)
                vib->leftLFOPhase -= 2.0f * M_PI;
            vib->rightLFOPhase += lfoInc;
            if (vib->rightLFOPhase >= 2.0f * M_PI)
                vib->rightLFOPhase -= 2.0f * M_PI;
        }
        return 0;
    }

    // Control thread: allows real-time adjustment of vibrato parameters.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomVibrato] Enter parameters: rate (Hz), depth (ms), baseDelay (ms), mix (0-1)\n"
                    << "e.g., \"5 2 5 1\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newRate, newDepth, newBaseDelay, newMix;
            if (!(iss >> newRate >> newDepth >> newBaseDelay >> newMix)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomVibrato] Invalid input. Please try again." << std::endl;
                continue;
            }
            // Basic sanity checks.
            if (newRate < 0.0f) newRate = 0.0f;
            if (newDepth < 0.0f) newDepth = 0.0f;
            if (newBaseDelay < newDepth) newBaseDelay = newDepth;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            rate.store(newRate);
            depth.store(newDepth);
            baseDelay.store(newBaseDelay);
            mix.store(newMix);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomVibrato] Updated parameters: rate = " << newRate
                    << " Hz, depth = " << newDepth << " ms, baseDelay = " << newBaseDelay
                    << " ms, mix = " << newMix << std::endl;
            }
        }
    }

public:
    PhantomVibrato(const char* client_name = "PhantomVibrato")
        : client(nullptr), leftWriteIndex(0), rightWriteIndex(0),
        leftLFOPhase(0.0f), rightLFOPhase(0.0f), running(true)
    {
        // Set default vibrato parameters.
        rate.store(5.0f);       // 5 Hz LFO
        depth.store(2.0f);      // 2 ms modulation depth
        baseDelay.store(5.0f);  // 5 ms base delay
        mix.store(1.0f);        // Fully wet (vibrato effect applied)

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomVibrato: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_left = jack_port_register(client, "in_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        in_right = jack_port_register(client, "in_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_left || !in_right || !out_left || !out_right) {
            jack_client_close(client);
            throw std::runtime_error("PhantomVibrato: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomVibrato: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomVibrato: Failed to activate JACK client");
        }

        // Compute delay buffer size (in samples) from MAX_DELAY_MS.
        delayBufferSize = static_cast<size_t>((MAX_DELAY_MS * sample_rate) / 1000.0f) + 1;
        leftDelayBuffer.resize(delayBufferSize, 0.0f);
        rightDelayBuffer.resize(delayBufferSize, 0.0f);

        leftWriteIndex = 0;
        rightWriteIndex = 0;

        control_thread = std::thread(&PhantomVibrato::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomVibrato] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomVibrato] Default parameters: rate = " << rate.load()
                << " Hz, depth = " << depth.load() << " ms, baseDelay = " << baseDelay.load()
                << " ms, mix = " << mix.load() << std::endl;
        }
    }

    ~PhantomVibrato() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomVibrato] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomVibrato] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomVibrato vib;
        vib.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomVibrato] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

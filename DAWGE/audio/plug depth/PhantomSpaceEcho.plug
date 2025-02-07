// PhantomSpaceEcho.cpp
// A simple mono space echo effect using JACK.
// The effect uses a circular delay buffer and produces three echo taps at delays of 
// 1x, 2x, and 3x the base delay. Each tap is weighted by a decay factor.
// The feedback loop writes the sum back into the delay buffer, and the final output
// is a mix between the dry signal and the echo sum.
//
// Real-time adjustable parameters:
//   - Base Delay (ms): fundamental delay time.
//   - Feedback: feedback amount (typically between 0 and 0.9).
//   - Decay: multiplier for successive echoes (0.0 to 1.0).
//   - Mix: dry/wet mix (0.0 = dry, 1.0 = fully echoed).
//
// Compile with:
//   g++ -std=c++11 PhantomSpaceEcho.cpp -ljack -lpthread -o PhantomSpaceEcho

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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

class PhantomSpaceEcho {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Granular echo parameters (all atomic for real-time control):
    atomic<float> baseDelay_ms;   // Base delay in milliseconds.
    atomic<float> feedback;       // Feedback factor (e.g., 0.0 to 0.9).
    atomic<float> decay;          // Decay factor for successive taps (0.0 to 1.0).
    atomic<float> mix;            // Mix between dry and echo (0.0 = dry, 1.0 = full echo).

    // Derived parameter: base delay in samples.
    float baseDelaySamples;

    // Delay buffer and its parameters.
    vector<float> delayBuffer;
    size_t bufferSize;  // For example, 2 seconds worth of samples.
    size_t writeIndex;

    // For controlling grain spawning timing.
    // (Not grain-based here; we simply process every sample.)

    // Utility: Linear interpolation from the delay buffer.
    float readDelay(float delayInSamples) {
        // Compute the read index relative to writeIndex.
        float readIndex = static_cast<float>(writeIndex) - delayInSamples;
        // Wrap around.
        while (readIndex < 0)
            readIndex += bufferSize;
        // Use linear interpolation.
        size_t index0 = static_cast<size_t>(floor(readIndex)) % bufferSize;
        size_t index1 = (index0 + 1) % bufferSize;
        float frac = readIndex - floor(readIndex);
        return (1.0f - frac) * delayBuffer[index0] + frac * delayBuffer[index1];
    }

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomSpaceEcho* echo = static_cast<PhantomSpaceEcho*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(echo->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(echo->out_port, nframes));

        // Compute derived base delay in samples.
        float currentBaseDelay_ms = echo->baseDelay_ms.load();
        echo->baseDelaySamples = currentBaseDelay_ms * echo->sample_rate / 1000.0f;

        float currentFeedback = echo->feedback.load();
        float currentDecay = echo->decay.load();
        float currentMix = echo->mix.load();

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];

            // For a space echo, we'll compute three taps:
            // Tap 1: delay = baseDelaySamples
            // Tap 2: delay = 2 * baseDelaySamples, scaled by decay.
            // Tap 3: delay = 3 * baseDelaySamples, scaled by decay^2.
            float tap1 = echo->readDelay(echo->baseDelaySamples);
            float tap2 = echo->readDelay(2.0f * echo->baseDelaySamples);
            float tap3 = echo->readDelay(3.0f * echo->baseDelaySamples);

            float echoSum = tap1 + currentDecay * tap2 + currentDecay * currentDecay * tap3;

            // Write the new sample into the delay buffer with feedback.
            echo->delayBuffer[echo->writeIndex] = dry + currentFeedback * echoSum;
            // Increment write index.
            echo->writeIndex = (echo->writeIndex + 1) % echo->bufferSize;

            // Output is a mix of dry and echo.
            out[i] = (1.0f - currentMix) * dry + currentMix * echoSum;
        }
        return 0;
    }

    // Control thread: allow real-time parameter updates.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomSpaceEcho] Enter parameters: base delay (ms), feedback (0-0.9), decay (0-1), mix (0-1)" << endl;
                cout << "e.g., \"300 0.7 0.5 0.8\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newBaseDelay, newFeedback, newDecay, newMix;
            if (!(iss >> newBaseDelay >> newFeedback >> newDecay >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomSpaceEcho] Invalid input. Please try again." << endl;
                continue;
            }
            if (newBaseDelay < 1.0f) newBaseDelay = 1.0f;
            if (newFeedback < 0.0f) newFeedback = 0.0f;
            if (newFeedback > 0.9f) newFeedback = 0.9f;
            if (newDecay < 0.0f) newDecay = 0.0f;
            if (newDecay > 1.0f) newDecay = 1.0f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            baseDelay_ms.store(newBaseDelay);
            feedback.store(newFeedback);
            decay.store(newDecay);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomSpaceEcho] Updated parameters:" << endl;
                cout << "  Base Delay = " << newBaseDelay << " ms" << endl;
                cout << "  Feedback = " << newFeedback << endl;
                cout << "  Decay = " << newDecay << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomSpaceEcho(const char* client_name = "PhantomSpaceEcho")
        : client(nullptr), running(true), writeIndex(0)
    {
        // Set default parameters.
        baseDelay_ms.store(300.0f);  // 300 ms base delay.
        feedback.store(0.7f);        // 70% feedback.
        decay.store(0.5f);           // 50% decay per tap.
        mix.store(0.8f);             // 80% wet mix.

        // Create a delay buffer of 2 seconds length.
        sample_rate = 44100;  // Default; will be updated by JACK.
        bufferSize = 44100 * 2;
        delayBuffer.resize(bufferSize, 0.0f);
        writeIndex = 0;

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomSpaceEcho: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);
        // Resize buffer in case sample_rate changed.
        bufferSize = sample_rate * 2;
        delayBuffer.resize(bufferSize, 0.0f);
        writeIndex = 0;

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomSpaceEcho: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomSpaceEcho: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomSpaceEcho: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomSpaceEcho::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomSpaceEcho] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomSpaceEcho] Default parameters:" << endl;
        cout << "  Base Delay = " << baseDelay_ms.load() << " ms" << endl;
        cout << "  Feedback = " << feedback.load() << endl;
        cout << "  Decay = " << decay.load() << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomSpaceEcho() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomSpaceEcho] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomSpaceEcho] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomSpaceEcho spaceEcho;
        spaceEcho.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomSpaceEcho] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

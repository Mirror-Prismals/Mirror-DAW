// PhantomGranular.cpp
// A simple mono granular delay/synthesis plugin using JACK.
// It records incoming audio into a delay buffer and spawns grains that are pitch‐shifted,
// windowed, and mixed to create granular textures.
// Real-time adjustable parameters:
//   - Grain Size (ms)
//   - Grain Density (grains per second)
//   - Pitch Shift (multiplier)
//   - Mix (0.0 = dry, 1.0 = fully granular)
//   - Randomness (0.0 to 1.0 for grain start position variation)
//
// Compile with:
//   g++ -std=c++11 PhantomGranular.cpp -ljack -lpthread -o PhantomGranular

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
#include <cstdlib>
#include <ctime>

using namespace std;

// Structure representing a grain.
struct Grain {
    float startPos;  // Starting index in the delay buffer (as a float)
    float pos;       // Current position within the grain (in samples)
    int length;      // Length of the grain in samples
    float speed;     // Playback speed factor (pitch shift)
    int remaining;   // Samples remaining to play
};

class PhantomGranular {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Granular parameters (atomic for real-time updates).
    atomic<float> grainSize_ms;   // Grain duration in milliseconds.
    atomic<float> grainDensity;   // Grains per second.
    atomic<float> pitchShift;     // Playback speed factor for grains.
    atomic<float> mix;            // Mix between dry and granular output.
    atomic<float> randomness;     // 0.0 to 1.0, controls start-position variation.

    // Derived parameters.
    int grainSize_samples;        // Grain size in samples.
    int grainTriggerInterval;     // Interval in samples between grain spawns.
    int sampleCounter;            // Counts samples to trigger new grain.

    // Delay buffer: holds a few seconds of audio.
    vector<float> delayBuffer;
    size_t bufferSize;            // Number of samples (e.g., sample_rate * 2 for 2 seconds).
    size_t writeIndex;            // Current write pointer.

    // Active grains.
    vector<Grain> activeGrains;
    mutex grainMutex;

    // Utility: Hann window function.
    inline float hannWindow(float n, int N) {
        if (N <= 1) return 1.0f;
        return 0.5f * (1.0f - cosf(2.0f * M_PI * n / (N - 1)));
    }

    // Spawn a new grain.
    void spawnGrain() {
        Grain grain;
        grain.length = grainSize_samples;
        grain.remaining = grain.length;
        grain.pos = 0.0f;
        // Set playback speed from pitchShift parameter.
        grain.speed = pitchShift.load();

        // Compute a baseline start position:
        // Use the current writeIndex minus the grainSize_samples (to get recent audio).
        float basePos = static_cast<float>(writeIndex) - grainSize_samples;
        if (basePos < 0) basePos += bufferSize;
        // Apply randomness: offset in range [-randomness * grainSize_samples, +randomness * grainSize_samples]
        float randFactor = ((static_cast<float>(rand()) / RAND_MAX) * 2.0f - 1.0f); // -1 to 1
        float offset = randFactor * randomness.load() * grainSize_samples;
        grain.startPos = basePos + offset;
        // Wrap-around:
        while (grain.startPos < 0) grain.startPos += bufferSize;
        while (grain.startPos >= bufferSize) grain.startPos -= bufferSize;

        lock_guard<mutex> lock(grainMutex);
        activeGrains.push_back(grain);
    }

    // Process all active grains and sum their output.
    float processGrains() {
        float sum = 0.0f;
        lock_guard<mutex> lock(grainMutex);
        // Iterate through grains; remove those that are finished.
        for (auto it = activeGrains.begin(); it != activeGrains.end(); ) {
            // Calculate the absolute read position in the delay buffer.
            float pos = it->startPos + it->pos;
            // Wrap-around.
            while (pos >= bufferSize) pos -= bufferSize;
            while (pos < 0) pos += bufferSize;
            int index0 = static_cast<int>(floor(pos));
            int index1 = (index0 + 1) % bufferSize;
            float frac = pos - index0;
            // Linear interpolation.
            float grainSample = (1.0f - frac) * delayBuffer[index0] + frac * delayBuffer[index1];
            // Apply a Hann window.
            float windowVal = hannWindow(it->pos, it->length);
            sum += grainSample * windowVal;

            // Update grain's position.
            it->pos += it->speed;
            it->remaining--;
            if (it->remaining <= 0) {
                it = activeGrains.erase(it);
            }
            else {
                ++it;
            }
        }
        return sum;
    }

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomGranular* pg = static_cast<PhantomGranular*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(pg->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(pg->out_port, nframes));

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Write the input sample into the delay buffer.
            pg->delayBuffer[pg->writeIndex] = dry;
            pg->writeIndex = (pg->writeIndex + 1) % pg->bufferSize;

            // Increment sample counter and check for grain spawn.
            pg->sampleCounter++;
            if (pg->sampleCounter >= pg->grainTriggerInterval) {
                pg->spawnGrain();
                pg->sampleCounter = 0;
            }

            // Sum active grains.
            float granularOutput = pg->processGrains();

            // Final output is mix of dry and granular outputs.
            float finalOutput = pg->mix.load() * granularOutput + (1.0f - pg->mix.load()) * dry;
            out[i] = finalOutput;
        }
        return 0;
    }

    // Control thread: update parameters via console.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomGranular] Enter parameters: grainSize (ms), grainDensity (grains/sec), pitchShift (multiplier), mix (0-1), randomness (0-1)" << endl;
                cout << "e.g., \"100 10 1.2 0.5 0.5\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newGrainSize, newDensity, newPitchShift, newMix, newRandomness;
            if (!(iss >> newGrainSize >> newDensity >> newPitchShift >> newMix >> newRandomness)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomGranular] Invalid input. Please try again." << endl;
                continue;
            }
            // Update parameters.
            grainSize_ms.store(newGrainSize);
            pitchShift.store(newPitchShift);
            mix.store(newMix);
            randomness.store(newRandomness);
            // Recompute derived parameters.
            grainSize_samples = static_cast<int>(newGrainSize * sample_rate / 1000.0f);
            if (newDensity <= 0) newDensity = 1;
            grainTriggerInterval = static_cast<int>(sample_rate / newDensity);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomGranular] Updated parameters:" << endl;
                cout << "  Grain Size = " << newGrainSize << " ms (" << grainSize_samples << " samples)" << endl;
                cout << "  Grain Density = " << newDensity << " grains/sec (interval = " << grainTriggerInterval << " samples)" << endl;
                cout << "  Pitch Shift = " << newPitchShift << endl;
                cout << "  Mix = " << newMix << endl;
                cout << "  Randomness = " << newRandomness << endl;
            }
        }
    }

public:
    PhantomGranular(const char* client_name = "PhantomGranular")
        : client(nullptr), running(true), writeIndex(0), sampleCounter(0)
    {
        // Initialize default parameters.
        grainSize_ms.store(100.0f);   // 100 ms grains.
        grainDensity.store(10.0f);    // 10 grains per second.
        pitchShift.store(1.0f);       // No pitch shift by default.
        mix.store(0.5f);              // 50% mix.
        randomness.store(0.5f);       // 50% randomness.

        // Initially, assume sample rate is 44100 Hz; will update after JACK client is opened.
        sample_rate = 44100;
        grainSize_samples = static_cast<int>(grainSize_ms.load() * sample_rate / 1000.0f);
        grainTriggerInterval = static_cast<int>(sample_rate / grainDensity.load());
        sampleCounter = 0;

        // Set delay buffer length to 2 seconds.
        bufferSize = sample_rate * 2;
        delayBuffer.resize(bufferSize, 0.0f);
        writeIndex = 0;

        // Clear active grains.
        activeGrains.clear();

        // Seed random number generator.
        srand(static_cast<unsigned int>(time(nullptr)));

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomGranular: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);
        // Recompute derived parameters with actual sample rate.
        grainSize_samples = static_cast<int>(grainSize_ms.load() * sample_rate / 1000.0f);
        grainTriggerInterval = static_cast<int>(sample_rate / grainDensity.load());
        bufferSize = sample_rate * 2;
        delayBuffer.resize(bufferSize, 0.0f);
        writeIndex = 0;

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomGranular: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomGranular: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomGranular: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomGranular::control_loop, this);

        {
            lock_guard<mutex> lock(print_mutex);
            cout << "[PhantomGranular] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
            cout << "[PhantomGranular] Default parameters:" << endl;
            cout << "  Grain Size = " << grainSize_ms.load() << " ms (" << grainSize_samples << " samples)" << endl;
            cout << "  Grain Density = " << grainDensity.load() << " grains/sec (interval = " << grainTriggerInterval << " samples)" << endl;
            cout << "  Pitch Shift = " << pitchShift.load() << endl;
            cout << "  Mix = " << mix.load() << endl;
            cout << "  Randomness = " << randomness.load() << endl;
        }
    }

    ~PhantomGranular() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomGranular] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomGranular] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomGranular granular;
        granular.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomGranular] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

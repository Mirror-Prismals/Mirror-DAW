// PhantomChorus.cpp
// A simple mono chorus plugin using JACK and standard C++.
// It delays the input signal by a modulated amount (using a circular delay buffer and an LFO)
// and blends the delayed (chorused) signal with the dry input.
// Real-time adjustable parameters:
//   - Base Delay (ms): The average delay time (e.g., 20 ms)
//   - Modulation Depth (ms): How much the delay is modulated (e.g., 5 ms)
//   - LFO Frequency (Hz): The modulation rate (e.g., 1-5 Hz)
//   - Mix: Blend between dry and chorus (0.0 = dry, 1.0 = fully chorused)
// 
// Compile with:
//   g++ -std=c++11 PhantomChorus.cpp -ljack -lpthread -o PhantomChorus

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
#include <string>
#include <vector>

using namespace std;

class PhantomChorus {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Chorus parameters (real-time adjustable).
    atomic<float> baseDelay_ms;      // Average delay in ms (e.g., 20 ms)
    atomic<float> modulationDepth_ms; // Modulation depth in ms (e.g., 5 ms)
    atomic<float> lfoFreq;           // LFO frequency in Hz (e.g., 2 Hz)
    atomic<float> mix;               // Mix between dry and chorused signal (0.0 to 1.0)

    // Internal state.
    vector<float> delayBuffer;
    size_t bufferSize;    // In samples (we'll allocate enough for 2 seconds of audio).
    size_t writeIndex;    // Current write position in the delay buffer.
    float lfoPhase;       // LFO phase accumulator (in radians).

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomChorus* chorus = static_cast<PhantomChorus*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(chorus->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(chorus->out_port, nframes));

        // Read current parameters.
        float currentBaseDelay_ms = chorus->baseDelay_ms.load();
        float currentDepth_ms = chorus->modulationDepth_ms.load();
        float currentLfoFreq = chorus->lfoFreq.load();
        float currentMix = chorus->mix.load();

        // Convert base delay and modulation depth from ms to samples.
        float baseDelaySamples = currentBaseDelay_ms * chorus->sample_rate / 1000.0f;
        float depthSamples = currentDepth_ms * chorus->sample_rate / 1000.0f;

        // Compute LFO phase increment per sample.
        float dt = 1.0f / chorus->sample_rate;
        float phaseInc = 2.0f * M_PI * currentLfoFreq * dt;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Always write the current sample into the delay buffer.
            chorus->delayBuffer[chorus->writeIndex] = dry;

            // Compute the modulated delay in samples.
            float modDelaySamples = baseDelaySamples + depthSamples * sinf(chorus->lfoPhase);
            // Compute the read pointer: it's writeIndex - modDelaySamples.
            float readPos = static_cast<float>(chorus->writeIndex) - modDelaySamples;
            // Wrap readPos if necessary.
            while (readPos < 0)
                readPos += chorus->bufferSize;
            while (readPos >= chorus->bufferSize)
                readPos -= chorus->bufferSize;
            // Linear interpolation:
            size_t index0 = static_cast<size_t>(floor(readPos));
            size_t index1 = (index0 + 1) % chorus->bufferSize;
            float frac = readPos - floor(readPos);
            float delayed = (1.0f - frac) * chorus->delayBuffer[index0] + frac * chorus->delayBuffer[index1];

            // Output is the blend of dry and delayed (chorused) signal.
            out[i] = (1.0f - currentMix) * dry + currentMix * delayed;

            // Increment write pointer.
            chorus->writeIndex = (chorus->writeIndex + 1) % chorus->bufferSize;
            // Advance LFO phase.
            chorus->lfoPhase += phaseInc;
            if (chorus->lfoPhase >= 2.0f * M_PI)
                chorus->lfoPhase -= 2.0f * M_PI;
        }
        return 0;
    }

    // Control thread: updates parameters in real time via console.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomChorus] Enter parameters:" << endl;
                cout << "Format: <BaseDelay_ms> <ModulationDepth_ms> <LFO_Frequency_Hz> <Mix (0.0-1.0)>" << endl;
                cout << "e.g., \"20 5 2 0.7\" (20 ms base, 5 ms depth, 2 Hz LFO, 70% wet) or 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newBaseDelay, newDepth, newLfoFreq, newMix;
            if (!(iss >> newBaseDelay >> newDepth >> newLfoFreq >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomChorus] Invalid input. Please try again." << endl;
                continue;
            }
            if (newBaseDelay < 0.0f) newBaseDelay = 0.0f;
            if (newDepth < 0.0f) newDepth = 0.0f;
            if (newLfoFreq < 0.0f) newLfoFreq = 0.0f;
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            baseDelay_ms.store(newBaseDelay);
            modulationDepth_ms.store(newDepth);
            lfoFreq.store(newLfoFreq);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomChorus] Updated parameters:" << endl;
                cout << "  Base Delay = " << newBaseDelay << " ms" << endl;
                cout << "  Modulation Depth = " << newDepth << " ms" << endl;
                cout << "  LFO Frequency = " << newLfoFreq << " Hz" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomChorus(const char* client_name = "PhantomChorus")
        : client(nullptr), running(true), baseDelay_ms(20.0f),
        modulationDepth_ms(5.0f), lfoFreq(2.0f), mix(0.7f),
        lfoPhase(0.0f), writeIndex(0)
    {
        // Open JACK client.
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomChorus: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);
        // Allocate a delay buffer large enough for 2 seconds of audio.
        bufferSize = sample_rate * 2;
        delayBuffer.resize(bufferSize, 0.0f);
        writeIndex = 0;
        lfoPhase = 0.0f;

        // Register mono input and output ports.
        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomChorus: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomChorus: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomChorus: Failed to activate JACK client");
        }

        // Launch control thread.
        control_thread = thread(&PhantomChorus::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomChorus] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomChorus] Default parameters:" << endl;
        cout << "  Base Delay = " << baseDelay_ms.load() << " ms" << endl;
        cout << "  Modulation Depth = " << modulationDepth_ms.load() << " ms" << endl;
        cout << "  LFO Frequency = " << lfoFreq.load() << " Hz" << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomChorus() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomChorus] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomChorus] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomChorus chorus;
        chorus.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomChorus] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

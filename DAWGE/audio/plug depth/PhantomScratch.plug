// PhantomScratch.cpp
// A simple mono vinyl scratch simulator using JACK and standard C++.
// It uses a circular delay buffer and a variable-speed read pointer to simulate the sound of a vinyl scratch.
// When scratch speed is zero, the plugin writes new samples to the buffer and outputs the dry signal.
// When scratch speed is nonzero, the delay buffer is "frozen" and the read pointer is advanced at the specified speed.
// A mix parameter blends the scratch-processed signal with the dry signal.
//
// Real-time adjustable parameters (via console):
//   - Scratch Speed (float): in samples per sample; 0.0 means normal mode,
//       positive values play forward (e.g., 1.0 is normal speed, 2.0 is double speed),
//       negative values play in reverse.
//   - Mix (0.0 = 100% dry, 1.0 = 100% scratch signal)
//
// Compile with:
//   g++ -std=c++11 PhantomScratch.cpp -ljack -lpthread -o PhantomScratch

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
#include <string>

using namespace std;

class PhantomScratch {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Parameters:
    // scratchSpeed: if 0.0, normal mode (no scratch effect, delay buffer updates continuously).
    // If nonzero, the delay buffer is frozen and the read pointer is advanced at this rate.
    atomic<float> scratchSpeed;  // in samples per sample. Typical values: 0.0 (normal), 1.0 (normal forward), -1.0 (normal reverse), >1 or < -1 for speed variations.
    atomic<float> mix;           // Mix between dry and scratch output (0.0 = dry, 1.0 = fully processed).

    // Delay buffer to store incoming audio.
    vector<float> delayBuffer;
    size_t bufferSize;  // e.g., 1 second of audio.
    size_t writeIndex;  // Position where new samples are written.
    float readPointer;  // Floating-point read pointer for scratch playback.

    // We also maintain a flag to know if we are currently in scratch mode.
    // When switching from normal to scratch mode, we "freeze" the delay buffer.
    // When switching back to normal, we resume updating the delay buffer.
    // (We simply check if scratchSpeed is nonzero.)

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomScratch* ps = static_cast<PhantomScratch*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(ps->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(ps->out_port, nframes));

        // Get current parameters.
        float currentScratchSpeed = ps->scratchSpeed.load();
        float currentMix = ps->mix.load();

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            float processed = 0.0f;
            if (fabs(currentScratchSpeed) < 1e-6f) {
                // Normal mode: update delay buffer with new sample.
                ps->delayBuffer[ps->writeIndex] = dry;
                // Set read pointer equal to writeIndex (or writeIndex - offset if desired).
                ps->readPointer = static_cast<float>(ps->writeIndex);
                // Output is just the dry signal (or could blend, but typically mix is 0 when not scratching).
                processed = dry;
            }
            else {
                // Scratch mode: do NOT update the delay buffer.
                // Instead, read from the buffer at the current read pointer.
                // Use linear interpolation.
                float readPos = ps->readPointer;
                // Wrap readPos into [0, bufferSize).
                while (readPos < 0) readPos += ps->bufferSize;
                while (readPos >= ps->bufferSize) readPos -= ps->bufferSize;
                size_t index0 = static_cast<size_t>(floor(readPos));
                size_t index1 = (index0 + 1) % ps->bufferSize;
                float frac = readPos - floor(readPos);
                processed = (1.0f - frac) * ps->delayBuffer[index0] + frac * ps->delayBuffer[index1];
                // Advance read pointer by scratch speed.
                ps->readPointer += currentScratchSpeed;
                // Wrap readPointer.
                if (ps->readPointer < 0)
                    ps->readPointer += ps->bufferSize;
                if (ps->readPointer >= ps->bufferSize)
                    ps->readPointer -= ps->bufferSize;
            }
            // In both modes, update the write index (only in normal mode we update the buffer).
            if (fabs(currentScratchSpeed) < 1e-6f) {
                ps->writeIndex = (ps->writeIndex + 1) % ps->bufferSize;
            }
            else {
                // In scratch mode, we might choose not to update the write pointer.
                // For a true scratch simulation, we "freeze" the buffer.
                // So we do nothing here.
            }
            // Blend the dry signal and the processed (scratch) signal.
            out[i] = (1.0f - currentMix) * dry + currentMix * processed;
        }
        return 0;
    }

    // Control thread: allows updating scratch parameters via console.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomScratch] Enter parameters: scratchSpeed (samples per sample, 0.0 = normal, e.g., 1.0, -1.0, 2.0, etc.) and mix (0.0-1.0)" << endl;
                cout << "e.g., \"1.0 1.0\" for normal forward scratch, \"-1.0 1.0\" for reverse scratch, or \"0 0\" to resume normal playback, or 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newSpeed, newMix;
            if (!(iss >> newSpeed >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomScratch] Invalid input. Please try again." << endl;
                continue;
            }
            // No clamping on newSpeed (it can be negative, positive, or zero).
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            scratchSpeed.store(newSpeed);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomScratch] Updated parameters:" << endl;
                cout << "  Scratch Speed = " << newSpeed << " samples per sample" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomScratch(const char* client_name = "PhantomScratch")
        : client(nullptr), running(true), scratchSpeed(0.0f), mix(0.0f)
    {
        // Default: normal mode (scratchSpeed 0 means no scratch), mix 0 (dry signal only).
        // Later, the user can set scratchSpeed to nonzero to engage the scratch effect.
        // Set up a delay buffer of, say, 2 seconds.
        sample_rate = 44100;  // We'll update after JACK client opens.
        bufferSize = sample_rate * 2;
        delayBuffer.resize(bufferSize, 0.0f);
        writeIndex = 0;
        readPointer = 0.0f;

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomScratch: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);
        // Resize delay buffer if necessary.
        bufferSize = sample_rate * 2;
        delayBuffer.resize(bufferSize, 0.0f);
        writeIndex = 0;
        readPointer = 0.0f;

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomScratch: Failed to register JACK ports");
        }
        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomScratch: Failed to set process callback");
        }
        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomScratch: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomScratch::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomScratch] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomScratch] Default parameters: scratchSpeed = " << scratchSpeed.load() << " (normal mode), mix = " << mix.load() << " (dry)" << endl;
    }

    ~PhantomScratch() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomScratch] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomScratch] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomScratch scratch;
        scratch.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomScratch] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

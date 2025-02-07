// PhantomSpectralFreezeNoFFTW.cpp
// A toy mono spectral-freeze style plugin using JACK and standard C++ only.
// Instead of using an FFT to freeze spectral content, this implementation captures
// a block of audio and then repeatedly outputs that block until updated.
// Real-time adjustable parameters:
//   - Freeze mode: toggled on/off.
//   - Mix: blending between dry input and the frozen block (0.0 = dry, 1.0 = fully frozen).
//
// Compile with:
//   g++ -std=c++11 PhantomSpectralFreezeNoFFTW.cpp -ljack -lpthread -o PhantomSpectralFreezeNoFFTW

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

const int BLOCK_SIZE = 1024;

class PhantomSpectralFreezeNoFFTW {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    atomic<bool> freeze;      // When true, hold the current frozen block.
    atomic<float> mix;        // Mix level (0.0 = dry, 1.0 = fully frozen).
    thread control_thread;
    mutex print_mutex;

    // Buffers for block processing.
    vector<float> inputBuffer;   // Accumulate BLOCK_SIZE samples.
    int inputIndex;              // Current write index into inputBuffer.
    vector<float> frozenBlock;   // Stored frozen block.
    int frozenOutIndex;          // Read index in frozenBlock for output.

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomSpectralFreezeNoFFTW* self = static_cast<PhantomSpectralFreezeNoFFTW*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(self->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(self->out_port, nframes));

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Append incoming sample to inputBuffer.
            if (self->inputIndex < BLOCK_SIZE) {
                self->inputBuffer[self->inputIndex] = dry;
                self->inputIndex++;
            }
            // When we have a full block...
            if (self->inputIndex >= BLOCK_SIZE) {
                // If freeze mode is not engaged, update the frozen block.
                if (!self->freeze.load()) {
                    // Copy the block into frozenBlock.
                    for (int j = 0; j < BLOCK_SIZE; j++) {
                        self->frozenBlock[j] = self->inputBuffer[j];
                    }
                }
                // Reset the input buffer index for the next block.
                self->inputIndex = 0;
                // Reset frozen output index.
                self->frozenOutIndex = 0;
            }
            // If a frozen block is available, read from it; otherwise, output dry.
            float frozenSample = 0.0f;
            if (self->frozenBlock.size() == BLOCK_SIZE) {
                frozenSample = self->frozenBlock[self->frozenOutIndex];
                self->frozenOutIndex++;
                if (self->frozenOutIndex >= BLOCK_SIZE)
                    self->frozenOutIndex = 0;
            }
            // Blend dry and frozen signals.
            out[i] = (1.0f - self->mix.load()) * dry + self->mix.load() * frozenSample;
        }
        return 0;
    }

    // Control thread: listens for commands.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomSpectralFreezeNoFFTW] Commands:" << endl;
                cout << "  'freeze on'  -> engage freeze mode" << endl;
                cout << "  'freeze off' -> update frozen block with new audio" << endl;
                cout << "  'mix X'      -> set mix level (0.0 to 1.0)" << endl;
                cout << "Type 'q' to quit." << endl;
                cout << "Enter command: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            string cmd;
            iss >> cmd;
            if (cmd == "freeze") {
                string state;
                iss >> state;
                if (state == "on") {
                    freeze.store(true);
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomSpectralFreezeNoFFTW] Freeze engaged." << endl;
                }
                else if (state == "off") {
                    freeze.store(false);
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomSpectralFreezeNoFFTW] Freeze disengaged." << endl;
                }
                else {
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomSpectralFreezeNoFFTW] Unknown freeze command. Use 'freeze on' or 'freeze off'." << endl;
                }
            }
            else if (cmd == "mix") {
                float newMix;
                if (!(iss >> newMix)) {
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomSpectralFreezeNoFFTW] Invalid mix value." << endl;
                    continue;
                }
                if (newMix < 0.0f) newMix = 0.0f;
                if (newMix > 1.0f) newMix = 1.0f;
                mix.store(newMix);
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomSpectralFreezeNoFFTW] Updated mix to " << newMix << endl;
            }
            else {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomSpectralFreezeNoFFTW] Unknown command." << endl;
            }
        }
    }

public:
    PhantomSpectralFreezeNoFFTW(const char* client_name = "PhantomSpectralFreezeNoFFTW")
        : client(nullptr), running(true), freeze(false), mix(1.0f), inputIndex(0), frozenOutIndex(0)
    {
        inputBuffer.resize(BLOCK_SIZE, 0.0f);
        frozenBlock.resize(BLOCK_SIZE, 0.0f);

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomSpectralFreezeNoFFTW: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomSpectralFreezeNoFFTW: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomSpectralFreezeNoFFTW: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomSpectralFreezeNoFFTW: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomSpectralFreezeNoFFTW::control_loop, this);

        {
            lock_guard<mutex> lock(print_mutex);
            cout << "[PhantomSpectralFreezeNoFFTW] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
            cout << "[PhantomSpectralFreezeNoFFTW] Default parameters: freeze off, mix = " << mix.load() << endl;
        }
    }

    ~PhantomSpectralFreezeNoFFTW() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomSpectralFreezeNoFFTW] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomSpectralFreezeNoFFTW] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomSpectralFreezeNoFFTW plugin;
        plugin.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomSpectralFreezeNoFFTW] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

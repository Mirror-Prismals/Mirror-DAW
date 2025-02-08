// PhantomPluckSynth.cpp
//
// A simple mono physical‑modeling synthesizer (plucked‑string) using the Karplus‑Strong algorithm.
// The synth runs as a JACK client and outputs synthesized audio to a mono port.
// 
// Real-time controllable parameters via the console:
//   - "freq <value>"   : Set the frequency in Hz (e.g., freq 440)
//   - "amp <value>"    : Set the amplitude (0.0 to 1.0) (e.g., amp 0.8)
//   - "damp <value>"   : Set the damping factor (e.g., damp 0.995)
//   - "trigger"        : Pluck the string (reinitialize the delay buffer with noise)
//   - "q"              : Quit
//
// Compile with:
//   g++ -std=c++11 PhantomPluckSynth.cpp -ljack -lpthread -o PhantomPluckSynth

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
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;

class PhantomPluckSynth {
private:
    jack_client_t* client;
    jack_port_t* out_port; // Mono output
    int sample_rate;

    atomic<bool> running;
    atomic<bool> trigger;  // When true, a new pluck will be triggered in the next process block
    atomic<float> frequency;  // Frequency in Hz
    atomic<float> amplitude;  // Amplitude (0.0 - 1.0)
    atomic<float> damping;    // Damping factor (typical values: 0.90 to 0.999)

    // Internal Karplus‑Strong state:
    vector<float> delayBuffer; // Delay line (the "string")
    int delayLength;           // Length of the delay line in samples (depends on frequency)
    int bufferIndex;           // Current read/write index in the delay buffer

    mutex print_mutex;
    thread control_thread;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomPluckSynth* synth = static_cast<PhantomPluckSynth*>(arg);

        // If a new pluck is triggered, reinitialize the delay buffer.
        if (synth->trigger.load()) {
            // Compute new delay length from frequency.
            float freq = synth->frequency.load();
            if (freq <= 0.0f)
                freq = 1.0f;
            synth->delayLength = static_cast<int>(synth->sample_rate / freq);
            if (synth->delayLength < 2)
                synth->delayLength = 2;
            // Resize the delay buffer.
            synth->delayBuffer.resize(synth->delayLength);
            // Fill the delay buffer with random noise in range [-1, 1] scaled by amplitude.
            float amp = synth->amplitude.load();
            for (int i = 0; i < synth->delayLength; i++) {
                float rnd = 2.0f * (static_cast<float>(rand()) / RAND_MAX) - 1.0f;
                synth->delayBuffer[i] = amp * rnd;
            }
            // Reset buffer index.
            synth->bufferIndex = 0;
            // Clear the trigger flag.
            synth->trigger.store(false);
        }

        float* out = static_cast<float*>(jack_port_get_buffer(synth->out_port, nframes));
        // If the delay buffer is not properly initialized, output silence.
        if (synth->delayBuffer.size() < 2) {
            for (jack_nframes_t i = 0; i < nframes; i++) {
                out[i] = 0.0f;
            }
            return 0;
        }

        float damp = synth->damping.load();
        // Process each sample.
        for (jack_nframes_t i = 0; i < nframes; i++) {
            // Output the current sample from the delay buffer.
            float y = synth->delayBuffer[synth->bufferIndex];
            out[i] = y;
            // Calculate next sample using Karplus‑Strong averaging.
            int nextIndex = (synth->bufferIndex + 1) % synth->delayLength;
            float nextSample = synth->delayBuffer[nextIndex];
            float newSample = damp * 0.5f * (y + nextSample);
            // Store the new sample in the delay buffer.
            synth->delayBuffer[synth->bufferIndex] = newSample;
            // Increment the buffer index.
            synth->bufferIndex = (synth->bufferIndex + 1) % synth->delayLength;
        }
        return 0;
    }

    // Control thread: listens for user commands.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomPluckSynth] Enter command:" << endl;
                cout << "Commands:" << endl;
                cout << "  freq <value>   - set frequency in Hz (e.g., freq 440)" << endl;
                cout << "  amp <value>    - set amplitude (0.0 to 1.0) (e.g., amp 0.8)" << endl;
                cout << "  damp <value>   - set damping factor (e.g., damp 0.995)" << endl;
                cout << "  trigger        - pluck the string" << endl;
                cout << "  q              - quit" << endl;
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
            if (cmd == "freq") {
                float newFreq;
                if (iss >> newFreq) {
                    frequency.store(newFreq);
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomPluckSynth] Frequency set to " << newFreq << " Hz" << endl;
                }
                else {
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomPluckSynth] Invalid frequency value." << endl;
                }
            }
            else if (cmd == "amp") {
                float newAmp;
                if (iss >> newAmp) {
                    if (newAmp < 0.0f) newAmp = 0.0f;
                    if (newAmp > 1.0f) newAmp = 1.0f;
                    amplitude.store(newAmp);
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomPluckSynth] Amplitude set to " << newAmp << endl;
                }
                else {
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomPluckSynth] Invalid amplitude value." << endl;
                }
            }
            else if (cmd == "damp") {
                float newDamp;
                if (iss >> newDamp) {
                    if (newDamp < 0.90f) newDamp = 0.90f;
                    if (newDamp > 0.999f) newDamp = 0.999f;
                    damping.store(newDamp);
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomPluckSynth] Damping set to " << newDamp << endl;
                }
                else {
                    lock_guard<mutex> lock(print_mutex);
                    cout << "[PhantomPluckSynth] Invalid damping value." << endl;
                }
            }
            else if (cmd == "trigger") {
                trigger.store(true);
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomPluckSynth] Triggered pluck." << endl;
            }
            else {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomPluckSynth] Unknown command." << endl;
            }
        }
    }

public:
    PhantomPluckSynth(const char* client_name = "PhantomPluckSynth")
        : client(nullptr), running(true), trigger(false),
        frequency(440.0f), amplitude(0.8f), damping(0.995f),
        delayLength(0), bufferIndex(0)
    {
        srand(static_cast<unsigned int>(time(nullptr)));

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomPluckSynth: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Register a mono output port.
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomPluckSynth: Failed to register JACK output port");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomPluckSynth: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomPluckSynth: Failed to activate JACK client");
        }

        // Initialize delay buffer (empty until triggered).
        delayBuffer.clear();
        delayLength = 0;
        bufferIndex = 0;

        control_thread = thread(&PhantomPluckSynth::control_loop, this);

        {
            lock_guard<mutex> lock(print_mutex);
            cout << "[PhantomPluckSynth] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
            cout << "[PhantomPluckSynth] Default parameters:" << endl;
            cout << "  Frequency = " << frequency.load() << " Hz" << endl;
            cout << "  Amplitude = " << amplitude.load() << endl;
            cout << "  Damping = " << damping.load() << endl;
            cout << "  (Type 'trigger' to pluck the string)" << endl;
        }
    }

    ~PhantomPluckSynth() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomPluckSynth] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomPluckSynth] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomPluckSynth synth;
        synth.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomPluckSynth] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

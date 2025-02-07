// PhantomSynth.cpp
// A simple mono synthesizer implemented as a JACK client in C++.
// It generates an oscillator whose parameters (frequency, amplitude, and waveform) can be controlled in real time.
// Supported waveforms: sine, square, saw, and triangle.
// Compile with:
//   g++ -std=c++11 PhantomSynth.cpp -ljack -lpthread -o PhantomSynth

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

using namespace std;

// Enumerate supported waveform types.
enum Waveform { SINE, SQUARE, SAW, TRIANGLE };

class PhantomSynth {
private:
    jack_client_t* client;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Synth parameters.
    atomic<float> frequency;  // in Hz.
    atomic<float> amplitude;  // 0.0 to 1.0.
    atomic<Waveform> waveform; // current waveform type.

    float phase;  // Current phase of the oscillator (in radians).

    // JACK process callback: generates synthesized audio.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomSynth* synth = static_cast<PhantomSynth*>(arg);
        float* out = static_cast<float*>(jack_port_get_buffer(synth->out_port, nframes));

        float freq = synth->frequency.load();
        float amp = synth->amplitude.load();
        Waveform wf = synth->waveform.load();

        // Compute phase increment per sample.
        float phaseInc = 2.0f * M_PI * freq / synth->sample_rate;
        float currentPhase = synth->phase;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float sample = 0.0f;
            // Generate waveform sample based on current waveform type.
            switch (wf) {
            case SINE:
                sample = sinf(currentPhase);
                break;
            case SQUARE:
                sample = (sinf(currentPhase) >= 0.0f) ? 1.0f : -1.0f;
                break;
            case SAW:
                // Sawtooth: map phase from 0 to 2π into -1 to 1.
                sample = 2.0f * (currentPhase / (2.0f * M_PI)) - 1.0f;
                break;
            case TRIANGLE:
                // Triangle: using asin(sin(x)) scaled to [-1, 1].
                sample = (2.0f / M_PI) * asinf(sinf(currentPhase));
                break;
            default:
                sample = sinf(currentPhase);
                break;
            }
            // Output the sample scaled by amplitude.
            out[i] = amp * sample;
            // Increment phase.
            currentPhase += phaseInc;
            if (currentPhase >= 2.0f * M_PI)
                currentPhase -= 2.0f * M_PI;
        }
        synth->phase = currentPhase;
        return 0;
    }

    // Control thread: allows real-time updating of synth parameters.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomSynth] Enter parameters:" << endl;
                cout << "Format: <frequency (Hz)> <amplitude (0.0-1.0)> <waveform (sine, square, saw, triangle)>" << endl;
                cout << "For example: \"440 0.8 sine\" or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newFreq, newAmp;
            string wfStr;
            if (!(iss >> newFreq >> newAmp >> wfStr)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomSynth] Invalid input. Please try again." << endl;
                continue;
            }
            if (newAmp < 0.0f) newAmp = 0.0f;
            if (newAmp > 1.0f) newAmp = 1.0f;
            frequency.store(newFreq);
            amplitude.store(newAmp);
            Waveform newWaveform = SINE;
            if (wfStr == "sine")
                newWaveform = SINE;
            else if (wfStr == "square")
                newWaveform = SQUARE;
            else if (wfStr == "saw")
                newWaveform = SAW;
            else if (wfStr == "triangle")
                newWaveform = TRIANGLE;
            else {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomSynth] Unknown waveform. Defaulting to sine." << endl;
                newWaveform = SINE;
            }
            waveform.store(newWaveform);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomSynth] Updated parameters:" << endl;
                cout << "  Frequency = " << newFreq << " Hz" << endl;
                cout << "  Amplitude = " << newAmp << endl;
                cout << "  Waveform = " << wfStr << endl;
            }
        }
    }

public:
    PhantomSynth(const char* client_name = "PhantomSynth")
        : client(nullptr), running(true), frequency(440.0f), amplitude(0.8f),
        waveform(SINE), phase(0.0f)
    {
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomSynth: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Register a mono output port.
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomSynth: Failed to register JACK output port");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomSynth: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomSynth: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomSynth::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomSynth] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomSynth] Default parameters: Frequency = " << frequency.load()
            << " Hz, Amplitude = " << amplitude.load() << ", Waveform = sine" << endl;
    }

    ~PhantomSynth() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomSynth] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomSynth] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomSynth synth;
        synth.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomSynth] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

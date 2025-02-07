// PhantomFreqShifter.cpp
// A simple mono frequency shifter using JACK and standard C++.
// It creates an approximate analytic signal via a 31-tap FIR Hilbert transformer,
// multiplies it by a complex exponential to shift the frequency,
// and outputs the real part of the result.
// Real-time adjustable parameters:
//   - Frequency Shift (Hz): the amount by which to shift the spectrum (can be positive or negative).
//   - Mix: blend between dry and frequency-shifted signals (0.0 = dry, 1.0 = fully shifted).
//
// Compile with:
//   g++ -std=c++11 PhantomFreqShifter.cpp -ljack -lpthread -o PhantomFreqShifter

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

using namespace std;

class PhantomFreqShifter {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Parameters:
    atomic<float> shiftHz;  // Frequency shift in Hz (can be negative)
    atomic<float> mix;      // Mix between dry and shifted (0.0 to 1.0)

    // Hilbert transform FIR parameters.
    static const int FIR_TAPS = 31;
    int firCenter;  // Center tap index (FIR_TAPS/2)
    vector<float> hilbertCoeffs;  // FIR coefficients for the Hilbert transformer.
    vector<float> firBuffer;      // Circular buffer to hold FIR_TAPS samples.
    int firIndex;                 // Current index in FIR buffer.

    // Phase accumulator for frequency shifting.
    float phase;

    // Precompute FIR Hilbert transformer coefficients using a Hamming window.
    void initHilbertCoeffs() {
        hilbertCoeffs.resize(FIR_TAPS, 0.0f);
        firBuffer.resize(FIR_TAPS, 0.0f);
        firCenter = FIR_TAPS / 2;
        // For an odd-length FIR Hilbert transformer, coefficients:
        // h[k] = (2 / (pi*(k - M))) * w[k] for k != M, h[M] = 0.
        for (int k = 0; k < FIR_TAPS; k++) {
            if (k == firCenter) {
                hilbertCoeffs[k] = 0.0f;
            }
            else {
                float n = float(k - firCenter);
                // Hamming window: w[k] = 0.54 - 0.46*cos(2*pi*k/(N-1))
                float window = 0.54f - 0.46f * cosf(2.0f * M_PI * k / (FIR_TAPS - 1));
                hilbertCoeffs[k] = (2.0f / (M_PI * n)) * window;
            }
        }
        // Initialize FIR buffer to zeros.
        firIndex = 0;
    }

    // Process one sample through the Hilbert transformer.
    float processHilbert(float sample) {
        // Insert the new sample into the FIR buffer at current index.
        firBuffer[firIndex] = sample;
        // Compute convolution using circular buffer.
        float imag = 0.0f;
        for (int k = 0; k < FIR_TAPS; k++) {
            int index = (firIndex - k + FIR_TAPS) % FIR_TAPS;
            imag += firBuffer[index] * hilbertCoeffs[k];
        }
        // Increment firIndex.
        firIndex = (firIndex + 1) % FIR_TAPS;
        return imag;
    }

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomFreqShifter* shifter = static_cast<PhantomFreqShifter*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(shifter->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(shifter->out_port, nframes));

        float currentShiftHz = shifter->shiftHz.load();
        float currentMix = shifter->mix.load();
        float dt = 1.0f / shifter->sample_rate;
        float phaseInc = 2.0f * M_PI * currentShiftHz * dt;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Get the imaginary part via Hilbert transform.
            float imag = shifter->processHilbert(dry);
            // Form the analytic signal: (dry + j*imag)
            // Multiply by the complex exponential: exp(j*phase) = cos(phase) + j*sin(phase)
            // The shifted signal is: dry*cos(phase) - imag*sin(phase)
            float shifted = dry * cosf(shifter->phase) - imag * sinf(shifter->phase);
            // Advance phase.
            shifter->phase += phaseInc;
            if (shifter->phase >= 2.0f * M_PI)
                shifter->phase -= 2.0f * M_PI;
            // Mix dry and shifted signals.
            out[i] = (1.0f - currentMix) * dry + currentMix * shifted;
        }
        return 0;
    }

    // Control thread: allows real-time parameter updates.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomFreqShifter] Enter parameters: frequency shift (Hz) and mix (0.0-1.0)" << endl;
                cout << "e.g., \"100 0.7\" (100 Hz shift, 70% shifted signal) or type 'q' to quit: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            float newShiftHz, newMix;
            if (!(iss >> newShiftHz >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomFreqShifter] Invalid input. Please try again." << endl;
                continue;
            }
            // No clamping needed for frequency shift; mix is clamped between 0 and 1.
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            shiftHz.store(newShiftHz);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomFreqShifter] Updated parameters:" << endl;
                cout << "  Frequency Shift = " << newShiftHz << " Hz" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomFreqShifter(const char* client_name = "PhantomFreqShifter")
        : client(nullptr), running(true), phase(0.0f)
    {
        // Set default parameters.
        shiftHz.store(100.0f);  // Default frequency shift: 100 Hz.
        mix.store(0.7f);        // 70% mix.

        // Initialize Hilbert transformer.
        initHilbertCoeffs();

        // Open JACK client.
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomFreqShifter: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Register mono input and output ports.
        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomFreqShifter: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomFreqShifter: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomFreqShifter: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomFreqShifter::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomFreqShifter] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomFreqShifter] Default parameters:" << endl;
        cout << "  Frequency Shift = " << shiftHz.load() << " Hz" << endl;
        cout << "  Mix = " << mix.load() << endl;
    }

    ~PhantomFreqShifter() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomFreqShifter] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomFreqShifter] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomFreqShifter shifter;
        shifter.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomFreqShifter] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

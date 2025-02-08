// PhantomDither.cpp
// A simple mono dither plugin using JACK and standard C++.
// This plugin simulates bit-depth reduction by quantizing the input to a target bit depth,
// but adds TPDF dither noise before quantization to decorrelate quantization errors.
// The output is a mix of the dry signal and the quantized (dithered) signal.
//
// Real-time adjustable parameters:
//   - Bit Depth: target bit depth (e.g., 16)
//   - Mix: blend between dry and dithered signal (0.0 = dry, 1.0 = fully quantized/dithered)
//
// Compile with:
//   g++ -std=c++11 PhantomDither.cpp -ljack -lpthread -o PhantomDither

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
#include <string>

using namespace std;

class PhantomDither {
private:
    jack_client_t* client;
    jack_port_t* in_port;
    jack_port_t* out_port;
    int sample_rate;

    atomic<bool> running;
    thread control_thread;
    mutex print_mutex;

    // Parameters:
    // Target bit depth (e.g., 16 bits).
    atomic<int> bitDepth;
    // Mix between dry and dithered output (0.0 = dry, 1.0 = fully dithered).
    atomic<float> mix;

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomDither* dither = static_cast<PhantomDither*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(dither->in_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(dither->out_port, nframes));

        int currentBitDepth = dither->bitDepth.load();
        float currentMix = dither->mix.load();
        // Compute quantization step.
        // For a signed signal in the range [-1, +1], we assume
        // step = 1 / (2^(bitDepth-1)). For example, for 16-bit, step ≈ 1/32768.
        float step = 1.0f / static_cast<float>(1 << (currentBitDepth - 1));

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Generate TPDF dither noise.
            // Generate two random floats in [0,1] and subtract them.
            float r1 = static_cast<float>(rand()) / RAND_MAX;
            float r2 = static_cast<float>(rand()) / RAND_MAX;
            // Scale by half the step to produce noise in [-step/2, step/2].
            float ditherNoise = (r1 - r2) * (step / 2.0f);
            // Add noise to the dry sample.
            float noisySample = dry + ditherNoise;
            // Quantize the result: round to nearest multiple of step.
            float quantized = round(noisySample / step) * step;
            // Blend the dry and quantized (dithered) signals.
            out[i] = (1.0f - currentMix) * dry + currentMix * quantized;
        }
        return 0;
    }

    // Control thread: allows updating parameters via console.
    void control_loop() {
        string line;
        while (running.load()) {
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "\n[PhantomDither] Enter parameters: bitDepth (e.g., 16) and mix (0.0-1.0)" << endl;
                cout << "For example: \"16 1.0\" for 16-bit dither, full effect; or \"24 0.0\" for 24-bit (effectively no dither), dry signal." << endl;
                cout << "Enter command: ";
            }
            if (!getline(cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            istringstream iss(line);
            int newBitDepth;
            float newMix;
            if (!(iss >> newBitDepth >> newMix)) {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomDither] Invalid input. Please try again." << endl;
                continue;
            }
            // Clamp mix to [0.0, 1.0].
            if (newMix < 0.0f) newMix = 0.0f;
            if (newMix > 1.0f) newMix = 1.0f;
            // Optionally, clamp bitDepth to a reasonable range (e.g., 8 to 24).
            if (newBitDepth < 8) newBitDepth = 8;
            if (newBitDepth > 24) newBitDepth = 24;
            bitDepth.store(newBitDepth);
            mix.store(newMix);
            {
                lock_guard<mutex> lock(print_mutex);
                cout << "[PhantomDither] Updated parameters:" << endl;
                cout << "  Bit Depth = " << newBitDepth << " bits" << endl;
                cout << "  Mix = " << newMix << endl;
            }
        }
    }

public:
    PhantomDither(const char* client_name = "PhantomDither")
        : client(nullptr), running(true), bitDepth(16), mix(1.0f)
    {
        // Seed random number generator.
        srand(static_cast<unsigned int>(time(nullptr)));

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw runtime_error("PhantomDither: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_port = jack_port_register(client, "in", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_port = jack_port_register(client, "out", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_port || !out_port) {
            jack_client_close(client);
            throw runtime_error("PhantomDither: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomDither: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw runtime_error("PhantomDither: Failed to activate JACK client");
        }

        control_thread = thread(&PhantomDither::control_loop, this);

        lock_guard<mutex> lock(print_mutex);
        cout << "[PhantomDither] Initialized. Sample rate: " << sample_rate << " Hz" << endl;
        cout << "[PhantomDither] Default parameters: Bit Depth = " << bitDepth.load()
            << " bits, Mix = " << mix.load() << " (fully dithered)" << endl;
    }

    ~PhantomDither() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        cout << "[PhantomDither] Running. Type 'q' in the control console to quit." << endl;
        while (running.load()) {
            this_thread::sleep_for(chrono::milliseconds(100));
        }
        cout << "[PhantomDither] Shutting down." << endl;
    }
};

int main() {
    try {
        PhantomDither dither;
        dither.run();
    }
    catch (const exception& e) {
        cerr << "[PhantomDither] Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}

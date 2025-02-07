// PhantomTape.cpp
// A simple real-time tape saturation (tape emulation) plugin using JACK
//
// Compile with:
//   g++ -std=c++11 PhantomTape.cpp -ljack -lpthread -o PhantomTape

#include <jack/jack.h>
#include <iostream>
#include <atomic>
#include <thread>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class PhantomTape {
private:
    jack_client_t* client;
    jack_port_t* input_port;
    jack_port_t* output_port;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Tape saturation parameters:
    // drive: Multiplier applied to the input before saturation (≥ 0).
    // mix: Wet/dry mix (0.0 = completely dry, 1.0 = fully processed).
    // cutoff: Cutoff frequency (Hz) for the low-pass filter to emulate high-frequency roll-off.
    // output_gain_dB: Overall output gain in dB (e.g., -20 dB up to +10 dB).
    std::atomic<float> drive;
    std::atomic<float> mix;
    std::atomic<float> cutoff;
    std::atomic<float> output_gain_dB;

    // Low-pass filter state variable (for a one-pole filter)
    float prev_sample;

    // JACK process callback: processes a block of audio samples.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomTape* tape = static_cast<PhantomTape*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(tape->input_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(tape->output_port, nframes));

        // Load parameters (atomically)
        float current_drive = tape->drive.load();
        float current_mix = tape->mix.load();
        float current_cutoff = tape->cutoff.load();
        float current_output_gain_dB = tape->output_gain_dB.load();
        // Convert output gain from dB to linear (linear_gain = 10^(dB/20))
        float linear_gain = powf(10.0f, current_output_gain_dB / 20.0f);

        // Compute the low-pass filter coefficient:
        // Using the one-pole filter: y[n] = alpha*x[n] + (1 - alpha)*y[n-1]
        float dt = 1.0f / tape->sample_rate;
        float RC = 1.0f / (2.0f * M_PI * current_cutoff);
        float alpha = dt / (RC + dt);

        // Process each sample in the current JACK frame:
        for (jack_nframes_t i = 0; i < nframes; i++) {
            float dry = in[i];
            // Apply drive and saturate the signal with tanh
            float saturated = tanhf(current_drive * dry);
            // Process the saturated signal through a low-pass filter
            float filtered = alpha * saturated + (1.0f - alpha) * tape->prev_sample;
            tape->prev_sample = filtered;
            // Blend the filtered (wet) signal with the dry signal
            float processed = current_mix * filtered + (1.0f - current_mix) * dry;
            // Apply the output gain (converted from dB)
            out[i] = processed * linear_gain;
        }
        return 0;
    }

    // Control thread: allows real-time parameter adjustments.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomTape] Enter new parameters: drive, mix (0-1), cutoff (Hz), output gain (dB)\n"
                    << "e.g., \"2.0 0.5 8000 0.0\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float new_drive, new_mix, new_cutoff, new_output_gain_dB;
            if (!(iss >> new_drive >> new_mix >> new_cutoff >> new_output_gain_dB)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomTape] Invalid input. Please try again." << std::endl;
                continue;
            }
            if (new_drive < 0.0f)
                new_drive = 0.0f;
            if (new_mix < 0.0f)
                new_mix = 0.0f;
            if (new_mix > 1.0f)
                new_mix = 1.0f;
            // Limit cutoff to a sensible range (e.g., 20 Hz to 15000 Hz)
            if (new_cutoff < 20.0f)
                new_cutoff = 20.0f;
            if (new_cutoff > 15000.0f)
                new_cutoff = 15000.0f;
            drive.store(new_drive);
            mix.store(new_mix);
            cutoff.store(new_cutoff);
            output_gain_dB.store(new_output_gain_dB);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomTape] Updated parameters:\n"
                    << "  drive = " << new_drive << "\n"
                    << "  mix = " << new_mix << "\n"
                    << "  cutoff = " << new_cutoff << " Hz\n"
                    << "  output gain = " << new_output_gain_dB << " dB" << std::endl;
            }
        }
    }

public:
    PhantomTape(const char* client_name = "PhantomTape")
        : client(nullptr), input_port(nullptr), output_port(nullptr),
        running(true), drive(2.0f), mix(0.5f), cutoff(8000.0f), output_gain_dB(0.0f),
        prev_sample(0.0f)
    {
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomTape: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        input_port = jack_port_register(client, "input", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        output_port = jack_port_register(client, "output", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!input_port || !output_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomTape: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomTape: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomTape: Failed to activate JACK client");
        }

        control_thread = std::thread(&PhantomTape::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomTape] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomTape] Default parameters:\n"
                << "  drive = " << drive.load() << "\n"
                << "  mix = " << mix.load() << "\n"
                << "  cutoff = " << cutoff.load() << " Hz\n"
                << "  output gain = " << output_gain_dB.load() << " dB" << std::endl;
        }
    }

    ~PhantomTape() {
        running.store(false);
        if (control_thread.joinable()) {
            control_thread.join();
        }
        if (client) {
            jack_client_close(client);
        }
    }

    void run() {
        std::cout << "[PhantomTape] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomTape] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomTape tape;
        tape.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomTape] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

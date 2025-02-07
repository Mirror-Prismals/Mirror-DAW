// PhantomComp.cpp
// A simple real-time compressor using JACK
//
// Compile with:
// g++ -std=c++11 PhantomComp.cpp -ljack -lpthread -o PhantomComp
//

#include <jack/jack.h>
#include <iostream>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>
#include <chrono>
#include <sstream>
#include <stdexcept>
#include <cmath>

// Utility: Convert decibels to a linear scale.
float dBToLinear(float dB) {
    return pow(10.0f, dB / 20.0f);
}

class PhantomComp {
private:
    jack_client_t* client;
    jack_port_t* input_port;
    jack_port_t* output_port;
    int sample_rate;

    // Real-time running flag and control thread.
    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Compressor parameters (with default values)
    // threshold is in dB (e.g., -20 dB means signals above 0.1 in linear domain)
    std::atomic<float> threshold;    // Default: -20 dB
    std::atomic<float> ratio;        // Default: 4:1
    std::atomic<float> attack;       // Attack time in milliseconds (default: 10 ms)
    std::atomic<float> release;      // Release time in milliseconds (default: 100 ms)
    std::atomic<float> makeup_gain;  // Linear gain applied after compression (default: 1.0)

    // The envelope detector state.
    float envelope;

    // JACK process callback: applies compression sample-by-sample.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomComp* comp = static_cast<PhantomComp*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(comp->input_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(comp->output_port, nframes));

        // Compute smoothing coefficients from attack/release times.
        // Using the formula: coeff = exp(-1/(time_constant * sample_rate))
        // Multiply time_constant in seconds by sample_rate.
        float attack_coeff = expf(-1000.0f / (comp->sample_rate * comp->attack.load()));
        float release_coeff = expf(-1000.0f / (comp->sample_rate * comp->release.load()));
        // Convert threshold from dB to linear.
        float thresh_linear = dBToLinear(comp->threshold.load());

        for (jack_nframes_t i = 0; i < nframes; i++) {
            float input = in[i];
            float abs_input = fabs(input);

            // Update the envelope detector:
            if (abs_input > comp->envelope)
                comp->envelope = attack_coeff * comp->envelope + (1 - attack_coeff) * abs_input;
            else
                comp->envelope = release_coeff * comp->envelope + (1 - release_coeff) * abs_input;

            // Compute gain reduction.
            float gain = 1.0f;
            if (comp->envelope > thresh_linear) {
                // For a signal above the threshold, the compressor reduces gain.
                // A simplified approach: compute how much over threshold (in linear ratio)
                float over = comp->envelope / thresh_linear;
                // The desired gain reduction is such that the output level is compressed by the ratio.
                // One common formulation is: gain = (over)^(1/ratio - 1)
                gain = pow(over, (1.0f / comp->ratio.load()) - 1.0f);
            }
            // Apply makeup gain.
            gain *= comp->makeup_gain.load();

            // Process the sample.
            out[i] = input * gain;
        }
        return 0;
    }

    // Control thread to allow real-time parameter adjustments.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomComp] Enter new parameters: threshold (dB), ratio, attack (ms), release (ms), makeup gain (linear) (or type 'q' to quit): ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float new_threshold, new_ratio, new_attack, new_release, new_makeup;
            if (!(iss >> new_threshold >> new_ratio >> new_attack >> new_release >> new_makeup)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomComp] Invalid input. Please try again." << std::endl;
                continue;
            }
            threshold.store(new_threshold);
            ratio.store(new_ratio);
            attack.store(new_attack);
            release.store(new_release);
            makeup_gain.store(new_makeup);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomComp] Updated parameters: threshold = " << new_threshold
                    << " dB, ratio = " << new_ratio << ":1, attack = " << new_attack
                    << " ms, release = " << new_release << " ms, makeup gain = " << new_makeup << std::endl;
            }
        }
    }

public:
    PhantomComp(const char* client_name = "PhantomComp")
        : client(nullptr), input_port(nullptr), output_port(nullptr),
        running(true), envelope(0.0f) {
        // Set default compressor parameters.
        threshold.store(-20.0f);
        ratio.store(4.0f);
        attack.store(10.0f);
        release.store(100.0f);
        makeup_gain.store(1.0f);

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomComp: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        input_port = jack_port_register(client, "input", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        output_port = jack_port_register(client, "output", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!input_port || !output_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomComp: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomComp: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomComp: Failed to activate JACK client");
        }

        // Start the control thread.
        control_thread = std::thread(&PhantomComp::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomComp] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomComp] Default parameters: threshold = " << threshold.load() << " dB, ratio = "
                << ratio.load() << ":1, attack = " << attack.load() << " ms, release = "
                << release.load() << " ms, makeup gain = " << makeup_gain.load() << std::endl;
        }
    }

    ~PhantomComp() {
        running.store(false);
        if (control_thread.joinable()) {
            control_thread.join();
        }
        if (client) {
            jack_client_close(client);
        }
    }

    void run() {
        std::cout << "[PhantomComp] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomComp] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomComp comp;
        comp.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomComp] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

// PhantomDist.cpp
// A simple real-time distortion effect using JACK with output gain in dB
//
// Compile with:
//   g++ -std=c++11 PhantomDist.cpp -ljack -lpthread -o PhantomDist

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

class PhantomDist {
private:
    jack_client_t* client;
    jack_port_t* input_port;
    jack_port_t* output_port;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Distortion parameters:
    // drive: multiplier for input signal before nonlinear processing (default: 2.0)
    // mix: wet/dry mix where 0.0 = dry and 1.0 = fully distorted (default: 0.5)
    // output_gain_dB: output gain specified in decibels (default: 0.0 dB for unity gain)
    std::atomic<float> drive;
    std::atomic<float> mix;
    std::atomic<float> output_gain_dB;

    // JACK process callback: applies distortion to each sample.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomDist* dist = static_cast<PhantomDist*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(dist->input_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(dist->output_port, nframes));

        float current_drive = dist->drive.load();
        float current_mix = dist->mix.load();
        float current_output_gain_dB = dist->output_gain_dB.load();
        // Convert output gain in dB to a linear multiplier.
        float current_output_gain = powf(10.0f, current_output_gain_dB / 20.0f);

        for (jack_nframes_t i = 0; i < nframes; ++i) {
            float input = in[i];
            // Apply drive and then soft clipping using tanh.
            float distorted = std::tanh(current_drive * input);
            // Mix the distorted signal with the original signal.
            float processed = current_mix * distorted + (1.0f - current_mix) * input;
            // Apply output gain (now specified in dB).
            out[i] = processed * current_output_gain;
        }

        return 0;
    }

    // Control thread: allows real-time adjustment of drive, mix, and output gain in dB.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomDist] Enter new drive, mix, and output gain (in dB, e.g., \"2.0 0.5 0.0\") "
                    "(drive must be >= 0; mix between 0.0 and 1.0; output gain from -inf up to +10 dB), "
                    "or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float new_drive, new_mix, new_output_gain_dB;
            if (!(iss >> new_drive >> new_mix >> new_output_gain_dB)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomDist] Invalid input. Please try again." << std::endl;
                continue;
            }
            if (new_drive < 0.0f)
                new_drive = 0.0f;
            if (new_mix < 0.0f)
                new_mix = 0.0f;
            if (new_mix > 1.0f)
                new_mix = 1.0f;
            // Clamp output gain dB to a maximum of +10 dB.
            if (new_output_gain_dB > 10.0f)
                new_output_gain_dB = 10.0f;
            // (Allow negative values to represent attenuation; extremely low values represent -infinity.)
            drive.store(new_drive);
            mix.store(new_mix);
            output_gain_dB.store(new_output_gain_dB);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomDist] Updated parameters: drive = " << new_drive
                    << ", mix = " << new_mix
                    << ", output gain = " << new_output_gain_dB << " dB" << std::endl;
            }
        }
    }

public:
    PhantomDist(const char* client_name = "PhantomDist")
        : client(nullptr), input_port(nullptr), output_port(nullptr),
        running(true), drive(2.0f), mix(0.5f), output_gain_dB(0.0f) {
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomDist: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        input_port = jack_port_register(client, "input", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        output_port = jack_port_register(client, "output", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!input_port || !output_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomDist: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomDist: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomDist: Failed to activate JACK client");
        }

        control_thread = std::thread(&PhantomDist::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomDist] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomDist] Default parameters: drive = " << drive.load()
                << ", mix = " << mix.load()
                << ", output gain = " << output_gain_dB.load() << " dB" << std::endl;
        }
    }

    ~PhantomDist() {
        running.store(false);
        if (control_thread.joinable()) {
            control_thread.join();
        }
        if (client) {
            jack_client_close(client);
        }
    }

    void run() {
        std::cout << "[PhantomDist] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomDist] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomDist dist;
        dist.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomDist] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

// PhantomWide.cpp
// A simple real-time stereo widening effect using JACK
//
// Compile with:
//   g++ -std=c++11 PhantomWide.cpp -ljack -lpthread -o PhantomWide

#include <jack/jack.h>
#include <iostream>
#include <atomic>
#include <thread>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <cmath>

class PhantomWide {
private:
    jack_client_t* client;
    jack_port_t* input_port_left;
    jack_port_t* input_port_right;
    jack_port_t* output_port_left;
    jack_port_t* output_port_right;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Stereo widening parameter:
    // side_gain: multiplier applied to the side channel.
    // Default 1.0 means no change.
    std::atomic<float> side_gain;

    // JACK process callback: processes stereo audio.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomWide* wide = static_cast<PhantomWide*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(wide->input_port_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(wide->input_port_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(wide->output_port_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(wide->output_port_right, nframes));

        float current_side_gain = wide->side_gain.load();

        // Process each sample frame
        for (jack_nframes_t i = 0; i < nframes; i++) {
            float L = inL[i];
            float R = inR[i];
            // Convert to mid-side representation
            float mid = 0.5f * (L + R);
            float side = 0.5f * (L - R);
            // Apply widening factor to the side channel
            side *= current_side_gain;
            // Reconstruct stereo signal
            outL[i] = mid + side;
            outR[i] = mid - side;
        }
        return 0;
    }

    // Control thread: allows real-time adjustment of the side gain.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomWide] Enter new side gain (e.g., 1.0 for no change, 1.5 to widen, 0.8 to narrow), or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float new_side_gain;
            if (!(iss >> new_side_gain)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomWide] Invalid input. Please try again." << std::endl;
                continue;
            }
            side_gain.store(new_side_gain);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomWide] Updated side gain: " << new_side_gain << std::endl;
            }
        }
    }

public:
    PhantomWide(const char* client_name = "PhantomWide")
        : client(nullptr),
        input_port_left(nullptr), input_port_right(nullptr),
        output_port_left(nullptr), output_port_right(nullptr),
        running(true), side_gain(1.0f) {
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomWide: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Register stereo input ports
        input_port_left = jack_port_register(client, "input_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        input_port_right = jack_port_register(client, "input_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        // Register stereo output ports
        output_port_left = jack_port_register(client, "output_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        output_port_right = jack_port_register(client, "output_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!input_port_left || !input_port_right || !output_port_left || !output_port_right) {
            jack_client_close(client);
            throw std::runtime_error("PhantomWide: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomWide: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomWide: Failed to activate JACK client");
        }

        // Start the control thread for real-time parameter adjustment.
        control_thread = std::thread(&PhantomWide::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomWide] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomWide] Default side gain: " << side_gain.load() << std::endl;
        }
    }

    ~PhantomWide() {
        running.store(false);
        if (control_thread.joinable()) {
            control_thread.join();
        }
        if (client) {
            jack_client_close(client);
        }
    }

    void run() {
        std::cout << "[PhantomWide] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomWide] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomWide wide;
        wide.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomWide] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

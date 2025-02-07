// PhantomEcho.cpp
// A real-time delay/echo effect processor using JACK

#include <iostream>
#include <jack/jack.h>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <sstream>

// PhantomEcho applies a delay (echo) effect with real-time control over delay time and feedback.
// It uses a circular buffer to store incoming samples and mixes delayed samples back into the output.

class PhantomEcho {
private:
    jack_client_t* client;
    jack_port_t* input_port;
    jack_port_t* output_port;

    // Circular buffer for delay (2 seconds max)
    std::vector<float> delay_buffer;
    std::atomic<size_t> write_index;
    size_t buffer_size; // in samples

    // Real-time adjustable parameters
    std::atomic<int> delay_time_ms;  // delay time in milliseconds
    std::atomic<float> feedback;     // feedback factor (0.0 - 1.0)

    int sample_rate; // obtained from JACK

    // Control thread for user parameter input
    std::atomic<bool> running;
    std::thread control_thread;

    // Mutex for printing (to avoid jumbled console output)
    std::mutex print_mutex;

    // JACK process callback: applies the delay effect sample-by-sample.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomEcho* echo = static_cast<PhantomEcho*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(echo->input_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(echo->output_port, nframes));

        // For each sample in the current JACK frame:
        for (jack_nframes_t i = 0; i < nframes; i++) {
            // Compute delay in samples from current delay_time_ms
            size_t delay_samples = static_cast<size_t>((echo->delay_time_ms.load() * echo->sample_rate) / 1000);
            // Calculate read index for the delayed sample
            size_t read_index = (echo->write_index.load() + echo->buffer_size - delay_samples) % echo->buffer_size;
            float delayed_sample = echo->delay_buffer[read_index];

            // Mix input and delayed signal
            float input_sample = in[i];
            float output_sample = input_sample + delayed_sample;

            // Write processed sample to output
            out[i] = output_sample;

            // Store new sample into the delay buffer with feedback applied
            echo->delay_buffer[echo->write_index.load()] = input_sample + delayed_sample * echo->feedback.load();

            // Increment write index circularly
            echo->write_index = (echo->write_index.load() + 1) % echo->buffer_size;
        }
        return 0;
    }

    // Control loop: runs on a separate thread to allow real-time parameter adjustments.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomEcho] Enter new delay time (ms) and feedback (0.0-1.0), separated by space (or type 'q' to quit): ";
            }
            if (!std::getline(std::cin, line))
                break;

            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }

            std::istringstream iss(line);
            int new_delay;
            float new_feedback;
            if (!(iss >> new_delay >> new_feedback)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomEcho] Invalid input. Please try again." << std::endl;
                continue;
            }

            // Clamp feedback between 0.0 and 1.0
            if (new_feedback < 0.0f) new_feedback = 0.0f;
            if (new_feedback > 1.0f) new_feedback = 1.0f;

            delay_time_ms.store(new_delay);
            feedback.store(new_feedback);

            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomEcho] Updated parameters: delay_time = " << new_delay
                          << " ms, feedback = " << new_feedback << std::endl;
            }
        }
    }

public:
    PhantomEcho(const char* client_name = "PhantomEcho")
        : client(nullptr), input_port(nullptr), output_port(nullptr),
          write_index(0), delay_time_ms(500), feedback(0.5f), running(true) {

        // Open JACK client
        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomEcho: Failed to open JACK client");
        }

        // Get the sample rate
        sample_rate = jack_get_sample_rate(client);

        // Allocate a 2-second delay buffer
        buffer_size = sample_rate * 2;
        delay_buffer.resize(buffer_size, 0.0f);

        // Register input and output ports
        input_port = jack_port_register(client, "input", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        output_port = jack_port_register(client, "output", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!input_port || !output_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomEcho: Failed to register JACK ports");
        }

        // Set the process callback
        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomEcho: Failed to set process callback");
        }

        // Activate the JACK client
        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomEcho: Failed to activate JACK client");
        }

        // Launch the control thread for real-time parameter tuning
        control_thread = std::thread(&PhantomEcho::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomEcho] Initialized. Sample rate: " << sample_rate << " Hz, "
                      << "Buffer size: " << buffer_size << " samples." << std::endl;
            std::cout << "[PhantomEcho] Default parameters: delay_time = " << delay_time_ms.load()
                      << " ms, feedback = " << feedback.load() << std::endl;
        }
    }

    ~PhantomEcho() {
        running.store(false);
        if (control_thread.joinable()) {
            control_thread.join();
        }
        if (client) {
            jack_client_close(client);
        }
    }

    // Runs until the user indicates to quit.
    void run() {
        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomEcho] Running. Press Enter (with 'q' input in control) to exit." << std::endl;
        }
        // Wait until the control thread signals termination.
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::lock_guard<std::mutex> lock(print_mutex);
        std::cout << "[PhantomEcho] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomEcho echo;
        echo.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomEcho] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

// PhantomPhaser.cpp
// A simple real-time stereo phaser effect using JACK
//
// Compile with:
//   g++ -std=c++11 PhantomPhaser.cpp -ljack -lpthread -o PhantomPhaser

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
#include <vector>

#define NUM_STAGES 4

// Structure for a single all-pass filter stage used in the phaser.
struct AllPassStage {
    float x_prev;
    float y_prev;
    float lfo_offset; // Unique offset for this stage's modulation.
    AllPassStage(float offset) : x_prev(0.0f), y_prev(0.0f), lfo_offset(offset) {}
};

class PhantomPhaser {
private:
    jack_client_t* client;
    jack_port_t* in_left, * in_right, * out_left, * out_right;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Phaser parameters
    std::atomic<float> rate;     // LFO frequency in Hz (e.g., 0.5 Hz)
    std::atomic<float> depth;    // Modulation depth (0–1)
    std::atomic<float> feedback; // Feedback amount (-1 to 1)
    std::atomic<float> mix;      // Dry/wet mix (0–1)

    // LFO phase accumulator (in radians)
    float lfo_phase;

    // Separate chains of all-pass filter stages for left and right channels.
    std::vector<AllPassStage> left_stages;
    std::vector<AllPassStage> right_stages;

    // Feedback state for each channel.
    float left_fb, right_fb;

    // JACK process callback: processes audio block by block.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomPhaser* phaser = static_cast<PhantomPhaser*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(phaser->in_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(phaser->in_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(phaser->out_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(phaser->out_right, nframes));

        float current_rate = phaser->rate.load();
        float current_depth = phaser->depth.load();
        float current_feedback = phaser->feedback.load();
        float current_mix = phaser->mix.load();

        float dt = 1.0f / phaser->sample_rate;
        float lfo_increment = 2.0f * M_PI * current_rate * dt;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            // --- Process Left Channel ---
            float inputL = inL[i];
            // Incorporate feedback into the input.
            float xL = inputL + phaser->left_fb * current_feedback;
            // Process through each all-pass filter stage.
            for (int stage = 0; stage < NUM_STAGES; stage++) {
                // Compute modulated coefficient for this stage.
                // Here, we choose a base coefficient of 0.5 and modulate it by 0.3 * depth.
                float a = 0.5f + 0.3f * current_depth * sinf(phaser->lfo_phase + phaser->left_stages[stage].lfo_offset);
                float y = -a * xL + phaser->left_stages[stage].x_prev + a * phaser->left_stages[stage].y_prev;
                phaser->left_stages[stage].x_prev = xL;
                phaser->left_stages[stage].y_prev = y;
                xL = y;
            }
            // Store the processed output for feedback.
            phaser->left_fb = xL;
            // Mix processed signal with dry signal.
            float outputL = current_mix * xL + (1.0f - current_mix) * inputL;
            outL[i] = outputL;

            // --- Process Right Channel ---
            float inputR = inR[i];
            float xR = inputR + phaser->right_fb * current_feedback;
            for (int stage = 0; stage < NUM_STAGES; stage++) {
                float a = 0.5f + 0.3f * current_depth * sinf(phaser->lfo_phase + phaser->right_stages[stage].lfo_offset);
                float y = -a * xR + phaser->right_stages[stage].x_prev + a * phaser->right_stages[stage].y_prev;
                phaser->right_stages[stage].x_prev = xR;
                phaser->right_stages[stage].y_prev = y;
                xR = y;
            }
            phaser->right_fb = xR;
            float outputR = current_mix * xR + (1.0f - current_mix) * inputR;
            outR[i] = outputR;

            // Advance the global LFO phase.
            phaser->lfo_phase += lfo_increment;
            if (phaser->lfo_phase >= 2.0f * M_PI)
                phaser->lfo_phase -= 2.0f * M_PI;
        }
        return 0;
    }

    // Control thread: allows real-time parameter updates.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomPhaser] Enter parameters: rate (Hz), depth (0-1), feedback (-1 to 1), mix (0-1)\n"
                    << "e.g., \"0.5 0.5 0.3 0.5\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float new_rate, new_depth, new_feedback, new_mix;
            if (!(iss >> new_rate >> new_depth >> new_feedback >> new_mix)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomPhaser] Invalid input. Please try again." << std::endl;
                continue;
            }
            rate.store(new_rate);
            depth.store(new_depth);
            feedback.store(new_feedback);
            mix.store(new_mix);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomPhaser] Updated parameters: rate = " << new_rate
                    << " Hz, depth = " << new_depth
                    << ", feedback = " << new_feedback
                    << ", mix = " << new_mix << std::endl;
            }
        }
    }

public:
    PhantomPhaser(const char* client_name = "PhantomPhaser")
        : client(nullptr), running(true), lfo_phase(0.0f), left_fb(0.0f), right_fb(0.0f)
    {
        // Set default phaser parameters.
        rate.store(0.5f);     // 0.5 Hz LFO
        depth.store(0.5f);    // Moderate modulation depth
        feedback.store(0.3f); // Some feedback
        mix.store(0.5f);      // 50/50 dry/wet mix

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomPhaser: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_left = jack_port_register(client, "in_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        in_right = jack_port_register(client, "in_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_left || !in_right || !out_left || !out_right) {
            jack_client_close(client);
            throw std::runtime_error("PhantomPhaser: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomPhaser: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomPhaser: Failed to activate JACK client");
        }

        // Initialize filter chains for both channels.
        // Distribute phase offsets evenly for the stages.
        for (int i = 0; i < NUM_STAGES; i++) {
            float offset = (2.0f * M_PI / NUM_STAGES) * i;
            left_stages.push_back(AllPassStage(offset));
            right_stages.push_back(AllPassStage(offset));
        }

        control_thread = std::thread(&PhantomPhaser::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomPhaser] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomPhaser] Default parameters: rate = " << rate.load()
                << " Hz, depth = " << depth.load()
                << ", feedback = " << feedback.load()
                << ", mix = " << mix.load() << std::endl;
        }
    }

    ~PhantomPhaser() {
        running.store(false);
        if (control_thread.joinable()) {
            control_thread.join();
        }
        if (client) {
            jack_client_close(client);
        }
    }

    void run() {
        std::cout << "[PhantomPhaser] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomPhaser] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomPhaser phaser;
        phaser.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomPhaser] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

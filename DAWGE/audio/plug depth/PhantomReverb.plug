// PhantomReverb.cpp
// A real-time reverb effect processor using JACK

#include <iostream>
#include <jack/jack.h>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>
#include <sstream>
#include <stdexcept>

// PhantomReverb uses four parallel comb filters and one all-pass filter
// to simulate a reverb tail. User-adjustable parameters include the comb
// filter feedback (which influences the decay time) and the wet/dry mix.

class PhantomReverb {
private:
    jack_client_t* client;
    jack_port_t* input_port;
    jack_port_t* output_port;
    int sample_rate;

    // Reverb parameters (adjustable in real time)
    std::atomic<float> comb_feedback;  // Should be between 0.0 and 1.0 (default: 0.8)
    std::atomic<float> mix;            // Wet/dry mix: 0.0 (dry) to 1.0 (wet) (default: 0.5)

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Comb filter structure
    struct CombFilter {
        std::vector<float> buffer;
        size_t index;
        size_t delay;    // delay length in samples
    };
    std::vector<CombFilter> comb_filters;

    // All-pass filter structure
    struct AllpassFilter {
        std::vector<float> buffer;
        size_t index;
        size_t delay;
        float feedback;  // Typically around 0.7
    } allpass_filter;

    // JACK process callback: applies the reverb effect on each audio frame.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomReverb* reverb = static_cast<PhantomReverb*>(arg);
        float* in = static_cast<float*>(jack_port_get_buffer(reverb->input_port, nframes));
        float* out = static_cast<float*>(jack_port_get_buffer(reverb->output_port, nframes));

        // For each sample in the current block:
        for (jack_nframes_t i = 0; i < nframes; i++) {
            float input_sample = in[i];

            // Process parallel comb filters
            float comb_sum = 0.0f;
            for (auto& cf : reverb->comb_filters) {
                // Retrieve delayed sample from the comb filter buffer
                float delayed = cf.buffer[cf.index];
                // Update the comb filter: current input plus feedback * delayed sample
                cf.buffer[cf.index] = input_sample + delayed * reverb->comb_feedback.load();
                comb_sum += delayed;
                // Increment the circular buffer index
                cf.index = (cf.index + 1) % cf.delay;
            }
            // Average the output of the comb filters
            comb_sum /= static_cast<float>(reverb->comb_filters.size());

            // Process through the all-pass filter
            float ap_delayed = reverb->allpass_filter.buffer[reverb->allpass_filter.index];
            float allpass_out = -reverb->allpass_filter.feedback * comb_sum + ap_delayed;
            reverb->allpass_filter.buffer[reverb->allpass_filter.index] = comb_sum + reverb->allpass_filter.feedback * allpass_out;
            reverb->allpass_filter.index = (reverb->allpass_filter.index + 1) % reverb->allpass_filter.delay;

            // Mix the wet (reverb) and dry (original) signals according to mix parameter
            float wet = allpass_out;
            float dry = input_sample;
            float output_sample = (1.0f - reverb->mix.load()) * dry + reverb->mix.load() * wet;
            out[i] = output_sample;
        }
        return 0;
    }

    // Control thread loop: allows updating the comb feedback and mix parameters.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomReverb] Enter new comb_feedback (0.0-1.0) and mix (0.0-1.0), separated by space (or type 'q' to quit): ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float new_feedback, new_mix;
            if (!(iss >> new_feedback >> new_mix)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomReverb] Invalid input. Please try again." << std::endl;
                continue;
            }
            // Clamp the values between 0.0 and 1.0
            if (new_feedback < 0.0f) new_feedback = 0.0f;
            if (new_feedback > 1.0f) new_feedback = 1.0f;
            if (new_mix < 0.0f) new_mix = 0.0f;
            if (new_mix > 1.0f) new_mix = 1.0f;
            comb_feedback.store(new_feedback);
            mix.store(new_mix);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomReverb] Updated parameters: comb_feedback = " << new_feedback
                    << ", mix = " << new_mix << std::endl;
            }
        }
    }

public:
    PhantomReverb(const char* client_name = "PhantomReverb")
        : client(nullptr), input_port(nullptr), output_port(nullptr),
        comb_feedback(0.8f), mix(0.5f), running(true) {

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomReverb: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Register one input and one output port
        input_port = jack_port_register(client, "input", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        output_port = jack_port_register(client, "output", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!input_port || !output_port) {
            jack_client_close(client);
            throw std::runtime_error("PhantomReverb: Failed to register JACK ports");
        }

        // Initialize four comb filters with chosen delay lengths (in samples)
        // (These numbers are taken from classic reverb designs and work well at 44100 Hz)
        int comb_delays[4] = { 1116, 1188, 1277, 1356 };
        for (int d : comb_delays) {
            CombFilter cf;
            cf.delay = static_cast<size_t>(d);
            cf.buffer.resize(cf.delay, 0.0f);
            cf.index = 0;
            comb_filters.push_back(cf);
        }

        // Initialize one all-pass filter (using a delay of 225 samples and feedback ~0.7)
        allpass_filter.delay = 225;
        allpass_filter.buffer.resize(allpass_filter.delay, 0.0f);
        allpass_filter.index = 0;
        allpass_filter.feedback = 0.7f;

        // Set JACK process callback
        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomReverb: Failed to set process callback");
        }

        // Activate the JACK client
        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomReverb: Failed to activate JACK client");
        }

        // Launch the control thread
        control_thread = std::thread(&PhantomReverb::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomReverb] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomReverb] Default parameters: comb_feedback = " << comb_feedback.load()
                << ", mix = " << mix.load() << std::endl;
        }
    }

    ~PhantomReverb() {
        running.store(false);
        if (control_thread.joinable()) {
            control_thread.join();
        }
        if (client) {
            jack_client_close(client);
        }
    }

    // Run until the user quits (via the control thread)
    void run() {
        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomReverb] Running. Press Enter (with 'q' in control) to exit." << std::endl;
        }
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomReverb] Shutting down." << std::endl;
        }
    }
};

int main() {
    try {
        PhantomReverb reverb;
        reverb.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomReverb] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

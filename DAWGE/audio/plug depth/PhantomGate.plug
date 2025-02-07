// PhantomGate.cpp
// A simple real-time noise gate plugin using JACK.
// When the computed envelope of the input signal falls below a specified threshold,
// the gate closes (outputting zero); otherwise, the original signal passes through.
// Real-time adjustable parameters: threshold (in dB), attack (ms), and release (ms).

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

// PhantomGate class encapsulates the JACK client and processing.
class PhantomGate {
private:
    jack_client_t* client;
    jack_port_t* in_left, * in_right, * out_left, * out_right;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Gate parameters (atomic variables)
    // Threshold is specified in dB (e.g., -40 dB) and converted to a linear value.
    std::atomic<float> threshold_dB;
    std::atomic<float> attackTime;   // in milliseconds
    std::atomic<float> releaseTime;  // in milliseconds

    // Envelope detectors for left and right channels.
    float leftEnvelope;
    float rightEnvelope;

    // Helper: Convert dB value to a linear amplitude.
    inline float dBToLinear(float dB) {
        return powf(10.0f, dB / 20.0f);
    }

    // JACK process callback: updates the envelope for each sample and gates the signal.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomGate* gate = static_cast<PhantomGate*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(gate->in_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(gate->in_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(gate->out_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(gate->out_right, nframes));

        float dt = 1.0f / gate->sample_rate; // seconds per sample
        float dt_ms = dt * 1000.0f;

        // Compute smoothing coefficients based on attack/release times.
        float att_coeff = expf(-dt_ms / gate->attackTime.load());
        float rel_coeff = expf(-dt_ms / gate->releaseTime.load());

        // Convert threshold from dB to linear amplitude.
        float linThreshold = gate->dBToLinear(gate->threshold_dB.load());

        for (jack_nframes_t i = 0; i < nframes; i++) {
            // Process left channel:
            float sampleL = inL[i];
            float absL = fabs(sampleL);
            if (absL > gate->leftEnvelope)
                gate->leftEnvelope = att_coeff * gate->leftEnvelope + (1.0f - att_coeff) * absL;
            else
                gate->leftEnvelope = rel_coeff * gate->leftEnvelope + (1.0f - rel_coeff) * absL;

            // Process right channel:
            float sampleR = inR[i];
            float absR = fabs(sampleR);
            if (absR > gate->rightEnvelope)
                gate->rightEnvelope = att_coeff * gate->rightEnvelope + (1.0f - att_coeff) * absR;
            else
                gate->rightEnvelope = rel_coeff * gate->rightEnvelope + (1.0f - rel_coeff) * absR;

            // If the envelope is below threshold, output silence; otherwise, pass the input.
            outL[i] = (gate->leftEnvelope < linThreshold) ? 0.0f : sampleL;
            outR[i] = (gate->rightEnvelope < linThreshold) ? 0.0f : sampleR;
        }
        return 0;
    }

    // Control thread: allows real-time parameter updates via the console.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomGate] Enter parameters: threshold (dB), attack (ms), release (ms)\n"
                    << "e.g., \"-40 10 50\" or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newThreshold, newAttack, newRelease;
            if (!(iss >> newThreshold >> newAttack >> newRelease)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomGate] Invalid input. Please try again." << std::endl;
                continue;
            }
            threshold_dB.store(newThreshold);
            attackTime.store(newAttack);
            releaseTime.store(newRelease);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomGate] Updated parameters: threshold = " << newThreshold
                    << " dB, attack = " << newAttack << " ms, release = " << newRelease << " ms" << std::endl;
            }
        }
    }

public:
    PhantomGate(const char* client_name = "PhantomGate")
        : client(nullptr), running(true), leftEnvelope(0.0f), rightEnvelope(0.0f)
    {
        // Set default parameters.
        threshold_dB.store(-40.0f);
        attackTime.store(10.0f);
        releaseTime.store(50.0f);

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomGate: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_left = jack_port_register(client, "in_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        in_right = jack_port_register(client, "in_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_left || !in_right || !out_left || !out_right) {
            jack_client_close(client);
            throw std::runtime_error("PhantomGate: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomGate: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomGate: Failed to activate JACK client");
        }

        control_thread = std::thread(&PhantomGate::control_loop, this);

        std::lock_guard<std::mutex> lock(print_mutex);
        std::cout << "[PhantomGate] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
        std::cout << "[PhantomGate] Default parameters: threshold = " << threshold_dB.load()
            << " dB, attack = " << attackTime.load() << " ms, release = " << releaseTime.load() << " ms" << std::endl;
    }

    ~PhantomGate() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomGate] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomGate] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomGate gate;
        gate.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomGate] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

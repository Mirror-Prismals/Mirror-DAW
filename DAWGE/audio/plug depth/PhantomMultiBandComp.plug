// PhantomMultibandComp.cpp
// A simple 4‑band multiband compressor using JACK.
// Splits the stereo input into 4 bands and applies compression per band.
// Controls for each band (threshold in dB, ratio, attack ms, release ms, makeup) are adjustable in real time.
//
// Compile with:
//   g++ -std=c++11 PhantomMultibandComp.cpp -ljack -lpthread -o PhantomMultibandComp

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

// ---------------------------------------------------
// Simple first-order filter state structures
struct LPF {
    float y_prev;
    LPF() : y_prev(0.0f) {}
};

struct HPF {
    float x_prev;
    float y_prev;
    HPF() : x_prev(0.0f), y_prev(0.0f) {}
};

// ---------------------------------------------------
// First-order low-pass filter function.
// Uses: y[n] = α * x[n] + (1-α) * y[n-1], with α = dt/(RC+dt)
inline float processLPF(float x, LPF& state, float cutoff, float dt) {
    float RC = 1.0f / (2.0f * M_PI * cutoff);
    float alpha = dt / (RC + dt);
    float y = alpha * x + (1.0f - alpha) * state.y_prev;
    state.y_prev = y;
    return y;
}

// First-order high-pass filter function.
// Uses: y[n] = β * (y[n-1] + x[n] - x[n-1]), with β = RC/(RC+dt)
inline float processHPF(float x, HPF& state, float cutoff, float dt) {
    float RC = 1.0f / (2.0f * M_PI * cutoff);
    float beta = RC / (RC + dt);
    float y = beta * (state.y_prev + x - state.x_prev);
    state.x_prev = x;
    state.y_prev = y;
    return y;
}

// ---------------------------------------------------
// PhantomMultibandComp class
class PhantomMultibandComp {
private:
    jack_client_t* client;
    jack_port_t* in_left, * in_right, * out_left, * out_right;
    int sample_rate;
    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // Fixed crossover frequencies:
    const float F1 = 200.0f;   // between Band 1 and 2
    const float F2 = 1000.0f;  // between Band 2 and 3
    const float F3 = 5000.0f;  // between Band 3 and 4

    // Filter states for splitting (per channel, index 0=left, 1=right)
    // Band 1 (Low): LPF at F1.
    LPF lpf_low[2];

    // Band 2 (Low-Mid): HPF at F1, then LPF at F2.
    HPF hpf_band2[2];
    LPF lpf_band2[2];

    // Band 3 (High-Mid): HPF at F2, then LPF at F3.
    HPF hpf_band3[2];
    LPF lpf_band3[2];

    // Band 4 (High): HPF at F3.
    HPF hpf_band4[2];

    // Compressor parameters for each band (4 bands).
    // Each parameter is stored as an atomic float.
    std::atomic<float> compThreshold[4]; // in dB (e.g., -20 dB)
    std::atomic<float> compRatio[4];       // e.g., 2:1, 3:1, etc.
    std::atomic<float> compAttack[4];      // in ms
    std::atomic<float> compRelease[4];     // in ms
    std::atomic<float> compMakeup[4];      // linear gain

    // Envelope state for compressor (per band, per channel: [band][channel])
    float envelope[4][2];

    // Helper function to convert dB to linear.
    inline float dBToLinear(float dB) {
        return powf(10.0f, dB / 20.0f);
    }

    // Compressor function: processes a sample x for a given band and channel.
    float compressSample(float x, int band, int channel, float dt) {
        float threshold_dB = compThreshold[band].load();
        float threshold_lin = dBToLinear(threshold_dB); // typically < 1 if threshold_dB is negative.
        float ratio = compRatio[band].load();
        float attack_ms = compAttack[band].load();
        float release_ms = compRelease[band].load();
        float makeup = compMakeup[band].load();

        // Compute attack and release coefficients.
        float attack_coeff = expf(-dt * 1000.0f / attack_ms);
        float release_coeff = expf(-dt * 1000.0f / release_ms);

        float x_abs = fabs(x);
        float& env = envelope[band][channel];
        if (x_abs > env)
            env = attack_coeff * env + (1.0f - attack_coeff) * x_abs;
        else
            env = release_coeff * env + (1.0f - release_coeff) * x_abs;

        float gain = 1.0f;
        if (env > threshold_lin && threshold_lin > 0)
            gain = powf(env / threshold_lin, 1.0f / ratio - 1.0f);

        return x * gain * makeup;
    }

    // JACK process callback.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        PhantomMultibandComp* comp = static_cast<PhantomMultibandComp*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(comp->in_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(comp->in_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(comp->out_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(comp->out_right, nframes));

        float dt = 1.0f / comp->sample_rate;

        for (jack_nframes_t i = 0; i < nframes; i++) {
            // Process left channel:
            float xL = inL[i];
            // Split into bands:
            // Band 1 (Low):
            float band1 = processLPF(xL, comp->lpf_low[0], comp->F1, dt);
            // Band 2 (Low-Mid): high-pass at F1, then low-pass at F2.
            float temp2 = processHPF(xL, comp->hpf_band2[0], comp->F1, dt);
            float band2 = processLPF(temp2, comp->lpf_band2[0], comp->F2, dt);
            // Band 3 (High-Mid): high-pass at F2, then low-pass at F3.
            float temp3 = processHPF(xL, comp->hpf_band3[0], comp->F2, dt);
            float band3 = processLPF(temp3, comp->lpf_band3[0], comp->F3, dt);
            // Band 4 (High):
            float band4 = processHPF(xL, comp->hpf_band4[0], comp->F3, dt);

            // Apply compression per band.
            float comp1 = comp->compressSample(band1, 0, 0, dt);
            float comp2 = comp->compressSample(band2, 1, 0, dt);
            float comp3 = comp->compressSample(band3, 2, 0, dt);
            float comp4 = comp->compressSample(band4, 3, 0, dt);

            // Sum bands.
            outL[i] = comp1 + comp2 + comp3 + comp4;

            // Process right channel similarly:
            float xR = inR[i];
            float band1_r = processLPF(xR, comp->lpf_low[1], comp->F1, dt);
            float temp2_r = processHPF(xR, comp->hpf_band2[1], comp->F1, dt);
            float band2_r = processLPF(temp2_r, comp->lpf_band2[1], comp->F2, dt);
            float temp3_r = processHPF(xR, comp->hpf_band3[1], comp->F2, dt);
            float band3_r = processLPF(temp3_r, comp->lpf_band3[1], comp->F3, dt);
            float band4_r = processHPF(xR, comp->hpf_band4[1], comp->F3, dt);

            float comp1_r = comp->compressSample(band1_r, 0, 1, dt);
            float comp2_r = comp->compressSample(band2_r, 1, 1, dt);
            float comp3_r = comp->compressSample(band3_r, 2, 1, dt);
            float comp4_r = comp->compressSample(band4_r, 3, 1, dt);

            outR[i] = comp1_r + comp2_r + comp3_r + comp4_r;
        }
        return 0;
    }

    // Control thread: allows updating compressor parameters per band.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[PhantomMultibandComp] Enter band number (1-4) and new parameters:\n"
                    << "Threshold (dB), Ratio, Attack (ms), Release (ms), Makeup (linear)\n"
                    << "For example: \"2 -18 3.0 10 100 1.0\" to update band 2, or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            int band;
            float newThreshold, newRatio, newAttack, newRelease, newMakeup;
            if (!(iss >> band >> newThreshold >> newRatio >> newAttack >> newRelease >> newMakeup)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomMultibandComp] Invalid input. Please try again." << std::endl;
                continue;
            }
            if (band < 1 || band > 4) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomMultibandComp] Band number must be between 1 and 4." << std::endl;
                continue;
            }
            int idx = band - 1;
            compThreshold[idx].store(newThreshold);
            compRatio[idx].store(newRatio);
            compAttack[idx].store(newAttack);
            compRelease[idx].store(newRelease);
            compMakeup[idx].store(newMakeup);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[PhantomMultibandComp] Updated band " << band << " parameters: "
                    << "Threshold = " << newThreshold << " dB, "
                    << "Ratio = " << newRatio << ", "
                    << "Attack = " << newAttack << " ms, "
                    << "Release = " << newRelease << " ms, "
                    << "Makeup = " << newMakeup << std::endl;
            }
        }
    }

public:
    PhantomMultibandComp(const char* client_name = "PhantomMultibandComp")
        : client(nullptr), running(true)
    {
        // Initialize compressor envelope states to zero.
        for (int b = 0; b < 4; b++) {
            for (int ch = 0; ch < 2; ch++) {
                envelope[b][ch] = 0.0f;
            }
        }
        // Set default compressor parameters.
        // Band 1 (Low)
        compThreshold[0].store(-20.0f);
        compRatio[0].store(2.0f);
        compAttack[0].store(10.0f);
        compRelease[0].store(100.0f);
        compMakeup[0].store(1.0f);
        // Band 2 (Low-Mid)
        compThreshold[1].store(-20.0f);
        compRatio[1].store(3.0f);
        compAttack[1].store(10.0f);
        compRelease[1].store(100.0f);
        compMakeup[1].store(1.0f);
        // Band 3 (High-Mid)
        compThreshold[2].store(-20.0f);
        compRatio[2].store(4.0f);
        compAttack[2].store(10.0f);
        compRelease[2].store(100.0f);
        compMakeup[2].store(1.0f);
        // Band 4 (High)
        compThreshold[3].store(-20.0f);
        compRatio[3].store(2.0f);
        compAttack[3].store(10.0f);
        compRelease[3].store(100.0f);
        compMakeup[3].store(1.0f);

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("PhantomMultibandComp: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        in_left = jack_port_register(client, "in_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        in_right = jack_port_register(client, "in_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_left || !in_right || !out_left || !out_right) {
            jack_client_close(client);
            throw std::runtime_error("PhantomMultibandComp: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomMultibandComp: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("PhantomMultibandComp: Failed to activate JACK client");
        }

        // Start the control thread.
        control_thread = std::thread(&PhantomMultibandComp::control_loop, this);

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[PhantomMultibandComp] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[PhantomMultibandComp] Default compressor parameters for bands:" << std::endl;
            for (int b = 0; b < 4; b++) {
                std::cout << "  Band " << (b + 1) << ": Threshold = " << compThreshold[b].load() << " dB, "
                    << "Ratio = " << compRatio[b].load() << ", "
                    << "Attack = " << compAttack[b].load() << " ms, "
                    << "Release = " << compRelease[b].load() << " ms, "
                    << "Makeup = " << compMakeup[b].load() << std::endl;
            }
        }
    }

    ~PhantomMultibandComp() {
        running.store(false);
        if (control_thread.joinable())
            control_thread.join();
        if (client)
            jack_client_close(client);
    }

    void run() {
        std::cout << "[PhantomMultibandComp] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[PhantomMultibandComp] Shutting down." << std::endl;
    }
};

int main() {
    try {
        PhantomMultibandComp multibandComp;
        multibandComp.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[PhantomMultibandComp] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

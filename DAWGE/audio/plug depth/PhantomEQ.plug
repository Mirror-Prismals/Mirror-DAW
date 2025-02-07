// OliveEQ.cpp
// A simple stereo three-band mastering EQ using JACK.
// Low Shelf at 200 Hz, Peaking at 1000 Hz, and High Shelf at 5000 Hz.
// Gains for each band are specified in dB.
// Compile with:
//   g++ -std=c++11 OliveEQ.cpp -ljack -lpthread -o OliveEQ

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <jack/jack.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <atomic>
#include <thread>
#include <mutex>
#include <chrono>
#include <stdexcept>

// ------------------ Biquad Filter Class ------------------
struct Biquad {
    // Coefficients (a0 normalized to 1)
    float b0, b1, b2, a1, a2;
    // State variables
    float x1, x2, y1, y2;
    Biquad() : b0(1), b1(0), b2(0), a1(0), a2(0),
        x1(0), x2(0), y1(0), y2(0) {
    }
    inline float process(float x) {
        float y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1;
        x1 = x;
        y2 = y1;
        y1 = y;
        return y;
    }
    inline void reset() { x1 = x2 = y1 = y2 = 0; }
};

// ------------------ Filter Coefficient Calculations ------------------
// Based on Robert Bristow-Johnson's Audio EQ Cookbook formulas

// Low-Shelf Filter Update
static void updateLowShelf(Biquad& bq, float fs, float f0, float dBgain, float Q) {
    float A = powf(10.0f, dBgain / 40.0f);
    float w0 = 2.0f * M_PI * f0 / fs;
    float cosw0 = cosf(w0);
    float sinw0 = sinf(w0);
    float S = 1.0f; // Shelf slope (fixed)
    float alpha = sinw0 / 2.0f * sqrtf((A + 1.0f / A) * (1.0f / S - 1.0f) + 2.0f);

    float b0 = A * ((A + 1.0f) - (A - 1.0f) * cosw0 + 2.0f * sqrtf(A) * alpha);
    float b1 = 2.0f * A * ((A - 1.0f) - (A + 1.0f) * cosw0);
    float b2 = A * ((A + 1.0f) - (A - 1.0f) * cosw0 - 2.0f * sqrtf(A) * alpha);
    float a0 = (A + 1.0f) + (A - 1.0f) * cosw0 + 2.0f * sqrtf(A) * alpha;
    float a1 = -2.0f * ((A - 1.0f) + (A + 1.0f) * cosw0);
    float a2 = (A + 1.0f) + (A - 1.0f) * cosw0 - 2.0f * sqrtf(A) * alpha;

    bq.b0 = b0 / a0;
    bq.b1 = b1 / a0;
    bq.b2 = b2 / a0;
    bq.a1 = a1 / a0;
    bq.a2 = a2 / a0;
}

// Peaking EQ Filter Update
static void updatePeaking(Biquad& bq, float fs, float f0, float dBgain, float Q) {
    float A = powf(10.0f, dBgain / 40.0f);
    float w0 = 2.0f * M_PI * f0 / fs;
    float cosw0 = cosf(w0);
    float sinw0 = sinf(w0);
    float alpha = sinw0 / (2.0f * Q);

    float b0 = 1.0f + alpha * A;
    float b1 = -2.0f * cosw0;
    float b2 = 1.0f - alpha * A;
    float a0 = 1.0f + alpha / A;
    float a1 = -2.0f * cosw0;
    float a2 = 1.0f - alpha / A;

    bq.b0 = b0 / a0;
    bq.b1 = b1 / a0;
    bq.b2 = b2 / a0;
    bq.a1 = a1 / a0;
    bq.a2 = a2 / a0;
}

// High-Shelf Filter Update
static void updateHighShelf(Biquad& bq, float fs, float f0, float dBgain, float Q) {
    float A = powf(10.0f, dBgain / 40.0f);
    float w0 = 2.0f * M_PI * f0 / fs;
    float cosw0 = cosf(w0);
    float sinw0 = sinf(w0);
    float S = 1.0f;
    float alpha = sinw0 / 2.0f * sqrtf((A + 1.0f / A) * (1.0f / S - 1.0f) + 2.0f);

    float b0 = A * ((A + 1.0f) + (A - 1.0f) * cosw0 + 2.0f * sqrtf(A) * alpha);
    float b1 = -2.0f * A * ((A - 1.0f) + (A + 1.0f) * cosw0);
    float b2 = A * ((A + 1.0f) + (A - 1.0f) * cosw0 - 2.0f * sqrtf(A) * alpha);
    float a0 = (A + 1.0f) - (A - 1.0f) * cosw0 + 2.0f * sqrtf(A) * alpha;
    float a1 = 2.0f * ((A - 1.0f) - (A + 1.0f) * cosw0);
    float a2 = (A + 1.0f) - (A - 1.0f) * cosw0 - 2.0f * sqrtf(A) * alpha;

    bq.b0 = b0 / a0;
    bq.b1 = b1 / a0;
    bq.b2 = b2 / a0;
    bq.a1 = a1 / a0;
    bq.a2 = a2 / a0;
}

// ------------------ OliveEQ Class (Mastering EQ) ------------------
class OliveEQ {
private:
    jack_client_t* client;
    // Stereo ports
    jack_port_t* in_left, * in_right, * out_left, * out_right;
    int sample_rate;

    std::atomic<bool> running;
    std::thread control_thread;
    std::mutex print_mutex;

    // EQ parameters (in dB) for each band
    std::atomic<float> lowGain;  // Low-shelf gain (default 0 dB)
    std::atomic<float> midGain;  // Mid peak gain (default 0 dB)
    std::atomic<float> highGain; // High-shelf gain (default 0 dB)

    // Fixed frequencies and Q factors (you can tweak these as desired)
    const float lowFreq = 200.0f;
    const float midFreq = 1000.0f;
    const float highFreq = 5000.0f;
    const float lowQ = 0.707f;
    const float midQ = 1.0f;
    const float highQ = 0.707f;

    // For each channel, we have three Biquad filters in series.
    // Left channel filters:
    Biquad leftLow, leftMid, leftHigh;
    // Right channel filters:
    Biquad rightLow, rightMid, rightHigh;

    // Update filter coefficients for a block (for both channels)
    void updateFilters() {
        float fs = static_cast<float>(sample_rate);
        // Load current gain settings from atomics
        float low_dB = lowGain.load();
        float mid_dB = midGain.load();
        float high_dB = highGain.load();
        // Update left channel filters
        updateLowShelf(leftLow, fs, lowFreq, low_dB, lowQ);
        updatePeaking(leftMid, fs, midFreq, mid_dB, midQ);
        updateHighShelf(leftHigh, fs, highFreq, high_dB, highQ);
        // Update right channel filters
        updateLowShelf(rightLow, fs, lowFreq, low_dB, lowQ);
        updatePeaking(rightMid, fs, midFreq, mid_dB, midQ);
        updateHighShelf(rightHigh, fs, highFreq, high_dB, highQ);
    }

    // JACK process callback: applies the EQ to stereo audio.
    static int process_callback(jack_nframes_t nframes, void* arg) {
        OliveEQ* eq = static_cast<OliveEQ*>(arg);
        float* inL = static_cast<float*>(jack_port_get_buffer(eq->in_left, nframes));
        float* inR = static_cast<float*>(jack_port_get_buffer(eq->in_right, nframes));
        float* outL = static_cast<float*>(jack_port_get_buffer(eq->out_left, nframes));
        float* outR = static_cast<float*>(jack_port_get_buffer(eq->out_right, nframes));

        // Update the filter coefficients at the beginning of the block
        eq->updateFilters();

        // Process each sample for left and right channels separately
        for (jack_nframes_t i = 0; i < nframes; i++) {
            // Left channel processing through low, mid, high filters in series
            float sampleL = inL[i];
            sampleL = eq->leftLow.process(sampleL);
            sampleL = eq->leftMid.process(sampleL);
            sampleL = eq->leftHigh.process(sampleL);
            outL[i] = sampleL;

            // Right channel processing
            float sampleR = inR[i];
            sampleR = eq->rightLow.process(sampleR);
            sampleR = eq->rightMid.process(sampleR);
            sampleR = eq->rightHigh.process(sampleR);
            outR[i] = sampleR;
        }
        return 0;
    }

    // Control thread: allows real-time adjustment of low, mid, and high gains.
    void control_loop() {
        std::string line;
        while (running.load()) {
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "\n[OliveEQ] Enter new gains for low, mid, and high bands (in dB), e.g., \"3.0 -2.0 4.0\", or type 'q' to quit: ";
            }
            if (!std::getline(std::cin, line))
                break;
            if (line == "q" || line == "Q") {
                running.store(false);
                break;
            }
            std::istringstream iss(line);
            float newLow, newMid, newHigh;
            if (!(iss >> newLow >> newMid >> newHigh)) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[OliveEQ] Invalid input. Please try again." << std::endl;
                continue;
            }
            lowGain.store(newLow);
            midGain.store(newMid);
            highGain.store(newHigh);
            {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[OliveEQ] Updated gains: low = " << newLow << " dB, mid = " << newMid << " dB, high = " << newHigh << " dB" << std::endl;
            }
        }
    }

public:
    OliveEQ(const char* client_name = "OliveEQ")
        : client(nullptr), running(true)
    {
        // Initialize EQ gains to 0 dB (unity gain)
        lowGain.store(0.0f);
        midGain.store(0.0f);
        highGain.store(0.0f);

        jack_status_t status;
        client = jack_client_open(client_name, JackNullOption, &status);
        if (!client) {
            throw std::runtime_error("OliveEQ: Failed to open JACK client");
        }
        sample_rate = jack_get_sample_rate(client);

        // Register stereo ports
        in_left = jack_port_register(client, "in_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        in_right = jack_port_register(client, "in_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        out_left = jack_port_register(client, "out_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        out_right = jack_port_register(client, "out_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        if (!in_left || !in_right || !out_left || !out_right) {
            jack_client_close(client);
            throw std::runtime_error("OliveEQ: Failed to register JACK ports");
        }

        if (jack_set_process_callback(client, process_callback, this) != 0) {
            jack_client_close(client);
            throw std::runtime_error("OliveEQ: Failed to set process callback");
        }

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("OliveEQ: Failed to activate JACK client");
        }

        // Start the control thread.
        control_thread = std::thread(&OliveEQ::control_loop, this);

        // (Optional) Reset filter states.
        leftLow.reset(); leftMid.reset(); leftHigh.reset();
        rightLow.reset(); rightMid.reset(); rightHigh.reset();

        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[OliveEQ] Initialized. Sample rate: " << sample_rate << " Hz" << std::endl;
            std::cout << "[OliveEQ] Default gains: low = " << lowGain.load() << " dB, mid = " << midGain.load() << " dB, high = " << highGain.load() << " dB" << std::endl;
        }
    }

    ~OliveEQ() {
        running.store(false);
        if (control_thread.joinable()) {
            control_thread.join();
        }
        if (client) {
            jack_client_close(client);
        }
    }

    void run() {
        std::cout << "[OliveEQ] Running. Type 'q' in the control console to quit." << std::endl;
        while (running.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        std::cout << "[OliveEQ] Shutting down." << std::endl;
    }
};

int main() {
    try {
        OliveEQ eq;
        eq.run();
    }
    catch (const std::exception& e) {
        std::cerr << "[OliveEQ] Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

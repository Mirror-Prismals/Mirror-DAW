#include <iostream>
#include <jack/jack.h>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <direct.h>
#include <algorithm>

struct Track {
    std::string filename;
    std::string role;
    std::vector<float> data;
};

class MirrorDaw {
private:
    std::vector<Track> tracks;
    jack_client_t* client;
    jack_port_t* input_port;
    jack_port_t* output_left;
    jack_port_t* output_right;
    bool is_recording;
    bool is_playing;
    std::vector<float> current_recording;
    std::map<std::string, std::vector<float>> playback_buffers;
    size_t playback_position;
    std::string tracks_dir;
    std::vector<std::string> valid_roles;

    static int process_callback(jack_nframes_t nframes, void* arg) {
        MirrorDaw* daw = static_cast<MirrorDaw*>(arg);

        // Handle recording
        if (daw->is_recording) {
            float* in = (float*)jack_port_get_buffer(daw->input_port, nframes);
            if (in != nullptr) {
                daw->current_recording.insert(
                    daw->current_recording.end(),
                    in,
                    in + nframes
                );
            }
        }

        // Handle playback
        if (daw->is_playing) {
            float* out_left = (float*)jack_port_get_buffer(daw->output_left, nframes);
            float* out_right = (float*)jack_port_get_buffer(daw->output_right, nframes);

            // Clear output buffers
            std::fill(out_left, out_left + nframes, 0.0f);
            std::fill(out_right, out_right + nframes, 0.0f);

            for (size_t i = 0; i < nframes; i++) {
                if (daw->playback_position >= daw->playback_buffers["left"].size()) {
                    daw->is_playing = false;
                    break;
                }

                // Left channel goes fully left
                out_left[i] += daw->playback_buffers["left"][daw->playback_position];

                // Right channel goes fully right
                out_right[i] += daw->playback_buffers["right"][daw->playback_position];

                // Front fill and sub are centered (split equally between L/R)
                float center_mix = (daw->playback_buffers["front_fill"][daw->playback_position] +
                    daw->playback_buffers["sub"][daw->playback_position]) * 0.5f;
                out_left[i] += center_mix;
                out_right[i] += center_mix;

                daw->playback_position++;
            }
        }

        return 0;
    }

    void save_track_to_wav(const std::string& filename, const std::vector<float>& data) {
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to create file: " + filename);
        }

        struct WavHeader {
            char riff[4] = { 'R', 'I', 'F', 'F' };
            uint32_t fileSize;
            char wave[4] = { 'W', 'A', 'V', 'E' };
            char fmt[4] = { 'f', 'm', 't', ' ' };
            uint32_t fmtSize = 16;
            uint16_t audioFormat = 3;
            uint16_t numChannels = 1;
            uint32_t sampleRate;
            uint32_t byteRate;
            uint16_t blockAlign;
            uint16_t bitsPerSample = 32;
            char data[4] = { 'd', 'a', 't', 'a' };
            uint32_t dataSize;
        } header;

        header.sampleRate = jack_get_sample_rate(client);
        header.byteRate = header.sampleRate * sizeof(float);
        header.blockAlign = sizeof(float);
        header.dataSize = data.size() * sizeof(float);
        header.fileSize = header.dataSize + sizeof(WavHeader) - 8;

        file.write(reinterpret_cast<char*>(&header), sizeof(WavHeader));
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
    }

public:
    MirrorDaw() : is_recording(false), is_playing(false), playback_position(0), tracks_dir("tracks") {
        valid_roles = { "left", "front_fill", "sub", "right" };
        _mkdir(tracks_dir.c_str());

        jack_status_t status;
        client = jack_client_open("MirrorDaw", JackNullOption, &status);
        if (client == nullptr) {
            throw std::runtime_error("Failed to connect to JACK server");
        }

        input_port = jack_port_register(client, "input", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
        output_left = jack_port_register(client, "output_left", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
        output_right = jack_port_register(client, "output_right", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);

        if (input_port == nullptr || output_left == nullptr || output_right == nullptr) {
            jack_client_close(client);
            throw std::runtime_error("Failed to create ports");
        }

        jack_set_process_callback(client, process_callback, this);

        if (jack_activate(client) != 0) {
            jack_client_close(client);
            throw std::runtime_error("Failed to activate JACK client");
        }
    }

    ~MirrorDaw() {
        if (client) {
            jack_client_close(client);
        }
    }

    void record_track(int track_number, const std::string& role) {
        if (std::find(valid_roles.begin(), valid_roles.end(), role) == valid_roles.end()) {
            throw std::runtime_error("Invalid role: " + role);
        }

        std::string filename = tracks_dir + "/track_" + std::to_string(track_number) + "_" + role + ".wav";

        current_recording.clear();
        is_recording = true;
        std::cout << "Recording Track " << track_number << " as role '" << role << "'..." << std::endl;
        std::cout << "Press Enter to stop recording..." << std::endl;
        std::cin.get();
        is_recording = false;

        Track new_track;
        new_track.filename = filename;
        new_track.role = role;
        new_track.data = current_recording;
        tracks.push_back(new_track);

        save_track_to_wav(filename, current_recording);
        std::cout << "Track " << track_number << " saved as " << filename << std::endl;
    }

    void mix_tracks_to_stems() {
        std::map<std::string, std::vector<float>> role_buffers;

        for (const auto& role : valid_roles) {
            role_buffers[role] = std::vector<float>();
        }

        size_t max_length = 0;
        for (const auto& track : tracks) {
            max_length = (track.data.size() > max_length) ? track.data.size() : max_length;
        }

        for (auto& buffer : role_buffers) {
            buffer.second.resize(max_length, 0.0f);
        }

        for (const auto& track : tracks) {
            auto& role_buffer = role_buffers[track.role];
            for (size_t i = 0; i < track.data.size(); i++) {
                role_buffer[i] += track.data[i];
            }
        }

        for (const auto& role : valid_roles) {
            std::string stem_filename = tracks_dir + "/stem_" + role + ".wav";
            save_track_to_wav(stem_filename, role_buffers[role]);
            std::cout << "Exported stem: " << stem_filename << std::endl;
        }
    }

    void play_master() {
        // Load all stems
        playback_buffers.clear();
        for (const auto& role : valid_roles) {
            std::string stem_filename = tracks_dir + "/stem_" + role + ".wav";
            std::vector<float> data;

            std::ifstream file(stem_filename, std::ios::binary);
            if (!file) {
                std::cout << "Warning: Could not find stem for " << role << std::endl;
                // Create empty buffer if stem doesn't exist
                playback_buffers[role] = std::vector<float>();
                continue;
            }

            // Skip WAV header
            file.seekg(44); // Standard WAV header size

            // Read the data
            while (file) {
                float sample;
                file.read(reinterpret_cast<char*>(&sample), sizeof(float));
                if (file) {
                    data.push_back(sample);
                }
            }
            playback_buffers[role] = data;
        }

        // Ensure all buffers are the same length
        size_t max_length = 0;
        for (const auto& buffer : playback_buffers) {
            max_length = (buffer.second.size() > max_length) ? buffer.second.size() : max_length;
        }
        for (auto& buffer : playback_buffers) {
            buffer.second.resize(max_length, 0.0f);
        }

        // Start playback
        playback_position = 0;
        is_playing = true;
        std::cout << "Playing master mix... Press Enter to stop." << std::endl;
        std::cin.get();
        is_playing = false;
    }
};

int main() {
    try {
        MirrorDaw daw;
        std::cout << "MirrorDaw initialized. Available commands:" << std::endl;
        std::cout << "1. Record track" << std::endl;
        std::cout << "2. Mix to stems" << std::endl;
        std::cout << "3. Exit" << std::endl;
        std::cout << "4. Play master mix" << std::endl;

        while (true) {
            std::cout << "\nEnter command (1-4): ";
            int choice;
            std::cin >> choice;
            std::cin.ignore();

            if (choice == 1) {
                int track_number;
                std::string role;
                std::cout << "Enter track number: ";
                std::cin >> track_number;
                std::cin.ignore();
                std::cout << "Enter role (left/front_fill/sub/right): ";
                std::getline(std::cin, role);
                daw.record_track(track_number, role);
            }
            else if (choice == 2) {
                daw.mix_tracks_to_stems();
            }
            else if (choice == 3) {
                break;
            }
            else if (choice == 4) {
                daw.play_master();
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

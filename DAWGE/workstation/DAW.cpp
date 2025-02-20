#include <iostream>
#include <jack/jack.h>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <direct.h>
#include <algorithm>
#include <stdexcept>
#include <thread>
#include <atomic>
#include <chrono>

// OpenGL / GLFW includes
#include <glad/glad.h>
#include <GLFW/glfw3.h>

//-----------------------------------------------------
// Data Structures for DAW
//-----------------------------------------------------

// Modified Track now includes a track number.
struct Track {
    int trackNumber;          // Track number relative to the lane.
    std::string role;
    std::vector<float> data;
    std::string filename;
};

// A Lane holds a lane number and a collection of clips.
struct Lane {
    int laneNumber;
    std::vector<Track> clips;
};

//-----------------------------------------------------
// Backend DAW Code (with lanes)
//-----------------------------------------------------

class MirrorDaw {
private:
    // Instead of a flat vector of tracks, we now have lanes.
    std::vector<Lane> lanes;
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

    // JACK process callback for recording and playback
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

        // Handle playback (unchanged)
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

                // Left channel
                out_left[i] += daw->playback_buffers["left"][daw->playback_position];
                // Right channel
                out_right[i] += daw->playback_buffers["right"][daw->playback_position];

                // Center channels (front_fill and sub, split equally)
                float center_mix = (daw->playback_buffers["front_fill"][daw->playback_position] +
                    daw->playback_buffers["sub"][daw->playback_position]) * 0.5f;
                out_left[i] += center_mix;
                out_right[i] += center_mix;

                daw->playback_position++;
            }
        }

        return 0;
    }

    // Save a vector of floats as a WAV file.
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
            uint16_t audioFormat = 3; // float format
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

    // Create a new lane (CLI option 5).
    void create_lane(int laneNumber) {
        auto it = std::find_if(lanes.begin(), lanes.end(), [laneNumber](const Lane& lane) {
            return lane.laneNumber == laneNumber;
            });
        if (it != lanes.end()) {
            throw std::runtime_error("Lane " + std::to_string(laneNumber) + " already exists.");
        }
        Lane newLane;
        newLane.laneNumber = laneNumber;
        lanes.push_back(newLane);
    }

    // Record a new clip into a specific lane with a given track number and role.
    void record_track(int laneNumber, int trackNumber, const std::string& role) {
        if (std::find(valid_roles.begin(), valid_roles.end(), role) == valid_roles.end()) {
            throw std::runtime_error("Invalid role: " + role);
        }

        // Find the target lane.
        auto laneIt = std::find_if(lanes.begin(), lanes.end(), [laneNumber](const Lane& lane) {
            return lane.laneNumber == laneNumber;
            });
        if (laneIt == lanes.end()) {
            throw std::runtime_error("Lane " + std::to_string(laneNumber) + " does not exist. Create it first.");
        }

        std::string filename = tracks_dir + "/lane_" + std::to_string(laneNumber) +
            "_track_" + std::to_string(trackNumber) + "_" + role + ".wav";

        current_recording.clear();
        is_recording = true;
        std::cout << "Recording Lane " << laneNumber << " Track " << trackNumber
            << " as role '" << role << "'..." << std::endl;
        std::cout << "Press Enter to stop recording..." << std::endl;
        std::cin.get();
        is_recording = false;

        // Look for an existing clip with the same track number and role in this lane.
        auto clipIt = std::find_if(laneIt->clips.begin(), laneIt->clips.end(),
            [trackNumber, &role](const Track& clip) {
                return clip.trackNumber == trackNumber && clip.role == role;
            });
        if (clipIt != laneIt->clips.end()) {
            // Overwrite existing clip.
            clipIt->data = current_recording;
            clipIt->filename = filename;
        }
        else {
            // Add new clip.
            Track newClip;
            newClip.trackNumber = trackNumber;
            newClip.role = role;
            newClip.data = current_recording;
            newClip.filename = filename;
            laneIt->clips.push_back(newClip);
        }
        save_track_to_wav(filename, current_recording);
        std::cout << "Lane " << laneNumber << " Track " << trackNumber
            << " saved as " << filename << std::endl;
    }

    // Mix all clips (from all lanes) into stems per role.
    void mix_tracks_to_stems() {
        std::map<std::string, std::vector<float>> role_buffers;
        for (const auto& role : valid_roles) {
            role_buffers[role] = std::vector<float>();
        }

        // First, find the maximum length among all clips (across lanes).
        size_t max_length = 0;
        for (const auto& lane : lanes) {
            for (const auto& clip : lane.clips) {
                if (clip.data.size() > max_length)
                    max_length = clip.data.size();
            }
        }
        // Resize buffers to the max length.
        for (auto& buffer : role_buffers) {
            buffer.second.resize(max_length, 0.0f);
        }
        // Sum clips from all lanes into the appropriate role buffer.
        for (const auto& lane : lanes) {
            for (const auto& clip : lane.clips) {
                auto& buffer = role_buffers[clip.role];
                for (size_t i = 0; i < clip.data.size(); i++) {
                    buffer[i] += clip.data[i];
                }
            }
        }

        // Save stems.
        for (const auto& role : valid_roles) {
            std::string stem_filename = tracks_dir + "/stem_" + role + ".wav";
            save_track_to_wav(stem_filename, role_buffers[role]);
            std::cout << "Exported stem: " << stem_filename << std::endl;
        }
    }

    // Play the master mix from the stored stems.
    void play_master() {
        // Load all stems.
        playback_buffers.clear();
        for (const auto& role : valid_roles) {
            std::string stem_filename = tracks_dir + "/stem_" + role + ".wav";
            std::vector<float> data;
            std::ifstream file(stem_filename, std::ios::binary);
            if (!file) {
                std::cout << "Warning: Could not find stem for " << role << std::endl;
                playback_buffers[role] = std::vector<float>();
                continue;
            }
            // Skip WAV header (44 bytes)
            file.seekg(44);
            while (file) {
                float sample;
                file.read(reinterpret_cast<char*>(&sample), sizeof(float));
                if (file)
                    data.push_back(sample);
            }
            playback_buffers[role] = data;
        }

        // Make sure all buffers are of equal length.
        size_t max_length = 0;
        for (const auto& buffer : playback_buffers) {
            if (buffer.second.size() > max_length)
                max_length = buffer.second.size();
        }
        for (auto& buffer : playback_buffers) {
            buffer.second.resize(max_length, 0.0f);
        }

        // Start playback.
        playback_position = 0;
        is_playing = true;
        std::cout << "Playing master mix... Press Enter to stop." << std::endl;
        std::cin.get();
        is_playing = false;
    }

    // PUBLIC GETTER for UI access (return lanes)
    const std::vector<Lane>& getLanes() const {
        return lanes;
    }

    // Returns the sample rate (retrieved from JACK)
    uint32_t getSampleRate() const {
        return jack_get_sample_rate(client);
    }

    // Returns the current playback position (in samples)
    size_t getPlaybackPosition() const {
        return playback_position;
    }
};

//-----------------------------------------------------
// UI Code using GLFW and OpenGL
//-----------------------------------------------------

class MirrorDawUI {
public:
    // The UI takes a reference to our backend DAW instance and a shared exit flag.
    MirrorDawUI(MirrorDaw& daw, std::atomic<bool>& exitFlag)
        : daw(daw), exitFlag(exitFlag) {
        if (!glfwInit()) {
            throw std::runtime_error("Failed to initialize GLFW");
        }
        // Remove window decorations (for windowed fullscreen)
        glfwWindowHint(GLFW_DECORATED, GLFW_FALSE);

        // Retrieve primary monitor's resolution.
        GLFWmonitor* monitor = glfwGetPrimaryMonitor();
        const GLFWvidmode* mode = glfwGetVideoMode(monitor);

        // Create a window the size of the primary monitor.
        window = glfwCreateWindow(mode->width, mode->height, "MirrorDaw - Timeline", nullptr, nullptr);
        glfwSetWindowPos(window, 0, 0);

        if (!window) {
            glfwTerminate();
            throw std::runtime_error("Failed to create GLFW window");
        }
        glfwMakeContextCurrent(window);
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            throw std::runtime_error("Failed to initialize GLAD");
        }
        // Set up a 2D orthographic projection (we use normalized device coordinates)
        glViewport(0, 0, mode->width, mode->height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1, 1, -1, 1, -1, 1);
        glMatrixMode(GL_MODELVIEW);
    }

    ~MirrorDawUI() {
        glfwDestroyWindow(window);
        glfwTerminate();
    }

    // Run the render loop.
    void run() {
        while (!glfwWindowShouldClose(window)) {
            if (exitFlag.load()) {
                glfwSetWindowShouldClose(window, true);
            }
            glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            drawTimeline();
            glfwSwapBuffers(window);
            glfwPollEvents();
            std::this_thread::sleep_for(std::chrono::milliseconds(16));
        }
    }

private:
    GLFWwindow* window;
    MirrorDaw& daw;
    std::atomic<bool>& exitFlag;

    // Draw a timeline showing each lane’s summed waveform along with a playhead.
    void drawTimeline() {
        const std::vector<Lane>& lanes = daw.getLanes();
        size_t numLanes = lanes.size();
        if (numLanes == 0)
            return;

        uint32_t sampleRate = daw.getSampleRate();

        // Divide vertical space equally among lanes (NDC from -1 to 1)
        float laneHeight = 2.0f / numLanes;

        for (size_t i = 0; i < numLanes; i++) {
            const Lane& lane = lanes[i];
            // Calculate vertical bounds for this lane.
            float top = 1.0f - i * laneHeight;
            float bottom = 1.0f - (i + 1) * laneHeight;
            float midY = (top + bottom) / 2.0f;
            float halfHeight = laneHeight / 2.0f;

            // Determine maximum length among clips in this lane.
            size_t laneMaxSamples = 0;
            for (const auto& clip : lane.clips) {
                if (clip.data.size() > laneMaxSamples)
                    laneMaxSamples = clip.data.size();
            }
            // Prepare a summed waveform.
            std::vector<float> summed(laneMaxSamples, 0.0f);
            for (const auto& clip : lane.clips) {
                for (size_t j = 0; j < clip.data.size(); j++) {
                    summed[j] += clip.data[j];
                }
            }

            // Compute lane duration (in seconds) for mapping time.
            float laneDuration = static_cast<float>(laneMaxSamples) / sampleRate;
            if (laneDuration <= 0.0f)
                laneDuration = 1.0f;

            // Draw the summed waveform.
            glBegin(GL_LINE_STRIP);
            glColor3f(0.0f, 1.0f, 0.0f); // green waveform
            size_t step = laneMaxSamples / 1000;
            if (step < 1) step = 1;
            for (size_t j = 0; j < laneMaxSamples; j += step) {
                float timeSec = static_cast<float>(j) / sampleRate;
                float x = -1.0f + 2.0f * (timeSec / laneDuration);
                float y = midY + (summed[j] * halfHeight);
                glVertex2f(x, y);
            }
            glEnd();

            // Draw border around lane area.
            glBegin(GL_LINE_LOOP);
            glColor3f(0.7f, 0.7f, 0.7f);
            glVertex2f(-1.0f, top);
            glVertex2f(1.0f, top);
            glVertex2f(1.0f, bottom);
            glVertex2f(-1.0f, bottom);
            glEnd();

            // If lane is empty, draw an empty blue rectangle.
            if (lane.clips.empty()) {
                glColor4f(0.0f, 0.0f, 1.0f, 0.3f);
                glBegin(GL_QUADS);
                glVertex2f(-1.0f, top);
                glVertex2f(1.0f, top);
                glVertex2f(1.0f, bottom);
                glVertex2f(-1.0f, bottom);
                glEnd();
            }
        }

        // Draw playhead (global, based on master mix duration).
        size_t playPos = daw.getPlaybackPosition();
        float playTimeSec = static_cast<float>(playPos) / sampleRate;
        // Determine global maximum duration across all lanes.
        float globalDuration = 0.0f;
        for (const auto& lane : lanes) {
            size_t laneMax = 0;
            for (const auto& clip : lane.clips)
                if (clip.data.size() > laneMax)
                    laneMax = clip.data.size();
            float dur = static_cast<float>(laneMax) / sampleRate;
            if (dur > globalDuration)
                globalDuration = dur;
        }
        if (globalDuration <= 0.0f)
            globalDuration = 1.0f;
        float playX = -1.0f + 2.0f * (playTimeSec / globalDuration);
        glBegin(GL_LINES);
        glColor3f(1.0f, 0.0f, 0.0f); // red playhead
        glVertex2f(playX, -1.0f);
        glVertex2f(playX, 1.0f);
        glEnd();
    }
};

//-----------------------------------------------------
// Main: CLI Controller and UI Thread
//-----------------------------------------------------

std::atomic<bool> g_exitFlag(false);

int main() {
    try {
        MirrorDaw daw;

        // Start the CLI in a separate thread.
        std::thread cliThread([&daw]() {
            std::cout << "MirrorDaw CLI Controller initialized. Available commands:" << std::endl;
            std::cout << "1. Record clip" << std::endl;
            std::cout << "2. Mix to stems" << std::endl;
            std::cout << "3. Exit" << std::endl;
            std::cout << "4. Play master mix" << std::endl;
            std::cout << "5. Create lane" << std::endl;
            while (true) {
                std::cout << "\nEnter command (1-5): ";
                int choice;
                std::cin >> choice;
                std::cin.ignore();
                if (choice == 5) {
                    int laneNumber;
                    std::cout << "Enter new lane number: ";
                    std::cin >> laneNumber;
                    std::cin.ignore();
                    try {
                        daw.create_lane(laneNumber);
                        std::cout << "Lane " << laneNumber << " created." << std::endl;
                    }
                    catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << std::endl;
                    }
                }
                else if (choice == 1) {
                    int laneNumber, trackNumber;
                    std::string role;
                    std::cout << "Enter lane number: ";
                    std::cin >> laneNumber;
                    std::cin.ignore();
                    std::cout << "Enter track number (within lane): ";
                    std::cin >> trackNumber;
                    std::cin.ignore();
                    std::cout << "Enter role (left/front_fill/sub/right): ";
                    std::getline(std::cin, role);
                    try {
                        daw.record_track(laneNumber, trackNumber, role);
                    }
                    catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << std::endl;
                    }
                }
                else if (choice == 2) {
                    try {
                        daw.mix_tracks_to_stems();
                    }
                    catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << std::endl;
                    }
                }
                else if (choice == 3) {
                    g_exitFlag.store(true);
                    break;
                }
                else if (choice == 4) {
                    try {
                        daw.play_master();
                    }
                    catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << std::endl;
                    }
                }
            }
            });

        // Run the UI in the main thread.
        MirrorDawUI ui(daw, g_exitFlag);
        ui.run();

        cliThread.join();
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

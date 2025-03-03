// SoundGalay.cpp
// Note: Removed GLM_ENABLE_EXPERIMENTAL and the experimental header.

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <jack/jack.h>
#include <jack/types.h>

#include <filesystem>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <cfloat>
#include <cmath>

// --- Minimal WAV Parser Structures ---
struct WAVHeader {
    char chunkId[4];       // "RIFF"
    uint32_t chunkSize;
    char format[4];        // "WAVE"
};

struct WAVFormat {
    char subchunk1Id[4];   // "fmt "
    uint32_t subchunk1Size;
    uint16_t audioFormat;  // PCM = 1
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
};

struct WAVData {
    char subchunk2Id[4];   // "data"
    uint32_t subchunk2Size;
};

// --- Minimal WAV Loader ---
// Reads a 16-bit PCM WAV file and converts samples to floats in [-1, 1].
bool load_wav_manual(const std::string& filepath, std::vector<float>& data, int& sampleRate) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        std::cerr << "Unable to open WAV file: " << filepath << std::endl;
        return false;
    }

    WAVHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));
    if (std::strncmp(header.chunkId, "RIFF", 4) != 0 ||
        std::strncmp(header.format, "WAVE", 4) != 0) {
        std::cerr << "Invalid WAV file format: " << filepath << std::endl;
        return false;
    }

    WAVFormat fmt;
    file.read(reinterpret_cast<char*>(&fmt), sizeof(WAVFormat));
    if (std::strncmp(fmt.subchunk1Id, "fmt ", 4) != 0) {
        std::cerr << "Missing fmt chunk in WAV file: " << filepath << std::endl;
        return false;
    }
    if (fmt.audioFormat != 1) { // Only PCM supported
        std::cerr << "Unsupported audio format (only PCM supported): " << filepath << std::endl;
        return false;
    }
    sampleRate = fmt.sampleRate;

    // Skip any extra bytes in the fmt chunk.
    if (fmt.subchunk1Size > 16) {
        file.seekg(fmt.subchunk1Size - 16, std::ios::cur);
    }

    // Read chunks until "data" chunk is found.
    WAVData dataHeader;
    while (true) {
        file.read(reinterpret_cast<char*>(&dataHeader), sizeof(WAVData));
        if (std::strncmp(dataHeader.subchunk2Id, "data", 4) == 0)
            break;
        file.seekg(dataHeader.subchunk2Size, std::ios::cur);
    }

    int bytesPerSample = fmt.bitsPerSample / 8;
    int numSamples = dataHeader.subchunk2Size / bytesPerSample;
    std::vector<int16_t> rawData(numSamples);
    file.read(reinterpret_cast<char*>(rawData.data()), dataHeader.subchunk2Size);

    data.resize(numSamples);
    for (int i = 0; i < numSamples; i++) {
        data[i] = rawData[i] / 32768.0f;
    }
    return true;
}

// --- Data Structures for Samples and Audio ---
struct SoundSample {
    std::string filename;
    std::vector<float> audioData; // Loaded on demand
    int sampleRate = 0;
    glm::vec3 avgColor;   // Computed from 6 random colors
    glm::vec2 position;   // Assigned via clustering or manual dragging
    float animTimer = 0.0f;  // Animation timer for reactive effect (in seconds)
    glm::vec2 velocity = glm::vec2(0.0f, 0.0f); // For gravity simulation
};

const float EFFECT_DURATION = 0.5f; // Duration of the click effect

// Global variables for controlling sample rate and debug mode.
int g_customSampleRate = 441000; // Playback sample rate override (in Hz)
bool g_debug = false;            // Toggle debug mode

// --- Gradient Control Points ---
// Each control point maps a sample rate to a color.
struct ControlPoint {
    int rate;         // Sample rate in Hz
    glm::vec3 color;  // Color (RGB components in 0-1)
};

std::vector<ControlPoint> controlPoints = {
    {4000,   glm::vec3(0.0f,       1.0f,       0.0f)},           // #00FF00
    {96000,  glm::vec3(128 / 255.0f, 0.0f,       255 / 255.0f)},     // #8000FF
    {220500, glm::vec3(0.0f,       0.0f,       255 / 255.0f)},       // #0000FF
    {384000, glm::vec3(128 / 255.0f, 128 / 255.0f, 255 / 255.0f)},     // #8080FF
    {441000, glm::vec3(161 / 255.0f, 0.0f,       161 / 255.0f)}        // #A100A1
};

// Global debug square color (will be computed from the gradient)
glm::vec3 g_debugColor;

// Compute the debug square color by interpolating between control points.
void update_debug_color() {
    // Clamp g_customSampleRate within the range of our control points.
    if (g_customSampleRate < controlPoints.front().rate)
        g_customSampleRate = controlPoints.front().rate;
    if (g_customSampleRate > controlPoints.back().rate)
        g_customSampleRate = controlPoints.back().rate;

    // Find the segment in which the current sample rate falls.
    ControlPoint cp1 = controlPoints.front(), cp2 = controlPoints.back();
    for (size_t i = 0; i < controlPoints.size() - 1; i++) {
        if (g_customSampleRate >= controlPoints[i].rate && g_customSampleRate <= controlPoints[i + 1].rate) {
            cp1 = controlPoints[i];
            cp2 = controlPoints[i + 1];
            break;
        }
    }
    // Compute interpolation factor t between cp1 and cp2.
    float t = float(g_customSampleRate - cp1.rate) / float(cp2.rate - cp1.rate);
    // Linearly interpolate between cp1.color and cp2.color.
    g_debugColor = cp1.color * (1.0f - t) + cp2.color * t;
}

// Global audio structure now uses a floating point position and stores original sample rate.
struct GlobalAudio {
    std::vector<float> buffer;
    double pos = 0.0;             // Floating-point playback position
    bool playing = false;
    int originalSampleRate = 44100;  // The sample rate of the loaded sample
};

GlobalAudio g_audio;

// --- JACK Globals and Process Callback ---
jack_client_t* jack_client = nullptr;
jack_port_t* output_port = nullptr;

// Updated process callback to use a fractional position.
int jack_process(jack_nframes_t nframes, void* arg) {
    float* out = (float*)jack_port_get_buffer(output_port, nframes);
    if (g_audio.playing) {
        // Compute the playback step: how much to advance per JACK frame.
        double step = static_cast<double>(g_customSampleRate) / g_audio.originalSampleRate;
        for (jack_nframes_t i = 0; i < nframes; i++) {
            int index = static_cast<int>(g_audio.pos);
            if (index < g_audio.buffer.size()) {
                out[i] = g_audio.buffer[index];
            }
            else {
                out[i] = 0.0f;
                g_audio.playing = false;
                break;
            }
            g_audio.pos += step;
        }
    }
    else {
        for (jack_nframes_t i = 0; i < nframes; i++)
            out[i] = 0.0f;
    }
    return 0;
}

bool init_jack() {
    jack_client = jack_client_open("SoundGalaxy", JackNullOption, nullptr);
    if (!jack_client) {
        std::cerr << "Failed to connect to JACK." << std::endl;
        return false;
    }
    jack_set_process_callback(jack_client, jack_process, nullptr);
    output_port = jack_port_register(jack_client, "output", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
    if (!output_port) {
        std::cerr << "Failed to register JACK output port." << std::endl;
        return false;
    }
    if (jack_activate(jack_client)) {
        std::cerr << "Cannot activate JACK client." << std::endl;
        return false;
    }
    return true;
}

// --- KMeans Clustering ---
void kmeans_cluster(std::vector<SoundSample>& samples, int k = 6, int iterations = 10) {
    int n = samples.size();
    if (n == 0) return;
    std::vector<glm::vec3> centroids(k);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, n - 1);
    for (int i = 0; i < k; i++) {
        centroids[i] = samples[dis(gen)].avgColor;
    }
    std::vector<int> labels(n, 0);
    for (int iter = 0; iter < iterations; iter++) {
        for (int i = 0; i < n; i++) {
            float minDist = FLT_MAX;
            int best = 0;
            for (int j = 0; j < k; j++) {
                float dist = glm::distance(samples[i].avgColor, centroids[j]);
                if (dist < minDist) {
                    minDist = dist;
                    best = j;
                }
            }
            labels[i] = best;
        }
        std::vector<glm::vec3> newCentroids(k, glm::vec3(0.0f));
        std::vector<int> counts(k, 0);
        for (int i = 0; i < n; i++) {
            newCentroids[labels[i]] += samples[i].avgColor;
            counts[labels[i]]++;
        }
        for (int j = 0; j < k; j++) {
            if (counts[j] > 0)
                centroids[j] = newCentroids[j] / static_cast<float>(counts[j]);
        }
    }
    std::vector<glm::vec2> clusterCenters(k);
    std::uniform_int_distribution<> xdis(300, 700);
    std::uniform_int_distribution<> ydis(300, 500);
    for (int j = 0; j < k; j++) {
        clusterCenters[j] = glm::vec2(xdis(gen), ydis(gen));
    }
    std::uniform_int_distribution<> offset(-100, 100);
    for (int i = 0; i < n; i++) {
        int cluster = labels[i];
        samples[i].position = clusterCenters[cluster] + glm::vec2(offset(gen), offset(gen));
    }
}

// --- Color Generation ---
void generate_colors(std::vector<SoundSample>& samples) {
    std::mt19937 gen(42); // fixed seed for consistency
    std::uniform_int_distribution<> dis(0, 255);
    for (auto& sample : samples) {
        glm::vec3 sum(0.0f);
        for (int i = 0; i < 6; i++) {
            glm::vec3 color(dis(gen), dis(gen), dis(gen));
            sum += color;
        }
        sample.avgColor = sum / 6.0f;
    }
}

// --- Sample Loading ---
void load_samples(const std::string& directory, std::vector<SoundSample>& samples) {
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.path().extension() == ".wav") {
            SoundSample sample;
            sample.filename = entry.path().string();
            samples.push_back(sample);
        }
    }
}

// --- OpenGL Utility: Draw a Circle ---
void draw_circle(const glm::vec2& center, float radius, const glm::vec3& color, int segments = 32) {
    glColor3f(color.r / 255.0f, color.g / 255.0f, color.b / 255.0f);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(center.x, center.y);
    for (int i = 0; i <= segments; i++) {
        float angle = i * 2.0f * 3.1415926f / segments;
        float x = center.x + cos(angle) * radius;
        float y = center.y + sin(angle) * radius;
        glVertex2f(x, y);
    }
    glEnd();
}

// --- Global Variables for View Transformation ---
float g_scale = 1.0f;
glm::vec2 g_offset(0.0f, 0.0f);
bool g_dragging = false;
double g_lastX = 0, g_lastY = 0;
std::vector<SoundSample> g_samples;

// Global variable to track which sample is currently being dragged
// -1 indicates no sample is being dragged.
int g_draggedSampleIndex = -1;

// --- GLFW Callbacks ---
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if (action == GLFW_PRESS) {
            g_dragging = true;
            glfwGetCursorPos(window, &g_lastX, &g_lastY);
        }
        else if (action == GLFW_RELEASE) {
            g_dragging = false;
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
        if (action == GLFW_PRESS) {
            double x, y;
            glfwGetCursorPos(window, &x, &y);
            float worldX = static_cast<float>(x) / g_scale - g_offset.x;
            float worldY = (700.0f - static_cast<float>(y)) / g_scale - g_offset.y;
            // Find a sample close enough to the cursor
            for (size_t i = 0; i < g_samples.size(); i++) {
                float dx = worldX - g_samples[i].position.x;
                float dy = worldY - g_samples[i].position.y;
                if (sqrt(dx * dx + dy * dy) <= 10.0f) {
                    g_draggedSampleIndex = i;
                    break;
                }
            }
        }
        else if (action == GLFW_RELEASE) {
            // Stop dragging the sample.
            g_draggedSampleIndex = -1;
        }
    }
    // Left-click for playing samples and triggering the reactive effect.
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        double x, y;
        glfwGetCursorPos(window, &x, &y);
        float worldX = static_cast<float>(x) / g_scale - g_offset.x;
        float worldY = (700.0f - static_cast<float>(y)) / g_scale - g_offset.y;
        for (auto& sample : g_samples) {
            float dx = worldX - sample.position.x;
            float dy = worldY - sample.position.y;
            if (sqrt(dx * dx + dy * dy) <= 10.0f) {
                // Load and play sample as before.
                if (sample.audioData.empty()) {
                    if (!load_wav_manual(sample.filename, sample.audioData, sample.sampleRate)) {
                        std::cerr << "Error loading audio: " << sample.filename << std::endl;
                        break;
                    }
                }
                g_audio.buffer = sample.audioData;
                g_audio.pos = 0.0;
                g_audio.playing = true;
                g_audio.originalSampleRate = sample.sampleRate;
                std::cout << "Now playing: " << sample.filename << std::endl;

                // Start the click effect: set the animation timer.
                sample.animTimer = EFFECT_DURATION;
                break;
            }
        }
    }
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    // Update pan if right button is dragging.
    if (g_dragging) {
        double dx = xpos - g_lastX;
        double dy = ypos - g_lastY;
        g_offset.x += static_cast<float>(dx) / g_scale;
        g_offset.y -= static_cast<float>(dy) / g_scale; // invert y
        g_lastX = xpos;
        g_lastY = ypos;
    }

    // If a sample is being dragged, update its position.
    if (g_draggedSampleIndex != -1) {
        // Convert the current cursor position to world coordinates.
        float worldX = static_cast<float>(xpos) / g_scale - g_offset.x;
        float worldY = (700.0f - static_cast<float>(ypos)) / g_scale - g_offset.y;
        g_samples[g_draggedSampleIndex].position = glm::vec2(worldX, worldY);
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    float scaleFactor = (yoffset > 0) ? 1.1f : 0.9f;
    g_scale *= scaleFactor;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    // Allow adjustments on key press and repeat events.
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
            g_audio.playing = false;
            std::cout << "Playback stopped." << std::endl;
        }
        else if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
            // Toggle debug mode on Q press.
            g_debug = !g_debug;
        }
        else if (key == GLFW_KEY_UP) {
            g_customSampleRate += 1000;
            if (g_customSampleRate > controlPoints.back().rate)
                g_customSampleRate = controlPoints.back().rate;
            update_debug_color();
        }
        else if (key == GLFW_KEY_DOWN) {
            g_customSampleRate -= 1000;
            if (g_customSampleRate < controlPoints.front().rate)
                g_customSampleRate = controlPoints.front().rate;
            update_debug_color();
        }
    }
}

int main() {
    // Initialize debug color based on the default sample rate.
    update_debug_color();

    // Initialize JACK.
    if (!init_jack())
        return -1;

    // Update this directory to point to your folder containing WAV files.
    std::string directory = R"(M:\port mira\audio_outputs\pcm)"; // <-- Change accordingly.
    load_samples(directory, g_samples);
    if (g_samples.empty()) {
        std::cerr << "No .wav files found in directory." << std::endl;
        return -1;
    }

    // Generate colors and cluster samples.
    generate_colors(g_samples);
    kmeans_cluster(g_samples);

    // Initialize GLFW and create an OpenGL window.
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW." << std::endl;
        return -1;
    }
    GLFWwindow* window = glfwCreateWindow(1000, 700, "Sound Galaxy Explorer", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window." << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD." << std::endl;
        return -1;
    }
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetKeyCallback(window, key_callback);

    // Set up an orthographic projection for 2D rendering.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1000, 0, 700, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    double lastFrameTime = glfwGetTime();
    // Main rendering loop.
    while (!glfwWindowShouldClose(window)) {
        double currentTime = glfwGetTime();
        float deltaTime = static_cast<float>(currentTime - lastFrameTime);
        lastFrameTime = currentTime;

        // Update animation timers for each sample.
        for (auto& sample : g_samples) {
            if (sample.animTimer > 0.0f) {
                sample.animTimer -= deltaTime;
                if (sample.animTimer < 0.0f)
                    sample.animTimer = 0.0f;
            }
        }

        // --- Gravity Simulation ---
        // For every pair of samples, compute a force based on their color similarity.
        // Similar colors (small normalized color difference) attract; dissimilar colors repel.
        for (size_t i = 0; i < g_samples.size(); i++) {
            glm::vec2 forceAccum(0.0f, 0.0f);
            for (size_t j = 0; j < g_samples.size(); j++) {
                if (i == j)
                    continue;
                glm::vec2 dir = g_samples[j].position - g_samples[i].position;
                float dist = glm::length(dir);
                if (dist < 0.1f)
                    continue; // avoid singularities

                // Normalize the average colors (they are in [0,255]) to [0,1]
                glm::vec3 colA = g_samples[i].avgColor / 255.0f;
                glm::vec3 colB = g_samples[j].avgColor / 255.0f;
                float normalizedColorDiff = glm::length(colA - colB); // in [0, sqrt(3)]
                // Map to a force factor: similar colors (diff ~ 0) yield positive force (attraction),
                // while very dissimilar colors (diff ~1.73) yield negative force (repulsion).
                float threshold = 0.15f; // adjust this threshold as needed
                float maxDiff = 1.732f; // maximum possible normalized difference (sqrt(3))
                float forceFactor;

                if (normalizedColorDiff < threshold) {
                    // Attraction remains the same.
                    forceFactor = (threshold - normalizedColorDiff) / threshold;
                }
                else {
                    // Repulsion is now 2x stronger.
                    forceFactor = -1.0f * (normalizedColorDiff - threshold) / (maxDiff - threshold);
                }




                // Calculate a simple force: proportional to forceFactor, inversely to distance.
                float constantForce = 10.0f; // tweak this constant for stronger/weaker interaction
                glm::vec2 force = (dir / dist) * constantForce * forceFactor / (dist + 1.0f);
                forceAccum += force;
            }
            // Update velocity and apply damping.
            g_samples[i].velocity += forceAccum * deltaTime;
            g_samples[i].velocity *= 0.99f; // damping to prevent oscillations
        }

        // Update positions based on velocities.
        for (auto& sample : g_samples) {
            sample.position += sample.velocity * deltaTime;
        }

        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT);

        // Render scene with pan and zoom.
        glPushMatrix();
        glTranslatef(g_offset.x * g_scale, g_offset.y * g_scale, 0);
        glScalef(g_scale, g_scale, 1.0f);
        for (auto& sample : g_samples) {
            float radius = 10.0f;
            glm::vec3 drawColor = sample.avgColor;

            if (sample.animTimer > 0.0f) {
                // Compute progress: 1 means just clicked, 0 means effect ended.
                float progress = sample.animTimer / EFFECT_DURATION;
                // At progress == 1, radius is 5.0; it grows to 10.0 as progress goes to 0.
                radius = 10.0f - 5.0f * progress;
                // Interpolate color from blue (when clicked) back to original color.
                drawColor = (glm::vec3(0.0f, 0.0f, 1.0f) * progress) + (sample.avgColor * (1.0f - progress));
            }
            draw_circle(sample.position, radius, drawColor);
        }
        glPopMatrix();

        // Render debug square if debug mode is enabled.
        if (g_debug) {
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glColor3f(g_debugColor.r, g_debugColor.g, g_debugColor.b);
            // Draw a square in the top left corner.
            glBegin(GL_QUADS);
            glVertex2f(10.0f, 640.0f);
            glVertex2f(60.0f, 640.0f);
            glVertex2f(60.0f, 690.0f);
            glVertex2f(10.0f, 690.0f);
            glEnd();
            glPopMatrix();
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    jack_client_close(jack_client);
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

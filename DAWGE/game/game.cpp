#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>
#include <numeric>

const unsigned int WINDOW_WIDTH = 800;
const unsigned int WINDOW_HEIGHT = 600;
const float RENDER_DISTANCE = 16.0f;
const int CHUNK_SIZE = 16;

// Camera variables
glm::vec3 cameraPos = glm::vec3(0.0f, 10.0f, 3.0f);  // Raised camera height
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

// Mouse variables
float yaw = -90.0f;
float pitch = 0.0f;
float lastX = WINDOW_WIDTH / 2.0f;
float lastY = WINDOW_HEIGHT / 2.0f;
bool firstMouse = true;

// Timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

class PerlinNoise {
private:
    std::vector<int> p;

    static double fade(double t) {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }

    static double lerp(double t, double a, double b) {
        return a + t * (b - a);
    }

    static double grad(int hash, double x, double y, double z) {
        int h = hash & 15;
        double u = h < 8 ? x : y;
        double v = h < 4 ? y : h == 12 || h == 14 ? x : z;
        return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
    }

public:
    PerlinNoise(int seed = 0) {
        p.resize(512);
        std::iota(p.begin(), p.begin() + 256, 0);
        std::mt19937 gen(seed);
        std::shuffle(p.begin(), p.begin() + 256, gen);
        for (int i = 0; i < 256; i++) {
            p[256 + i] = p[i];
        }
    }

    double noise(double x, double y, double z) {
        int X = static_cast<int>(std::floor(x)) & 255;
        int Y = static_cast<int>(std::floor(y)) & 255;
        int Z = static_cast<int>(std::floor(z)) & 255;

        x -= std::floor(x);
        y -= std::floor(y);
        z -= std::floor(z);

        double u = fade(x);
        double v = fade(y);
        double w = fade(z);

        int A = p[X] + Y;
        int AA = p[A] + Z;
        int AB = p[A + 1] + Z;
        int B = p[X + 1] + Y;
        int BA = p[B] + Z;
        int BB = p[B + 1] + Z;

        return lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z),
            grad(p[BA], x - 1, y, z)),
            lerp(u, grad(p[AB], x, y - 1, z),
                grad(p[BB], x - 1, y - 1, z))),
            lerp(v, lerp(u, grad(p[AA + 1], x, y, z - 1),
                grad(p[BA + 1], x - 1, y, z - 1)),
                lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                    grad(p[BB + 1], x - 1, y - 1, z - 1))));
    }

    double ridgeNoise(double x, double y, double z) {
        double n = noise(x, y, z);
        n = 1.0 - std::abs(n);
        n = n * n;
        return n;
    }
};

// Initialize noise generators with different seeds
PerlinNoise continentalNoise(1);
PerlinNoise elevationNoise(2);
PerlinNoise ridgeNoise(3);

struct TerrainPoint {
    double height;
    bool isLand;
};

TerrainPoint getTerrainHeight(double x, double z) {
    const double CONTINENTAL_SCALE = 100.0;
    const double ELEVATION_SCALE = 50.0;
    const double RIDGE_SCALE = 25.0;

    double continental = continentalNoise.noise(x / CONTINENTAL_SCALE, 0, z / CONTINENTAL_SCALE);
    continental = (continental + 1.0) / 2.0;

    bool isLand = continental > 0.48;  // Slightly lower threshold for more land

    if (!isLand) {
        return { -4.0, false };
    }

    double elevation = elevationNoise.noise(x / ELEVATION_SCALE, 0, z / ELEVATION_SCALE);
    elevation = (elevation + 1.0) / 2.0;

    double ridge = ridgeNoise.ridgeNoise(x / RIDGE_SCALE, 0, z / RIDGE_SCALE);

    double height = elevation * 8.0;  // Base height
    height += ridge * 12.0;  // Add mountains

    return { height, true };
}

struct ChunkPos {
    int x, z;
    bool operator==(const ChunkPos& other) const {
        return x == other.x && z == other.z;
    }
};

namespace std {
    template<>
    struct hash<ChunkPos> {
        size_t operator()(const ChunkPos& k) const {
            return hash<int>()(k.x) ^ (hash<int>()(k.z) << 1);
        }
    };
}

struct Chunk {
    std::vector<glm::vec3> waterPositions;
    std::vector<glm::vec3> stonePositions;
    bool needsMeshUpdate;

    Chunk() : needsMeshUpdate(true) {}
};

std::unordered_map<ChunkPos, Chunk> chunks;

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    const float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;

    pitch = std::min(pitch, 89.0f);
    pitch = std::max(pitch, -89.0f);

    glm::vec3 front;
    front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front.y = sin(glm::radians(pitch));
    front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraFront = glm::normalize(front);
} // || BREAK |>
void processInput(GLFWwindow* window) {
    const float cameraSpeed = 20.0f * deltaTime;  // Increased for better terrain viewing
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos += cameraSpeed * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos -= cameraSpeed * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        cameraPos += cameraSpeed * cameraUp;
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        cameraPos -= cameraSpeed * cameraUp;
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void generateChunkMesh(Chunk& chunk, int chunkX, int chunkZ) {
    if (!chunk.needsMeshUpdate) return;

    chunk.waterPositions.clear();
    chunk.stonePositions.clear();

    for (int x = 0; x < CHUNK_SIZE; x++) {
        for (int z = 0; z < CHUNK_SIZE; z++) {
            double worldX = chunkX * CHUNK_SIZE + x;
            double worldZ = chunkZ * CHUNK_SIZE + z;

            TerrainPoint terrain = getTerrainHeight(worldX, worldZ);

            // Water at y=0 if not land
            if (!terrain.isLand) {
                chunk.waterPositions.push_back(glm::vec3(worldX, 0.0f, worldZ));
            }

            // Generate stone columns
            int height = static_cast<int>(std::floor(terrain.height));
            for (int y = height; y >= -4; y--) {
                chunk.stonePositions.push_back(glm::vec3(worldX, y, worldZ));
            }
        }
    }

    chunk.needsMeshUpdate = false;
}

void updateChunks() {
    int playerChunkX = static_cast<int>(floor(cameraPos.x / CHUNK_SIZE));
    int playerChunkZ = static_cast<int>(floor(cameraPos.z / CHUNK_SIZE));
    int renderDistance = static_cast<int>(RENDER_DISTANCE);

    // Remove far chunks
    for (auto it = chunks.begin(); it != chunks.end();) {
        int dx = abs(it->first.x - playerChunkX);
        int dz = abs(it->first.z - playerChunkZ);
        if (dx > renderDistance || dz > renderDistance) {
            it = chunks.erase(it);
        }
        else {
            ++it;
        }
    }

    // Add or update nearby chunks
    for (int x = playerChunkX - renderDistance; x <= playerChunkX + renderDistance; x++) {
        for (int z = playerChunkZ - renderDistance; z <= playerChunkZ + renderDistance; z++) {
            ChunkPos pos{ x, z };
            if (chunks.find(pos) == chunks.end()) {
                chunks[pos] = Chunk();
            }
            generateChunkMesh(chunks[pos], x, z);
        }
    }
}

int main() {
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW\n";
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Minecraft Clone", nullptr, nullptr);
    if (!window) {
        std::cout << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD\n";
        return -1;
    } // || BREAK |>
    const char* vertexShaderSource = R"(
       #version 330 core
       layout (location = 0) in vec3 aPos;
       layout (location = 1) in vec3 aColor;
       layout (location = 2) in vec3 aOffset;
       out vec3 ourColor;
       uniform mat4 model;
       uniform mat4 view;
       uniform mat4 projection;
       uniform bool isWater;
       void main() {
           vec3 pos = aPos;
           if (gl_InstanceID > 0) {
               pos += aOffset;
           }
           gl_Position = projection * view * model * vec4(pos, 1.0);
           ourColor = isWater ? vec3(0.0, 0.0, 1.0) : aColor;
       }
   )";

    const char* fragmentShaderSource = R"(
       #version 330 core
       in vec3 ourColor;
       out vec4 FragColor;
       uniform bool isWater;
       void main() {
           float alpha = isWater ? 0.5 : 1.0;
           FragColor = vec4(ourColor, alpha);
       }
   )";

    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    float vertices[] = {
        -0.5f, -0.5f,  0.5f,  0.5f, 0.5f, 0.5f,
         0.5f, -0.5f,  0.5f,  0.5f, 0.5f, 0.5f,
         0.5f,  0.5f,  0.5f,  0.5f, 0.5f, 0.5f,
        -0.5f,  0.5f,  0.5f,  0.5f, 0.5f, 0.5f,
        -0.5f, -0.5f, -0.5f,  0.5f, 0.5f, 0.5f,
         0.5f, -0.5f, -0.5f,  0.5f, 0.5f, 0.5f,
         0.5f,  0.5f, -0.5f,  0.5f, 0.5f, 0.5f,
        -0.5f,  0.5f, -0.5f,  0.5f, 0.5f, 0.5f,
    };

    // Red vertices for origin cube
    float redVertices[] = {
        -0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
        -0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
        -0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 0.0f,
         0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 0.0f,
         0.5f,  0.5f, -0.5f,  1.0f, 0.0f, 0.0f,
        -0.5f,  0.5f, -0.5f,  1.0f, 0.0f, 0.0f,
    }; // || BREAK |>
    unsigned int indices[] = {
       0, 1, 2,  2, 3, 0,
       1, 5, 6,  6, 2, 1,
       5, 4, 7,  7, 6, 5,
       4, 0, 3,  3, 7, 4,
       3, 2, 6,  6, 7, 3,
       4, 5, 1,  1, 0, 4
    };

    unsigned int VAO, redVAO, VBO, redVBO, waterVAO, stoneVAO, EBO, instanceVBO;
    glGenVertexArrays(1, &VAO);
    glGenVertexArrays(1, &redVAO);
    glGenVertexArrays(1, &waterVAO);
    glGenVertexArrays(1, &stoneVAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &redVBO);
    glGenBuffers(1, &EBO);
    glGenBuffers(1, &instanceVBO);

    // Setup regular VAO
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Setup red cube VAO
    glBindVertexArray(redVAO);
    glBindBuffer(GL_ARRAY_BUFFER, redVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(redVertices), redVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Setup water and stone VAOs with instancing
    glBindVertexArray(waterVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    glBindVertexArray(stoneVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // || BREAK |>
    while (!glfwWindowShouldClose(window)) {
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);
        updateChunks();

        glClearColor(0.53f, 0.81f, 0.92f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(shaderProgram);

        glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
        glm::mat4 projection = glm::perspective(glm::radians(45.0f),
            static_cast<float>(WINDOW_WIDTH) / static_cast<float>(WINDOW_HEIGHT), 0.1f, 1000.0f);

        // Draw stone first
        glUniform1i(glGetUniformLocation(shaderProgram, "isWater"), 0);
        glm::mat4 model = glm::mat4(1.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));

        glBindVertexArray(stoneVAO);
        for (const auto& chunk : chunks) {
            if (!chunk.second.stonePositions.empty()) {
                glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
                glBufferData(GL_ARRAY_BUFFER, chunk.second.stonePositions.size() * sizeof(glm::vec3),
                    chunk.second.stonePositions.data(), GL_DYNAMIC_DRAW);
                glDrawElementsInstanced(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0,
                    static_cast<GLsizei>(chunk.second.stonePositions.size()));
            }
        }

        // Draw red origin cube at y=1
        glUniform1i(glGetUniformLocation(shaderProgram, "isWater"), 0);
        model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        glBindVertexArray(redVAO);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);

        // Draw water last for transparency
        glUniform1i(glGetUniformLocation(shaderProgram, "isWater"), 1);
        model = glm::mat4(1.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        glBindVertexArray(waterVAO);
        for (const auto& chunk : chunks) {
            if (!chunk.second.waterPositions.empty()) {
                glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
                glBufferData(GL_ARRAY_BUFFER, chunk.second.waterPositions.size() * sizeof(glm::vec3),
                    chunk.second.waterPositions.data(), GL_DYNAMIC_DRAW);
                glDrawElementsInstanced(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0,
                    static_cast<GLsizei>(chunk.second.waterPositions.size()));
            }
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteVertexArrays(1, &redVAO);
    glDeleteVertexArrays(1, &waterVAO);
    glDeleteVertexArrays(1, &stoneVAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &redVBO);
    glDeleteBuffers(1, &EBO);
    glDeleteBuffers(1, &instanceVBO);
    glDeleteProgram(shaderProgram);

    glfwTerminate();
    return 0;
} // ##!!

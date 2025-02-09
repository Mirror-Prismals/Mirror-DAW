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
#include <algorithm>

// ----------------------------------------------------------------------
// Configuration constants
// ----------------------------------------------------------------------
const unsigned int WINDOW_WIDTH = 800;
const unsigned int WINDOW_HEIGHT = 600;
// Increase render distance to 100 chunks (each chunk is 16 blocks wide)
const float RENDER_DISTANCE = 100.0f;
const int CHUNK_SIZE = 16;

// ----------------------------------------------------------------------
// Camera variables
// ----------------------------------------------------------------------
glm::vec3 cameraPos   = glm::vec3(0.0f, 10.0f, 3.0f);  // Raised for terrain view.
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

// Mouse state
float yaw = -90.0f;
float pitch = 0.0f;
float lastX = WINDOW_WIDTH / 2.0f;
float lastY = WINDOW_HEIGHT / 2.0f;
bool firstMouse = true;

// Timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

// ----------------------------------------------------------------------
// Perlin Noise (for terrain generation)
// ----------------------------------------------------------------------
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
        double v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
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

// Initialize three different noise generators.
PerlinNoise continentalNoise(1);
PerlinNoise elevationNoise(2);
PerlinNoise ridgeNoise(3);

// ----------------------------------------------------------------------
// Terrain generation
// ----------------------------------------------------------------------
struct TerrainPoint {
    double height;
    bool isLand;
};

TerrainPoint getTerrainHeight(double x, double z) {
    const double CONTINENTAL_SCALE = 100.0;
    const double ELEVATION_SCALE   = 50.0;
    const double RIDGE_SCALE       = 25.0;

    double continental = continentalNoise.noise(x / CONTINENTAL_SCALE, 0, z / CONTINENTAL_SCALE);
    continental = (continental + 1.0) / 2.0;
    bool isLand = continental > 0.48; // Slightly lower threshold for more land.

    if (!isLand)
        return { -4.0, false };

    double elevation = elevationNoise.noise(x / ELEVATION_SCALE, 0, z / ELEVATION_SCALE);
    elevation = (elevation + 1.0) / 2.0;
    double ridge = ridgeNoise.ridgeNoise(x / RIDGE_SCALE, 0, z / RIDGE_SCALE);
    double height = elevation * 8.0;
    height += ridge * 12.0;
    return { height, true };
}

// ----------------------------------------------------------------------
// Chunk and Chunk Position Structures
// ----------------------------------------------------------------------
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
    std::vector<glm::vec3> treeTrunkPositions;
    std::vector<glm::vec3> treeLeafPositions;
    bool needsMeshUpdate;
    Chunk() : needsMeshUpdate(true) {}
};

std::unordered_map<ChunkPos, Chunk> chunks;

// ----------------------------------------------------------------------
// Frustum Culling
// ----------------------------------------------------------------------
struct Plane {
    glm::vec3 normal;
    float d;
};

std::vector<Plane> extractFrustumPlanes(const glm::mat4& VP) {
    std::vector<Plane> planes(6);
    // Left
    planes[0].normal.x = VP[0][3] + VP[0][0];
    planes[0].normal.y = VP[1][3] + VP[1][0];
    planes[0].normal.z = VP[2][3] + VP[2][0];
    planes[0].d = VP[3][3] + VP[3][0];
    // Right
    planes[1].normal.x = VP[0][3] - VP[0][0];
    planes[1].normal.y = VP[1][3] - VP[1][0];
    planes[1].normal.z = VP[2][3] - VP[2][0];
    planes[1].d = VP[3][3] - VP[3][0];
    // Bottom
    planes[2].normal.x = VP[0][3] + VP[0][1];
    planes[2].normal.y = VP[1][3] + VP[1][1];
    planes[2].normal.z = VP[2][3] + VP[2][1];
    planes[2].d = VP[3][3] + VP[3][1];
    // Top
    planes[3].normal.x = VP[0][3] - VP[0][1];
    planes[3].normal.y = VP[1][3] - VP[1][1];
    planes[3].normal.z = VP[2][3] - VP[2][1];
    planes[3].d = VP[3][3] - VP[3][1];
    // Near
    planes[4].normal.x = VP[0][3] + VP[0][2];
    planes[4].normal.y = VP[1][3] + VP[1][2];
    planes[4].normal.z = VP[2][3] + VP[2][2];
    planes[4].d = VP[3][3] + VP[3][2];
    // Far
    planes[5].normal.x = VP[0][3] - VP[0][2];
    planes[5].normal.y = VP[1][3] - VP[1][2];
    planes[5].normal.z = VP[2][3] - VP[2][2];
    planes[5].d = VP[3][3] - VP[3][2];
    
    // Normalize all planes.
    for (int i = 0; i < 6; i++) {
        float length = glm::length(planes[i].normal);
        planes[i].normal /= length;
        planes[i].d /= length;
    }
    return planes;
}

bool aabbInFrustum(const std::vector<Plane>& planes, const glm::vec3& min, const glm::vec3& max) {
    for (int i = 0; i < 6; i++) {
        glm::vec3 p;
        p.x = (planes[i].normal.x >= 0) ? max.x : min.x;
        p.y = (planes[i].normal.y >= 0) ? max.y : min.y;
        p.z = (planes[i].normal.z >= 0) ? max.z : min.z;
        if (glm::dot(planes[i].normal, p) + planes[i].d < 0)
            return false;
    }
    return true;
}

// ----------------------------------------------------------------------
// Pine Tree Canopy Generator for Tall Trees
// ----------------------------------------------------------------------
//
// This function generates a conical pine tree canopy using horizontal rings.
// With a 60-block trunk, the canopy now extends 4x lower than before.
// For example, with trunkHeight = 60 and canopyOffset = 50, the canopy starts at
// groundHeight + 60 - 50 = groundHeight + 10 and spans canopyLayers blocks upward.
std::vector<glm::vec3> generatePineCanopy(int groundHeight, int trunkHeight, int trunkThickness, double worldX, double worldZ) {
    std::vector<glm::vec3> leafPositions;
    
    // Adjusted parameters (you mentioned using 50 works perfectly):
    int canopyOffset = 50;      // The canopy now starts 50 blocks below the trunk top.
    int canopyLayers = 80;      // Number of canopy layers; you can tweak this as desired.
    int canopyBase = groundHeight + trunkHeight - canopyOffset;
    
    // Define ring radii.
    float bottomRadius = 8.0f;  // Largest ring at the bottom.
    float topRadius = 2.0f;     // Smallest ring at the top.
    float ringThickness = 1.0f; // Thickness of each ring.
    
    // For a trunk of size trunkThickness x trunkThickness, assume the trunk's center:
    float centerOffset = (trunkThickness - 1) / 2.0f;
    
    // Build the canopy in layers.
    for (int layer = 0; layer < canopyLayers; layer++) {
        // Linearly interpolate the ring radius.
        float currentRadius = bottomRadius - layer * ((bottomRadius - topRadius) / (canopyLayers - 1));
        int yPos = canopyBase + layer;
        int range = static_cast<int>(std::ceil(currentRadius));
        for (int dx = -range; dx <= range; dx++) {
            for (int dz = -range; dz <= range; dz++) {
                float dist = std::sqrt(dx * dx + dz * dz);
                if (std::abs(dist - currentRadius) < ringThickness) {
                    leafPositions.push_back(glm::vec3(
                        worldX + centerOffset + dx,
                        yPos,
                        worldZ + centerOffset + dz
                    ));
                }
            }
        }
    }
    return leafPositions;
}

// ----------------------------------------------------------------------
// Chunk Mesh Generation: terrain, stone columns, water, and trees.
// ----------------------------------------------------------------------
void generateChunkMesh(Chunk& chunk, int chunkX, int chunkZ) {
    if (!chunk.needsMeshUpdate)
        return;
    
    chunk.waterPositions.clear();
    chunk.stonePositions.clear();
    chunk.treeTrunkPositions.clear();
    chunk.treeLeafPositions.clear();
    
    for (int x = 0; x < CHUNK_SIZE; x++) {
        for (int z = 0; z < CHUNK_SIZE; z++) {
            double worldX = chunkX * CHUNK_SIZE + x;
            double worldZ = chunkZ * CHUNK_SIZE + z;
            TerrainPoint terrain = getTerrainHeight(worldX, worldZ);
            
            if (!terrain.isLand) {
                chunk.waterPositions.push_back(glm::vec3(worldX, 0.0f, worldZ));
            } else {
                int groundHeight = static_cast<int>(std::floor(terrain.height));
                for (int y = groundHeight; y >= -4; y--) {
                    chunk.stonePositions.push_back(glm::vec3(worldX, y, worldZ));
                }
                if (terrain.height > 2.0) {
                    int intWorldX = static_cast<int>(worldX);
                    int intWorldZ = static_cast<int>(worldZ);
                    int hashVal = std::abs((intWorldX * 73856093) ^ (intWorldZ * 19349663));
                    if (hashVal % 10000 < 1) {  // 1 in 10,000 chance
                        int trunkHeight = 60;       // 60-block tall trunk.
                        int trunkThickness = 4;     // 4x4 trunk.
                        for (int i = 1; i <= trunkHeight; i++) {
                            for (int tx = 0; tx < trunkThickness; tx++) {
                                for (int tz = 0; tz < trunkThickness; tz++) {
                                    chunk.treeTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                }
                            }
                        }
                        std::vector<glm::vec3> canopy = generatePineCanopy(groundHeight, trunkHeight, trunkThickness, worldX, worldZ);
                        chunk.treeLeafPositions.insert(chunk.treeLeafPositions.end(), canopy.begin(), canopy.end());
                    }
                }
            }
        }
    }
    
    chunk.needsMeshUpdate = false;
}

// ----------------------------------------------------------------------
// Update chunks based on player position.
// ----------------------------------------------------------------------
void updateChunks() {
    int playerChunkX = static_cast<int>(floor(cameraPos.x / CHUNK_SIZE));
    int playerChunkZ = static_cast<int>(floor(cameraPos.z / CHUNK_SIZE));
    int renderDistance = static_cast<int>(RENDER_DISTANCE);
    
    for (auto it = chunks.begin(); it != chunks.end();) {
        int dx = abs(it->first.x - playerChunkX);
        int dz = abs(it->first.z - playerChunkZ);
        if (dx > renderDistance || dz > renderDistance)
            it = chunks.erase(it);
        else
            ++it;
    }
    
    for (int x = playerChunkX - renderDistance; x <= playerChunkX + renderDistance; x++) {
        for (int z = playerChunkZ - renderDistance; z <= playerChunkZ + renderDistance; z++) {
            ChunkPos pos { x, z };
            if (chunks.find(pos) == chunks.end())
                chunks[pos] = Chunk();
            generateChunkMesh(chunks[pos], x, z);
        }
    }
}

// ----------------------------------------------------------------------
// Callback functions for window and input
// ----------------------------------------------------------------------
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
    yaw   += xoffset;
    pitch += yoffset;
    pitch = std::min(pitch, 89.0f);
    pitch = std::max(pitch, -89.0f);
    glm::vec3 front;
    front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front.y = sin(glm::radians(pitch));
    front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraFront = glm::normalize(front);
}

void processInput(GLFWwindow* window) {
    const float cameraSpeed = 20.0f * deltaTime;
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

// ----------------------------------------------------------------------
// Main function
// ----------------------------------------------------------------------
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
    }
    
    // ------------------------------------------------------------------
    // Shader setup
    // ------------------------------------------------------------------
    // Shader uses a "blockType" uniform:
    //   0 = stone/red (default)
    //   1 = water (blue, semi-transparent)
    //   2 = tree trunk (brown)
    //   3 = tree leaves (green)
    const char* vertexShaderSource = R"(
       #version 330 core
       layout (location = 0) in vec3 aPos;
       layout (location = 1) in vec3 aColor;
       layout (location = 2) in vec3 aOffset;
       out vec3 ourColor;
       uniform mat4 model;
       uniform mat4 view;
       uniform mat4 projection;
       uniform int blockType;
       void main() {
           vec3 pos = aPos;
           if (gl_InstanceID > 0)
               pos += aOffset;
           gl_Position = projection * view * model * vec4(pos, 1.0);
           if (blockType == 1)
               ourColor = vec3(0.0, 0.0, 1.0);
           else if (blockType == 2)
               ourColor = vec3(0.55, 0.27, 0.07);
           else if (blockType == 3)
               ourColor = vec3(0.0, 0.8, 0.0);
           else
               ourColor = aColor;
       }
    )";
    
    const char* fragmentShaderSource = R"(
       #version 330 core
       in vec3 ourColor;
       out vec4 FragColor;
       uniform int blockType;
       void main() {
           float alpha = (blockType == 1) ? 0.5 : 1.0;
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
    
    // ------------------------------------------------------------------
    // Set up vertex data for a cube.
    // ------------------------------------------------------------------
    float vertices[] = {
        // positions          // colors (unused for most blocks)
        -0.5f, -0.5f,  0.5f,   0.5f, 0.5f, 0.5f,
         0.5f, -0.5f,  0.5f,   0.5f, 0.5f, 0.5f,
         0.5f,  0.5f,  0.5f,   0.5f, 0.5f, 0.5f,
        -0.5f,  0.5f,  0.5f,   0.5f, 0.5f, 0.5f,
        -0.5f, -0.5f, -0.5f,   0.5f, 0.5f, 0.5f,
         0.5f, -0.5f, -0.5f,   0.5f, 0.5f, 0.5f,
         0.5f,  0.5f, -0.5f,   0.5f, 0.5f, 0.5f,
        -0.5f,  0.5f, -0.5f,   0.5f, 0.5f, 0.5f,
    };
    
    float redVertices[] = {
        -0.5f, -0.5f,  0.5f,   1.0f, 0.0f, 0.0f,
         0.5f, -0.5f,  0.5f,   1.0f, 0.0f, 0.0f,
         0.5f,  0.5f,  0.5f,   1.0f, 0.0f, 0.0f,
        -0.5f,  0.5f,  0.5f,   1.0f, 0.0f, 0.0f,
        -0.5f, -0.5f, -0.5f,   1.0f, 0.0f, 0.0f,
         0.5f, -0.5f, -0.5f,   1.0f, 0.0f, 0.0f,
         0.5f,  0.5f, -0.5f,   1.0f, 0.0f, 0.0f,
        -0.5f,  0.5f, -0.5f,   1.0f, 0.0f, 0.0f,
    };
    
    unsigned int indices[] = {
       0, 1, 2,  2, 3, 0,
       1, 5, 6,  6, 2, 1,
       5, 4, 7,  7, 6, 5,
       4, 0, 3,  3, 7, 4,
       3, 2, 6,  6, 7, 3,
       4, 5, 1,  1, 0, 4
    };
    
    // ------------------------------------------------------------------
    // Generate and set up VAOs/VBOs for cube and instancing.
    // ------------------------------------------------------------------
    unsigned int VAO, redVAO, VBO, redVBO, waterVAO, stoneVAO, EBO, instanceVBO;
    unsigned int treeTrunkVAO, treeLeafVAO;
    glGenVertexArrays(1, &VAO);
    glGenVertexArrays(1, &redVAO);
    glGenVertexArrays(1, &waterVAO);
    glGenVertexArrays(1, &stoneVAO);
    glGenVertexArrays(1, &treeTrunkVAO);
    glGenVertexArrays(1, &treeLeafVAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &redVBO);
    glGenBuffers(1, &EBO);
    glGenBuffers(1, &instanceVBO);
    
    // Regular cube VAO.
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    
    // Red cube VAO.
    glBindVertexArray(redVAO);
    glBindBuffer(GL_ARRAY_BUFFER, redVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(redVertices), redVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    
    // Water VAO with instancing.
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
    
    // Stone VAO with instancing.
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
    
    // Tree trunk VAO with instancing.
    glBindVertexArray(treeTrunkVAO);
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
    
    // Tree leaves VAO with instancing.
    glBindVertexArray(treeLeafVAO);
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
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // ------------------------------------------------------------------
    // Main render loop
    // ------------------------------------------------------------------
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
        // Increase far clip to 10,000.0 to see far.
        glm::mat4 projection = glm::perspective(glm::radians(45.0f),
                                  static_cast<float>(WINDOW_WIDTH) / static_cast<float>(WINDOW_HEIGHT),
                                  0.1f, 10000.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        
        // ------------------------------------------------------------------
        // Frustum culling and global batching of instance positions.
        // ------------------------------------------------------------------
        glm::mat4 VP = projection * view;
        std::vector<Plane> frustumPlanes = extractFrustumPlanes(VP);
        
        // Global instance arrays for each block type.
        std::vector<glm::vec3> globalStoneInstances;
        std::vector<glm::vec3> globalWaterInstances;
        std::vector<glm::vec3> globalTreeTrunkInstances;
        std::vector<glm::vec3> globalTreeLeafInstances;
        
        // For each chunk, compute its AABB and if visible, append its instance positions.
        for (const auto& entry : chunks) {
            const ChunkPos& pos = entry.first;
            const Chunk& chunk = entry.second;
            // Define the AABB of the chunk:
            glm::vec3 aabbMin(pos.x * CHUNK_SIZE, -4.0f, pos.z * CHUNK_SIZE);
            glm::vec3 aabbMax((pos.x + 1) * CHUNK_SIZE, 150.0f, (pos.z + 1) * CHUNK_SIZE);
            if (!aabbInFrustum(frustumPlanes, aabbMin, aabbMax))
                continue; // Skip chunks outside the view.
            
            globalStoneInstances.insert(globalStoneInstances.end(), chunk.stonePositions.begin(), chunk.stonePositions.end());
            globalWaterInstances.insert(globalWaterInstances.end(), chunk.waterPositions.begin(), chunk.waterPositions.end());
            globalTreeTrunkInstances.insert(globalTreeTrunkInstances.end(), chunk.treeTrunkPositions.begin(), chunk.treeTrunkPositions.end());
            globalTreeLeafInstances.insert(globalTreeLeafInstances.end(), chunk.treeLeafPositions.begin(), chunk.treeLeafPositions.end());
        }
        
        // ------------------------------------------------------------------
        // Render each block type with one instanced draw call per type.
        // ------------------------------------------------------------------
        // Helper lambda to update instance VBO and draw.
        auto drawInstances = [&](unsigned int vao, int blockType, const std::vector<glm::vec3>& instances) {
            if (instances.empty())
                return;
            glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), blockType);
            glm::mat4 model = glm::mat4(1.0f);
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(vao);
            glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(glm::vec3), instances.data(), GL_DYNAMIC_DRAW);
            glDrawElementsInstanced(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0, static_cast<GLsizei>(instances.size()));
        };
        
        // Draw stone blocks (blockType 0).
        drawInstances(stoneVAO, 0, globalStoneInstances);
        // Draw red origin cube (blockType 0).
        glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), 0);
        {
            glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(redVAO);
            glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        }
        // Draw tree trunks (blockType 2).
        drawInstances(treeTrunkVAO, 2, globalTreeTrunkInstances);
        // Draw tree leaves (blockType 3).
        drawInstances(treeLeafVAO, 3, globalTreeLeafInstances);
        // Draw water blocks last for transparency (blockType 1).
        drawInstances(waterVAO, 1, globalWaterInstances);
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    // ------------------------------------------------------------------
    // Cleanup
    // ------------------------------------------------------------------
    glDeleteVertexArrays(1, &VAO);
    glDeleteVertexArrays(1, &redVAO);
    glDeleteVertexArrays(1, &waterVAO);
    glDeleteVertexArrays(1, &stoneVAO);
    glDeleteVertexArrays(1, &treeTrunkVAO);
    glDeleteVertexArrays(1, &treeLeafVAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &redVBO);
    glDeleteBuffers(1, &EBO);
    glDeleteBuffers(1, &instanceVBO);
    glDeleteProgram(shaderProgram);
    
    glfwTerminate();
    return 0;
}

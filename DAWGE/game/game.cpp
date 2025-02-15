// ======================================================================
// VoxelGame.cpp
// A single–file version of a Minecraft–style clone with multiple biomes,
// dynamic skybox, sun/moon and now water/swimming/prone/paraglide mechanics.
// ======================================================================

#define GLM_ENABLE_EXPERIMENTAL

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <ctime>

// ---------------------- ChunkPos Definition and Hash Specialization ----------------------
struct ChunkPos {
    int x, z;
    ChunkPos() : x(0), z(0) {}
    ChunkPos(int _x, int _z) : x(_x), z(_z) {}
    bool operator==(const ChunkPos& other) const { return x == other.x && z == other.z; }
};

namespace std {
    template <>
    struct hash<ChunkPos> {
        std::size_t operator()(const ChunkPos& cp) const {
            return std::hash<int>()(cp.x) ^ (std::hash<int>()(cp.z) << 1);
        }
    };
}

// ---------------------- Global Constants ----------------------
const unsigned int WINDOW_WIDTH = 1920;
const unsigned int WINDOW_HEIGHT = 1080;
const float RENDER_DISTANCE = 10.0f;
const int CHUNK_SIZE = 16;
const int MIN_Y = -1;
const float WATER_SURFACE = 0.0f; // water is drawn at y=0

// ---------------------- Global Variables ----------------------

// Player/camera and movement variables
glm::vec3 waterLilyLaunchMomentum = glm::vec3(0.0f);
// Global variable for the sprint tint effect (extra blue tint)
// Track sprint timer for the incremental speed boost.
static float sprintTimer = 0.0f;
// For clarity, also track whether sprint is currently active.
static bool sprintActive = false;


glm::vec3 cameraPos = glm::vec3(0.0f, 10.0f, 3.0f);
float cameraYaw = -90.0f;
float pitch = 0.0f;
glm::vec3 velocity = glm::vec3(0.0f);
float jumpVelocity = 0.0f;
// Instead of walkingMode, we now have playerMode:
// 0 = Standing (normal walking), 1 = Prone (crawling on land),
// 2 = Swimming (in water), 3 = Paragliding (in air)
int playerMode = 0;
const float boostImpulse = 40.0f;
const float baseAcceleration = 20.0f;
const float dragFactor = 0.995f;
const float gravityForce = 9.81f * 0.1f;
const float walkSpeed = 10.0f;
const float walkGravity = 9.81f;
const float walkJumpImpulse = 4.8f;
float eyeLevelOffset = 1.6f;  // used when in standing/paragliding mode
float lastX = WINDOW_WIDTH / 2.0f;
float lastY = WINDOW_HEIGHT / 2.0f;
bool firstMouse = true;
float deltaTime = 0.0f;
float lastFrame = 0.0f;
bool fullscreenMap = false;
std::unordered_set<ChunkPos> visitedChunks;
bool bigMapDirty = true;
std::vector<float> bigMapInterleaved;
float bigMapPanX = 0.0f;
float bigMapPanZ = 0.0f;
bool minimapEnabled = true;

// ---------------------- Utility Functions for Player Bounding Box ----------------------
glm::vec3 getPlayerBoxMin() {
    // For all modes, horizontal extents remain the same.
    return glm::vec3(-0.3f, 0.0f, -0.3f);
}
glm::vec3 getPlayerBoxMax() {
    // Standing and paraglide: 2 blocks tall; prone and swimming: 1 block tall.
    if (playerMode == 1 || playerMode == 2)
        return glm::vec3(0.3f, 1.0f, 0.3f);
    else
        return glm::vec3(0.3f, 2.0f, 0.3f);
}

// ---------------------- Skybox Data and Shaders ----------------------
float skyboxQuadVertices[] = {
    // positions for full-screen quad (in clip space)
    -1.0f,  1.0f,
    -1.0f, -1.0f,
     1.0f, -1.0f,
    -1.0f,  1.0f,
     1.0f, -1.0f,
     1.0f,  1.0f
};

const char* skyboxVertexShaderSource = R"(
#version 330 core
layout(location = 0) in vec2 aPos;
out vec2 TexCoord;
void main(){
    TexCoord = aPos * 0.5 + 0.5;
    gl_Position = vec4(aPos, 0.0, 1.0);
}
)";

const char* skyboxFragmentShaderSource = R"(
#version 330 core
in vec2 TexCoord;
out vec4 FragColor;

uniform vec3 skyTop;
uniform vec3 skyBottom;
uniform float time;        // Global time
uniform vec2 starOffset;   // Derived from camera rotation

// Simple random function based on input coordinates.
float rand(vec2 co) {
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

void main(){
    // Base sky gradient.
    vec3 color = mix(skyBottom, skyTop, TexCoord.y);

    // Use the starOffset (derived from camera yaw/pitch) to adjust the star field.
    vec2 starUV = (TexCoord + starOffset) * 1000.0;
    vec2 cellId = floor(starUV);
    
    // Only a few cells become stars.
    float starValue = rand(cellId);
    
    // Lower threshold to increase star density.
    if(starValue > 0.995) {
        // Two independent random values: one for frequency and one for phase.
        float freqRand  = rand(cellId + vec2(1.0, 1.0));
        float phaseRand = rand(cellId + vec2(2.0, 3.0));
        
        // Each star gets a random frequency between 3.0 and 5.0; slow twinkle with time * 0.2.
        float frequency = 3.0 + freqRand * 2.0;
        float twinkle = 0.5 + 0.5 * sin(time * frequency * 0.2 + phaseRand * 6.2831);

        // Create a round star shape.
        vec2 f = fract(starUV) - 0.5;
        float dist = length(f);
        float starAlpha = smoothstep(0.3, 0.2, dist); // Adjust for star size/softness

        // Add the star (additively) to the sky color.
        color += vec3(twinkle) * starAlpha;
    }
    
    FragColor = vec4(color, 1.0);
}
)";

GLuint skyboxVAO, skyboxVBO;
void setupSkyboxQuad() {
    glGenVertexArrays(1, &skyboxVAO);
    glGenBuffers(1, &skyboxVBO);
    glBindVertexArray(skyboxVAO);
    glBindBuffer(GL_ARRAY_BUFFER, skyboxVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(skyboxQuadVertices), skyboxQuadVertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

// ---------------------- Sun/Moon Data and Shader Sources ----------------------
// We will render the sun and moon as billboards using a simple quad.
const char* sunMoonVertexShaderSource = R"(
#version 330 core
layout(location = 0) in vec2 aPos;
out vec2 TexCoord;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main(){
    TexCoord = aPos * 0.5 + 0.5;
    gl_Position = projection * view * model * vec4(aPos, 0.0, 1.0);
}
)";

const char* sunMoonFragmentShaderSource = R"(
#version 330 core
in vec2 TexCoord;
out vec4 FragColor;
uniform vec3 color;       // Sun or moon base color
uniform float brightness; // Intensity factor
void main(){
    // Calculate distance from center of the quad.
    float d = distance(TexCoord, vec2(0.5));
    // Create a solid disk with a soft edge.
    float diskAlpha = smoothstep(0.5, 0.45, d);
    // Add an outer glow effect.
    float glowAlpha = 1.0 - smoothstep(0.45, 0.5, d);

    // Combine disk and glow.
    float finalAlpha = clamp(diskAlpha + 0.3 * glowAlpha, 0.0, 1.0);
    FragColor = vec4(color * brightness, finalAlpha * brightness);
}
)";

GLuint sunMoonShaderProgram;
GLuint sunMoonVAO, sunMoonVBO;

// Setup a quad for sun/moon (in local object space, centered at origin)
void setupSunMoonQuad() {
    float quadVertices[] = {
        -1.0f,  1.0f,
        -1.0f, -1.0f,
         1.0f, -1.0f,
        -1.0f,  1.0f,
         1.0f, -1.0f,
         1.0f,  1.0f
    };
    glGenVertexArrays(1, &sunMoonVAO);
    glGenBuffers(1, &sunMoonVBO);
    glBindVertexArray(sunMoonVAO);
    glBindBuffer(GL_ARRAY_BUFFER, sunMoonVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

// ---------------------- Sky Color Interpolation ----------------------
struct SkyColorKey {
    float time; // Fraction of day (0.0–1.0)
    glm::vec3 top;
    glm::vec3 bottom;
};

// Key colors (converted from hex to [0,1] floats)
SkyColorKey skyKeys[5] = {
    { 0.0f,  glm::vec3(16 / 255.0f, 16 / 255.0f, 48 / 255.0f),  glm::vec3(0.0f, 0.0f, 0.0f) },     // Night
    { 0.25f, glm::vec3(0.0f, 0.0f, 1.0f),                 glm::vec3(128 / 255.0f, 128 / 255.0f, 1.0f) }, // Morning
    { 0.5f,  glm::vec3(135 / 255.0f, 206 / 255.0f, 235 / 255.0f), glm::vec3(254 / 255.0f, 254 / 255.0f, 254 / 255.0f) }, // Afternoon
    { 0.75f, glm::vec3(0.0f, 128 / 255.0f, 128 / 255.0f),       glm::vec3(1.0f, 71 / 255.0f, 0.0f) },          // Evening
    { 1.0f,  glm::vec3(16 / 255.0f, 16 / 255.0f, 48 / 255.0f),    glm::vec3(0.0f, 0.0f, 0.0f) }      // Night
};

void getCurrentSkyColors(float dayFraction, glm::vec3& currentTop, glm::vec3& currentBottom) {
    int lowerIndex = 0, upperIndex = 0;
    for (int i = 0; i < 4; i++) {
        if (dayFraction >= skyKeys[i].time && dayFraction <= skyKeys[i + 1].time) {
            lowerIndex = i;
            upperIndex = i + 1;
            break;
        }
    }
    float t = (dayFraction - skyKeys[lowerIndex].time) / (skyKeys[upperIndex].time - skyKeys[lowerIndex].time);
    currentTop = glm::mix(skyKeys[lowerIndex].top, skyKeys[upperIndex].top, t);
    currentBottom = glm::mix(skyKeys[lowerIndex].bottom, skyKeys[upperIndex].bottom, t);
}

// ---------------------- Forward Declarations ----------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn);
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void processInput(GLFWwindow* window);
glm::ivec3 raycastForBlock(bool place);

// ---------------------- Perlin Noise Class ----------------------
class PerlinNoise {
private:
    std::vector<int> p;
    static double fade(double t) { return t * t * t * (t * (t * 6 - 15) + 10); }
    static double lerp(double t, double a, double b) { return a + t * (b - a); }
    static double grad(int hash, double x, double y, double z) {
        int h = hash & 15;
        double u = (h < 8) ? x : y;
        double v = (h < 4) ? y : ((h == 12 || h == 14) ? x : z);
        return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
    }
public:
    PerlinNoise(int seed = 0) {
        p.resize(512);
        std::iota(p.begin(), p.begin() + 256, 0);
        std::mt19937 gen(seed);
        std::shuffle(p.begin(), p.begin() + 256, gen);
        for (int i = 0; i < 256; i++)
            p[256 + i] = p[i];
    }
    double noise(double x, double y, double z) {
        int X = static_cast<int>(std::floor(x)) & 255;
        int Y = static_cast<int>(std::floor(y)) & 255;
        int Z = static_cast<int>(std::floor(z)) & 255;
        x -= std::floor(x); y -= std::floor(y); z -= std::floor(z);
        double u = fade(x), v = fade(y), w = fade(z);
        int A = p[X] + Y, AA = p[A] + Z, AB = p[A + 1] + Z;
        int B = p[X + 1] + Y, BA = p[B] + Z, BB = p[B + 1] + Z;
        return lerp(w,
            lerp(v,
                lerp(u, grad(p[AA], x, y, z),
                    grad(p[BA], x - 1, y, z)),
                lerp(u, grad(p[AB], x, y - 1, z),
                    grad(p[BB], x - 1, y - 1, z))),
            lerp(v,
                lerp(u, grad(p[AA + 1], x, y, z - 1),
                    grad(p[BA + 1], x - 1, y, z - 1)),
                lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                    grad(p[BB + 1], x - 1, y - 1, z - 1))));
    }
};

PerlinNoise continentalNoise(1);
PerlinNoise elevationNoise(2);
PerlinNoise ridgeNoise(3);
PerlinNoise caveNoise(5);
PerlinNoise auroraNoise(4);
PerlinNoise lavaCaveNoise(6);
GLuint compileShaderProgram(const char* vShaderSrc, const char* fShaderSrc) {
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vShaderSrc, NULL);
    glCompileShader(vertexShader);
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "Vertex Shader Compilation Error:\n" << infoLog << "\n";
    }
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fShaderSrc, NULL);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "Fragment Shader Compilation Error:\n" << infoLog << "\n";
    }
    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(program, 512, NULL, infoLog);
        std::cout << "Shader Program Linking Error:\n" << infoLog << "\n";
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return program;
}

// ---------------------- Terrain Generation ----------------------
struct TerrainPoint { double height; bool isLand; };

TerrainPoint getTerrainHeight(double x, double z) {
    const double CONTINENTAL_SCALE = 100.0;
    const double ELEVATION_SCALE = 50.0;
    const double RIDGE_SCALE = 25.0;
    double continental = continentalNoise.noise(x / CONTINENTAL_SCALE, 0, z / CONTINENTAL_SCALE);
    continental = (continental + 1.0) / 2.0;
    bool isLand = continental > 0.48;
    if (!isLand)
        return { -4.0, false };
    double elevation = elevationNoise.noise(x / ELEVATION_SCALE, 0, z / ELEVATION_SCALE);
    elevation = (elevation + 1.0) / 2.0;
    double ridge = ridgeNoise.noise(x / RIDGE_SCALE, 0, z / RIDGE_SCALE);
    double height = elevation * 8.0 + ridge * 12.0;
    int chunkX = static_cast<int>(std::floor(x / (double)CHUNK_SIZE));
    int chunkZ = static_cast<int>(std::floor(z / (double)CHUNK_SIZE));
    if (chunkX < -20 && chunkX >= -40)
        height = elevation * 128.0 + ridge * 96.0;
    if (chunkZ >= 290 && chunkZ < 1024)
        return { -4.0, false };
    if (chunkZ <= -200 && chunkZ > -256)
        return { -4.0, false };
    return { height, true };
}

int getChunkTopBlock(int cx, int cz) {
    double samples[5][2];
    samples[0][0] = cx * CHUNK_SIZE + CHUNK_SIZE / 2.0;
    samples[0][1] = cz * CHUNK_SIZE + CHUNK_SIZE / 2.0;
    samples[1][0] = cx * CHUNK_SIZE;
    samples[1][1] = cz * CHUNK_SIZE;
    samples[2][0] = cx * CHUNK_SIZE + CHUNK_SIZE;
    samples[2][1] = cz * CHUNK_SIZE;
    samples[3][0] = cx * CHUNK_SIZE;
    samples[3][1] = cz * CHUNK_SIZE + CHUNK_SIZE;
    samples[4][0] = cx * CHUNK_SIZE + CHUNK_SIZE;
    samples[4][1] = cz * CHUNK_SIZE + CHUNK_SIZE;
    int waterSamples = 0;
    for (int i = 0; i < 5; i++) {
        TerrainPoint tp = getTerrainHeight(samples[i][0], samples[i][1]);
        if (!tp.isLand)
            waterSamples++;
    }
    if (waterSamples > 0)
        return 1;
    if (cx >= 200)
        return 22;
    if (cz <= -256)
        return 23;
    return 0;
}

// ---------------------- Chunk Structure ----------------------
struct Chunk {
    std::vector<glm::vec3> waterPositions;
    std::vector<glm::vec3> grassPositions;
    std::vector<glm::vec3> sandPositions;
    std::vector<glm::vec3> snowPositions;
    std::vector<glm::vec3> dirtPositions;
    std::vector<glm::vec3> deepStonePositions;
    std::vector<glm::vec3> lavaPositions;
    std::vector<glm::vec3> treeTrunkPositions;
    std::vector<glm::vec3> treeLeafPositions;
    std::vector<glm::vec3> firLeafPositions;
    std::vector<glm::vec3> waterLilyPositions;
    std::vector<glm::vec3> fallenTreeTrunkPositions;
    std::vector<glm::vec3> oakTrunkPositions;
    std::vector<glm::vec3> oakLeafPositions;
    std::vector<glm::vec3> leafPilePositions;
    std::vector<glm::vec3> bushSmallPositions;
    std::vector<glm::vec3> bushMediumPositions;
    std::vector<glm::vec3> bushLargePositions;
    std::vector<glm::vec4> branchPositions;
    std::vector<glm::vec3> ancientTrunkPositions;
    std::vector<glm::vec3> ancientLeafPositions;
    std::vector<glm::vec3> ancientBranchPositions;
    std::vector<glm::vec3> auroraPositions;
    std::vector<glm::vec3> icePositions;
    bool needsMeshUpdate;
    Chunk() : needsMeshUpdate(true) {}
};

std::unordered_map<ChunkPos, Chunk> chunks;

// ---------------------- Quadtree Structures ----------------------
struct Plane { glm::vec3 normal; float d; };

std::vector<Plane> extractFrustumPlanes(const glm::mat4& VP) {
    std::vector<Plane> planes(6);
    planes[0].normal.x = VP[0][3] + VP[0][0];
    planes[0].normal.y = VP[1][3] + VP[1][0];
    planes[0].normal.z = VP[2][3] + VP[2][0];
    planes[0].d = VP[3][3] + VP[3][0];
    planes[1].normal.x = VP[0][3] - VP[0][0];
    planes[1].normal.y = VP[1][3] - VP[1][0];
    planes[1].normal.z = VP[2][3] - VP[2][0];
    planes[1].d = VP[3][3] - VP[3][0];
    planes[2].normal.x = VP[0][3] + VP[0][1];
    planes[2].normal.y = VP[1][3] + VP[1][1];
    planes[2].normal.z = VP[2][3] + VP[2][1];
    planes[2].d = VP[3][3] + VP[3][1];
    planes[3].normal.x = VP[0][3] - VP[0][1];
    planes[3].normal.y = VP[1][3] - VP[1][1];
    planes[3].normal.z = VP[2][3] - VP[2][1];
    planes[3].d = VP[3][3] - VP[3][1];
    planes[4].normal.x = VP[0][3] + VP[0][2];
    planes[4].normal.y = VP[1][3] + VP[1][2];
    planes[4].normal.z = VP[2][3] + VP[2][2];
    planes[4].d = VP[3][3] + VP[3][2];
    planes[5].normal.x = VP[0][3] - VP[0][2];
    planes[5].normal.y = VP[1][3] - VP[1][2];
    planes[5].normal.z = VP[2][3] - VP[2][2];
    planes[5].d = VP[3][3] - VP[3][2];
    for (int i = 0; i < 6; i++) {
        float len = glm::length(planes[i].normal);
        planes[i].normal /= len;
        planes[i].d /= len;
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

// ---------------------- Helper: Tree Collision Check ----------------------
bool treeCollision(const std::vector<glm::vec3>& trunkArray, const glm::vec3& base) {
    for (const auto& p : trunkArray)
        if (glm::distance(p, base) < 3.0f)
            return true;
    return false;
}

// ---------------------- Tree Generation Functions ----------------------
// ---------------------- Pine Tree Generation Functions ----------------------

// --- Pine Canopy Generation Function ---
// Uses an "effective" trunk height (base trunk + extra trunk logs) to position the canopy.
std::vector<glm::vec3> generatePineCanopy(int groundHeight, int effectiveTrunkHeight, int trunkThickness, double worldX, double worldZ) {
    std::vector<glm::vec3> leafPositions;
    int canopyOffset = 70;       // How far below the effective trunk top the canopy starts
    int canopyLayers = 80;       // Number of canopy layers
    int canopyBase = groundHeight + effectiveTrunkHeight - canopyOffset;
    float bottomRadius = 8.0f;     // Radius at the bottom of the canopy
    float topRadius = 2.0f;     // Radius at the top of the canopy
    float centerOffset = (trunkThickness - 1) / 2.0f;  // Centers the canopy relative to the trunk

    for (int layer = 0; layer < canopyLayers; layer++) {
        // Linearly interpolate the radius from bottomRadius to topRadius.
        float currentRadius = bottomRadius - layer * ((bottomRadius - topRadius) / (canopyLayers - 1));
        int yPos = canopyBase + layer;
        int range = static_cast<int>(std::ceil(currentRadius));
        for (int dx = -range; dx <= range; dx++) {
            for (int dz = -range; dz <= range; dz++) {
                float dist = std::sqrt(dx * dx + dz * dz);
                // Fill the entire circle for this layer.
                if (dist <= currentRadius) {
                    leafPositions.push_back(glm::vec3(worldX + centerOffset + dx, yPos, worldZ + centerOffset + dz));
                }
            }
        }
    }
    return leafPositions;
}

// --- Pine Tree Generation Function (Trunk + Extra Logs + Canopy) ---
// "trunkHeight" is the base trunk height, and "extraTrunkLogs" is how many extra trunk blocks to add on top.
void generatePineTree(Chunk& chunk, int groundHeight, int trunkHeight, int trunkThickness, double worldX, double worldZ, int extraTrunkLogs) {
    // --- Generate Main Trunk Logs (base trunk) ---
    for (int i = 1; i <= trunkHeight; i++) {
        for (int tx = 0; tx < trunkThickness; tx++) {
            for (int tz = 0; tz < trunkThickness; tz++) {
                glm::vec3 pos(worldX + tx, groundHeight + i, worldZ + tz);
                chunk.treeTrunkPositions.push_back(pos);
            }
        }
    }

    // --- Generate Extra Trunk Logs (to extend the trunk upward) ---
    for (int i = trunkHeight + 1; i <= trunkHeight + extraTrunkLogs; i++) {
        for (int tx = 0; tx < trunkThickness; tx++) {
            for (int tz = 0; tz < trunkThickness; tz++) {
                glm::vec3 pos(worldX + tx, groundHeight + i, worldZ + tz);
                chunk.treeTrunkPositions.push_back(pos);
            }
        }
    }

    // --- Generate Canopy Using the Effective Trunk Height ---
    int effectiveTrunkHeight = trunkHeight + extraTrunkLogs;
    std::vector<glm::vec3> canopy = generatePineCanopy(groundHeight, effectiveTrunkHeight, trunkThickness, worldX, worldZ);
    chunk.treeLeafPositions.insert(chunk.treeLeafPositions.end(), canopy.begin(), canopy.end());
}




std::vector<glm::vec3> generateFirCanopy(int groundHeight, int trunkHeight, int trunkThickness, double worldX, double worldZ) {
    std::vector<glm::vec3> leafPositions;
    int centerY = groundHeight + trunkHeight;
    float radius = 7.0f;
    for (int dy = -static_cast<int>(radius); dy <= static_cast<int>(radius); dy++) {
        for (int dx = -static_cast<int>(radius); dx <= static_cast<int>(radius); dx++) {
            for (int dz = -static_cast<int>(radius); dz <= static_cast<int>(radius); dz++) {
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (dist < radius) {
                    leafPositions.push_back(glm::vec3(worldX + trunkThickness / 2.0f + dx, centerY + dy, worldZ + trunkThickness / 2.0f + dz));
                }
            }
        }
    }
    return leafPositions;
}

std::vector<glm::vec3> generateOakCanopy(int groundHeight, int trunkHeight, int trunkThickness, double worldX, double worldZ) {
    std::vector<glm::vec3> leaves;
    int centerY = groundHeight + trunkHeight + 2;
    float radius = 4.0f;
    float centerOffset = trunkThickness / 2.0f;
    for (int dy = -static_cast<int>(radius); dy <= static_cast<int>(radius); dy++) {
        for (int dx = -static_cast<int>(radius); dx <= static_cast<int>(radius); dx++) {
            for (int dz = -static_cast<int>(radius); dz <= static_cast<int>(radius); dz++) {
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (dist < radius) {
                    leaves.push_back(glm::vec3(worldX + centerOffset + dx, centerY + dy, worldZ + centerOffset + dz));
                }
            }
        }
    }
    return leaves;
}

// ---------------------- Raycasting ----------------------
glm::ivec3 raycastForBlock(bool place) {
    float t = 0.0f;
    while (t < 5.0f) {
        glm::vec3 front;
        front.x = cos(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        front = glm::normalize(front);
        glm::vec3 p = cameraPos + t * front;
        glm::ivec3 candidate = glm::ivec3(std::round(p.x), std::round(p.y), std::round(p.z));
        bool exists = false;
        int chunkX = static_cast<int>(std::floor(candidate.x / static_cast<float>(CHUNK_SIZE)));
        int chunkZ = static_cast<int>(std::floor(candidate.z / static_cast<float>(CHUNK_SIZE)));
        ChunkPos cp2{ chunkX, chunkZ };
        if (chunks.find(cp2) != chunks.end()) {
            Chunk& ch = chunks[cp2];
            for (const auto& vec : { ch.waterPositions, ch.grassPositions, ch.dirtPositions,
                                      ch.deepStonePositions, ch.lavaPositions, ch.treeTrunkPositions,
                                      ch.treeLeafPositions, ch.firLeafPositions, ch.waterLilyPositions,
                                      ch.fallenTreeTrunkPositions, ch.oakTrunkPositions, ch.oakLeafPositions,
                                      ch.leafPilePositions, ch.bushSmallPositions, ch.bushMediumPositions,
                                      ch.bushLargePositions, ch.sandPositions, ch.snowPositions })
            {
                for (const auto& pos : vec) {
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), 0.5f))) {
                        exists = true;
                        break;
                    }
                }
                if (exists)
                    break;
            }
        }
        TerrainPoint terrain = getTerrainHeight(p.x, p.z);
        if (!exists && terrain.isLand && candidate.y <= static_cast<int>(std::floor(terrain.height)))
            exists = true;
        if (exists) {
            if (!place)
                return candidate;
            else {
                glm::vec3 center = glm::vec3(candidate) + glm::vec3(0.5f);
                glm::vec3 diff = p - center;
                if (std::abs(diff.x) > std::abs(diff.y) && std::abs(diff.x) > std::abs(diff.z))
                    return candidate + glm::ivec3((diff.x > 0) ? 1 : -1, 0, 0);
                else if (std::abs(diff.y) > std::abs(diff.x) && std::abs(diff.y) > std::abs(diff.z))
                    return candidate + glm::ivec3(0, (diff.y > 0) ? 1 : -1, 0);
                else
                    return candidate + glm::ivec3(0, 0, (diff.z > 0) ? 1 : -1);
            }
        }
        t += 0.1f;
    }
    return glm::ivec3(-10000, -10000, -10000);
}

// ---------------------- Mouse Button Callback ----------------------
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            glm::ivec3 pos = raycastForBlock(false);
            if (pos.x != -10000) {
                // Removal logic (not implemented)
            }
        }
        else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            glm::ivec3 hitBlock = raycastForBlock(false);
            if (hitBlock.x != -10000) {
                int blockTypeToPlace = 0;
                // Determine block type...
                glm::ivec3 pos = raycastForBlock(true);
                if (pos.x != -10000) {
                    // Place block (not implemented)
                }
            }
        }
    }
}

// ---------------------- Quadtree Structures ----------------------
struct QuadtreeItem {
    ChunkPos pos;
    Chunk* chunk;
};

struct QuadtreeNode {
    int minX, minZ, maxX, maxZ;
    std::vector<QuadtreeItem> items;
    bool subdivided;
    QuadtreeNode* children[4];
    static const int capacity = 10;
    QuadtreeNode(int minX, int minZ, int maxX, int maxZ)
        : minX(minX), minZ(minZ), maxX(maxX), maxZ(maxZ), subdivided(false) {
        for (int i = 0; i < 4; i++) children[i] = nullptr;
    }
    ~QuadtreeNode() {
        for (int i = 0; i < 4; i++)
            if (children[i])
                delete children[i];
    }
    glm::vec3 getMinWorld() const { return glm::vec3(minX * CHUNK_SIZE, MIN_Y, minZ * CHUNK_SIZE); }
    glm::vec3 getMaxWorld() const { return glm::vec3((maxX + 1) * CHUNK_SIZE, 150.0f, (maxZ + 1) * CHUNK_SIZE); }
    bool contains(const ChunkPos& pos) const { return pos.x >= minX && pos.x <= maxX && pos.z >= minZ && pos.z <= maxZ; }
    void subdivide() {
        int midX = (minX + maxX) / 2, midZ = (minZ + maxZ) / 2;
        children[0] = new QuadtreeNode(minX, minZ, midX, midZ);
        children[1] = new QuadtreeNode(midX + 1, minZ, maxX, midZ);
        children[2] = new QuadtreeNode(minX, midZ + 1, midX, maxZ);
        children[3] = new QuadtreeNode(midX + 1, midZ + 1, maxX, maxZ);
        subdivided = true;
        for (const auto& item : items) {
            for (int i = 0; i < 4; i++) {
                if (children[i]->contains(item.pos)) {
                    children[i]->items.push_back(item);
                    break;
                }
            }
        }
        items.clear();
    }
    void insert(const QuadtreeItem& item) {
        if (!contains(item.pos))
            return;
        if (!subdivided && items.size() < capacity)
            items.push_back(item);
        else {
            if (!subdivided)
                subdivide();
            for (int i = 0; i < 4; i++) {
                if (children[i]->contains(item.pos)) {
                    children[i]->insert(item);
                    return;
                }
            }
        }
    }
    void query(const std::vector<Plane>& frustum, std::vector<Chunk*>& out) {
        glm::vec3 nodeMin = getMinWorld(), nodeMax = getMaxWorld();
        if (!aabbInFrustum(frustum, nodeMin, nodeMax))
            return;
        if (!subdivided) {
            for (const auto& item : items)
                out.push_back(item.chunk);
        }
        else {
            for (int i = 0; i < 4; i++)
                children[i]->query(frustum, out);
        }
    }
};

struct Quadtree {
    QuadtreeNode* root;
    Quadtree(int minX, int minZ, int maxX, int maxZ) { root = new QuadtreeNode(minX, minZ, maxX, maxZ); }
    ~Quadtree() { delete root; }
    void insert(const ChunkPos& pos, Chunk* chunk) {
        QuadtreeItem item{ pos, chunk };
        root->insert(item);
    }
    std::vector<Chunk*> query(const std::vector<Plane>& frustum) {
        std::vector<Chunk*> result;
        root->query(frustum, result);
        return result;
    }
};

// ---------------------- Collision Handling ----------------------
void handleCollision() {
    // Get the player's current AABB.
    glm::vec3 playerMin = cameraPos + getPlayerBoxMin();
    glm::vec3 playerMax = cameraPos + getPlayerBoxMax();

    // For each nearby chunk:
    int chunkX = static_cast<int>(std::floor(cameraPos.x / CHUNK_SIZE));
    int chunkZ = static_cast<int>(std::floor(cameraPos.z / CHUNK_SIZE));
    for (int cx = chunkX - 1; cx <= chunkX + 1; cx++) {
        for (int cz = chunkZ - 1; cz <= chunkZ + 1; cz++) {
            ChunkPos cp(cx, cz);
            if (chunks.find(cp) != chunks.end()) {
                Chunk& ch = chunks[cp];

                auto resolveAABB = [&](const glm::vec3& boxMin, const glm::vec3& boxMax) {
                    // Check if overlapping
                    if (playerMax.x > boxMin.x && playerMin.x < boxMax.x &&
                        playerMax.y > boxMin.y && playerMin.y < boxMax.y &&
                        playerMax.z > boxMin.z && playerMin.z < boxMax.z) {

                        // Calculate penetration depth along each axis
                        float penX = std::min(playerMax.x - boxMin.x, boxMax.x - playerMin.x);
                        float penY = std::min(playerMax.y - boxMin.y, boxMax.y - playerMin.y);
                        float penZ = std::min(playerMax.z - boxMin.z, boxMax.z - playerMin.z);

                        // Only resolve if penetration is significant
                        const float threshold = 0.01f;
                        if (penX < threshold && penY < threshold && penZ < threshold)
                            return;

                        // Resolve along the axis with the smallest penetration.
                        if (penX <= penY && penX <= penZ)
                        {
                            if (cameraPos.x < boxMin.x)
                                cameraPos.x -= penX;
                            else
                                cameraPos.x += penX;
                        }
                        else if (penY <= penX && penY <= penZ)
                        {
                            if (cameraPos.y < boxMin.y)
                                cameraPos.y -= penY;
                            else
                                cameraPos.y += penY;
                            // Also zero out vertical velocity if moving downward.
                            if (velocity.y < 0)
                                velocity.y = 0;
                        }
                        else
                        {
                            if (cameraPos.z < boxMin.z)
                                cameraPos.z -= penZ;
                            else
                                cameraPos.z += penZ;
                        }

                        // Recalculate player's AABB after adjustment.
                        playerMin = cameraPos + getPlayerBoxMin();
                        playerMax = cameraPos + getPlayerBoxMax();
                    }
                    };

                // Check against each solid block in the chunk.
                auto checkBlocks = [&](const std::vector<glm::vec3>& blocks) {
                    for (const auto& pos : blocks) {
                        glm::vec3 blockMin = pos;
                        glm::vec3 blockMax = pos + glm::vec3(1.0f);
                        resolveAABB(blockMin, blockMax);
                    }
                    };

                checkBlocks(ch.grassPositions);
                checkBlocks(ch.sandPositions);
                checkBlocks(ch.snowPositions);
                checkBlocks(ch.dirtPositions);
                checkBlocks(ch.treeTrunkPositions);
                checkBlocks(ch.oakTrunkPositions);
                checkBlocks(ch.ancientTrunkPositions);
                checkBlocks(ch.waterLilyPositions);
                checkBlocks(ch.deepStonePositions);
                // ... add other block arrays as needed.
            }
        }
    }

    // For water areas, only apply minimal sinking correction.
    TerrainPoint tp = getTerrainHeight(cameraPos.x, cameraPos.z);
    if (!tp.isLand && cameraPos.y < -1.0f)
        cameraPos.y = -1.0f;
}
// ---------------------- Chunk Mesh Generation ----------------------
void generateChunkMesh(Chunk& chunk, int chunkX, int chunkZ) {
    if (!chunk.needsMeshUpdate)
        return;
    chunk.waterPositions.clear();
    chunk.grassPositions.clear();
    chunk.sandPositions.clear();
    chunk.snowPositions.clear();
    chunk.dirtPositions.clear();
    chunk.deepStonePositions.clear();
    chunk.lavaPositions.clear();
    chunk.treeTrunkPositions.clear();
    chunk.treeLeafPositions.clear();
    chunk.firLeafPositions.clear();
    chunk.waterLilyPositions.clear();
    chunk.fallenTreeTrunkPositions.clear();
    chunk.oakTrunkPositions.clear();
    chunk.oakLeafPositions.clear();
    chunk.leafPilePositions.clear();
    chunk.bushSmallPositions.clear();
    chunk.bushMediumPositions.clear();
    chunk.bushLargePositions.clear();
    chunk.branchPositions.clear();
    chunk.ancientTrunkPositions.clear();
    chunk.ancientLeafPositions.clear();
    chunk.ancientBranchPositions.clear();
    chunk.auroraPositions.clear();
    chunk.icePositions.clear();
    for (int x = 0; x < CHUNK_SIZE; x++) {
        for (int z = 0; z < CHUNK_SIZE; z++) {
            double worldX = chunkX * CHUNK_SIZE + x;
            double worldZ = chunkZ * CHUNK_SIZE + z;
            TerrainPoint terrain = getTerrainHeight(worldX, worldZ);
            if (!terrain.isLand) {
                chunk.waterPositions.push_back(glm::vec3(worldX, 0.0f, worldZ));
                if (x > 3 && x < CHUNK_SIZE - 3 && z > 3 && z < CHUNK_SIZE - 3 &&
                    (x % 7 == 3) && (z % 7 == 3)) {
                    bool canPlaceLily = true;
                    for (int dx = -3; dx <= 3; dx++) {
                        for (int dz = -3; dz <= 3; dz++) {
                            TerrainPoint neighbor = getTerrainHeight(worldX + dx, worldZ + dz);
                            if (neighbor.isLand) { canPlaceLily = false; break; }
                        }
                        if (!canPlaceLily) break;
                    }
                    if (canPlaceLily) {
                        int hashVal = std::abs((static_cast<int>(worldX) * 91321) ^ (static_cast<int>(worldZ) * 7817));
                        if (hashVal % 100 < 1) {
                            for (int dx = -6; dx < 6; dx++) {
                                for (int dz = -6; dz < 6; dz++) {
                                    // Skip the 3-block corners
                                    if ((dx <= -5 || dx >= 4) && (dz <= -5 || dz >= 4))
                                        continue;
                                    chunk.waterLilyPositions.push_back(glm::vec3(worldX + dx, 0.2f, worldZ + dz));
                                }
                            }
                        }
                    }
                }
                for (int y = -1; y >= MIN_Y; y--) {
                    double caveVal = caveNoise.noise(worldX * 0.04, y * 0.04, worldZ * 0.04);
                    if (caveVal < 0.6)
                        chunk.deepStonePositions.push_back(glm::vec3(worldX, y, worldZ));
                    else {
                        if (y < 0)
                            chunk.waterPositions.push_back(glm::vec3(worldX, y, worldZ));
                    }
                }
            }
            else {
                int groundHeight = static_cast<int>(std::floor(terrain.height));
                if (chunkX >= 160)
                    chunk.sandPositions.push_back(glm::vec3(worldX, groundHeight, worldZ));
                else if (chunkZ <= -160)
                    chunk.snowPositions.push_back(glm::vec3(worldX, groundHeight, worldZ));
                else
                    chunk.grassPositions.push_back(glm::vec3(worldX, groundHeight, worldZ));
                chunk.dirtPositions.push_back(glm::vec3(worldX, groundHeight - 1, worldZ));
                for (int y = groundHeight - 2; y >= MIN_Y; y--) {
                    if (y >= 0) {
                        chunk.deepStonePositions.push_back(glm::vec3(worldX, y, worldZ));
                    }
                    else {
                        double caveVal = caveNoise.noise(worldX * 0.1, y * 0.1, worldZ * 0.1);
                        if (caveVal < -0.8)
                            chunk.deepStonePositions.push_back(glm::vec3(worldX, y, worldZ));
                        else {
                            double liquidVal = lavaCaveNoise.noise(worldX * 0.02, y * 0.02, worldZ * 0.02);
                            if (liquidVal < 0.3) {
                                if (worldZ / CHUNK_SIZE <= -20)
                                    chunk.icePositions.push_back(glm::vec3(worldX, y, worldZ));
                                else
                                    chunk.waterPositions.push_back(glm::vec3(worldX, y, worldZ));
                            }
                            else
                                chunk.lavaPositions.push_back(glm::vec3(worldX, y, worldZ));
                        }
                    }
                }
                int currentChunkZ = static_cast<int>(std::floor(worldZ / (double)CHUNK_SIZE));
                if (terrain.height > 2.0) {
                    int intWorldX = static_cast<int>(worldX);
                    int intWorldZ = static_cast<int>(worldZ);
                    // Parameters for pine trees
                    int trunkHeight = 80, trunkThickness = 4;
                    int extraBottom = 15;  // extra trunk blocks below ground (only for the trunk)
                    int extraHeight = 90;  // extra trunk blocks on top (if desired)
                    // ----- Pine Tree Branch for far chunks -----
                    if (currentChunkZ <= -40) {
                        int hashValPine = std::abs((intWorldX * 73856093) ^ (intWorldZ * 19349663));
                        // For pine trees, we want the canopy to remain at groundHeight + trunkHeight.
                        // So we add trunk blocks from groundHeight + 1 up to groundHeight + trunkHeight
                        // and then add extra trunk blocks below ground (without affecting canopy).
                        glm::vec3 pineBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                        if (hashValPine % 2000 < 1 && !treeCollision(chunk.treeTrunkPositions, pineBase)) {
                            // Draw main trunk (above ground)
                            for (int i = 1; i <= trunkHeight; i++) {
                                for (int tx = 0; tx < trunkThickness; tx++) {
                                    for (int tz = 0; tz < trunkThickness; tz++) {
                                        chunk.treeTrunkPositions.push_back(
                                            glm::vec3(worldX + tx, groundHeight + i, worldZ + tz)
                                        );
                                    }
                                }
                            }
                            // Draw extra trunk blocks below ground
                            for (int i = 0; i < extraBottom; i++) {
                                for (int tx = 0; tx < trunkThickness; tx++) {
                                    for (int tz = 0; tz < trunkThickness; tz++) {
                                        chunk.treeTrunkPositions.push_back(
                                            glm::vec3(worldX + tx, groundHeight - i, worldZ + tz)
                                        );
                                    }
                                }
                            }
                            // Optionally, if you still want extra height on top, you could add:
                            for (int i = trunkHeight + 1; i <= trunkHeight + extraHeight; i++) {
                                for (int tx = 0; tx < trunkThickness; tx++) {
                                    for (int tz = 0; tz < trunkThickness; tz++) {
                                        chunk.treeTrunkPositions.push_back(
                                            glm::vec3(worldX + tx, groundHeight + i, worldZ + tz)
                                        );
                                    }
                                }
                            }
                            // Generate canopy at the original trunk height (unchanged)
                            std::vector<glm::vec3> pineCanopy = generatePineCanopy(groundHeight, trunkHeight, trunkThickness, worldX, worldZ);
                            chunk.treeLeafPositions.insert(chunk.treeLeafPositions.end(), pineCanopy.begin(), pineCanopy.end());
                        }
                    }
                    // ----- Pine Tree Branch for other areas -----
                    else if ((chunkX < 20) || (currentChunkZ >= 40)) {
                        if (currentChunkZ < 40) {
                            int hashValPine = std::abs((intWorldX * 73856093) ^ (intWorldZ * 19349663));
                            glm::vec3 pineBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                            if (hashValPine % 2000 < 1 && !treeCollision(chunk.treeTrunkPositions, pineBase)) {
                                // Draw main trunk above ground
                                for (int i = 1; i <= trunkHeight; i++) {
                                    for (int tx = 0; tx < trunkThickness; tx++) {
                                        for (int tz = 0; tz < trunkThickness; tz++) {
                                            chunk.treeTrunkPositions.push_back(
                                                glm::vec3(worldX + tx, groundHeight + i, worldZ + tz)
                                            );
                                        }
                                    }
                                }
                                // Extra blocks below ground
                                for (int i = 0; i < extraBottom; i++) {
                                    for (int tx = 0; tx < trunkThickness; tx++) {
                                        for (int tz = 0; tz < trunkThickness; tz++) {
                                            chunk.treeTrunkPositions.push_back(
                                                glm::vec3(worldX + tx, groundHeight - i, worldZ + tz)
                                            );
                                        }
                                    }
                                }
                                // Original (only converting the very top block):
                                for (int i = 1; i <= trunkHeight; i++) {
                                    for (int tx = 0; tx < trunkThickness; tx++) {
                                        for (int tz = 0; tz < trunkThickness; tz++) {
                                            glm::vec3 pos = glm::vec3(worldX + tx, groundHeight + i, worldZ + tz);
                                            // If this is one of the top 8 blocks, turn it into a leaf block.
                                            if (i > trunkHeight - 1)
                                                chunk.treeLeafPositions.push_back(pos);
                                            else
                                                chunk.treeTrunkPositions.push_back(pos);
                                        }
                                    }
                                }

                                std::vector<glm::vec3> pineCanopy = generatePineCanopy(groundHeight, trunkHeight, trunkThickness, worldX, worldZ);
                                chunk.treeLeafPositions.insert(chunk.treeLeafPositions.end(), pineCanopy.begin(), pineCanopy.end());
                            }
                        } {
                            int hashValFir = std::abs((intWorldX * 83492791) ^ (intWorldZ * 19349663));
                            glm::vec3 firBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                            if (hashValFir % 2000 < 1 && !treeCollision(chunk.treeTrunkPositions, firBase)) {
                                int trunkHeightFir = 40, trunkThicknessFir = 3;
                                for (int i = 1; i <= trunkHeightFir; i++) {
                                    for (int tx = 0; tx < trunkThicknessFir; tx++) {
                                        for (int tz = 0; tz < trunkThicknessFir; tz++) {
                                            chunk.treeTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                        }
                                    }
                                }
                                std::vector<glm::vec3> firCanopy = generateFirCanopy(groundHeight, trunkHeightFir, trunkThicknessFir, worldX, worldZ);
                                chunk.firLeafPositions.insert(chunk.firLeafPositions.end(), firCanopy.begin(), firCanopy.end());
                            }
                        } {
                            int hashValOak = std::abs((intWorldX * 92821) ^ (intWorldZ * 123457));
                            glm::vec3 oakBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                            if (hashValOak % 1000 < 1 && !treeCollision(chunk.oakTrunkPositions, oakBase)) {
                                int trunkHeightOak = 7, trunkThicknessOak = 2;
                                for (int i = 1; i <= trunkHeightOak; i++) {
                                    for (int tx = 0; tx < trunkThicknessOak; tx++) {
                                        for (int tz = 0; tz < trunkThicknessOak; tz++) {
                                            chunk.oakTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                        }
                                    }
                                }
                                std::vector<glm::vec3> oakCanopy = generateOakCanopy(groundHeight, trunkHeightOak, trunkThicknessOak, worldX, worldZ);
                                chunk.oakLeafPositions.insert(chunk.oakLeafPositions.end(), oakCanopy.begin(), oakCanopy.end());
                            }
                        } {
                            int hashValAncient = std::abs((intWorldX * 112233) ^ (intWorldZ * 445566));
                            glm::vec3 ancientBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                            if (hashValAncient % 3000 < 1 && !treeCollision(chunk.ancientTrunkPositions, ancientBase)) {
                                int trunkHeightAncient = 30, trunkThicknessAncient = 3;
                                for (int i = 1; i <= trunkHeightAncient; i++) {
                                    for (int tx = 0; tx < trunkThicknessAncient; tx++) {
                                        for (int tz = 0; tz < trunkThicknessAncient; tz++) {
                                            chunk.ancientTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                        }
                                    }
                                }
                                int centerY = groundHeight + trunkHeightAncient;
                                float canopyRadius = 5.0f;
                                for (int dy = -static_cast<int>(canopyRadius); dy <= static_cast<int>(canopyRadius); dy++) {
                                    for (int dx = -static_cast<int>(canopyRadius); dx <= static_cast<int>(canopyRadius); dx++) {
                                        for (int dz = -static_cast<int>(canopyRadius); dz <= static_cast<int>(canopyRadius); dz++) {
                                            float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                                            if (dist < canopyRadius)
                                                chunk.ancientLeafPositions.push_back(glm::vec3(worldX + trunkThicknessAncient / 2.0f + dx, centerY + dy, worldZ + trunkThicknessAncient / 2.0f + dz));
                                        }
                                    }
                                }
                                int branchBaseHeights[4] = { 7, 13, 19, 25 };
                                for (int b = 0; b < 4; b++) {
                                    int randomOffset = (rand() % 3) - 1;
                                    int branchStart = branchBaseHeights[b] + randomOffset;
                                    float branchRot = (b * 90.0f) * (3.14159f / 180.0f);
                                    glm::vec3 branchStartPos = glm::vec3(worldX + trunkThicknessAncient / 2.0f, groundHeight + branchStart, worldZ + trunkThicknessAncient / 2.0f);
                                    int branchLength = 10 + (rand() % 3);
                                    for (int i = 1; i <= branchLength; i++) {
                                        float bx = cos(branchRot) * i;
                                        float bz = sin(branchRot) * i;
                                        glm::vec3 branchBlockPos = branchStartPos + glm::vec3(bx, 0, bz);
                                        chunk.ancientBranchPositions.push_back(branchBlockPos);
                                    }
                                    glm::vec3 tip = branchStartPos + glm::vec3(cos(branchRot) * (branchLength + 1), 0, sin(branchRot) * (branchLength + 1));
                                    for (int dx = -1; dx <= 1; dx++) {
                                        for (int dy = -1; dy <= 1; dy++) {
                                            for (int dz = -1; dz <= 1; dz++) {
                                                if (glm::length(glm::vec3(dx, dy, dz)) < 1.5f)
                                                    chunk.ancientLeafPositions.push_back(tip + glm::vec3(dx, dy, dz));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    int hashValFallen = std::abs((intWorldX * 92821) ^ (intWorldZ * 68917));
                    bool nearWater = false;
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dz = -1; dz <= 1; dz++) {
                            TerrainPoint neighbor = getTerrainHeight(worldX + dx, worldZ + dz);
                            if (!neighbor.isLand) { nearWater = true; break; }
                        }
                        if (nearWater) break;
                    }
                    if (nearWater && hashValFallen % 500 < 1) {
                        int maxSearch = 20;
                        float angle = static_cast<float>(hashValFallen % 360);
                        float rad = glm::radians(angle);
                        int backLength = 0;
                        for (; backLength < maxSearch; backLength++) {
                            double sampleX = worldX - (backLength + 1) * cos(rad);
                            double sampleZ = worldZ - (backLength + 1) * sin(rad);
                            TerrainPoint sample = getTerrainHeight(sampleX, sampleZ);
                            if (!sample.isLand) break;
                        }
                        int forwardLength = 0;
                        for (; forwardLength < maxSearch; forwardLength++) {
                            double sampleX = worldX + (forwardLength + 1) * cos(rad);
                            double sampleZ = worldZ + (forwardLength + 1) * sin(rad);
                            TerrainPoint sample = getTerrainHeight(sampleX, sampleZ);
                            if (!sample.isLand) break;
                        }
                        int totalLength = backLength + forwardLength + 1;
                        if (totalLength >= 6) {
                            int thickness = 2;
                            for (int i = 0; i < totalLength; i++) {
                                float posX = worldX - backLength * cos(rad) + i * cos(rad);
                                float posZ = worldZ - backLength * sin(rad) + i * sin(rad);
                                for (int tx = 0; tx < thickness; tx++) {
                                    for (int tz = 0; tz < thickness; tz++) {
                                        float localX = posX + tx - thickness / 2.0f;
                                        float localZ = posZ + tz - thickness / 2.0f;
                                        chunk.fallenTreeTrunkPositions.push_back(glm::vec3(localX, groundHeight + 1, worldZ + tz));
                                    }
                                }
                            }
                        }
                    }
                    int hashValPile = std::abs((intWorldX * 412871) ^ (intWorldZ * 167591));
                    if (hashValPile % 300 < 1) {
                        int pileSize = (hashValPile % 4) + 3;
                        for (int i = 0; i < pileSize; i++) {
                            int px = (hashValPile + i * 13) % 3 - 1;
                            int pz = (hashValPile + i * 7) % 3 - 1;
                            float placeX = worldX + px;
                            float placeZ = worldZ + pz;
                            chunk.leafPilePositions.push_back(glm::vec3(placeX, groundHeight + 1, worldZ + pz));
                        }
                    }
                    {
                        int hashValBushSmall = std::abs((intWorldX * 17771) ^ (intWorldZ * 55117));
                        if (hashValBushSmall % 700 < 1) {
                            int centerY = groundHeight + 1;
                            float radius = 1.0f;
                            for (int dx = -1; dx <= 1; dx++) {
                                for (int dz = -1; dz <= 1; dz++) {
                                    if (glm::length(glm::vec2(dx, dz)) <= radius)
                                        chunk.bushSmallPositions.push_back(glm::vec3(worldX + dx, centerY, worldZ + dz));
                                }
                            }
                        }
                        int hashValBushMed = std::abs((intWorldX * 18323) ^ (intWorldZ * 51511));
                        if (hashValBushMed % 1000 < 2) {
                            int centerY = groundHeight + 1;
                            float radius = 2.0f;
                            for (int dx = -2; dx <= 2; dx++) {
                                for (int dz = -2; dz <= 2; dz++) {
                                    if (glm::length(glm::vec2(dx, dz)) <= radius)
                                        chunk.bushMediumPositions.push_back(glm::vec3(worldX + dx, centerY, worldZ + dz));
                                }
                            }
                        }
                        int hashValBushLarge = std::abs((intWorldX * 23719) ^ (intWorldZ * 41389));
                        if (hashValBushLarge % 1200 < 1) {
                            int centerY = groundHeight + 1;
                            float radius = 3.0f;
                            for (int dx = -3; dx <= 3; dx++) {
                                for (int dz = -3; dz <= 3; dz++) {
                                    if (glm::length(glm::vec2(dx, dz)) <= radius)
                                        chunk.bushLargePositions.push_back(glm::vec3(worldX + dx, centerY, worldZ + dz));
                                }
                            }
                        }
                    }
                    {
                        int hashValBranch = std::abs((intWorldX * 12345) ^ (intWorldZ * 6789));
                        if (hashValBranch % 1000 < 1) {
                            float rot = (hashValBranch % 360) * (3.14159f / 180.0f);
                            chunk.branchPositions.push_back(glm::vec4(worldX + 0.5f, groundHeight + 0.5f, worldZ + 0.5f, rot));
                        }
                    }
                }
            }
        }
    }
    for (int x = 0; x < CHUNK_SIZE; x++) {
        for (int z = 0; z < CHUNK_SIZE; z++) {
            double worldX = chunkX * CHUNK_SIZE + x;
            double worldZ = chunkZ * CHUNK_SIZE + z;
            for (int y = 135; y <= 136; y++) {
                double n = auroraNoise.noise(worldX * 0.1, y * 0.1, worldZ * 0.1);
                if (n > 0.44)
                    chunk.auroraPositions.push_back(glm::vec3(worldX, y, worldZ));
            }
        }
    }
    chunk.needsMeshUpdate = false;
    visitedChunks.insert(ChunkPos(chunkX, chunkZ));
    bigMapDirty = true;
}


// ---------------------- Chunk Update ----------------------
void updateChunks() {
    int playerChunkX = static_cast<int>(std::floor(cameraPos.x / CHUNK_SIZE));
    int playerChunkZ = static_cast<int>(std::floor(cameraPos.z / CHUNK_SIZE));
    int renderDistance = static_cast<int>(RENDER_DISTANCE);
    for (auto it = chunks.begin(); it != chunks.end();) {
        int dx = std::abs(it->first.x - playerChunkX);
        int dz = std::abs(it->first.z - playerChunkZ);
        if (dx > renderDistance || dz > renderDistance)
            it = chunks.erase(it);
        else
            ++it;
    }
    for (int x = playerChunkX - renderDistance; x <= playerChunkX + renderDistance; x++) {
        for (int z = playerChunkZ - renderDistance; z <= playerChunkZ + renderDistance; z++) {
            ChunkPos pos{ x, z };
            if (chunks.find(pos) == chunks.end())
                chunks[pos] = Chunk();
            generateChunkMesh(chunks[pos], x, z);
        }
    }
}

// ---------------------- Input Handling ----------------------

void processInput(GLFWwindow* window) {
    static bool nWasPressed = false;
    if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) {
        if (!nWasPressed) {
            minimapEnabled = !minimapEnabled;
            nWasPressed = true;
        }
    }
    else {
        nWasPressed = false;
    }

    // Toggle mode (prone / swim / paraglide) with P.
    static bool pWasPressed = false;
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
        if (!pWasPressed) {
            TerrainPoint tp = getTerrainHeight(cameraPos.x, cameraPos.z);
            float groundY = static_cast<float>(std::floor(tp.height)) + 1.0f;
            bool onGround = (cameraPos.y <= groundY + 0.1f);
            bool actuallyInWater = (!tp.isLand) && (cameraPos.y <= WATER_SURFACE + 0.1f);

            bool inAir = !onGround;
            if (actuallyInWater) {
                // Toggle swimming mode: if not in swimming, switch to it; otherwise revert to standing.
                playerMode = (playerMode == 2 ? 0 : 2);
            }
            else if (inAir) {
                playerMode = (playerMode == 3 ? 0 : 3);
            }
            else {
                // On land
                playerMode = (playerMode == 1 ? 0 : 1);
            }
            pWasPressed = true;
        }
    }
    else {
        pWasPressed = false;
    }

    // If fullscreen map is active, handle panning and return.
    if (fullscreenMap) {
        const float panSpeed = 500.0f * deltaTime;
        bool panned = false;
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) { bigMapPanX -= panSpeed; panned = true; }
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) { bigMapPanX += panSpeed; panned = true; }
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) { bigMapPanZ -= panSpeed; panned = true; }
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) { bigMapPanZ += panSpeed; panned = true; }
        if (panned) { bigMapDirty = true; }
        return;
    }
    const float waterJumpImpulse = 6.0f; // adjust as desired

    // Handle movement based on current playerMode.
    switch (playerMode) {
 
    case 0: // Standing/Walking
    {
        // -----------------------------------------------
        // 1) Calculate intended horizontal velocity
        // -----------------------------------------------
        glm::vec3 forward = glm::vec3(cos(glm::radians(cameraYaw)), 0.0f, sin(glm::radians(cameraYaw)));
        glm::vec3 right = glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f));
        glm::vec3 walkDir(0.0f);
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            walkDir += forward;
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            walkDir -= forward;
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            walkDir -= right;
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            walkDir += right;
        if (glm::length(walkDir) > 0.001f)
            walkDir = glm::normalize(walkDir);

        float currentSpeed = walkSpeed;
        static bool sprintActive = false;
        if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
            sprintActive = true;
        else
            sprintActive = false;
        if (sprintActive)
            currentSpeed *= 2.0f;

        // Desired horizontal velocity.
        glm::vec3 desiredHorizVel = walkDir * currentSpeed;

        // -----------------------------------------------
        // 2) Smoothly blend the horizontal velocity
        // -----------------------------------------------
        glm::vec3 oldVel = velocity;
        glm::vec3 newVel = oldVel;
        float blendFactor = 10.0f * deltaTime; // Adjust as needed.
        glm::vec3 oldHoriz(oldVel.x, 0.0f, oldVel.z);
        glm::vec3 nextHoriz = glm::mix(oldHoriz, desiredHorizVel, blendFactor);
        newVel.x = nextHoriz.x;
        newVel.z = nextHoriz.z;
        velocity = newVel;

        // -----------------------------------------------
        // 3) Handle Jump Input and Double-Jump Mantle
        // -----------------------------------------------
        // Check if we're on the ground.
        TerrainPoint tp = getTerrainHeight(cameraPos.x, cameraPos.z);
        float groundY = floor(tp.height) + 1.0f;
        bool onGround = (cameraPos.y <= groundY + 0.1f);
        bool inWater = (!tp.isLand && cameraPos.y <= WATER_SURFACE + 0.1f);

        // Static state flags.
        static bool hasJumped = false;    // Set true after a normal jump.
        static bool mantleUsed = false;   // Prevents multiple mantles per airborne state.
        static bool spaceWasPressed = false; // Edge detector.
        static float lastMantleTime = 0.0f;  // Cooldown timer.
        float currentTime = static_cast<float>(glfwGetTime());
        const float mantleCooldown = 0.5f; // Cooldown in seconds.

        // Reset jump state on landing.
        if (onGround) {
            hasJumped = false;
            mantleUsed = false;
        }

        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
            if (!spaceWasPressed) {
                if (onGround || inWater) {
                    // Normal jump.
                    if (inWater)
                        velocity.y = waterJumpImpulse; // e.g. 6.0f
                    else
                        velocity.y = walkJumpImpulse;  // e.g. 4.8f
                    hasJumped = true;
                }
                else if (hasJumped && !mantleUsed && (currentTime - lastMantleTime >= mantleCooldown)) {
                    // Only attempt mantle if already jumped and falling slowly.
                    // Require that you're falling but not too fast.
                    if (velocity.y < -0.1f && velocity.y > -5.0f) {
                        // Use your raycast function to check for a candidate ledge.
                        glm::ivec3 candidate = raycastForBlock(false);
                        if (candidate.x != -10000) {
                            // Assume blocks are 1 unit tall.
                            float blockTop = candidate.y + 1.0f;
                            // Vertical thresholds.
                            const float mantleMinOffset = 0.2f;
                            const float mantleMaxOffset = 1.5f;
                            if (blockTop > cameraPos.y + mantleMinOffset && blockTop <= cameraPos.y + mantleMaxOffset) {
                                // Compute candidate block's center.
                                glm::vec3 candidateCenter = glm::vec3(candidate) + glm::vec3(0.5f);
                                // Horizontal distance.
                                float horizDist = glm::length(glm::vec2(candidateCenter.x - cameraPos.x,
                                    candidateCenter.z - cameraPos.z));
                                // Ensure candidate is in front: dot product between (candidateCenter - cameraPos) and forward.
                                float forwardDot = glm::dot(glm::normalize(candidateCenter - cameraPos), forward);
                                // Define thresholds:
                                const float minMantleHorizDistance = 0.5f;  // must not be too close (to avoid wall collisions)
                                const float maxMantleHorizDistance = 1.0f;  // candidate should be reasonably in front.
                                if (horizDist >= minMantleHorizDistance &&
                                    horizDist <= maxMantleHorizDistance &&
                                    forwardDot > 0.7f) // candidate is mostly in front
                                {
                                    // Apply mantle impulse.
                                    const float mantleUpImpulse = 6.0f;      // upward boost (tweak as needed)
                                    const float mantleForwardImpulse = 4.0f; // forward boost (tweak as needed)
                                    velocity.y = mantleUpImpulse;
                                    velocity += forward * mantleForwardImpulse;
                                    mantleUsed = true;
                                    lastMantleTime = currentTime;
                                }
                            }
                        }
                    }
                }
                spaceWasPressed = true;
            }
        }
        else {
            spaceWasPressed = false;
        }

        // -----------------------------------------------
        // 4) Apply Gravity
        // -----------------------------------------------
        velocity.y -= walkGravity * deltaTime;

        // -----------------------------------------------
        // 5) Update Position Separately (Horizontal & Vertical)
        // -----------------------------------------------
        // Horizontal update (X/Z):
        glm::vec3 horizVel(velocity.x, 0.0f, velocity.z);
        glm::vec3 newPosHoriz = cameraPos + horizVel * deltaTime;
        cameraPos.x = newPosHoriz.x;
        cameraPos.z = newPosHoriz.z;
        // Vertical update (Y):
        glm::vec3 vertVel(0.0f, velocity.y, 0.0f);
        glm::vec3 newPosVert = cameraPos + vertVel * deltaTime;
        cameraPos.y = newPosVert.y;

        break;
    }



    case 1: // Prone (crawling on land)
    {
        // Calculate horizontal movement only.
        glm::vec3 forward;
        forward.x = cos(glm::radians(cameraYaw));
        forward.y = 0.0f;
        forward.z = sin(glm::radians(cameraYaw));
        forward = glm::normalize(forward);

        glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
        glm::vec3 moveDir(0.0f);
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            moveDir += forward;
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            moveDir -= forward;
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            moveDir -= right;
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            moveDir += right;
        if (glm::length(moveDir) > 0.001f)
            moveDir = glm::normalize(moveDir);

        float proneSpeed = walkSpeed * 0.5f;
        cameraPos += moveDir * proneSpeed * deltaTime;

        // Immediately clamp vertical position to ground level.
        TerrainPoint tp = getTerrainHeight(cameraPos.x, cameraPos.z);
        float groundY = std::floor(tp.height) + 1.0f;
        cameraPos.y = groundY;
    }
    break;
    case 2: // Swimming mode
    {
        // Horizontal movement: use WASD based on horizontal view only.
        glm::vec3 forward;
        forward.x = cos(glm::radians(cameraYaw));
        forward.y = 0.0f; // ignore vertical component for horizontal movement
        forward.z = sin(glm::radians(cameraYaw));
        forward = glm::normalize(forward);

        glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
        glm::vec3 horizontalDir(0.0f);
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            horizontalDir += forward;
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            horizontalDir -= forward;
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            horizontalDir += right;
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            horizontalDir -= right;
        if (glm::length(horizontalDir) > 0.001f)
            horizontalDir = glm::normalize(horizontalDir);

        float swimSpeed = 6.0f;
        cameraPos += horizontalDir * swimSpeed * deltaTime;

        // Vertical movement: in swimming mode, SPACE and LEFT SHIFT control vertical motion linearly.
        float verticalSpeed = 0.0f;
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            verticalSpeed = 2.0f; // constant upward speed
        else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            verticalSpeed = -2.0f; // constant downward speed

        // Apply vertical movement directly.
        cameraPos.y += verticalSpeed * deltaTime;

        // Optional: if no vertical input, you might let the player slowly sink:
        if (verticalSpeed == 0.0f)
        {
            // Only sink if you're below the water surface.
            if (cameraPos.y > WATER_SURFACE)
                cameraPos.y -= 0.5f * deltaTime;
        }

        // Also, if the player is at or above the water surface,
        // you could add a gentle bobbing effect:
        if (cameraPos.y >= WATER_SURFACE && verticalSpeed >= 0.0f)
        {
            cameraPos.y = WATER_SURFACE + 0.1f * sin(glfwGetTime() * 2.0f);
        }
    }
    break;

    case 3: // Paraglider
    {
        // 1) Compute the look direction
        glm::vec3 viewDir;
        viewDir.x = cos(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        viewDir.y = sin(glm::radians(pitch));
        viewDir.z = sin(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        viewDir = glm::normalize(viewDir);

        // 2) Gently turn velocity horizontally toward the camera’s heading
        float steeringFactor = 0.02f;
        float speed = glm::length(velocity);
        if (speed > 0.01f) {
            // Only rotate the XZ plane portion
            glm::vec3 velXZ(velocity.x, 0.0f, velocity.z);
            float horizSpeed = glm::length(velXZ);
            if (horizSpeed > 0.001f) {
                glm::vec3 horizDir = glm::normalize(velXZ);
                glm::vec3 desiredDir = glm::normalize(glm::vec3(viewDir.x, 0.0f, viewDir.z));
                glm::vec3 newDir = glm::mix(horizDir, desiredDir, steeringFactor);
                newDir = glm::normalize(newDir) * horizSpeed;
                velocity.x = newDir.x;
                velocity.z = newDir.z;
            }
        }

        // 3) Adjust speed by pitch (similar to your existing logic)
        //    Negative pitch => accelerate forward, positive => slow
        float pitchFactor = -pitch / 90.0f; // pitch: -89..+89 => pitchFactor in [-, +]
        float accel = baseAcceleration * pitchFactor;
        velocity += viewDir * accel * deltaTime;

        // 4) Reduced gravity
        float gliderGravity = gravityForce * 0.3f;
        velocity.y -= gliderGravity * deltaTime;

        // 5) Apply drag
        velocity *= dragFactor; // e.g. ~0.995

        // 6) Remove rocket–boost code
        //    If you had something like:
        //        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {...}
        //    we remove it entirely. No boost.

        // 7) Move the player
        cameraPos += velocity * deltaTime;

        // 8) Optional: if you want to auto–exit paraglider mode on ground
        TerrainPoint tp = getTerrainHeight(cameraPos.x, cameraPos.z);
        float groundY = floor(tp.height) + 1.0f;
        bool landed = (cameraPos.y <= groundY + 0.1f);
        if (landed) {
            // Snap to ground and revert to standing
            cameraPos.y = groundY;
            velocity.y = 0.0f;
            playerMode = 0;
        }
    }
    break;


    }

    handleCollision();
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// ---------------------- Mouse Callback ----------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);
    if (firstMouse) { lastX = xpos; lastY = ypos; firstMouse = false; }
    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos; lastY = ypos;
    const float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;
    cameraYaw += xoffset;
    pitch += yoffset;
    if (pitch > 89.0f)
        pitch = 89.0f;
    if (pitch < -89.0f)
        pitch = -89.0f;
}

// ---------------------- Toggle Map Mode ----------------------
void toggleMapMode(GLFWwindow* window) {
    static bool mWasPressed = false;
    if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) {
        if (!mWasPressed) {
            fullscreenMap = !fullscreenMap;
            if (fullscreenMap) {
                bigMapPanX = 0.0f;
                bigMapPanZ = 0.0f;
            }
            mWasPressed = true;
        }
    }
    else {
        mWasPressed = false;
    }
}

// ---------------------- Shader Sources for Main Scene, Minimap, and Sun/Moon ----------------------

// --- Main Scene Vertex Shader (with normals, texcoords, instancing) ---
const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec3 aPos;       // Vertex position
layout (location = 1) in vec3 aNormal;      // Vertex normal
layout (location = 2) in vec2 aTexCoord;    // Texture coordinate
layout (location = 3) in vec3 aOffset;      // Instance offset (constant per instance)
layout (location = 4) in float aRotation;   // Instance rotation

out vec2 TexCoord;
out vec3 ourColor;
out float instanceDistance;
out vec3 Normal;
out vec3 WorldPos;  // World-space position

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform int blockType;
uniform vec3 blockColors[25];
uniform vec3 cameraPos;
uniform float time;

void main(){
    vec3 pos;
    vec3 normal = aNormal;
    
    // For blocks that are not type 14 or 18:
    if(blockType != 14 && blockType != 18) {
        // WATER BLOCK (blockType == 1): Scale up and juggle uniformly.
        if(blockType == 1){
            float scaleFactor = 1.2; // Increase water block size
            // Scale the base vertex positions.
            vec3 scaledPos = aPos * scaleFactor;
            // Compute a uniform displacement for the whole instance.
            vec3 waterDisplacement;
            waterDisplacement.x = sin(time * 0.5 + aOffset.x * 1.3 + aOffset.y * 0.7) * 0.2;
            waterDisplacement.y = sin(time * 0.5 + aOffset.y * 1.3 + aOffset.z * 0.7) * 0.2;
            waterDisplacement.z = sin(time * 0.5 + aOffset.z * 1.3 + aOffset.x * 0.7) * 0.2;
            // Add scaled position, instance offset, and uniform displacement.
            pos = scaledPos + aOffset + waterDisplacement;
        }
        // WATER LILY (blockType == 5): Scale up and juggle uniformly.
        else if(blockType == 5) {
            float scaleFactor = 1.2; // Increase water lily size
            // If your water lily model is already centered, scale it.
            vec3 scaledPos = aPos * scaleFactor;
            // Compute a uniform displacement.
            vec3 lilyDisplacement;
            lilyDisplacement.x = sin(time * 0.5 + aOffset.x) * 0.1;
            lilyDisplacement.y = sin(time * 0.5 + aOffset.y) * 0.1;
            lilyDisplacement.z = sin(time * 0.5 + aOffset.z) * 0.1;
            pos = scaledPos + aOffset + lilyDisplacement;
        }
        // LEAF BLOCKS: Apply a uniform displacement.
        else if(blockType == 3 || blockType == 7 || blockType == 9 || blockType == 17){
            vec3 leafDisplacement;
            leafDisplacement.x = sin((aOffset.x + time) * 0.3) * 0.05;
            leafDisplacement.y = cos((aOffset.y + time) * 0.3) * 0.05;
            leafDisplacement.z = sin((aOffset.z + time) * 0.3) * 0.05;
            pos = aPos + aOffset + leafDisplacement;
        }
        // All other blocks: simply add the instance offset.
        else {
            pos = aPos + aOffset;
        }
    } 
    else {
        // Special handling for block types 14 and 18 (rotate/scaled models)
        if(blockType == 14) {
            float angle = aRotation;
            mat3 rot = mat3(
                cos(angle), 0.0, sin(angle),
                0.0,        1.0, 0.0,
               -sin(angle), 0.0, cos(angle)
            );
            mat3 scaleMat = mat3(0.3, 0.0, 0.0,
                                 0.0, 0.8, 0.0,
                                 0.0, 0.0, 0.3);
            pos = rot * (scaleMat * aPos) + aOffset;
            normal = rot * aNormal;
        } 
        else {
            pos = aPos + aOffset;
        }
    }
    
    // Additional adjustment for blockType 19 (if needed)
    if(blockType == 19){
        pos.y += sin(time + aOffset.x * 0.1) * 0.5;
    }
    
    // Compute world-space position.
    vec4 worldPos4 = model * vec4(pos, 1.0);
    WorldPos = worldPos4.xyz;
    
    gl_Position = projection * view * worldPos4;
    ourColor = blockColors[blockType];
    TexCoord = aTexCoord;
    
    // Compute instance distance for effects.
    if(blockType != 14 && blockType != 18) {
        if(gl_InstanceID > 0)
            instanceDistance = length(aOffset - cameraPos);
        else
            instanceDistance = length(vec3(model[3]) - cameraPos);
    } else {
        instanceDistance = length(aOffset - cameraPos);
    }
    
    Normal = normalize(mat3(model) * normal);
}
)";
// --- Main Scene Fragment Shader ---
const char* fragmentShaderSource = R"(
#version 330 core
in vec2 TexCoord;
in vec3 ourColor;
in float instanceDistance;
in vec3 Normal;
in vec3 WorldPos;  // World-space position from the vertex shader

out vec4 FragColor;

uniform int blockType;
uniform vec3 blockColors[25];
uniform vec3 lightDir;
uniform vec3 ambientLight;
uniform vec3 diffuseLight;
uniform float time;

// Simple noise function that returns a pseudo-random value for a given cell
float generateNoise(vec2 cell) {
    return fract(sin(dot(cell, vec2(12.9898, 78.233))) * 43758.5453);
}

void main(){
    // Special handling for translucent aurora blocks (blockType 19)
    if(blockType == 19){
        FragColor = vec4(ourColor, 0.1);
        return;
    }
    
    // Define grid parameters: both grid overlay and noise sampling use 12x12 cells.
    float gridSize = 12.0;
    float lineWidth = 0.03;
    float leafNoiseSize = 12.0;
    
    // For leaf blocks, we want the grid overlay plus a randomized noise "crack" effect.
    if(blockType == 3 || blockType == 7 || blockType == 9 || blockType == 17){
        // 1. Determine grid lines using the unshifted TexCoord.
        vec2 f = fract(TexCoord * gridSize);
        bool isGridLine = (f.x < lineWidth || f.y < lineWidth);
        
        // 2. Compute a per-block seed using the block's integer world coordinates.
        // Using floor(WorldPos.xy) makes sure that all fragments for a given block share the same seed.
        vec2 blockCoord = floor(WorldPos.xy);
        vec2 seed = fract(blockCoord * 0.12345);
        
        // 3. Use the seed to offset the TexCoord when sampling noise, using leafNoiseSize.
        vec2 cell = floor((TexCoord + seed) * leafNoiseSize);
        float noiseVal = generateNoise(cell);
        
        // 4. Threshold the noise to determine the crack transparency.
        float crackThreshold = 0.8;
        float noiseAlpha = (noiseVal > crackThreshold) ? 0.0 : 1.0;
        
        // 5. For grid lines, force black color and full opacity.
        float finalAlpha = isGridLine ? 1.0 : noiseAlpha;
        vec3 finalColor = isGridLine ? vec3(0.0) : ourColor;
        
        // 6. Apply simple diffuse lighting.
        vec3 norm = normalize(Normal);
        float diff = max(dot(norm, normalize(lightDir)), 0.0);
        vec3 lighting = ambientLight + diffuseLight * diff;
        finalColor *= lighting;
        
        FragColor = vec4(finalColor, finalAlpha);
        return;
    }
    
    // WATER BLOCKS (blockType 1): apply a wave effect.
    if(blockType == 1){
        vec3 waterColor = blockColors[1];
        float wave1 = sin(WorldPos.x * 0.1 + time * 2.0);
        float wave2 = cos(WorldPos.z * 0.1 + time * 2.0);
        float wave = (wave1 + wave2) * 0.5;
        waterColor *= (1.0 + wave * 0.1);
        
        vec3 norm = normalize(Normal);
        float diff = max(dot(norm, normalize(lightDir)), 0.0);
        vec3 lighting = ambientLight + diffuseLight * diff;
        vec3 finalColor = waterColor * lighting;
        
        FragColor = vec4(finalColor, 0.3);
    }
    else {
        // NON-WATER, NON-LEAF BLOCKS: grid overlay effect.
        vec2 f = fract(TexCoord * gridSize);
        vec3 baseColor;
        if(f.x < lineWidth || f.y < lineWidth)
            baseColor = vec3(0.0, 0.0, 0.0);
        else {
            float factor = instanceDistance / 100.0;
            vec3 offset = vec3(0.03 * factor, 0.03 * factor, 0.05 * factor);
            baseColor = ourColor + offset;
            baseColor = clamp(baseColor, 0.0, 1.0);
        }
        
        vec3 norm = normalize(Normal);
        float diff = max(dot(norm, normalize(lightDir)), 0.0);
        vec3 lighting = ambientLight + diffuseLight * diff;
        vec3 finalColor = baseColor * lighting;
        
        FragColor = vec4(finalColor, 1.0);
    }
}



)";

// --- Minimap Vertex Shader ---
const char* minimapVertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec3 aColor;
out vec3 ourColor;
uniform mat4 ortho;
void main(){
    ourColor = aColor;
    gl_Position = ortho * vec4(aPos, 0.0, 1.0);
}
)";

// --- Minimap Fragment Shader ---
const char* minimapFragmentShaderSource = R"(
#version 330 core
in vec3 ourColor;
out vec4 FragColor;
void main(){
    FragColor = vec4(ourColor, 1.0);
}
)";

// ---------------------- Sun/Moon Shaders (compiled above) ----------------------

// ---------------------- Cube Vertex Data ----------------------
float cubeVertices[] = {
    // Each vertex: position (3), normal (3), texCoord (2)
    // Front face (normal 0,0,1)
   -0.5f, -0.5f,  0.5f,    0,0,1,    0.0f, 0.0f,
    0.5f, -0.5f,  0.5f,    0,0,1,    1.0f, 0.0f,
    0.5f,  0.5f,  0.5f,    0,0,1,    1.0f, 1.0f,
    0.5f,  0.5f,  0.5f,    0,0,1,    1.0f, 1.0f,
   -0.5f,  0.5f,  0.5f,    0,0,1,    0.0f, 1.0f,
   -0.5f, -0.5f,  0.5f,    0,0,1,    0.0f, 0.0f,
   // Right face (normal 1,0,0)
   0.5f, -0.5f,  0.5f,    1,0,0,    0.0f, 0.0f,
   0.5f, -0.5f, -0.5f,    1,0,0,    1.0f, 0.0f,
   0.5f,  0.5f, -0.5f,    1,0,0,    1.0f, 1.0f,
   0.5f,  0.5f, -0.5f,    1,0,0,    1.0f, 1.0f,
   0.5f,  0.5f,  0.5f,    1,0,0,    0.0f, 1.0f,
   0.5f, -0.5f,  0.5f,    1,0,0,    0.0f, 0.0f,
   // Back face (normal 0,0,-1)
   0.5f, -0.5f, -0.5f,    0,0,-1,   0.0f, 0.0f,
  -0.5f, -0.5f, -0.5f,    0,0,-1,   1.0f, 0.0f,
  -0.5f,  0.5f, -0.5f,    0,0,-1,   1.0f, 1.0f,
  -0.5f,  0.5f, -0.5f,    0,0,-1,   1.0f, 1.0f,
   0.5f,  0.5f, -0.5f,    0,0,-1,   0.0f, 1.0f,
   0.5f, -0.5f, -0.5f,    0,0,-1,   0.0f, 0.0f,
   // Left face (normal -1,0,0)
  -0.5f, -0.5f, -0.5f,   -1,0,0,    0.0f, 0.0f,
  -0.5f, -0.5f,  0.5f,   -1,0,0,    1.0f, 0.0f,
  -0.5f,  0.5f,  0.5f,   -1,0,0,    1.0f, 1.0f,
  -0.5f,  0.5f,  0.5f,   -1,0,0,    1.0f, 1.0f,
  -0.5f,  0.5f, -0.5f,   -1,0,0,    0.0f, 1.0f,
  -0.5f, -0.5f, -0.5f,   -1,0,0,    0.0f, 0.0f,
  // Top face (normal 0,1,0)
 -0.5f,  0.5f,  0.5f,    0,1,0,    0.0f, 0.0f,
  0.5f,  0.5f,  0.5f,    0,1,0,    1.0f, 0.0f,
  0.5f,  0.5f, -0.5f,    0,1,0,    1.0f, 1.0f,
  0.5f,  0.5f, -0.5f,    0,1,0,    1.0f, 1.0f,
 -0.5f,  0.5f, -0.5f,    0,1,0,    0.0f, 1.0f,
 -0.5f,  0.5f,  0.5f,    0,1,0,    0.0f, 0.0f,
 // Bottom face (normal 0,-1,0)
-0.5f, -0.5f, -0.5f,    0,-1,0,   0.0f, 0.0f,
 0.5f, -0.5f, -0.5f,    0,-1,0,   1.0f, 0.0f,
 0.5f, -0.5f,  0.5f,    0,-1,0,   1.0f, 1.0f,
 0.5f, -0.5f,  0.5f,    0,-1,0,   1.0f, 1.0f,
-0.5f, -0.5f,  0.5f,    0,-1,0,   0.0f, 1.0f,
-0.5f, -0.5f, -0.5f,    0,-1,0,   0.0f, 0.0f
};

// ---------------------- Minimap Rendering ----------------------
void renderMinimap(GLuint minimapShaderProgram, GLuint minimapVAO, GLuint minimapVBO) {
    if (!minimapEnabled)
        return;

    if (!fullscreenMap) {
        static double lastMapUpdateTime = 0.0;
        static std::vector<float> cachedInterleaved;
        double currentTime = glfwGetTime();
        const float region = 96.0f;
        if (currentTime - lastMapUpdateTime > 1.0 || cachedInterleaved.empty()) {
            lastMapUpdateTime = currentTime;
            std::vector<glm::vec2> vertices;
            std::vector<glm::vec3> colors;
            int startX = static_cast<int>(cameraPos.x - region);
            int endX = static_cast<int>(cameraPos.x + region);
            int startZ = static_cast<int>(cameraPos.z - region);
            int endZ = static_cast<int>(cameraPos.z + region);
            for (int z = startZ; z < endZ; z++) {
                for (int x = startX; x < endX; x++) {
                    TerrainPoint tp = getTerrainHeight(x + 0.5, z + 0.5);
                    int blockType;
                    if (!tp.isLand)
                        blockType = 1;
                    else {
                        int cx = x / CHUNK_SIZE;
                        int cz = z / CHUNK_SIZE;
                        if (cx >= 120)
                            blockType = 22;
                        else if (cz <= -40)
                            blockType = 23;
                        else
                            blockType = 0;
                    }
                    vertices.push_back(glm::vec2(x, z));
                    vertices.push_back(glm::vec2(x + CHUNK_SIZE, z));
                    vertices.push_back(glm::vec2(x + CHUNK_SIZE, z + CHUNK_SIZE));
                    vertices.push_back(glm::vec2(x, z));
                    vertices.push_back(glm::vec2(x + CHUNK_SIZE, z + CHUNK_SIZE));
                    vertices.push_back(glm::vec2(x, z + CHUNK_SIZE));
                    glm::vec3 col;
                    switch (blockType) {
                    case 0:  col = glm::vec3(0.19f, 0.66f, 0.32f); break;
                    case 1:  col = glm::vec3(0.0f, 0.5f, 0.5f); break;
                    case 22: col = glm::vec3(0.93f, 0.79f, 0.69f); break;
                    case 23: col = glm::vec3(0.95f, 0.95f, 1.0f); break;
                    default: col = glm::vec3(1.0f); break;
                    }
                    for (int i = 0; i < 6; i++)
                        colors.push_back(col);
                }
            }
            cachedInterleaved.clear();
            for (size_t i = 0; i < vertices.size(); i++) {
                cachedInterleaved.push_back(vertices[i].x);
                cachedInterleaved.push_back(vertices[i].y);
                cachedInterleaved.push_back(colors[i].r);
                cachedInterleaved.push_back(colors[i].g);
                cachedInterleaved.push_back(colors[i].b);
            }
        }
        glViewport(WINDOW_WIDTH - 200, WINDOW_HEIGHT - 200, 200, 200);
        glm::mat4 ortho = glm::ortho((float)(cameraPos.x - region), (float)(cameraPos.x + region),
            (float)(cameraPos.z - region), (float)(cameraPos.z + region),
            -1.0f, 1.0f);
        glUseProgram(minimapShaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(minimapShaderProgram, "ortho"), 1, GL_FALSE, glm::value_ptr(ortho));
        glBindBuffer(GL_ARRAY_BUFFER, minimapVBO);
        glBufferData(GL_ARRAY_BUFFER, cachedInterleaved.size() * sizeof(float), cachedInterleaved.data(), GL_DYNAMIC_DRAW);
        glBindVertexArray(minimapVAO);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);
        glDrawArrays(GL_TRIANGLES, 0, cachedInterleaved.size() / 5);
        glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
    }
    else {
        if (bigMapDirty) {
            bigMapInterleaved.clear();
            std::vector<glm::vec2> vertices;
            std::vector<glm::vec3> colors;
            int regionChunks = 400;
            int halfRegion = regionChunks / 2;
            int centerChunkX = static_cast<int>(round(bigMapPanX / CHUNK_SIZE));
            int centerChunkZ = static_cast<int>(round(bigMapPanZ / CHUNK_SIZE));
            int startChunkX = centerChunkX - halfRegion;
            int endChunkX = centerChunkX + halfRegion;
            int startChunkZ = centerChunkZ - halfRegion;
            int endChunkZ = centerChunkZ + halfRegion;
            for (int cz = startChunkZ; cz < endChunkZ; cz++) {
                for (int cx = startChunkX; cx < endChunkX; cx++) {
                    ChunkPos cp(cx, cz);
                    int blockType;
                    if (chunks.find(cp) != chunks.end()) {
                        Chunk& ch = chunks[cp];
                        int waterCount = 0;
                        for (const auto& pos : ch.waterPositions) {
                            if (std::abs(pos.y - 0.0f) < 0.1f)
                                waterCount++;
                        }
                        if (waterCount > 5)
                            blockType = 1;
                        else if (!ch.sandPositions.empty())
                            blockType = 22;
                        else if (!ch.snowPositions.empty())
                            blockType = 23;
                        else
                            blockType = 0;
                    }
                    else {
                        blockType = getChunkTopBlock(cx, cz);
                    }
                    bool visited = (visitedChunks.find(cp) != visitedChunks.end());
                    glm::vec3 col;
                    switch (blockType) {
                    case 0:  col = glm::vec3(0.19f, 0.66f, 0.32f); break;
                    case 1:  col = glm::vec3(0.0f, 0.5f, 0.5f); break;
                    case 22: col = glm::vec3(0.93f, 0.79f, 0.69f); break;
                    case 23: col = glm::vec3(0.95f, 0.95f, 1.0f); break;
                    default: col = glm::vec3(1.0f); break;
                    }
                    if (!visited && blockType != 1)
                        col *= 0.5f;
                    int startX = cx * CHUNK_SIZE;
                    int startZ = cz * CHUNK_SIZE;
                    vertices.push_back(glm::vec2(startX, startZ));
                    vertices.push_back(glm::vec2(startX + CHUNK_SIZE, startZ));
                    vertices.push_back(glm::vec2(startX + CHUNK_SIZE, startZ + CHUNK_SIZE));
                    vertices.push_back(glm::vec2(startX, startZ));
                    vertices.push_back(glm::vec2(startX + CHUNK_SIZE, startZ + CHUNK_SIZE));
                    vertices.push_back(glm::vec2(startX, startZ + CHUNK_SIZE));
                    for (int i = 0; i < 6; i++)
                        colors.push_back(col);
                }
            }
            bigMapInterleaved.clear();
            for (size_t i = 0; i < vertices.size(); i++) {
                bigMapInterleaved.push_back(vertices[i].x);
                bigMapInterleaved.push_back(vertices[i].y);
                bigMapInterleaved.push_back(colors[i].r);
                bigMapInterleaved.push_back(colors[i].g);
                bigMapInterleaved.push_back(colors[i].b);
            }
            bigMapDirty = false;
        }
        glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
        int regionChunks = 400;
        int mapWidth = regionChunks * CHUNK_SIZE;
        int mapHeight = regionChunks * CHUNK_SIZE;
        float halfW = mapWidth / 2.0f;
        float halfH = mapHeight / 2.0f;
        glm::mat4 ortho = glm::ortho(bigMapPanX - halfW, bigMapPanX + halfW, bigMapPanZ - halfH, bigMapPanZ + halfH, -1.0f, 1.0f);
        glUseProgram(minimapShaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(minimapShaderProgram, "ortho"), 1, GL_FALSE, glm::value_ptr(ortho));
        glBindBuffer(GL_ARRAY_BUFFER, minimapVBO);
        glBufferData(GL_ARRAY_BUFFER, bigMapInterleaved.size() * sizeof(float), bigMapInterleaved.data(), GL_DYNAMIC_DRAW);
        glBindVertexArray(minimapVAO);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);
        glDrawArrays(GL_TRIANGLES, 0, bigMapInterleaved.size() / 5);
        std::vector<glm::vec2> gridVerts;
        for (int x = static_cast<int>(bigMapPanX - halfW); x <= static_cast<int>(bigMapPanX + halfW); x += CHUNK_SIZE) {
            gridVerts.push_back(glm::vec2(x, bigMapPanZ - halfH));
            gridVerts.push_back(glm::vec2(x, bigMapPanZ + halfH));
        }
        for (int z = static_cast<int>(bigMapPanZ - halfH); z <= static_cast<int>(bigMapPanZ + halfH); z += CHUNK_SIZE) {
            gridVerts.push_back(glm::vec2(bigMapPanX - halfW, z));
            gridVerts.push_back(glm::vec2(bigMapPanX + halfW, z));
        }
        glUniform3f(glGetUniformLocation(minimapShaderProgram, "ourColor"), 0.3f, 0.3f, 0.3f);
        glBufferData(GL_ARRAY_BUFFER, gridVerts.size() * sizeof(glm::vec2), gridVerts.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
        glEnableVertexAttribArray(0);
        glDrawArrays(GL_LINES, 0, gridVerts.size());
        std::vector<glm::vec2> arrowVerts = {
            glm::vec2(0.0f, 8.0f),
            glm::vec2(-4.0f, -4.0f),
            glm::vec2(-4.0f, -4.0f),
            glm::vec2(0.0f, -1.0f),
            glm::vec2(0.0f, -1.0f),
            glm::vec2(4.0f, -4.0f),
            glm::vec2(4.0f, -4.0f),
            glm::vec2(0.0f, 8.0f)
        };
        float arrowAngle = glm::radians(-cameraYaw - 90.0f);
        for (auto& v : arrowVerts) {
            float tx = v.x, ty = v.y;
            v.x = tx * cos(arrowAngle) - ty * sin(arrowAngle);
            v.y = tx * sin(arrowAngle) + ty * cos(arrowAngle);
            v += glm::vec2(cameraPos.x, cameraPos.z);
        }
        glBufferData(GL_ARRAY_BUFFER, arrowVerts.size() * sizeof(glm::vec2), arrowVerts.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
        glEnableVertexAttribArray(0);
        glUniform3f(glGetUniformLocation(minimapShaderProgram, "ourColor"), 1.0f, 0.0f, 0.0f);
        glDrawArrays(GL_LINE_STRIP, 0, arrowVerts.size());
        std::vector<glm::vec2> spawnVerts = {
            glm::vec2(-5.0f, 0.0f), glm::vec2(5.0f, 0.0f),
            glm::vec2(0.0f, -5.0f), glm::vec2(0.0f, 5.0f)
        };
        glBufferData(GL_ARRAY_BUFFER, spawnVerts.size() * sizeof(glm::vec2), spawnVerts.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
        glEnableVertexAttribArray(0);
        glUniform3f(glGetUniformLocation(minimapShaderProgram, "ourColor"), 1.0f, 1.0f, 1.0f);
        glDrawArrays(GL_LINES, 0, spawnVerts.size());
    }
}

// ---------------------- Main Function ----------------------
int main() {
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW\n";
        return -1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWmonitor* primaryMonitor = glfwGetPrimaryMonitor();
    const GLFWvidmode* mode = glfwGetVideoMode(primaryMonitor);
    GLFWwindow* window = glfwCreateWindow(mode->width, mode->height, "Minecraft Clone", primaryMonitor, nullptr);

    if (!window) {
        std::cout << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD\n";
        return -1;
    }
    // Compile shader programs
    GLuint shaderProgram = compileShaderProgram(vertexShaderSource, fragmentShaderSource);
    GLuint minimapShaderProgram = compileShaderProgram(minimapVertexShaderSource, minimapFragmentShaderSource);
    GLuint skyboxShaderProgram = compileShaderProgram(skyboxVertexShaderSource, skyboxFragmentShaderSource);
    sunMoonShaderProgram = compileShaderProgram(sunMoonVertexShaderSource, sunMoonFragmentShaderSource);
    setupSkyboxQuad();
    setupSunMoonQuad();

    // Setup VAOs/VBOs for scene geometry (a lot of them)...
    GLuint VAO, redVAO, waterVAO, grassVAO, treeTrunkVAO, treeLeafVAO, waterLilyVAO, fallenTreeVAO, firLeafVAO;
    GLuint oakTrunkVAO, oakLeafVAO, leafPileVAO, bushSmallVAO, bushMediumVAO, bushLargeVAO;
    GLuint ancientTrunkVAO, ancientLeafVAO;
    GLuint branchVAO, ancientBranchVAO;
    GLuint dirtVAO, deepStoneVAO, lavaVAO, sandVAO, snowVAO, iceVAO;
    GLuint minimapVAO, minimapVBO;
    glGenVertexArrays(1, &VAO);
    glGenVertexArrays(1, &redVAO);
    glGenVertexArrays(1, &waterVAO);
    glGenVertexArrays(1, &grassVAO);
    glGenVertexArrays(1, &treeTrunkVAO);
    glGenVertexArrays(1, &treeLeafVAO);
    glGenVertexArrays(1, &waterLilyVAO);
    glGenVertexArrays(1, &fallenTreeVAO);
    glGenVertexArrays(1, &firLeafVAO);
    glGenVertexArrays(1, &oakTrunkVAO);
    glGenVertexArrays(1, &oakLeafVAO);
    glGenVertexArrays(1, &leafPileVAO);
    glGenVertexArrays(1, &bushSmallVAO);
    glGenVertexArrays(1, &bushMediumVAO);
    glGenVertexArrays(1, &bushLargeVAO);
    glGenVertexArrays(1, &ancientTrunkVAO);
    glGenVertexArrays(1, &ancientLeafVAO);
    glGenVertexArrays(1, &branchVAO);
    glGenVertexArrays(1, &ancientBranchVAO);
    glGenVertexArrays(1, &dirtVAO);
    glGenVertexArrays(1, &deepStoneVAO);
    glGenVertexArrays(1, &lavaVAO);
    glGenVertexArrays(1, &sandVAO);
    glGenVertexArrays(1, &snowVAO);
    glGenVertexArrays(1, &iceVAO);
    glGenVertexArrays(1, &minimapVAO);
    glGenBuffers(1, &minimapVBO);
    GLuint VBO, instanceVBO;
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &instanceVBO);
    GLuint branchInstanceVBO, ancientBranchInstanceVBO;
    glGenBuffers(1, &branchInstanceVBO);
    glGenBuffers(1, &ancientBranchInstanceVBO);

    auto setupVAOFunc = [&](GLuint vao, bool instanced) {
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
        glEnableVertexAttribArray(2);
        if (instanced) {
            glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
            glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
            glEnableVertexAttribArray(3);
            glVertexAttribDivisor(3, 1);
        }
        };
    setupVAOFunc(VAO, false);
    setupVAOFunc(redVAO, false);
    setupVAOFunc(waterVAO, true);
    setupVAOFunc(grassVAO, true);
    setupVAOFunc(treeTrunkVAO, true);
    setupVAOFunc(treeLeafVAO, true);
    setupVAOFunc(waterLilyVAO, true);
    setupVAOFunc(fallenTreeVAO, true);
    setupVAOFunc(firLeafVAO, true);
    setupVAOFunc(oakTrunkVAO, true);
    setupVAOFunc(oakLeafVAO, true);
    setupVAOFunc(leafPileVAO, true);
    setupVAOFunc(bushSmallVAO, true);
    setupVAOFunc(bushMediumVAO, true);
    setupVAOFunc(bushLargeVAO, true);
    setupVAOFunc(ancientTrunkVAO, true);
    setupVAOFunc(ancientLeafVAO, true);
    setupVAOFunc(dirtVAO, true);
    setupVAOFunc(deepStoneVAO, true);
    setupVAOFunc(lavaVAO, true);
    setupVAOFunc(sandVAO, true);
    setupVAOFunc(snowVAO, true);
    setupVAOFunc(iceVAO, true);

    glBindVertexArray(branchVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, branchInstanceVBO);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    glEnableVertexAttribArray(3);
    glVertexAttribDivisor(3, 1);
    glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)(sizeof(glm::vec3)));
    glEnableVertexAttribArray(4);
    glVertexAttribDivisor(4, 1);

    glBindVertexArray(ancientBranchVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, ancientBranchInstanceVBO);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(3);
    glVertexAttribDivisor(3, 1);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // ---------------------- Main Render Loop ----------------------
    while (!glfwWindowShouldClose(window)) {
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        processInput(window);
        toggleMapMode(window);
        updateChunks();
        glClearColor(0.53f, 0.81f, 0.92f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Compute view and projection matrices.
        glm::vec3 front;
        front.x = cos(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        front = glm::normalize(front);
        // Set eye position based on mode: lower if prone or swimming.
        glm::vec3 eyePos;
        if (playerMode == 1 || playerMode == 2)
            eyePos = cameraPos + glm::vec3(0.0f, 0.5f, 0.0f);
        else
            eyePos = cameraPos + glm::vec3(0.0f, eyeLevelOffset, 0.0f);
        glm::mat4 view = glm::lookAt(eyePos, eyePos + front, glm::vec3(0.0f, 1.0f, 0.0f));
        glm::mat4 projection = glm::perspective(glm::radians(103.0f),
            static_cast<float>(WINDOW_WIDTH) / static_cast<float>(WINDOW_HEIGHT),
            0.1f, 10000.0f);

        // Get Eastern local time.
        time_t currentTimeT = time(0);
        tm localTimeInfo;
        localtime_s(&localTimeInfo, &currentTimeT);
        int secondsSinceMidnight = localTimeInfo.tm_hour * 3600 + localTimeInfo.tm_min * 60 + localTimeInfo.tm_sec;
        float dayFraction = secondsSinceMidnight / 86400.0f;
        float angle = dayFraction * 2.0f * 3.14159f;
        float rawHour = localTimeInfo.tm_hour + localTimeInfo.tm_min / 60.0f + localTimeInfo.tm_sec / 3600.0f;
        float hour = fmod(rawHour, 24.0f);

        glm::vec3 sunDir, moonDir;
        float brightnessMain = 0.0f;

        if (hour >= 6.0f && hour < 18.0f) {
            float u = (hour - 6.0f) / 12.0f;
            sunDir = glm::vec3(cos(u * 3.14159f), sin(u * 3.14159f), 0.0f);
            moonDir = -sunDir;
            brightnessMain = sin(u * 3.14159f);
        }
        else {
            float adjustedHour = (hour < 6.0f) ? hour + 24.0f : hour;
            float u = (adjustedHour - 18.0f) / 12.0f;
            moonDir = glm::vec3(cos(u * 3.14159f), sin(u * 3.14159f), 0.0f);
            sunDir = -moonDir;
            brightnessMain = 0.0f;
        }

        glm::vec3 ambientLightMain = glm::vec3(0.2f + brightnessMain * 0.3f);
        glm::vec3 diffuseLightMain = glm::vec3(0.3f + brightnessMain * 0.7f);

        glm::vec3 sunWorldPos = cameraPos + sunDir * 1000.0f;
        glm::vec3 moonWorldPos = cameraPos + moonDir * 1000.0f;

        glm::vec3 skyTop, skyBottom;
        getCurrentSkyColors(dayFraction, skyTop, skyBottom);
        float offsetX = cameraYaw / 360.0f;
        float offsetY = pitch / 180.0f;
        glUniform2f(glGetUniformLocation(skyboxShaderProgram, "starOffset"), offsetX, offsetY);

        // --- Draw Skybox ---
        glDepthMask(GL_FALSE);
        glUseProgram(skyboxShaderProgram);
        glUniform1f(glGetUniformLocation(skyboxShaderProgram, "time"), currentFrame);
        glUniform3fv(glGetUniformLocation(skyboxShaderProgram, "skyTop"), 1, glm::value_ptr(skyTop));
        glUniform3fv(glGetUniformLocation(skyboxShaderProgram, "skyBottom"), 1, glm::value_ptr(skyBottom));
        glBindVertexArray(skyboxVAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
        glDepthMask(GL_TRUE);

        // --- Compute Sun and Moon positions based on real time ---
        glm::vec3 sunDirBillboard(0.0f), moonDirBillboard(0.0f);
        float sunBrightness = 0.0f, moonBrightness = 0.0f;
        if (hour >= 6.0f && hour < 18.0f) {
            float u = (hour - 6.0f) / 12.0f;
            float elev = sin(u * 3.14159f);
            float azimuth = glm::radians(90.0f + u * 180.0f);
            sunDirBillboard.x = cos(glm::radians(elev * 90.0f)) * sin(azimuth);
            sunDirBillboard.y = elev;
            sunDirBillboard.z = cos(glm::radians(elev * 90.0f)) * cos(azimuth);
            sunBrightness = elev;
            moonBrightness = 0.0f;
        }
        else {
            float v = 0.0f;
            if (hour < 6.0f)
                v = hour / 6.0f;
            else
                v = (hour - 18.0f) / 6.0f;
            float elev = sin(v * 3.14159f);
            float azimuth = glm::radians(90.0f + v * 180.0f);
            moonDirBillboard.x = cos(glm::radians(elev * 90.0f)) * sin(azimuth);
            moonDirBillboard.y = elev;
            moonDirBillboard.z = cos(glm::radians(elev * 90.0f)) * cos(azimuth);
            moonBrightness = elev;
            sunBrightness = 0.0f;
        }

        // --- Render Sun/Moon ---
        glDepthMask(GL_FALSE);
        glUseProgram(sunMoonShaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(sunMoonShaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(sunMoonShaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        if (sunBrightness > 0.01f) {
            glm::mat4 model = glm::translate(glm::mat4(1.0f), sunWorldPos);
            glm::mat3 camRotation = glm::mat3(view);
            glm::mat3 billboard = glm::inverse(camRotation);
            model *= glm::mat4(billboard);
            model = glm::scale(model, glm::vec3(50.0f));
            glUniformMatrix4fv(glGetUniformLocation(sunMoonShaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glUniform3f(glGetUniformLocation(sunMoonShaderProgram, "color"), 1.0f, 1.0f, 0.0f);
            glUniform1f(glGetUniformLocation(sunMoonShaderProgram, "brightness"), sunBrightness);
            glBindVertexArray(sunMoonVAO);
            glDrawArrays(GL_TRIANGLES, 0, 6);
        }
        if (moonBrightness > 0.01f) {
            glm::mat4 moonModel = glm::translate(glm::mat4(1.0f), moonWorldPos);
            glm::mat3 camRotation = glm::mat3(view);
            glm::mat3 billboard = glm::inverse(camRotation);
            moonModel *= glm::mat4(billboard);
            moonModel = glm::scale(moonModel, glm::vec3(60.0f));
            glUniformMatrix4fv(glGetUniformLocation(sunMoonShaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(moonModel));
            glUniform3f(glGetUniformLocation(sunMoonShaderProgram, "color"), 0.8f, 0.8f, 1.0f);
            glUniform1f(glGetUniformLocation(sunMoonShaderProgram, "brightness"), moonBrightness * 2.0f);
            glBindVertexArray(sunMoonVAO);
            glDrawArrays(GL_TRIANGLES, 0, 6);
        }
        glDepthMask(GL_TRUE);

        // --- Set up main scene shader ---
        glUseProgram(shaderProgram);
        glUniform1f(glGetUniformLocation(shaderProgram, "time"), currentFrame);
        glUniform3fv(glGetUniformLocation(shaderProgram, "lightDir"), 1, glm::value_ptr(sunDir));
        glUniform3fv(glGetUniformLocation(shaderProgram, "ambientLight"), 1, glm::value_ptr(ambientLightMain));
        glUniform3fv(glGetUniformLocation(shaderProgram, "diffuseLight"), 1, glm::value_ptr(diffuseLightMain));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniform3fv(glGetUniformLocation(shaderProgram, "cameraPos"), 1, glm::value_ptr(cameraPos));

        glm::vec3 blockColors[25];
        blockColors[0] = glm::vec3(0.19f, 0.66f, 0.32f);
        blockColors[1] = glm::vec3(0.0f, 0.5f, 1.0f);
        blockColors[2] = glm::vec3(0.29f, 0.21f, 0.13f);
        blockColors[3] = glm::vec3(0.07f, 0.46f, 0.34f);
        blockColors[4] = glm::vec3(1.0f, 0.0f, 0.0f);
        blockColors[5] = glm::vec3(0.2f, 0.7f, 0.2f);
        blockColors[6] = glm::vec3(0.45f, 0.22f, 0.07f);
        blockColors[7] = glm::vec3(0.13f, 0.54f, 0.13f);
        blockColors[8] = glm::vec3(0.55f, 0.27f, 0.07f);
        blockColors[9] = glm::vec3(0.36f, 0.6f, 0.33f);
        blockColors[10] = glm::vec3(0.44f, 0.39f, 0.32f);
        blockColors[11] = glm::vec3(0.35f, 0.43f, 0.30f);
        blockColors[12] = glm::vec3(0.52f, 0.54f, 0.35f);
        blockColors[13] = glm::vec3(0.6f, 0.61f, 0.35f);
        blockColors[14] = glm::vec3(0.4f, 0.3f, 0.2f);
        blockColors[15] = glm::vec3(0.43f, 0.39f, 0.34f);
        blockColors[16] = glm::vec3(0.4f, 0.25f, 0.1f);
        blockColors[17] = glm::vec3(0.2f, 0.5f, 0.2f);
        blockColors[18] = glm::vec3(0.3f, 0.2f, 0.1f);
        blockColors[19] = glm::vec3(1.0f, 1.0f, 1.0f);
        blockColors[20] = glm::vec3(0.5f, 0.5f, 0.5f);
        blockColors[21] = glm::vec3(1.0f, 0.5f, 0.0f);
        blockColors[22] = glm::vec3(0.93f, 0.79f, 0.69f);
        blockColors[23] = glm::vec3(0.95f, 0.95f, 1.0f);
        blockColors[24] = glm::vec3(0.8f, 0.9f, 1.0f);
        glUniform3fv(glGetUniformLocation(shaderProgram, "blockColors"), 25, glm::value_ptr(blockColors[0]));

        int playerChunkX = static_cast<int>(std::floor(cameraPos.x / CHUNK_SIZE));
        int playerChunkZ = static_cast<int>(std::floor(cameraPos.z / CHUNK_SIZE));
        int qtMinX = playerChunkX - static_cast<int>(RENDER_DISTANCE);
        int qtMaxX = playerChunkX + static_cast<int>(RENDER_DISTANCE);
        int qtMinZ = playerChunkZ - static_cast<int>(RENDER_DISTANCE);
        int qtMaxZ = playerChunkZ + static_cast<int>(RENDER_DISTANCE);
        Quadtree qt(qtMinX, qtMinZ, qtMaxX, qtMaxZ);
        for (auto& entry : chunks) {
            const ChunkPos& pos = entry.first;
            if (pos.x >= qtMinX && pos.x <= qtMaxX && pos.z >= qtMinZ && pos.z <= qtMaxZ)
                qt.insert(pos, &entry.second);
        }
        std::vector<Chunk*> visibleChunks = qt.query(extractFrustumPlanes(projection * view));

        std::vector<glm::vec3> globalGrassInstances;
        std::vector<glm::vec3> globalSandInstances;
        std::vector<glm::vec3> globalSnowInstances;
        std::vector<glm::vec3> globalDirtInstances;
        std::vector<glm::vec3> globalDeepStoneInstances;
        std::vector<glm::vec3> globalWaterInstances;
        std::vector<glm::vec3> globalIceInstances;
        std::vector<glm::vec3> globalLavaInstances;
        std::vector<glm::vec3> globalTreeTrunkInstances;
        std::vector<glm::vec3> globalPineLeafInstances;
        std::vector<glm::vec3> globalFirLeafInstances;
        std::vector<glm::vec3> globalWaterLilyInstances;
        std::vector<glm::vec3> globalFallenTreeTrunkInstances;
        std::vector<glm::vec3> globalOakTrunkInstances;
        std::vector<glm::vec3> globalOakLeafInstances;
        std::vector<glm::vec3> globalLeafPileInstances;
        std::vector<glm::vec3> globalBushSmallInstances;
        std::vector<glm::vec3> globalBushMediumInstances;
        std::vector<glm::vec3> globalBushLargeInstances;
        std::vector<glm::vec3> globalAncientTrunkInstances;
        std::vector<glm::vec3> globalAncientLeafInstances;
        std::vector<glm::vec3> globalAncientBranchInstances;
        std::vector<glm::vec4> globalBranchInstances;
        std::vector<glm::vec3> globalAuroraInstances;
        for (Chunk* chunk : visibleChunks) {
            globalGrassInstances.insert(globalGrassInstances.end(), chunk->grassPositions.begin(), chunk->grassPositions.end());
            globalSandInstances.insert(globalSandInstances.end(), chunk->sandPositions.begin(), chunk->sandPositions.end());
            globalSnowInstances.insert(globalSnowInstances.end(), chunk->snowPositions.begin(), chunk->snowPositions.end());
            globalDirtInstances.insert(globalDirtInstances.end(), chunk->dirtPositions.begin(), chunk->dirtPositions.end());
            globalDeepStoneInstances.insert(globalDeepStoneInstances.end(), chunk->deepStonePositions.begin(), chunk->deepStonePositions.end());
            globalWaterInstances.insert(globalWaterInstances.end(), chunk->waterPositions.begin(), chunk->waterPositions.end());
            globalIceInstances.insert(globalIceInstances.end(), chunk->icePositions.begin(), chunk->icePositions.end());
            globalLavaInstances.insert(globalLavaInstances.end(), chunk->lavaPositions.begin(), chunk->lavaPositions.end());
            globalTreeTrunkInstances.insert(globalTreeTrunkInstances.end(), chunk->treeTrunkPositions.begin(), chunk->treeTrunkPositions.end());
            globalPineLeafInstances.insert(globalPineLeafInstances.end(), chunk->treeLeafPositions.begin(), chunk->treeLeafPositions.end());
            globalFirLeafInstances.insert(globalFirLeafInstances.end(), chunk->firLeafPositions.begin(), chunk->firLeafPositions.end());
            globalWaterLilyInstances.insert(globalWaterLilyInstances.end(), chunk->waterLilyPositions.begin(), chunk->waterLilyPositions.end());
            globalFallenTreeTrunkInstances.insert(globalFallenTreeTrunkInstances.end(), chunk->fallenTreeTrunkPositions.begin(), chunk->fallenTreeTrunkPositions.end());
            globalOakTrunkInstances.insert(globalOakTrunkInstances.end(), chunk->oakTrunkPositions.begin(), chunk->oakTrunkPositions.end());
            globalOakLeafInstances.insert(globalOakLeafInstances.end(), chunk->oakLeafPositions.begin(), chunk->oakLeafPositions.end());
            globalLeafPileInstances.insert(globalLeafPileInstances.end(), chunk->leafPilePositions.begin(), chunk->leafPilePositions.end());
            globalBushSmallInstances.insert(globalBushSmallInstances.end(), chunk->bushSmallPositions.begin(), chunk->bushSmallPositions.end());
            globalBushMediumInstances.insert(globalBushMediumInstances.end(), chunk->bushMediumPositions.begin(), chunk->bushMediumPositions.end());
            globalBushLargeInstances.insert(globalBushLargeInstances.end(), chunk->bushLargePositions.begin(), chunk->bushLargePositions.end());
            globalAncientTrunkInstances.insert(globalAncientTrunkInstances.end(), chunk->ancientTrunkPositions.begin(), chunk->ancientTrunkPositions.end());
            globalAncientLeafInstances.insert(globalAncientLeafInstances.end(), chunk->ancientLeafPositions.begin(), chunk->ancientLeafPositions.end());
            globalAncientBranchInstances.insert(globalAncientBranchInstances.end(), chunk->ancientBranchPositions.begin(), chunk->ancientBranchPositions.end());
            globalBranchInstances.insert(globalBranchInstances.end(), chunk->branchPositions.begin(), chunk->branchPositions.end());
            globalAuroraInstances.insert(globalAuroraInstances.end(), chunk->auroraPositions.begin(), chunk->auroraPositions.end());
        }
        auto drawInstances = [&](GLuint vao, int blockType, const std::vector<glm::vec3>& instances) {
            if (instances.empty()) return;
            glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), blockType);
            glm::mat4 model = glm::mat4(1.0f);
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(vao);
            glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(glm::vec3), instances.data(), GL_DYNAMIC_DRAW);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(instances.size()));
            };
        auto drawBranchInstances = [&](GLuint vao, GLuint branchVBO, int blockType, const std::vector<glm::vec4>& instances) {
            if (instances.empty()) return;
            glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), blockType);
            glm::mat4 model = glm::mat4(1.0f);
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(vao);
            glBindBuffer(GL_ARRAY_BUFFER, branchVBO);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(glm::vec4), instances.data(), GL_DYNAMIC_DRAW);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(instances.size()));
            };
        auto drawAncientBranchInstances = [&](GLuint vao, int blockType, const std::vector<glm::vec3>& instances) {
            if (instances.empty()) return;
            glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), blockType);
            glm::mat4 model = glm::mat4(1.0f);
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(vao);
            glBindBuffer(GL_ARRAY_BUFFER, ancientBranchInstanceVBO);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(glm::vec3), instances.data(), GL_DYNAMIC_DRAW);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(instances.size()));
            };
        drawInstances(grassVAO, 0, globalGrassInstances);
        drawInstances(sandVAO, 22, globalSandInstances);
        drawInstances(snowVAO, 23, globalSnowInstances);
        drawInstances(dirtVAO, 15, globalDirtInstances);
        drawInstances(deepStoneVAO, 20, globalDeepStoneInstances);
        drawInstances(waterVAO, 1, globalWaterInstances);
        drawInstances(iceVAO, 24, globalIceInstances);
        drawInstances(lavaVAO, 21, globalLavaInstances);
        drawInstances(treeTrunkVAO, 2, globalTreeTrunkInstances);
        if (playerChunkZ < 40)
            drawInstances(treeLeafVAO, 3, globalPineLeafInstances);
        drawInstances(firLeafVAO, 7, globalFirLeafInstances);
        drawInstances(waterLilyVAO, 5, globalWaterLilyInstances);
        drawInstances(fallenTreeVAO, 6, globalFallenTreeTrunkInstances);
        drawInstances(oakTrunkVAO, 8, globalOakTrunkInstances);
        drawInstances(oakLeafVAO, 9, globalOakLeafInstances);
        drawInstances(leafPileVAO, 10, globalLeafPileInstances);
        drawInstances(bushSmallVAO, 11, globalBushSmallInstances);
        drawInstances(bushMediumVAO, 12, globalBushMediumInstances);
        drawInstances(bushLargeVAO, 13, globalBushLargeInstances);
        drawInstances(ancientTrunkVAO, 16, globalAncientTrunkInstances);
        drawInstances(ancientLeafVAO, 17, globalAncientLeafInstances);
        drawAncientBranchInstances(ancientBranchVAO, 18, globalAncientBranchInstances);
        drawBranchInstances(branchVAO, branchInstanceVBO, 14, globalBranchInstances);
        drawInstances(waterVAO, 19, globalAuroraInstances);
        glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), 4);
        {
            glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(VAO);
            glDrawArrays(GL_TRIANGLES, 0, 36);
        }
        glm::ivec3 selectedBlock = raycastForBlock(false);
        if (selectedBlock.x != -10000) {
            glm::mat4 outlineModel = glm::translate(glm::mat4(1.0f), glm::vec3(selectedBlock));
            outlineModel = glm::scale(outlineModel, glm::vec3(1.05f));
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(outlineModel));
            glm::vec3 outlineColor = glm::vec3(1.0f, 1.0f, 1.0f);
            glm::vec3 oldBlockColor0 = blockColors[0];
            blockColors[0] = outlineColor;
            glUniform3fv(glGetUniformLocation(shaderProgram, "blockColors"), 25, glm::value_ptr(blockColors[0]));
            glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), 0);
            glDisable(GL_DEPTH_TEST);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glLineWidth(2.0f);
            glBindVertexArray(redVAO);
            glDrawArrays(GL_TRIANGLES, 0, 36);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_DEPTH_TEST);
            blockColors[0] = oldBlockColor0;
            glUniform3fv(glGetUniformLocation(shaderProgram, "blockColors"), 25, glm::value_ptr(blockColors[0]));
        }
        renderMinimap(minimapShaderProgram, minimapVAO, minimapVBO);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glDeleteVertexArrays(1, &VAO);
    glDeleteVertexArrays(1, &redVAO);
    glDeleteVertexArrays(1, &waterVAO);
    glDeleteVertexArrays(1, &grassVAO);
    glDeleteVertexArrays(1, &treeTrunkVAO);
    glDeleteVertexArrays(1, &treeLeafVAO);
    glDeleteVertexArrays(1, &waterLilyVAO);
    glDeleteVertexArrays(1, &fallenTreeVAO);
    glDeleteVertexArrays(1, &firLeafVAO);
    glDeleteVertexArrays(1, &oakTrunkVAO);
    glDeleteVertexArrays(1, &oakLeafVAO);
    glDeleteVertexArrays(1, &leafPileVAO);
    glDeleteVertexArrays(1, &bushSmallVAO);
    glDeleteVertexArrays(1, &bushMediumVAO);
    glDeleteVertexArrays(1, &bushLargeVAO);
    glDeleteVertexArrays(1, &ancientTrunkVAO);
    glDeleteVertexArrays(1, &ancientLeafVAO);
    glDeleteVertexArrays(1, &branchVAO);
    glDeleteVertexArrays(1, &ancientBranchVAO);
    glDeleteVertexArrays(1, &dirtVAO);
    glDeleteVertexArrays(1, &deepStoneVAO);
    glDeleteVertexArrays(1, &lavaVAO);
    glDeleteVertexArrays(1, &sandVAO);
    glDeleteVertexArrays(1, &snowVAO);
    glDeleteVertexArrays(1, &iceVAO);
    glDeleteVertexArrays(1, &minimapVAO);
    glDeleteVertexArrays(1, &skyboxVAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &instanceVBO);
    glDeleteBuffers(1, &branchInstanceVBO);
    glDeleteBuffers(1, &ancientBranchInstanceVBO);
    glDeleteBuffers(1, &minimapVBO);
    glDeleteBuffers(1, &skyboxVBO);
    glDeleteBuffers(1, &sunMoonVBO);
    glDeleteVertexArrays(1, &sunMoonVAO);
    glDeleteProgram(shaderProgram);
    glDeleteProgram(minimapShaderProgram);
    glDeleteProgram(skyboxShaderProgram);
    glDeleteProgram(sunMoonShaderProgram);
    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

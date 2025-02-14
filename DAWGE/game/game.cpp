// ======================================================================
// VoxelGame.cpp
// A single–file version of a Minecraft–style clone with multiple biomes.
//
// Mapping:
//  0  = Grass,
//  1  = Water,
//  2  = Pine/Fir trunk,
//  3  = Pine leaves,
//  4  = Origin debug,
//  5  = Water lily,
//  6  = Fallen log,
//  7  = Fir leaves,
//  8  = Oak trunk,
//  9  = Oak leaves,
// 10  = Leaf pile,
// 11  = Bush (small),
// 12  = Bush (medium),
// 13  = Bush (large),
// 14  = Ground branch,
// 15  = Dirt,
// 16  = Ancient trunk,
// 17  = Ancient leaves,
// 18  = Ancient branch,
// 19  = Aurora block,
// 20  = Deep Stone,
// 21  = Lava,
// 22  = Sand,         // Desert top
// 23  = Snow,         // North pole top
// 24  = Ice           // North pole ice
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
    template<>
    struct hash<ChunkPos> {
        std::size_t operator()(const ChunkPos& cp) const {
            return std::hash<int>()(cp.x) ^ (std::hash<int>()(cp.z) << 1);
        }
    };
}

// ---------------------- Global Constants ----------------------
const unsigned int WINDOW_WIDTH = 1206;
const unsigned int WINDOW_HEIGHT = 832;
const float RENDER_DISTANCE = 5.0f;
const int CHUNK_SIZE = 16;
const int MIN_Y = -1;

// ---------------------- Global Variables ----------------------

// We treat cameraPos as the player's feet.
glm::vec3 cameraPos = glm::vec3(0.0f, 10.0f, 3.0f);

// Free–look variables (mouse controls view)
float cameraYaw = -90.0f;
float pitch = 0.0f;

// For flight mode, the velocity vector
glm::vec3 velocity = glm::vec3(0.0f);

// For walking mode, a jump velocity
float jumpVelocity = 0.0f;

// Toggle between flight (elytra) mode and walking mode.
// Walking mode is enabled by default.
bool walkingMode = true;

// In flight mode, pressing SPACE gives a boost.
const float boostImpulse = 40.0f;

// Constants for flight (elytra) mode
const float baseAcceleration = 20.0f;    // Acceleration when diving (depends on pitch)
const float dragFactor = 0.995f;           // Air drag per frame
const float gravityForce = 9.81f * 0.1f;   // Reduced gravity for flight

// Constants for walking mode
const float walkSpeed = 10.0f;           // Walking speed in blocks per second
const float walkGravity = 9.81f;         // Normal gravity for walking
// Set jump impulse so that maximum jump height is ~1.2 blocks (impulse^2/(2*g) ~ 1.2)
const float walkJumpImpulse = 4.8f;

// Player collision box: 0.6 blocks wide, 2.0 blocks tall.
// Here cameraPos is the player's feet. (Collision uses these offsets.)
const glm::vec3 playerBoxOffsetMin = glm::vec3(-0.3f, 0.0f, -0.3f);
const glm::vec3 playerBoxOffsetMax = glm::vec3(0.3f, 2.0f, 0.3f);

// When rendering, we want the camera at eye-level. Adjust the eye-level offset here.
const float eyeLevelOffset = 1.6f;

// Mouse state for free–look
float lastX = WINDOW_WIDTH / 2.0f;
float lastY = WINDOW_HEIGHT / 2.0f;
bool firstMouse = true;

// Timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

// Fullscreen (big) map mode flag
bool fullscreenMap = false;

// Set of visited chunks (for the fullscreen map)
std::unordered_set<ChunkPos> visitedChunks;

// For fullscreen map caching
bool bigMapDirty = true;
std::vector<float> bigMapInterleaved;

// Big map panning globals
float bigMapPanX = 0.0f;
float bigMapPanZ = 0.0f;

// ---------------------- Forward Declarations ----------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn);
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void processInput(GLFWwindow* window);
glm::ivec3 raycastForBlock(bool place);
void renderMinimap(GLuint minimapShaderProgram, GLuint minimapVAO, GLuint minimapVBO);
void toggleMapMode(GLFWwindow* window);
void handleCollision();

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

// ---------------------- Global Chunk Storage ----------------------
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
std::vector<glm::vec3> generatePineCanopy(int groundHeight, int trunkHeight, int trunkThickness, double worldX, double worldZ) {
    std::vector<glm::vec3> leafPositions;
    int canopyOffset = 50;
    int canopyLayers = 80;
    int canopyBase = groundHeight + trunkHeight - canopyOffset;
    float bottomRadius = 8.0f;
    float topRadius = 2.0f;
    float ringThickness = 1.0f;
    float centerOffset = (trunkThickness - 1) / 2.0f;
    for (int layer = 0; layer < canopyLayers; layer++) {
        float currentRadius = bottomRadius - layer * ((bottomRadius - topRadius) / (canopyLayers - 1));
        int yPos = canopyBase + layer;
        int range = static_cast<int>(std::ceil(currentRadius));
        for (int dx = -range; dx <= range; dx++) {
            for (int dz = -range; dz <= range; dz++) {
                float dist = std::sqrt(dx * dx + dz * dz);
                if (std::abs(dist - currentRadius) < ringThickness) {
                    leafPositions.push_back(glm::vec3(worldX + centerOffset + dx, yPos, worldZ + centerOffset + dz));
                }
            }
        }
    }
    return leafPositions;
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
// Use AABB resolution based on player's collision box. The player's box is:
//    min = cameraPos + (-0.3, 0, -0.3)
//    max = cameraPos + (0.3, 2.0, 0.3)
void handleCollision() {
    glm::vec3 playerMin = cameraPos + glm::vec3(-0.3f, 0.0f, -0.3f);
    glm::vec3 playerMax = cameraPos + glm::vec3(0.3f, 2.0f, 0.3f);

    int chunkX = static_cast<int>(std::floor(cameraPos.x / CHUNK_SIZE));
    int chunkZ = static_cast<int>(std::floor(cameraPos.z / CHUNK_SIZE));
    for (int cx = chunkX - 1; cx <= chunkX + 1; cx++) {
        for (int cz = chunkZ - 1; cz <= chunkZ + 1; cz++) {
            ChunkPos cp(cx, cz);
            if (chunks.find(cp) != chunks.end()) {
                Chunk& ch = chunks[cp];
                auto resolveAABB = [&](const glm::vec3& boxMin, const glm::vec3& boxMax) {
                    if (playerMax.x > boxMin.x && playerMin.x < boxMax.x &&
                        playerMax.y > boxMin.y && playerMin.y < boxMax.y &&
                        playerMax.z > boxMin.z && playerMin.z < boxMax.z) {
                        float overlapX = std::min(playerMax.x - boxMin.x, boxMax.x - playerMin.x);
                        float overlapY = std::min(playerMax.y - boxMin.y, boxMax.y - playerMin.y);
                        float overlapZ = std::min(playerMax.z - boxMin.z, boxMax.z - playerMin.z);
                        if (overlapX < overlapY && overlapX < overlapZ) {
                            if (cameraPos.x < boxMin.x)
                                cameraPos.x -= overlapX;
                            else
                                cameraPos.x += overlapX;
                        }
                        else if (overlapY < overlapX && overlapY < overlapZ) {
                            if (cameraPos.y < boxMin.y)
                                cameraPos.y -= overlapY;
                            else
                                cameraPos.y += overlapY;
                        }
                        else {
                            if (cameraPos.z < boxMin.z)
                                cameraPos.z -= overlapZ;
                            else
                                cameraPos.z += overlapZ;
                        }
                        playerMin = cameraPos + glm::vec3(-0.3f, 0.0f, -0.3f);
                        playerMax = cameraPos + glm::vec3(0.3f, 2.0f, 0.3f);
                    }
                    };

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
            }
        }
    }
    TerrainPoint tp = getTerrainHeight(cameraPos.x, cameraPos.z);
    float terrainY = static_cast<float>(std::floor(tp.height)) + 1.0f;
    if (cameraPos.y < terrainY)
        cameraPos.y = terrainY;
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
                            for (int dx = -7; dx < 7; dx++) {
                                for (int dz = -7; dz < 7; dz++) {
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
                    if (currentChunkZ <= -40) {
                        int hashValPine = std::abs((intWorldX * 73856093) ^ (intWorldZ * 19349663));
                        glm::vec3 pineBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                        if (hashValPine % 2000 < 1 && !treeCollision(chunk.treeTrunkPositions, pineBase)) {
                            int trunkHeight = 60, trunkThickness = 4;
                            for (int i = 1; i <= trunkHeight; i++) {
                                for (int tx = 0; tx < trunkThickness; tx++) {
                                    for (int tz = 0; tz < trunkThickness; tz++) {
                                        chunk.treeTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                    }
                                }
                            }
                            std::vector<glm::vec3> pineCanopy = generatePineCanopy(groundHeight, trunkHeight, trunkThickness, worldX, worldZ);
                            chunk.treeLeafPositions.insert(chunk.treeLeafPositions.end(), pineCanopy.begin(), pineCanopy.end());
                        }
                    }
                    else if ((chunkX < 20) || (currentChunkZ >= 40)) {
                        if (currentChunkZ < 40) {
                            int hashValPine = std::abs((intWorldX * 73856093) ^ (intWorldZ * 19349663));
                            glm::vec3 pineBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                            if (hashValPine % 2000 < 1 && !treeCollision(chunk.treeTrunkPositions, pineBase)) {
                                int trunkHeight = 60, trunkThickness = 4;
                                for (int i = 1; i <= trunkHeight; i++) {
                                    for (int tx = 0; tx < trunkThickness; tx++) {
                                        for (int tz = 0; tz < trunkThickness; tz++) {
                                            chunk.treeTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                        }
                                    }
                                }
                                std::vector<glm::vec3> pineCanopy = generatePineCanopy(groundHeight, trunkHeight, trunkThickness, worldX, worldZ);
                                chunk.treeLeafPositions.insert(chunk.treeLeafPositions.end(), pineCanopy.begin(), pineCanopy.end());
                            }
                        }
                        {
                            int hashValFir = std::abs((intWorldX * 83492791) ^ (intWorldZ * 19349663));
                            glm::vec3 firBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                            if (hashValFir % 2000 < 1 && !treeCollision(chunk.treeTrunkPositions, firBase)) {
                                int trunkHeight = 40, trunkThickness = 3;
                                for (int i = 1; i <= trunkHeight; i++) {
                                    for (int tx = 0; tx < trunkThickness; tx++) {
                                        for (int tz = 0; tz < trunkThickness; tz++) {
                                            chunk.treeTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                        }
                                    }
                                }
                                std::vector<glm::vec3> firCanopy = generateFirCanopy(groundHeight, trunkHeight, trunkThickness, worldX, worldZ);
                                chunk.firLeafPositions.insert(chunk.firLeafPositions.end(), firCanopy.begin(), firCanopy.end());
                            }
                        }
                        {
                            int hashValOak = std::abs((intWorldX * 92821) ^ (intWorldZ * 123457));
                            glm::vec3 oakBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                            if (hashValOak % 1000 < 1 && !treeCollision(chunk.oakTrunkPositions, oakBase)) {
                                int trunkHeight = 7, trunkThickness = 2;
                                for (int i = 1; i <= trunkHeight; i++) {
                                    for (int tx = 0; tx < trunkThickness; tx++) {
                                        for (int tz = 0; tz < trunkThickness; tz++) {
                                            chunk.oakTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                        }
                                    }
                                }
                                std::vector<glm::vec3> oakCanopy = generateOakCanopy(groundHeight, trunkHeight, trunkThickness, worldX, worldZ);
                                chunk.oakLeafPositions.insert(chunk.oakLeafPositions.end(), oakCanopy.begin(), oakCanopy.end());
                            }
                        }
                        {
                            int hashValAncient = std::abs((intWorldX * 112233) ^ (intWorldZ * 445566));
                            glm::vec3 ancientBase = glm::vec3(worldX, groundHeight + 1, worldZ);
                            if (hashValAncient % 3000 < 1 && !treeCollision(chunk.ancientTrunkPositions, ancientBase)) {
                                int trunkHeight = 30, trunkThickness = 3;
                                for (int i = 1; i <= trunkHeight; i++) {
                                    for (int tx = 0; tx < trunkThickness; tx++) {
                                        for (int tz = 0; tz < trunkThickness; tz++) {
                                            chunk.ancientTrunkPositions.push_back(glm::vec3(worldX + tx, groundHeight + i, worldZ + tz));
                                        }
                                    }
                                }
                                int centerY = groundHeight + trunkHeight;
                                float canopyRadius = 5.0f;
                                for (int dy = -static_cast<int>(canopyRadius); dy <= static_cast<int>(canopyRadius); dy++) {
                                    for (int dx = -static_cast<int>(canopyRadius); dx <= static_cast<int>(canopyRadius); dx++) {
                                        for (int dz = -static_cast<int>(canopyRadius); dz <= static_cast<int>(canopyRadius); dz++) {
                                            float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                                            if (dist < canopyRadius)
                                                chunk.ancientLeafPositions.push_back(glm::vec3(worldX + trunkThickness / 2.0f + dx, centerY + dy, worldZ + trunkThickness / 2.0f + dz));
                                        }
                                    }
                                }
                                int branchBaseHeights[4] = { 7, 13, 19, 25 };
                                for (int b = 0; b < 4; b++) {
                                    int randomOffset = (rand() % 3) - 1;
                                    int branchStart = branchBaseHeights[b] + randomOffset;
                                    float branchRot = (b * 90.0f) * (3.14159f / 180.0f);
                                    glm::vec3 branchStartPos = glm::vec3(worldX + trunkThickness / 2.0f, groundHeight + branchStart, worldZ + trunkThickness / 2.0f);
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
// Supports two modes: Flight (elytra) and Walking.
// Press 'P' to toggle modes. Walking is the default.
// In Walking mode, WASD moves relative to the horizontal view,
// SPACE causes a jump with an impulse for ~1.2 block height and normal gravity.
// In Flight mode, velocity is steered toward the view with reduced gravity and drag.
// Collision handling is applied in both modes.
void processInput(GLFWwindow* window) {
    static bool pWasPressed = false;
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
        if (!pWasPressed) {
            walkingMode = !walkingMode;
            pWasPressed = true;
            velocity = glm::vec3(0.0f);
            jumpVelocity = 0.0f;
        }
    }
    else {
        pWasPressed = false;
    }

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

    if (!walkingMode) {
        // Flight mode:
        glm::vec3 viewDir;
        viewDir.x = cos(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        viewDir.y = sin(glm::radians(pitch));
        viewDir.z = sin(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        viewDir = glm::normalize(viewDir);
        float steeringFactor = 0.02f;
        velocity = glm::mix(velocity, viewDir * glm::length(velocity), steeringFactor);
        float pitchFactor = -pitch / 90.0f;
        float accel = baseAcceleration * pitchFactor;
        velocity += viewDir * accel * deltaTime;
        velocity.y -= gravityForce * deltaTime;
        velocity *= dragFactor;
        cameraPos += velocity * deltaTime;
        static bool spaceWasPressed = false;
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
            if (!spaceWasPressed) {
                velocity += viewDir * boostImpulse;
                spaceWasPressed = true;
            }
        }
        else {
            spaceWasPressed = false;
        }
    }
    else {
        // Walking mode:
        glm::vec3 forward;
        forward.x = cos(glm::radians(cameraYaw));
        forward.y = 0.0f;
        forward.z = sin(glm::radians(cameraYaw));
        forward = glm::normalize(forward);
        glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
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
        cameraPos += walkDir * walkSpeed * deltaTime;
        TerrainPoint tp = getTerrainHeight(cameraPos.x, cameraPos.z);
        float groundY = static_cast<float>(std::floor(tp.height)) + 1.0f;
        bool onGround = (cameraPos.y <= groundY + 0.1f);
        if (onGround) {
            jumpVelocity = 0.0f;
            cameraPos.y = groundY;
            if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
                jumpVelocity = walkJumpImpulse;
        }
        jumpVelocity -= walkGravity * deltaTime;
        cameraPos.y += jumpVelocity * deltaTime;
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

// ---------------------- Shader Sources ----------------------
const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexCoord;
layout (location = 2) in vec3 aOffset;
layout (location = 3) in float aRotation;
out vec2 TexCoord;
out vec3 ourColor;
out float instanceDistance;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform int blockType;
uniform vec3 blockColors[25];
uniform vec3 cameraPos;
uniform float time;
void main(){
    vec3 pos = aPos;
    if(blockType != 14 && blockType != 18)
        pos += aOffset;
    else {
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
            pos = rot * (scaleMat * pos) + aOffset;
        } else {
            pos += aOffset;
        }
    }
    if(blockType == 19){
        pos.y += sin(time + aOffset.x * 0.1) * 0.5;
    }
    gl_Position = projection * view * model * vec4(pos, 1.0);
    ourColor = blockColors[blockType];
    TexCoord = aTexCoord;
    if(blockType != 14 && blockType != 18){
        if(gl_InstanceID > 0)
            instanceDistance = length(aOffset - cameraPos);
        else
            instanceDistance = length(vec3(model[3]) - cameraPos);
    } else {
        instanceDistance = length(aOffset - cameraPos);
    }
}
)";

const char* fragmentShaderSource = R"(
#version 330 core
in vec2 TexCoord;
in vec3 ourColor;
in float instanceDistance;
out vec4 FragColor;
uniform int blockType;
uniform vec3 blockColors[25];
void main(){
    if(blockType == 19){
        FragColor = vec4(ourColor, 0.1);
        return;
    }
    float gridSize = 24.0;
    float lineWidth = 0.03;
    vec2 f = fract(TexCoord * gridSize);
    if(f.x < lineWidth || f.y < lineWidth)
        FragColor = vec4(0.0, 0.0, 0.0, 1.0);
    else {
        float factor = instanceDistance / 100.0;
        vec3 offset = vec3(0.03 * factor, 0.03 * factor, 0.05 * factor);
        vec3 finalColor = ourColor + offset;
        finalColor = clamp(finalColor, 0.0, 1.0);
        FragColor = vec4(finalColor, 1.0);
    }
}
)";

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

const char* minimapFragmentShaderSource = R"(
#version 330 core
in vec3 ourColor;
out vec4 FragColor;
void main(){
    FragColor = vec4(ourColor, 1.0);
}
)";

// ---------------------- Cube Vertex Data ----------------------
float cubeVertices[] = {
    // positions          // texture Coords
    // Front face
   -0.5f, -0.5f,  0.5f,   0.0f, 0.0f,
    0.5f, -0.5f,  0.5f,   1.0f, 0.0f,
    0.5f,  0.5f,  0.5f,   1.0f, 1.0f,
    0.5f,  0.5f,  0.5f,   1.0f, 1.0f,
   -0.5f,  0.5f,  0.5f,   0.0f, 1.0f,
   -0.5f, -0.5f,  0.5f,   0.0f, 0.0f,
   // Right face
   0.5f, -0.5f,  0.5f,   0.0f, 0.0f,
   0.5f, -0.5f, -0.5f,   1.0f, 0.0f,
   0.5f,  0.5f, -0.5f,   1.0f, 1.0f,
   0.5f,  0.5f, -0.5f,   1.0f, 1.0f,
   0.5f,  0.5f,  0.5f,   0.0f, 1.0f,
   0.5f, -0.5f,  0.5f,   0.0f, 0.0f,
   // Back face
   0.5f, -0.5f, -0.5f,   0.0f, 0.0f,
  -0.5f, -0.5f, -0.5f,   1.0f, 0.0f,
  -0.5f,  0.5f, -0.5f,   1.0f, 1.0f,
  -0.5f,  0.5f, -0.5f,   1.0f, 1.0f,
   0.5f,  0.5f, -0.5f,   0.0f, 1.0f,
   0.5f, -0.5f, -0.5f,   0.0f, 0.0f,
   // Left face
  -0.5f, -0.5f, -0.5f,   0.0f, 0.0f,
  -0.5f, -0.5f,  0.5f,   1.0f, 0.0f,
  -0.5f,  0.5f,  0.5f,   1.0f, 1.0f,
  -0.5f,  0.5f,  0.5f,   1.0f, 1.0f,
  -0.5f,  0.5f, -0.5f,   0.0f, 1.0f,
  -0.5f, -0.5f, -0.5f,   0.0f, 0.0f,
  // Top face
-0.5f,  0.5f,  0.5f,   0.0f, 0.0f,
 0.5f,  0.5f,  0.5f,   1.0f, 0.0f,
 0.5f,  0.5f, -0.5f,   1.0f, 1.0f,
 0.5f,  0.5f, -0.5f,   1.0f, 1.0f,
-0.5f,  0.5f, -0.5f,   0.0f, 1.0f,
-0.5f,  0.5f,  0.5f,   0.0f, 0.0f,
// Bottom face
-0.5f, -0.5f, -0.5f,   0.0f, 0.0f,
 0.5f, -0.5f, -0.5f,   1.0f, 0.0f,
 0.5f, -0.5f,  0.5f,   1.0f, 1.0f,
 0.5f, -0.5f,  0.5f,   1.0f, 1.0f,
-0.5f, -0.5f,  0.5f,   0.0f, 1.0f,
-0.5f, -0.5f, -0.5f,   0.0f, 0.0f
};

// ---------------------- Minimap Rendering ----------------------
void renderMinimap(GLuint minimapShaderProgram, GLuint minimapVAO, GLuint minimapVBO) {
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
                    vertices.push_back(glm::vec2(x + 1, z));
                    vertices.push_back(glm::vec2(x + 1, z + 1));
                    vertices.push_back(glm::vec2(x, z));
                    vertices.push_back(glm::vec2(x + 1, z + 1));
                    vertices.push_back(glm::vec2(x, z + 1));
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
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD\n";
        return -1;
    }

    // ---------------------- Shader Compilation ----------------------
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    {
        int success;
        glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
            std::cout << "Vertex Shader Compilation Error:\n" << infoLog << "\n";
        }
    }
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    {
        int success;
        glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
            std::cout << "Fragment Shader Compilation Error:\n" << infoLog << "\n";
        }
    }
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    {
        int success;
        glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
            std::cout << "Shader Program Linking Error:\n" << infoLog << "\n";
        }
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // ---------------------- Minimap Shader Compilation ----------------------
    unsigned int minimapVertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(minimapVertexShader, 1, &minimapVertexShaderSource, NULL);
    glCompileShader(minimapVertexShader);
    {
        int success;
        glGetShaderiv(minimapVertexShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetShaderInfoLog(minimapVertexShader, 512, NULL, infoLog);
            std::cout << "Minimap Vertex Shader Compilation Error:\n" << infoLog << "\n";
        }
    }
    unsigned int minimapFragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(minimapFragmentShader, 1, &minimapFragmentShaderSource, NULL);
    glCompileShader(minimapFragmentShader);
    {
        int success;
        glGetShaderiv(minimapFragmentShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetShaderInfoLog(minimapFragmentShader, 512, NULL, infoLog);
            std::cout << "Minimap Fragment Shader Compilation Error:\n" << infoLog << "\n";
        }
    }
    unsigned int minimapShaderProgram = glCreateProgram();
    glAttachShader(minimapShaderProgram, minimapVertexShader);
    glAttachShader(minimapShaderProgram, minimapFragmentShader);
    glLinkProgram(minimapShaderProgram);
    {
        int success;
        glGetProgramiv(minimapShaderProgram, GL_LINK_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetProgramInfoLog(minimapShaderProgram, 512, NULL, infoLog);
            std::cout << "Minimap Shader Program Linking Error:\n" << infoLog << "\n";
        }
    }
    glDeleteShader(minimapVertexShader);
    glDeleteShader(minimapFragmentShader);

    // ---------------------- Setup VAOs and VBOs ----------------------
    unsigned int VAO, redVAO, waterVAO, grassVAO, treeTrunkVAO, treeLeafVAO, waterLilyVAO, fallenTreeVAO, firLeafVAO;
    unsigned int oakTrunkVAO, oakLeafVAO, leafPileVAO, bushSmallVAO, bushMediumVAO, bushLargeVAO;
    unsigned int ancientTrunkVAO, ancientLeafVAO;
    unsigned int branchVAO, ancientBranchVAO;
    unsigned int dirtVAO, deepStoneVAO, lavaVAO, sandVAO, snowVAO, iceVAO;
    unsigned int minimapVAO, minimapVBO;
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
    unsigned int VBO, instanceVBO;
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &instanceVBO);
    GLuint branchInstanceVBO, ancientBranchInstanceVBO;
    glGenBuffers(1, &branchInstanceVBO);
    glGenBuffers(1, &ancientBranchInstanceVBO);

    auto setupVAOFunc = [&](unsigned int vao, bool instanced) {
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);
        if (instanced) {
            glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
            glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
            glEnableVertexAttribArray(2);
            glVertexAttribDivisor(2, 1);
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
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, branchInstanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)(sizeof(glm::vec3)));
    glEnableVertexAttribArray(3);
    glVertexAttribDivisor(3, 1);

    glBindVertexArray(ancientBranchVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, ancientBranchInstanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

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
        glUseProgram(shaderProgram);
        glUniform1f(glGetUniformLocation(shaderProgram, "time"), currentFrame);

        // --- Compute view matrix with eye-level offset ---
        glm::vec3 front;
        front.x = cos(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(cameraYaw)) * cos(glm::radians(pitch));
        front = glm::normalize(front);
        // Here we add the eye level offset (e.g. 1.6 blocks above cameraPos)
        glm::vec3 eyePos = cameraPos + glm::vec3(0.0f, eyeLevelOffset, 0.0f);
        glm::mat4 view = glm::lookAt(eyePos, eyePos + front, glm::vec3(0.0f, 1.0f, 0.0f));

        glm::mat4 projection = glm::perspective(glm::radians(45.0f),
            static_cast<float>(WINDOW_WIDTH) / static_cast<float>(WINDOW_HEIGHT),
            0.1f, 10000.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniform3fv(glGetUniformLocation(shaderProgram, "cameraPos"), 1, glm::value_ptr(cameraPos));

        glm::vec3 blockColors[25];
        blockColors[0] = glm::vec3(0.19f, 0.66f, 0.32f);
        blockColors[1] = glm::vec3(0.0f, 0.5f, 0.5f);
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

        auto drawInstances = [&](unsigned int vao, int blockType, const std::vector<glm::vec3>& instances) {
            if (instances.empty()) return;
            glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), blockType);
            glm::mat4 model = glm::mat4(1.0f);
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(vao);
            glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(glm::vec3), instances.data(), GL_DYNAMIC_DRAW);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(instances.size()));
            };
        auto drawBranchInstances = [&](unsigned int vao, unsigned int branchVBO, int blockType, const std::vector<glm::vec4>& instances) {
            if (instances.empty()) return;
            glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), blockType);
            glm::mat4 model = glm::mat4(1.0f);
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(vao);
            glBindBuffer(GL_ARRAY_BUFFER, branchVBO);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(glm::vec4), instances.data(), GL_DYNAMIC_DRAW);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(instances.size()));
            };
        auto drawAncientBranchInstances = [&](unsigned int vao, int blockType, const std::vector<glm::vec3>& instances) {
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

        // Render player selection outline.
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
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &instanceVBO);
    glDeleteBuffers(1, &branchInstanceVBO);
    glDeleteBuffers(1, &ancientBranchInstanceVBO);
    glDeleteBuffers(1, &minimapVBO);
    glDeleteProgram(shaderProgram);
    glDeleteProgram(minimapShaderProgram);
    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

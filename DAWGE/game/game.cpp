// ======================================================================
// VoxelGame.cpp
// A single-file version with several features:
//   1) Block outline (mining) feature with depth-testing disabled for outline.
//   2) Larger epsilon tolerance for mining tree logs and leaves.
//   3) Per-face UV coordinates and a fragment shader that “pixelates” each face
//      into a 24×24 grid with grid lines.
//   4) A third type of tree added (Oak), distinct from Pine and Fir.
//   5) Procedurally generated leaf piles on the ground.
//   6) Procedurally generated bushes (small, medium, large) with larger, distinct clusters.
//   7) Independent color controls for all new items (oak trunks, oak leaves,
//      leaf piles, and the 3 bush sizes).
//   8) Atmospheric / aerial perspective: blocks far from the player appear lighter
//      and bluer. Now the effect blends smoothly with a gradient. For every 100
//      blocks, red and green are increased by 0.05 and blue by 0.07, blended
//      continuously with distance.
//   9) NEW: Ground branches – thin procedural branches on the ground.
//  10) NEW: Ancient tree (fourth tree type) with trunk, canopy, and procedural
//      branches. Ancient tree branches are now full blocks. Four branches are
//      generated at fixed relative heights ("low", "low-mid", "high-mid", "high")
//      with a small random fluctuation, extended longer, and with a mini canopy
//      of leaves at each tip.
//  11) NEW: A simple collision check for trees so that trees of the same type
//      do not spawn too close to each other.
// ======================================================================

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

// ==================== Global Configuration ====================
const unsigned int WINDOW_WIDTH = 1206;
const unsigned int WINDOW_HEIGHT = 832;
const float RENDER_DISTANCE = 64.0f; // For testing.
const int CHUNK_SIZE = 16;

// ==================== Global Camera Variables ====================
glm::vec3 cameraPos = glm::vec3(0.0f, 10.0f, 3.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

float yaw = -90.0f;
float pitch = 0.0f;
float lastX = WINDOW_WIDTH / 2.0f;
float lastY = WINDOW_HEIGHT / 2.0f;
bool  firstMouse = true;

float deltaTime = 0.0f;
float lastFrame = 0.0f;

// ==================== Forward Declarations ====================
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn);
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void processInput(GLFWwindow* window);
glm::ivec3 raycastForBlock(bool place);

// ==================== Chunk and Block Modifications ====================

struct ChunkPos {
    int x, z;
    bool operator==(const ChunkPos& other) const { return x == other.x && z == other.z; }
};

namespace std {
    template<>
    struct hash<ChunkPos> {
        size_t operator()(const ChunkPos& k) const {
            return hash<int>()(k.x) ^ (hash<int>()(k.z) << 1);
        }
    };
}

// The Chunk now also stores vectors for ground branches and ancient trees.
struct Chunk {
    std::vector<glm::vec3> waterPositions;
    std::vector<glm::vec3> stonePositions;
    std::vector<glm::vec3> treeTrunkPositions;     // Pine/Fir trunk
    std::vector<glm::vec3> treeLeafPositions;      // Pine leaves
    std::vector<glm::vec3> firLeafPositions;       // Fir leaves
    std::vector<glm::vec3> waterLilyPositions;
    std::vector<glm::vec3> fallenTreeTrunkPositions;

    // Oak tree arrays:
    std::vector<glm::vec3> oakTrunkPositions;
    std::vector<glm::vec3> oakLeafPositions;

    // Leaf piles and bushes:
    std::vector<glm::vec3> leafPilePositions;
    std::vector<glm::vec3> bushSmallPositions;
    std::vector<glm::vec3> bushMediumPositions;
    std::vector<glm::vec3> bushLargePositions;

    // NEW: Ground branches (drawn as thin blocks)
    // Stored as vec4: xyz = position, w = rotation.
    std::vector<glm::vec4> branchPositions;

    // NEW: Ancient tree arrays (fourth tree type)
    std::vector<glm::vec3> ancientTrunkPositions; // Block type 16
    std::vector<glm::vec3> ancientLeafPositions;  // Block type 17
    // Ancient branches are now full blocks.
    std::vector<glm::vec3> ancientBranchPositions; // Block type 18

    bool needsMeshUpdate;
    Chunk() : needsMeshUpdate(true) {}
};

struct ivec3_hash {
    std::size_t operator()(const glm::ivec3& v) const {
        return ((std::hash<int>()(v.x) ^ (std::hash<int>()(v.y) << 1)) >> 1)
            ^ (std::hash<int>()(v.z) << 1);
    }
};

// Global block modifications: key is a block position, value indicates modification:
//   -1 = removal; nonnegative values indicate block type:
//   0 = Stone, 1 = Water, 2 = Pine/Fir trunk, 3 = Pine leaves, 4 = Origin debug,
//   5 = Water lily, 6 = Fallen log, 7 = Fir leaves, 8 = Oak trunk, 9 = Oak leaves,
//  10 = Leaf pile, 11 = Bush (small), 12 = Bush (medium), 13 = Bush (large),
//  14 = Ground branch, 16 = Ancient trunk, 17 = Ancient leaves, 18 = Ancient branch.
std::unordered_map<glm::ivec3, int, ivec3_hash> blockModifications;

// Global chunk storage.
std::unordered_map<ChunkPos, Chunk> chunks;

// ==================== Utility: Perlin Noise ====================
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
        return lerp(w,
            lerp(v, lerp(u, grad(p[AA], x, y, z),
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
        return n * n;
    }
};

PerlinNoise continentalNoise(1);
PerlinNoise elevationNoise(2);
PerlinNoise ridgeNoise(3);

// ==================== Terrain Generation ====================
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
    double ridge = ridgeNoise.ridgeNoise(x / RIDGE_SCALE, 0, z / RIDGE_SCALE);
    double height = elevation * 8.0 + ridge * 12.0;
    return { height, true };
}

// ==================== Frustum Culling ====================
struct Plane { glm::vec3 normal; float d; };

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

// ==================== Helper: Tree Collision Check ====================
// Returns true if any block in trunkArray is within 3.0 units of base.
bool treeCollision(const std::vector<glm::vec3>& trunkArray, const glm::vec3& base) {
    for (const auto& p : trunkArray) {
        if (glm::distance(p, base) < 3.0f)
            return true;
    }
    return false;
}

// ==================== Tree Generation ====================
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

std::vector<glm::vec3> generateFirCanopy(int groundHeight, int trunkHeight, int trunkThickness, double worldX, double worldZ) {
    std::vector<glm::vec3> leafPositions;
    int centerY = groundHeight + trunkHeight;
    float radius = 7.0f;
    for (int dy = -static_cast<int>(radius); dy <= static_cast<int>(radius); dy++) {
        for (int dx = -static_cast<int>(radius); dx <= static_cast<int>(radius); dx++) {
            for (int dz = -static_cast<int>(radius); dz <= static_cast<int>(radius); dz++) {
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (dist < radius) {
                    leafPositions.push_back(glm::vec3(
                        worldX + trunkThickness / 2.0f + dx,
                        centerY + dy,
                        worldZ + trunkThickness / 2.0f + dz
                    ));
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
                    leaves.push_back(glm::vec3(
                        worldX + centerOffset + dx,
                        centerY + dy,
                        worldZ + centerOffset + dz
                    ));
                }
            }
        }
    }
    return leaves;
}

// ==================== Raycasting ====================
glm::ivec3 raycastForBlock(bool place) {
    float t = 0.0f;
    float stoneTol = 0.5f, trunkTol = 0.5f, leafTol = 0.6f, waterTol = 0.5f;
    float lilyTol = 0.5f, fallenTol = 0.5f, oakTrunkTol = 0.5f, oakLeafTol = 0.6f;
    float leafPileTol = 0.5f, bushTol = 0.5f;
    while (t < 5.0f) {
        glm::vec3 p = cameraPos + t * cameraFront;
        glm::ivec3 candidate = glm::ivec3(std::round(p.x), std::round(p.y), std::round(p.z));
        if (blockModifications.find(candidate) != blockModifications.end() && blockModifications[candidate] == -1) {
            t += 0.1f;
            continue;
        }
        bool exists = false;
        int chunkX = static_cast<int>(std::floor(candidate.x / static_cast<float>(CHUNK_SIZE)));
        int chunkZ = static_cast<int>(std::floor(candidate.z / static_cast<float>(CHUNK_SIZE)));
        ChunkPos cp{ chunkX, chunkZ };
        if (chunks.find(cp) != chunks.end()) {
            Chunk& ch = chunks[cp];
            for (const auto& pos : ch.stonePositions)
                if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), stoneTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.treeTrunkPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), trunkTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.treeLeafPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), leafTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.firLeafPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), leafTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.waterPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), waterTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.waterLilyPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), lilyTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.fallenTreeTrunkPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), fallenTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.oakTrunkPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), oakTrunkTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.oakLeafPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), oakLeafTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.leafPilePositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), leafPileTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.bushSmallPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), bushTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.bushMediumPositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), bushTol))) { exists = true; break; }
            if (!exists)
                for (const auto& pos : ch.bushLargePositions)
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), bushTol))) { exists = true; break; }
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

// ==================== Mouse Button Callback ====================
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            glm::ivec3 pos = raycastForBlock(false);
            if (pos.x != -10000) {
                blockModifications[pos] = -1;
                int chunkX = static_cast<int>(std::floor(pos.x / static_cast<float>(CHUNK_SIZE)));
                int chunkZ = static_cast<int>(std::floor(pos.z / static_cast<float>(CHUNK_SIZE)));
                ChunkPos cp{ chunkX, chunkZ };
                if (chunks.find(cp) != chunks.end())
                    chunks[cp].needsMeshUpdate = true;
            }
        }
        else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            glm::ivec3 hitBlock = raycastForBlock(false);
            if (hitBlock.x != -10000) {
                int chunkX = static_cast<int>(std::floor(hitBlock.x / static_cast<float>(CHUNK_SIZE)));
                int chunkZ = static_cast<int>(std::floor(hitBlock.z / static_cast<float>(CHUNK_SIZE)));
                ChunkPos cp{ chunkX, chunkZ };
                int blockTypeToPlace = 0;
                if (chunks.find(cp) != chunks.end()) {
                    Chunk& ch = chunks[cp];
                    auto within = [&](const std::vector<glm::vec3>& vec, float tol) -> bool {
                        for (const auto& v : vec) {
                            if (glm::all(glm::epsilonEqual(v, glm::vec3(hitBlock), tol)))
                                return true;
                        }
                        return false;
                        };
                    if (within(ch.treeTrunkPositions, 0.5f))
                        blockTypeToPlace = 2;
                    else if (within(ch.treeLeafPositions, 0.6f))
                        blockTypeToPlace = 3;
                    else if (within(ch.firLeafPositions, 0.6f))
                        blockTypeToPlace = 7;
                    else if (within(ch.waterPositions, 0.5f))
                        blockTypeToPlace = 1;
                    else if (within(ch.waterLilyPositions, 0.5f))
                        blockTypeToPlace = 5;
                    else if (within(ch.fallenTreeTrunkPositions, 0.5f))
                        blockTypeToPlace = 6;
                    else if (within(ch.oakTrunkPositions, 0.5f))
                        blockTypeToPlace = 8;
                    else if (within(ch.oakLeafPositions, 0.6f))
                        blockTypeToPlace = 9;
                    else if (within(ch.leafPilePositions, 0.5f))
                        blockTypeToPlace = 10;
                    else if (within(ch.bushSmallPositions, 0.5f))
                        blockTypeToPlace = 11;
                    else if (within(ch.bushMediumPositions, 0.5f))
                        blockTypeToPlace = 12;
                    else if (within(ch.bushLargePositions, 0.5f))
                        blockTypeToPlace = 13;
                    else
                        blockTypeToPlace = 0;
                }
                glm::ivec3 pos = raycastForBlock(true);
                if (pos.x != -10000) {
                    blockModifications[pos] = blockTypeToPlace;
                    int placeChunkX = static_cast<int>(std::floor(pos.x / static_cast<float>(CHUNK_SIZE)));
                    int placeChunkZ = static_cast<int>(std::floor(pos.z / static_cast<float>(CHUNK_SIZE)));
                    ChunkPos cp2{ placeChunkX, placeChunkZ };
                    if (chunks.find(cp2) != chunks.end())
                        chunks[cp2].needsMeshUpdate = true;
                }
            }
        }
    }
}

// ==================== Quadtree for Spatial Partitioning ====================
struct QuadtreeItem { ChunkPos pos; Chunk* chunk; };

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
    ~QuadtreeNode() { for (int i = 0; i < 4; i++) if (children[i]) delete children[i]; }
    glm::vec3 getMinWorld() const { return glm::vec3(minX * CHUNK_SIZE, -4.0f, minZ * CHUNK_SIZE); }
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
                if (children[i]->contains(item.pos)) { children[i]->items.push_back(item); break; }
            }
        }
        items.clear();
    }
    void insert(const QuadtreeItem& item) {
        if (!contains(item.pos)) return;
        if (!subdivided && items.size() < capacity)
            items.push_back(item);
        else {
            if (!subdivided) subdivide();
            for (int i = 0; i < 4; i++) {
                if (children[i]->contains(item.pos)) { children[i]->insert(item); return; }
            }
        }
    }
    void query(const std::vector<Plane>& frustum, std::vector<Chunk*>& out) {
        glm::vec3 nodeMin = getMinWorld(), nodeMax = getMaxWorld();
        if (!aabbInFrustum(frustum, nodeMin, nodeMax)) return;
        if (!subdivided) {
            for (const auto& item : items) out.push_back(item.chunk);
        }
        else {
            for (int i = 0; i < 4; i++) children[i]->query(frustum, out);
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

// ==================== Chunk Mesh Generation ====================
void generateChunkMesh(Chunk& chunk, int chunkX, int chunkZ) {
    if (!chunk.needsMeshUpdate)
        return;
    // Clear old data
    chunk.waterPositions.clear();
    chunk.stonePositions.clear();
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
            }
            else {
                int groundHeight = static_cast<int>(std::floor(terrain.height));
                for (int y = groundHeight; y >= -4; y--) {
                    chunk.stonePositions.push_back(glm::vec3(worldX, y, worldZ));
                }
                if (terrain.height > 2.0) {
                    int intWorldX = static_cast<int>(worldX);
                    int intWorldZ = static_cast<int>(worldZ);

                    // --- Pine Tree Generation ---
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

                    // --- Fir Tree Generation ---
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

                    // --- Oak Tree Generation ---
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

                    // --- Ancient Tree Generation ---
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
                        // Generate ancient canopy (spherical pattern of leaves)
                        int centerY = groundHeight + trunkHeight;
                        float canopyRadius = 5.0f;
                        for (int dy = -static_cast<int>(canopyRadius); dy <= static_cast<int>(canopyRadius); dy++) {
                            for (int dx = -static_cast<int>(canopyRadius); dx <= static_cast<int>(canopyRadius); dx++) {
                                for (int dz = -static_cast<int>(canopyRadius); dz <= static_cast<int>(canopyRadius); dz++) {
                                    float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                                    if (dist < canopyRadius) {
                                        chunk.ancientLeafPositions.push_back(glm::vec3(worldX + trunkThickness / 2.0f + dx, centerY + dy, worldZ + trunkThickness / 2.0f + dz));
                                    }
                                }
                            }
                        }
                        // Generate FOUR ancient branches at fixed relative heights:
                        // low, low-mid, high-mid, high.
                        int branchBaseHeights[4] = { 7, 13, 19, 25 }; // relative to trunk (for trunkHeight=30)
                        for (int b = 0; b < 4; b++) {
                            // Add a small random fluctuation of -1, 0, or +1.
                            int randomOffset = (rand() % 3) - 1;
                            int branchStart = branchBaseHeights[b] + randomOffset;
                            float branchRot = (b * 90.0f) * (3.14159f / 180.0f); // fixed cardinal directions.
                            glm::vec3 branchStartPos = glm::vec3(worldX + trunkThickness / 2.0f, groundHeight + branchStart, worldZ + trunkThickness / 2.0f);
                            int branchLength = 10 + (rand() % 3); // 10-12 blocks long.
                            for (int i = 1; i <= branchLength; i++) {
                                float bx = cos(branchRot) * i;
                                float bz = sin(branchRot) * i;
                                glm::vec3 branchBlockPos = branchStartPos + glm::vec3(bx, 0, bz);
                                // Store as full block positions.
                                chunk.ancientBranchPositions.push_back(branchBlockPos);
                            }
                            // At the tip, add a mini canopy of leaves.
                            glm::vec3 tip = branchStartPos + glm::vec3(cos(branchRot) * (branchLength + 1), 0, sin(branchRot) * (branchLength + 1));
                            for (int dx = -1; dx <= 1; dx++) {
                                for (int dy = -1; dy <= 1; dy++) {
                                    for (int dz = -1; dz <= 1; dz++) {
                                        if (glm::length(glm::vec3(dx, dy, dz)) < 1.5f) {
                                            chunk.ancientLeafPositions.push_back(tip + glm::vec3(dx, dy, dz));
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // --- Fallen Log Generation (existing) ---
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

                    // --- Leaf Piles Generation (existing) ---
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

                    // --- Enhanced Bush Generation ---
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

                    // --- Ground Branches Generation (existing) ---
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
    // Apply block modifications (existing)
    for (const auto& mod : blockModifications) {
        glm::ivec3 pos = mod.first;
        int modType = mod.second;
        if (pos.x >= chunkX * CHUNK_SIZE && pos.x < (chunkX + 1) * CHUNK_SIZE &&
            pos.z >= chunkZ * CHUNK_SIZE && pos.z < (chunkZ + 1) * CHUNK_SIZE) {
            auto removeBlock = [&](std::vector<glm::vec3>& vec) {
                vec.erase(std::remove_if(vec.begin(), vec.end(), [&](const glm::vec3& v) {
                    return std::abs(v.x - pos.x) < 0.1f &&
                        std::abs(v.y - pos.y) < 0.1f &&
                        std::abs(v.z - pos.z) < 0.1f;
                    }), vec.end());
                };
            if (modType == -1) {
                removeBlock(chunk.stonePositions);
                removeBlock(chunk.treeTrunkPositions);
                removeBlock(chunk.treeLeafPositions);
                removeBlock(chunk.firLeafPositions);
                removeBlock(chunk.waterPositions);
                removeBlock(chunk.waterLilyPositions);
                removeBlock(chunk.fallenTreeTrunkPositions);
                removeBlock(chunk.oakTrunkPositions);
                removeBlock(chunk.oakLeafPositions);
                removeBlock(chunk.leafPilePositions);
                removeBlock(chunk.bushSmallPositions);
                removeBlock(chunk.bushMediumPositions);
                removeBlock(chunk.bushLargePositions);
            }
            else {
                auto existsIn = [&](const std::vector<glm::vec3>& vec) -> bool {
                    for (const auto& v : vec) {
                        if (std::abs(v.x - pos.x) < 0.1f &&
                            std::abs(v.y - pos.y) < 0.1f &&
                            std::abs(v.z - pos.z) < 0.1f)
                            return true;
                    }
                    return false;
                    };
                if (!existsIn(chunk.stonePositions) &&
                    !existsIn(chunk.treeTrunkPositions) &&
                    !existsIn(chunk.treeLeafPositions) &&
                    !existsIn(chunk.firLeafPositions) &&
                    !existsIn(chunk.waterPositions) &&
                    !existsIn(chunk.waterLilyPositions) &&
                    !existsIn(chunk.fallenTreeTrunkPositions) &&
                    !existsIn(chunk.oakTrunkPositions) &&
                    !existsIn(chunk.oakLeafPositions) &&
                    !existsIn(chunk.leafPilePositions) &&
                    !existsIn(chunk.bushSmallPositions) &&
                    !existsIn(chunk.bushMediumPositions) &&
                    !existsIn(chunk.bushLargePositions)) {
                    switch (modType) {
                    case 0:  chunk.stonePositions.push_back(glm::vec3(pos));  break;
                    case 1:  chunk.waterPositions.push_back(glm::vec3(pos));  break;
                    case 2:  chunk.treeTrunkPositions.push_back(glm::vec3(pos));  break;
                    case 3:  chunk.treeLeafPositions.push_back(glm::vec3(pos));   break;
                    case 5:  chunk.waterLilyPositions.push_back(glm::vec3(pos));  break;
                    case 6:  chunk.fallenTreeTrunkPositions.push_back(glm::vec3(pos)); break;
                    case 7:  chunk.firLeafPositions.push_back(glm::vec3(pos));   break;
                    case 8:  chunk.oakTrunkPositions.push_back(glm::vec3(pos));   break;
                    case 9:  chunk.oakLeafPositions.push_back(glm::vec3(pos));    break;
                    case 10: chunk.leafPilePositions.push_back(glm::vec3(pos));   break;
                    case 11: chunk.bushSmallPositions.push_back(glm::vec3(pos));  break;
                    case 12: chunk.bushMediumPositions.push_back(glm::vec3(pos)); break;
                    case 13: chunk.bushLargePositions.push_back(glm::vec3(pos));  break;
                    default:
                        chunk.stonePositions.push_back(glm::vec3(pos));
                        break;
                    }
                }
            }
        }
    }
    chunk.needsMeshUpdate = false;
}

// ==================== Chunk Update ====================
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
            if (chunks.find(pos) == chunks.end()) {
                chunks[pos] = Chunk();
            }
            generateChunkMesh(chunks[pos], x, z);
        }
    }
}

// ==================== CONTROLLER (Input Handling) ====================
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);
    if (firstMouse) { lastX = xpos; lastY = ypos; firstMouse = false; }
    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos; lastY = ypos;
    const float sensitivity = 0.1f;
    xoffset *= sensitivity; yoffset *= sensitivity;
    yaw += xoffset; pitch += yoffset;
    pitch = std::min(pitch, 89.0f);
    pitch = std::max(pitch, -89.0f);
    glm::vec3 front;
    front.x = std::cos(glm::radians(yaw)) * std::cos(glm::radians(pitch));
    front.y = std::sin(glm::radians(pitch));
    front.z = std::sin(glm::radians(yaw)) * std::cos(glm::radians(pitch));
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

// ==================== VIEW (Rendering) ====================
//
// We define a 24×24 “pixel art” look on each block face using per-face UV coordinates and a fragment shader
// that quantizes the UVs into 24 steps, drawing grid lines at each cell boundary.
// Note: The blockColors uniform array size has been increased to 19.
// For branch–type blocks (blockType 14 and 18), the vertex shader reads an extra instance attribute (rotation) at location 3.
// For ancient branches (blockType 18) we now draw full blocks.
// ------------------------- Shader Sources -------------------------
const char* vertexShaderSource = R"(
   #version 330 core
   layout (location = 0) in vec3 aPos;
   layout (location = 1) in vec2 aTexCoord;
   // For normal (non-branch) instancing, location 2 holds a vec3 offset.
   // For branch types, if used, location 2 holds a vec3 offset and location 3 holds a float rotation.
   out vec2 TexCoord;
   out vec3 ourColor;
   out float instanceDistance;
   uniform mat4 model;
   uniform mat4 view;
   uniform mat4 projection;
   uniform int blockType;
   // Array of 19 block colors.
   uniform vec3 blockColors[19];
   // Atmospheric perspective uniform.
   uniform vec3 cameraPos;
   layout (location = 2) in vec3 aOffset;
   layout (location = 3) in float aRotation;
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
   void main(){
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

// ------------------------- Cube Vertex Data -------------------------
float cubeVertices[] = {
    // positions            // texture Coords
    // Front face
   -0.5f, -0.5f,  0.5f,     0.0f, 0.0f,
    0.5f, -0.5f,  0.5f,     1.0f, 0.0f,
    0.5f,  0.5f,  0.5f,     1.0f, 1.0f,
    0.5f,  0.5f,  0.5f,     1.0f, 1.0f,
   -0.5f,  0.5f,  0.5f,     0.0f, 1.0f,
   -0.5f, -0.5f,  0.5f,     0.0f, 0.0f,
   // Right face
   0.5f, -0.5f,  0.5f,     0.0f, 0.0f,
   0.5f, -0.5f, -0.5f,     1.0f, 0.0f,
   0.5f,  0.5f, -0.5f,     1.0f, 1.0f,
   0.5f,  0.5f, -0.5f,     1.0f, 1.0f,
   0.5f,  0.5f,  0.5f,     0.0f, 1.0f,
   0.5f, -0.5f,  0.5f,     0.0f, 0.0f,
   // Back face
   0.5f, -0.5f, -0.5f,     0.0f, 0.0f,
  -0.5f, -0.5f, -0.5f,     1.0f, 0.0f,
  -0.5f,  0.5f, -0.5f,     1.0f, 1.0f,
  -0.5f,  0.5f, -0.5f,     1.0f, 1.0f,
   0.5f,  0.5f, -0.5f,     0.0f, 1.0f,
   0.5f, -0.5f, -0.5f,     0.0f, 0.0f,
   // Left face
  -0.5f, -0.5f, -0.5f,     0.0f, 0.0f,
  -0.5f, -0.5f,  0.5f,     1.0f, 0.0f,
  -0.5f,  0.5f,  0.5f,     1.0f, 1.0f,
  -0.5f,  0.5f,  0.5f,     1.0f, 1.0f,
  -0.5f,  0.5f, -0.5f,     0.0f, 1.0f,
  -0.5f, -0.5f, -0.5f,     0.0f, 0.0f,
  // Top face
 -0.5f,  0.5f,  0.5f,     0.0f, 0.0f,
  0.5f,  0.5f,  0.5f,     1.0f, 0.0f,
  0.5f,  0.5f, -0.5f,     1.0f, 1.0f,
  0.5f,  0.5f, -0.5f,     1.0f, 1.0f,
 -0.5f,  0.5f, -0.5f,     0.0f, 1.0f,
 -0.5f,  0.5f,  0.5f,     0.0f, 0.0f,
 // Bottom face
-0.5f, -0.5f, -0.5f,     0.0f, 0.0f,
 0.5f, -0.5f, -0.5f,     1.0f, 0.0f,
 0.5f, -0.5f,  0.5f,     1.0f, 1.0f,
 0.5f, -0.5f,  0.5f,     1.0f, 1.0f,
-0.5f, -0.5f,  0.5f,     0.0f, 1.0f,
-0.5f, -0.5f, -0.5f,     0.0f, 0.0f
};

// ------------------------- MAIN FUNCTION -------------------------
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

    // ------------------------- Shader Compilation -------------------------
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

    // ------------------------- Setup VAOs and VBOs -------------------------
    unsigned int VAO, redVAO, waterVAO, stoneVAO, treeTrunkVAO, treeLeafVAO, waterLilyVAO, fallenTreeVAO, firLeafVAO;
    unsigned int oakTrunkVAO, oakLeafVAO, leafPileVAO, bushSmallVAO, bushMediumVAO, bushLargeVAO;
    unsigned int ancientTrunkVAO, ancientLeafVAO;
    // For branch-type blocks:
    unsigned int branchVAO, branchInstanceVBO;
    unsigned int ancientBranchVAO, ancientBranchInstanceVBO;

    unsigned int VBO, instanceVBO;
    glGenVertexArrays(1, &VAO);
    glGenVertexArrays(1, &redVAO);
    glGenVertexArrays(1, &waterVAO);
    glGenVertexArrays(1, &stoneVAO);
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

    glGenBuffers(1, &VBO);
    glGenBuffers(1, &instanceVBO);
    glGenBuffers(1, &branchInstanceVBO);
    glGenBuffers(1, &ancientBranchInstanceVBO);

    // Helper: setup a generic cube VAO for non-branch blocks.
    auto setupVAO = [&](unsigned int vao, bool instanced) {
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
    setupVAO(VAO, false);
    setupVAO(redVAO, false);
    setupVAO(waterVAO, true);
    setupVAO(stoneVAO, true);
    setupVAO(treeTrunkVAO, true);
    setupVAO(treeLeafVAO, true);
    setupVAO(waterLilyVAO, true);
    setupVAO(fallenTreeVAO, true);
    setupVAO(firLeafVAO, true);
    setupVAO(oakTrunkVAO, true);
    setupVAO(oakLeafVAO, true);
    setupVAO(leafPileVAO, true);
    setupVAO(bushSmallVAO, true);
    setupVAO(bushMediumVAO, true);
    setupVAO(bushLargeVAO, true);
    setupVAO(ancientTrunkVAO, true);
    setupVAO(ancientLeafVAO, true);

    // Setup VAO for ground branches (type 14)
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

    // Setup VAO for ancient branches (type 18) as full blocks.
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

    // ------------------------- Main Render Loop -------------------------
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
            static_cast<float>(WINDOW_WIDTH) / static_cast<float>(WINDOW_HEIGHT),
            0.1f, 10000.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniform3fv(glGetUniformLocation(shaderProgram, "cameraPos"), 1, glm::value_ptr(cameraPos));

        // Define 19 block colors.
        glm::vec3 blockColors[19];
        blockColors[0] = glm::vec3(0.19f, 0.66f, 0.32f); // Stone
        blockColors[1] = glm::vec3(0.0f, 0.5f, 0.5f);     // Water
        blockColors[2] = glm::vec3(0.29f, 0.21f, 0.13f);  // Pine/Fir trunk
        blockColors[3] = glm::vec3(0.07f, 0.46f, 0.34f);  // Pine leaves
        blockColors[4] = glm::vec3(1.0f, 0.0f, 0.0f);     // Origin debug
        blockColors[5] = glm::vec3(0.2f, 0.7f, 0.2f);     // Water lily
        blockColors[6] = glm::vec3(0.45f, 0.22f, 0.07f);  // Fallen log
        blockColors[7] = glm::vec3(0.13f, 0.54f, 0.13f);  // Fir leaves
        blockColors[8] = glm::vec3(0.55f, 0.27f, 0.07f);  // Oak trunk
        blockColors[9] = glm::vec3(0.36f, 0.6f, 0.33f);   // Oak leaves
        blockColors[10] = glm::vec3(0.44f, 0.39f, 0.32f);  // Leaf pile
        blockColors[11] = glm::vec3(0.35f, 0.43f, 0.30f);  // Bush small
        blockColors[12] = glm::vec3(0.52f, 0.54f, 0.35f);  // Bush medium
        blockColors[13] = glm::vec3(0.6f, 0.61f, 0.35f);   // Bush large
        blockColors[14] = glm::vec3(0.4f, 0.3f, 0.2f);     // Ground branch
        // Index 15 unused.
        blockColors[16] = glm::vec3(0.4f, 0.25f, 0.1f);    // Ancient trunk
        blockColors[17] = glm::vec3(0.2f, 0.5f, 0.2f);     // Ancient leaves
        blockColors[18] = glm::vec3(0.3f, 0.2f, 0.1f);     // Ancient branch
        glUniform3fv(glGetUniformLocation(shaderProgram, "blockColors"), 19, glm::value_ptr(blockColors[0]));

        // Build quadtree over chunks.
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

        // Collect instance data.
        std::vector<glm::vec3> globalStoneInstances;
        std::vector<glm::vec3> globalWaterInstances;
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

        for (Chunk* chunk : visibleChunks) {
            globalStoneInstances.insert(globalStoneInstances.end(), chunk->stonePositions.begin(), chunk->stonePositions.end());
            globalWaterInstances.insert(globalWaterInstances.end(), chunk->waterPositions.begin(), chunk->waterPositions.end());
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
        }

        // Helper lambda for instanced drawing (non-branch blocks)
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
        // Helper lambda for branch drawing (for ground branches, instance data as vec4)
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
        // For ancient branches (type 18)
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

        // Draw regular block types:
        drawInstances(stoneVAO, 0, globalStoneInstances);
        drawInstances(waterVAO, 1, globalWaterInstances);
        drawInstances(treeTrunkVAO, 2, globalTreeTrunkInstances);
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

        // Draw the origin cube (blockType = 4)
        glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), 4);
        {
            glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(VAO);
            glDrawArrays(GL_TRIANGLES, 0, 36);
        }

        // ==================== Block Outline (existing) ====================
        glm::ivec3 selectedBlock = raycastForBlock(false);
        if (selectedBlock.x != -10000) {
            glm::mat4 outlineModel = glm::translate(glm::mat4(1.0f), glm::vec3(selectedBlock));
            outlineModel = glm::scale(outlineModel, glm::vec3(1.05f));
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(outlineModel));
            glm::vec3 outlineColor = glm::vec3(1.0f, 1.0f, 1.0f);
            glm::vec3 oldBlockColor0 = blockColors[0];
            blockColors[0] = outlineColor;
            glUniform3fv(glGetUniformLocation(shaderProgram, "blockColors"), 19, glm::value_ptr(blockColors[0]));
            glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), 0);
            glDisable(GL_DEPTH_TEST);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glLineWidth(2.0f);
            glBindVertexArray(redVAO);
            glDrawArrays(GL_TRIANGLES, 0, 36);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_DEPTH_TEST);
            blockColors[0] = oldBlockColor0;
            glUniform3fv(glGetUniformLocation(shaderProgram, "blockColors"), 19, glm::value_ptr(blockColors[0]));
        }
        // ====================================================================

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // ------------------------- Cleanup -------------------------
    glDeleteVertexArrays(1, &VAO);
    glDeleteVertexArrays(1, &redVAO);
    glDeleteVertexArrays(1, &waterVAO);
    glDeleteVertexArrays(1, &stoneVAO);
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

    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &instanceVBO);
    glDeleteBuffers(1, &branchInstanceVBO);
    glDeleteBuffers(1, &ancientBranchInstanceVBO);
    glDeleteProgram(shaderProgram);
    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

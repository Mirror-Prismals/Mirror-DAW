// ======================================================================
// VoxelGame.cpp
// A single-file version with the block outline (mining) feature added,
// modified so that raycasting checks all block types and forces the outline
// to be visible (by temporarily disabling depth testing).
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
const unsigned int WINDOW_WIDTH = 800;
const unsigned int WINDOW_HEIGHT = 600;
const float RENDER_DISTANCE = 10.0f;
const int CHUNK_SIZE = 16;

// ==================== Global Camera Variables ====================
glm::vec3 cameraPos = glm::vec3(0.0f, 10.0f, 3.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

float yaw = -90.0f;
float pitch = 0.0f;
float lastX = WINDOW_WIDTH / 2.0f;
float lastY = WINDOW_HEIGHT / 2.0f;
bool firstMouse = true;

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

struct Chunk {
    std::vector<glm::vec3> waterPositions;
    std::vector<glm::vec3> stonePositions;
    std::vector<glm::vec3> treeTrunkPositions;
    std::vector<glm::vec3> treeLeafPositions; // Pine leaves
    std::vector<glm::vec3> firLeafPositions;  // Fir leaves
    std::vector<glm::vec3> waterLilyPositions;
    std::vector<glm::vec3> fallenTreeTrunkPositions;
    bool needsMeshUpdate;
    Chunk() : needsMeshUpdate(true) {}
};

struct ivec3_hash {
    std::size_t operator()(const glm::ivec3& v) const {
        return ((std::hash<int>()(v.x) ^ (std::hash<int>()(v.y) << 1)) >> 1)
            ^ (std::hash<int>()(v.z) << 1);
    }
};

// Global block modifications: key is a block position, value indicates modification type (-1 = removed, 0 = placed)
std::unordered_map<glm::ivec3, int, ivec3_hash> blockModifications;

// Global chunk storage
std::unordered_map<ChunkPos, Chunk> chunks;

// ==================== Utility: Perlin Noise ====================
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

// ==================== Raycasting ====================
// Modified: now checks if the candidate block exists in any block array in its chunk,
// in addition to the terrain test.
glm::ivec3 raycastForBlock(bool place) {
    float t = 0.0f;
    // You may adjust the step increment if needed (0.1f used here)
    while (t < 5.0f) {
        glm::vec3 p = cameraPos + t * cameraFront;
        glm::ivec3 candidate = glm::ivec3(std::round(p.x), std::round(p.y), std::round(p.z));
        // Skip if this candidate has been removed
        if (blockModifications.find(candidate) != blockModifications.end() &&
            blockModifications[candidate] == -1) {
            t += 0.1f;
            continue;
        }
        bool exists = false;
        // Determine which chunk this candidate is in.
        int chunkX = static_cast<int>(std::floor(candidate.x / static_cast<float>(CHUNK_SIZE)));
        int chunkZ = static_cast<int>(std::floor(candidate.z / static_cast<float>(CHUNK_SIZE)));
        ChunkPos cp{ chunkX, chunkZ };
        if (chunks.find(cp) != chunks.end()) {
            Chunk& ch = chunks[cp];
            // Check in stone positions.
            for (const auto& pos : ch.stonePositions) {
                if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), 0.1f))) {
                    exists = true; break;
                }
            }
            // Check in tree trunk positions.
            if (!exists) {
                for (const auto& pos : ch.treeTrunkPositions) {
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), 0.1f))) {
                        exists = true; break;
                    }
                }
            }
            // Check in pine tree leaf positions.
            if (!exists) {
                for (const auto& pos : ch.treeLeafPositions) {
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), 0.1f))) {
                        exists = true; break;
                    }
                }
            }
            // Check in fir tree leaf positions.
            if (!exists) {
                for (const auto& pos : ch.firLeafPositions) {
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), 0.1f))) {
                        exists = true; break;
                    }
                }
            }
            // Check in water positions.
            if (!exists) {
                for (const auto& pos : ch.waterPositions) {
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), 0.1f))) {
                        exists = true; break;
                    }
                }
            }
            // Check in water lily positions.
            if (!exists) {
                for (const auto& pos : ch.waterLilyPositions) {
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), 0.1f))) {
                        exists = true; break;
                    }
                }
            }
            // Check in fallen tree trunk positions.
            if (!exists) {
                for (const auto& pos : ch.fallenTreeTrunkPositions) {
                    if (glm::all(glm::epsilonEqual(pos, glm::vec3(candidate), 0.1f))) {
                        exists = true; break;
                    }
                }
            }
        }
        // Fallback: if the terrain test passes (for ground-level blocks)
        TerrainPoint terrain = getTerrainHeight(p.x, p.z);
        if (terrain.isLand && candidate.y <= static_cast<int>(std::floor(terrain.height))) {
            exists = true;
        }
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
// Modified: after modifying a block, mark the corresponding chunk as needing an update.
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
            glm::ivec3 pos = raycastForBlock(true);
            if (pos.x != -10000) {
                blockModifications[pos] = 0;
                int chunkX = static_cast<int>(std::floor(pos.x / static_cast<float>(CHUNK_SIZE)));
                int chunkZ = static_cast<int>(std::floor(pos.z / static_cast<float>(CHUNK_SIZE)));
                ChunkPos cp{ chunkX, chunkZ };
                if (chunks.find(cp) != chunks.end())
                    chunks[cp].needsMeshUpdate = true;
            }
        }
    }
}

// ==================== Quadtree for Spatial Partitioning ====================
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
            if (children[i]) delete children[i];
    }
    glm::vec3 getMinWorld() const { return glm::vec3(minX * CHUNK_SIZE, -4.0f, minZ * CHUNK_SIZE); }
    glm::vec3 getMaxWorld() const { return glm::vec3((maxX + 1) * CHUNK_SIZE, 150.0f, (maxZ + 1) * CHUNK_SIZE); }
    bool contains(const ChunkPos& pos) const {
        return pos.x >= minX && pos.x <= maxX && pos.z >= minZ && pos.z <= maxZ;
    }
    void subdivide() {
        int midX = (minX + maxX) / 2;
        int midZ = (minZ + maxZ) / 2;
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
        if (!contains(item.pos)) return;
        if (!subdivided && items.size() < capacity)
            items.push_back(item);
        else {
            if (!subdivided) subdivide();
            for (int i = 0; i < 4; i++) {
                if (children[i]->contains(item.pos)) {
                    children[i]->insert(item);
                    return;
                }
            }
        }
    }
    void query(const std::vector<Plane>& frustum, std::vector<Chunk*>& out) {
        glm::vec3 nodeMin = getMinWorld();
        glm::vec3 nodeMax = getMaxWorld();
        if (!aabbInFrustum(frustum, nodeMin, nodeMax)) return;
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
    Quadtree(int minX, int minZ, int maxX, int maxZ) {
        root = new QuadtreeNode(minX, minZ, maxX, maxZ);
    }
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
    chunk.waterPositions.clear();
    chunk.stonePositions.clear();
    chunk.treeTrunkPositions.clear();
    chunk.treeLeafPositions.clear();
    chunk.firLeafPositions.clear();
    chunk.waterLilyPositions.clear();
    chunk.fallenTreeTrunkPositions.clear();
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
                        if (hashVal % 300 < 1) {
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
                    int hashValPine = std::abs((intWorldX * 73856093) ^ (intWorldZ * 19349663));
                    if (hashValPine % 10000 < 1) {
                        int trunkHeight = 60;
                        int trunkThickness = 4;
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
                    int hashValFir = std::abs((intWorldX * 83492791) ^ (intWorldZ * 19349663));
                    if (hashValFir % 12000 < 1) {
                        int trunkHeight = 40;
                        int trunkThickness = 3;
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
                }
            }
        }
    }
    // Apply block modifications for stone blocks.
    for (const auto& mod : blockModifications) {
        glm::ivec3 pos = mod.first;
        int modType = mod.second;
        if (pos.x >= chunkX * CHUNK_SIZE && pos.x < (chunkX + 1) * CHUNK_SIZE &&
            pos.z >= chunkZ * CHUNK_SIZE && pos.z < (chunkZ + 1) * CHUNK_SIZE) {
            if (modType == -1) {
                auto& stones = chunk.stonePositions;
                stones.erase(std::remove_if(stones.begin(), stones.end(), [&](const glm::vec3& v) {
                    return std::abs(v.x - pos.x) < 0.1f &&
                        std::abs(v.y - pos.y) < 0.1f &&
                        std::abs(v.z - pos.z) < 0.1f;
                    }), stones.end());
            }
            else {
                bool exists = false;
                for (const auto& v : chunk.stonePositions) {
                    if (std::abs(v.x - pos.x) < 0.1f &&
                        std::abs(v.y - pos.y) < 0.1f &&
                        std::abs(v.z - pos.z) < 0.1f) {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    chunk.stonePositions.push_back(glm::vec3(pos));
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
            if (chunks.find(pos) == chunks.end())
                chunks[pos] = Chunk();
            generateChunkMesh(chunks[pos], x, z);
        }
    }
}

// ==================== CONTROLLER (Input Handling) ====================
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
// Shader sources
const char* vertexShaderSource = R"(
   #version 330 core
   layout (location = 0) in vec3 aPos;
   layout (location = 2) in vec3 aOffset;
   out vec3 ourColor;
   uniform mat4 model;
   uniform mat4 view;
   uniform mat4 projection;
   uniform int blockType;
   uniform vec3 blockColors[8];
   void main() {
       vec3 pos = aPos;
       if (gl_InstanceID > 0)
           pos += aOffset;
       gl_Position = projection * view * model * vec4(pos, 1.0);
       ourColor = blockColors[blockType];
   }
)";
const char* fragmentShaderSource = R"(
   #version 330 core
   in vec3 ourColor;
   out vec4 FragColor;
   uniform int blockType;
   float random(vec2 st) {
       return fract(sin(dot(st.xy, vec2(12.9898,78.233))) * 43758.5453123);
   }
   void main() {
       if(blockType == 0) {
           float grain = random(gl_FragCoord.xy);
           FragColor = vec4(ourColor * (0.8 + 0.2 * grain), 1.0);
       } else {
           float alpha = (blockType == 1) ? 0.5 : 1.0;
           FragColor = vec4(ourColor, alpha);
       }
   }
)";

// ==================== MAIN FUNCTION ====================
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

    // Shader compilation and linking.
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
    // Set up vertex data for a cube (positions only).
    // ------------------------------------------------------------------
    float vertices[] = {
        -0.5f, -0.5f,  0.5f,
         0.5f, -0.5f,  0.5f,
         0.5f,  0.5f,  0.5f,
        -0.5f,  0.5f,  0.5f,
        -0.5f, -0.5f, -0.5f,
         0.5f, -0.5f, -0.5f,
         0.5f,  0.5f, -0.5f,
        -0.5f,  0.5f, -0.5f,
    };
    float redVertices[] = {
        -0.5f, -0.5f,  0.5f,
         0.5f, -0.5f,  0.5f,
         0.5f,  0.5f,  0.5f,
        -0.5f,  0.5f,  0.5f,
        -0.5f, -0.5f, -0.5f,
         0.5f, -0.5f, -0.5f,
         0.5f,  0.5f, -0.5f,
        -0.5f,  0.5f, -0.5f,
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
    // Generate and set up VAOs/VBOs.
    // (We create VAOs for multiple block types.)
    // ------------------------------------------------------------------
    unsigned int VAO, redVAO, VBO, redVBO, waterVAO, stoneVAO, EBO, instanceVBO;
    unsigned int treeTrunkVAO, treeLeafVAO, waterLilyVAO, fallenTreeVAO, firLeafVAO;
    glGenVertexArrays(1, &VAO);
    glGenVertexArrays(1, &redVAO);
    glGenVertexArrays(1, &waterVAO);
    glGenVertexArrays(1, &stoneVAO);
    glGenVertexArrays(1, &treeTrunkVAO);
    glGenVertexArrays(1, &treeLeafVAO);
    glGenVertexArrays(1, &waterLilyVAO);
    glGenVertexArrays(1, &fallenTreeVAO);
    glGenVertexArrays(1, &firLeafVAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &redVBO);
    glGenBuffers(1, &EBO);
    glGenBuffers(1, &instanceVBO);

    // Setup Regular Cube VAO.
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Setup Red Cube VAO.
    glBindVertexArray(redVAO);
    glBindBuffer(GL_ARRAY_BUFFER, redVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(redVertices), redVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Setup Water VAO.
    glBindVertexArray(waterVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    // Setup Stone VAO.
    glBindVertexArray(stoneVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    // Setup Tree Trunk VAO.
    glBindVertexArray(treeTrunkVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    // Setup Pine Tree Leaves VAO.
    glBindVertexArray(treeLeafVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    // Setup Water Lily VAO.
    glBindVertexArray(waterLilyVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    // Setup Fallen Tree (Log) VAO.
    glBindVertexArray(fallenTreeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    // Setup Fir Tree Leaves VAO.
    glBindVertexArray(firLeafVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // ------------------- Main Render Loop -------------------
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

        // Set uniform blockColors (8 entries)
        glm::vec3 blockColors[8];
        blockColors[0] = glm::vec3(0.35f, 0.5f, 0.39f);    // Stone (with noise grain)
        blockColors[1] = glm::vec3(0.0f, 0.5f, 0.5f);       // Water
        blockColors[2] = glm::vec3(0.55f, 0.27f, 0.07f);    // Tree trunk
        blockColors[3] = glm::vec3(0.5f, 0.5f, 0.0f);       // Pine tree leaves
        blockColors[4] = glm::vec3(1.0f, 0.0f, 0.0f);       // Origin cube
        blockColors[5] = glm::vec3(0.2f, 0.7f, 0.2f);       // Water lily
        blockColors[6] = glm::vec3(0.45f, 0.22f, 0.07f);    // Fallen log
        blockColors[7] = glm::vec3(0.0f, 0.8f, 0.0f);       // Fir tree leaves
        glUniform3fv(glGetUniformLocation(shaderProgram, "blockColors"), 8, glm::value_ptr(blockColors[0]));

        // Build a quadtree over loaded chunks.
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

        // Batch instance arrays from visible chunks.
        std::vector<glm::vec3> globalStoneInstances;
        std::vector<glm::vec3> globalWaterInstances;
        std::vector<glm::vec3> globalTreeTrunkInstances;
        std::vector<glm::vec3> globalPineLeafInstances;
        std::vector<glm::vec3> globalFirLeafInstances;
        std::vector<glm::vec3> globalWaterLilyInstances;
        std::vector<glm::vec3> globalFallenTreeTrunkInstances;
        for (Chunk* chunk : visibleChunks) {
            globalStoneInstances.insert(globalStoneInstances.end(), chunk->stonePositions.begin(), chunk->stonePositions.end());
            globalWaterInstances.insert(globalWaterInstances.end(), chunk->waterPositions.begin(), chunk->waterPositions.end());
            globalTreeTrunkInstances.insert(globalTreeTrunkInstances.end(), chunk->treeTrunkPositions.begin(), chunk->treeTrunkPositions.end());
            globalPineLeafInstances.insert(globalPineLeafInstances.end(), chunk->treeLeafPositions.begin(), chunk->treeLeafPositions.end());
            globalFirLeafInstances.insert(globalFirLeafInstances.end(), chunk->firLeafPositions.begin(), chunk->firLeafPositions.end());
            globalWaterLilyInstances.insert(globalWaterLilyInstances.end(), chunk->waterLilyPositions.begin(), chunk->waterLilyPositions.end());
            globalFallenTreeTrunkInstances.insert(globalFallenTreeTrunkInstances.end(), chunk->fallenTreeTrunkPositions.begin(), chunk->fallenTreeTrunkPositions.end());
        }

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

        drawInstances(stoneVAO, 0, globalStoneInstances);
        glUniform1i(glGetUniformLocation(shaderProgram, "blockType"), 4);
        {
            glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
            glBindVertexArray(redVAO);
            glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        }
        drawInstances(treeTrunkVAO, 2, globalTreeTrunkInstances);
        drawInstances(treeLeafVAO, 3, globalPineLeafInstances);
        drawInstances(firLeafVAO, 7, globalFirLeafInstances);
        drawInstances(waterVAO, 1, globalWaterInstances);
        drawInstances(waterLilyVAO, 5, globalWaterLilyInstances);
        drawInstances(fallenTreeVAO, 6, globalFallenTreeTrunkInstances);

        // ==================== NEW FEATURE: Block Outline ====================
        // Raycast to get the block you're looking at.
        glm::ivec3 selectedBlock = raycastForBlock(false);
        if (selectedBlock.x != -10000) {
            // Create a model matrix to position an outline cube at the block's location.
            glm::mat4 outlineModel = glm::translate(glm::mat4(1.0f), glm::vec3(selectedBlock));
            outlineModel = glm::scale(outlineModel, glm::vec3(1.05f)); // Scale up slightly for the outline
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(outlineModel));
            // Override the blockColors uniform with an outline color (white).
            glm::vec3 outlineColor = glm::vec3(1.0f, 1.0f, 1.0f);
            glUniform3fv(glGetUniformLocation(shaderProgram, "blockColors"), 1, glm::value_ptr(outlineColor));
            // Disable depth testing so the outline is drawn on top.
            glDisable(GL_DEPTH_TEST);
            // Switch to wireframe mode to render the outline.
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glLineWidth(2.0f);
            glBindVertexArray(redVAO);
            glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
            // Restore fill mode and re-enable depth testing.
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_DEPTH_TEST);
        }
        // ====================================================================

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Cleanup.
    glDeleteVertexArrays(1, &VAO);
    glDeleteVertexArrays(1, &redVAO);
    glDeleteVertexArrays(1, &waterVAO);
    glDeleteVertexArrays(1, &stoneVAO);
    glDeleteVertexArrays(1, &treeTrunkVAO);
    glDeleteVertexArrays(1, &treeLeafVAO);
    glDeleteVertexArrays(1, &waterLilyVAO);
    glDeleteVertexArrays(1, &fallenTreeVAO);
    glDeleteVertexArrays(1, &firLeafVAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &redVBO);
    glDeleteBuffers(1, &EBO);
    glDeleteBuffers(1, &instanceVBO);
    glDeleteProgram(shaderProgram);
    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

// Plate tectonics terrain with stress viz, optimized faster
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>

const int SCRW = 1280, SCRH = 720;
// reduce render radius to reduce block count
const int RADIUS = 32;
const int maxHeight = 40; // max height cap for blocks

struct Plate {
    glm::vec2 seed;
    glm::vec3 color;
    glm::vec2 velocity;
    Plate(glm::vec2 s, glm::vec3 c, glm::vec2 v) : seed(s), color(c), velocity(v) {}
};

std::vector<Plate> plates;

float hash(int x, int z) {
    int n = x * 73856093 ^ z * 19349663;
    n = (n << 13) ^ n;
    return 1.0f - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0f;
}
float noise(float x, float z) {
    int xi = floor(x), zi = floor(z);
    float xf = x - xi, zf = z - zi;
    float v00 = hash(xi, zi);
    float v10 = hash(xi + 1, zi);
    float v01 = hash(xi, zi + 1);
    float v11 = hash(xi + 1, zi + 1);
    float u = xf * xf * (3 - 2 * xf), v = zf * zf * (3 - 2 * zf);
    return (v00 * (1 - u) + v10 * u) * (1 - v) + (v01 * (1 - u) + v11 * u) * v;
}
float fbm(float x, float z) {
    float sum = 0, amp = 1, freq = 1, maxSum = 0;
    for (int i = 0; i < 4; i++) { // reduce octaves for some perf gain
        sum += noise(x * freq, z * freq) * amp;
        maxSum += amp;
        amp *= 0.5;
        freq *= 2;
    }
    return sum / maxSum;
}

#define NBLEND 2 // fewer blends == faster
void nearestPlates(glm::vec2 p, int* ids, float* ds) {
    for (int i = 0; i < NBLEND; i++) { ids[i] = -1; ds[i] = 1e9; }
    for (size_t i = 0; i < plates.size(); i++) {
        float d = glm::distance(p, plates[i].seed);
        for (int j = 0; j < NBLEND; j++) {
            if (d < ds[j]) {
                for (int k = NBLEND - 1; k > j; k--) { ids[k] = ids[k - 1]; ds[k] = ds[k - 1]; }
                ids[j] = int(i); ds[j] = d; break;
            }
        }
    }
}

float terrainHeightWithStress(float x, float z, float& stress, glm::vec3& col) {
    glm::vec2 p(x, z);
    int ids[NBLEND]; float ds[NBLEND];
    nearestPlates(p, ids, ds);
    float h = 0, wsum = 0;
    float maxAmp = 0;
    glm::vec3 colorSum(0);
    float n = fbm(x * 0.05, z * 0.05);
    for (int j = 0; j < NBLEND; j++) {
        if (ids[j] < 0) continue;
        Plate& P = plates[ids[j]];
        float d = ds[j];
        float w = 1.f / (0.1f + pow(d / 64.f, 2));
        wsum += w;

        float r = glm::length(p - P.seed);
        float dome = 5 * pow(1 - glm::clamp(r / 64.f, 0.f, 1.f), 2);

        float amp = 0;
        const int DX[8] = { -1,-1,0,1,1,1,0,-1 }, DZ[8] = { 0,1,1,1,0,-1,-1,-1 };
        for (int k = 0; k < 8; k++) {
            glm::vec2 q = p + glm::vec2(DX[k], DZ[k]);
            int nids[NBLEND]; float nds[NBLEND];
            nearestPlates(q, nids, nds);
            int nbid = nids[0]; if (nbid == ids[j]) continue;
            glm::vec2 dir = glm::normalize(q - p);
            Plate& N = plates[nbid];
            glm::vec2 rel = P.velocity - N.velocity;
            float proj = glm::dot(rel, dir);
            amp += (proj > 0) ? proj : proj * 0.5f;
        }
        maxAmp = (fabs(amp) > fabs(maxAmp)) ? amp : maxAmp;

        h += w * (10.f + dome + 15.f * amp * n);
        colorSum += w * P.color;
    }
    stress = maxAmp;
    col = colorSum / wsum;
    return h / wsum;
}

glm::vec3 camPos(0, 30, 0), vel(0);
float yaw = -90, pitch = 0, dt = 0, lastTime = 0;
bool onGround = false, firstMouse = true; float lastX = SCRW / 2, lastY = SCRH / 2;

void mouse_cb(GLFWwindow*, double xpos, double ypos) {
    if (firstMouse) { lastX = (float)xpos; lastY = (float)ypos; firstMouse = false; }
    float dx = (float)xpos - lastX, dy = lastY - (float)ypos;
    lastX = (float)xpos; lastY = (float)ypos;
    yaw += dx * 0.1f; pitch += dy * 0.1f;
    if (pitch > 89.f)pitch = 89.f; if (pitch < -89.f)pitch = -89.f;
}
void handle(GLFWwindow* w) {
    glm::vec3 front(glm::cos(glm::radians(yaw)) * glm::cos(glm::radians(pitch)), 0,
        glm::sin(glm::radians(yaw)) * glm::cos(glm::radians(pitch)));
    front = glm::normalize(front);
    glm::vec3 right = glm::normalize(glm::cross(front, { 0,1,0 }));
    glm::vec3 move(0);
    if (glfwGetKey(w, GLFW_KEY_W) == GLFW_PRESS)move += front;
    if (glfwGetKey(w, GLFW_KEY_S) == GLFW_PRESS)move -= front;
    if (glfwGetKey(w, GLFW_KEY_A) == GLFW_PRESS)move -= right;
    if (glfwGetKey(w, GLFW_KEY_D) == GLFW_PRESS)move += right;
    if (glm::length(move) > 0.01f)move = glm::normalize(move);
    vel.x = move.x * 8; vel.z = move.z * 8;
    if (glfwGetKey(w, GLFW_KEY_SPACE) == GLFW_PRESS && onGround) { vel.y = 7; onGround = false; }
}
void collision() {
    float dummyStress; glm::vec3 dummyCol;
    float h = terrainHeightWithStress(camPos.x, camPos.z, dummyStress, dummyCol);
    if (camPos.y - 1.6f < h) { camPos.y = h + 1.6f; vel.y = 0; onGround = true; }
}

float cube[] = { // same
-0.5f,-0.5f,-0.5f,0,0,-1, 0,0, 0.5f,-0.5f,-0.5f,0,0,-1,1,0, 0.5f, 0.5f,-0.5f,0,0,-1,1,1,
0.5f, 0.5f,-0.5f,0,0,-1,1,1,-0.5f, 0.5f,-0.5f,0,0,-1,0,1,-0.5f,-0.5f,-0.5f,0,0,-1,0,0,
// Front
-0.5f,-0.5f,0.5f,0,0,1,0,0,0.5f,-0.5f,0.5f,0,0,1,1,0,0.5f, 0.5f,0.5f,0,0,1,1,1,
0.5f, 0.5f,0.5f,0,0,1,1,1,-0.5f, 0.5f,0.5f,0,0,1,0,1,-0.5f,-0.5f,0.5f,0,0,1,0,0,
// Left
-0.5f,0.5f,0.5f,-1,0,0,1,0,-0.5f,0.5f,-0.5f,-1,0,0,1,1,-0.5f,-0.5f,-0.5f,-1,0,0,0,1,
-0.5f,-0.5f,-0.5f,-1,0,0,0,1,-0.5f,-0.5f,0.5f,-1,0,0,0,0,-0.5f,0.5f,0.5f,-1,0,0,1,0,
// Right
0.5f,0.5f,0.5f,1,0,0,1,0,0.5f,0.5f,-0.5f,1,0,0,1,1,0.5f,-0.5f,-0.5f,1,0,0,0,1,
0.5f,-0.5f,-0.5f,1,0,0,0,1,0.5f,-0.5f,0.5f,1,0,0,0,0,0.5f,0.5f,0.5f,1,0,0,1,0,
// Bottom
-0.5f,-0.5f,-0.5f,0,-1,0,0,1,0.5f,-0.5f,-0.5f,0,-1,0,1,1,0.5f,-0.5f,0.5f,0,-1,0,1,0,
0.5f,-0.5f,0.5f,0,-1,0,1,0,-0.5f,-0.5f,0.5f,0,-1,0,0,0,-0.5f,-0.5f,-0.5f,0,-1,0,0,1,
// Top
-0.5f,0.5f,-0.5f,0,1,0,0,1,0.5f,0.5f,-0.5f,0,1,0,1,1,0.5f,0.5f,0.5f,0,1,0,1,0,
0.5f,0.5f,0.5f,0,1,0,1,0,-0.5f,0.5f,0.5f,0,1,0,0,0,-0.5f,0.5f,-0.5f,0,1,0,0,1,
};

const char* vsh = "#version 330 core\n"
"layout(location=0) in vec3 pos;"
"layout(location=1) in vec3 normal;"
"layout(location=2) in vec2 uv;"
"layout(location=3) in vec3 offset;"
"layout(location=4) in vec3 color0;"
"uniform mat4 model,view,proj;"
"out vec3 vColor;"
"void main(){ gl_Position=proj*view*model*vec4(pos+offset,1); vColor=color0; }";
const char* fsh = "#version 330 core\n"
"in vec3 vColor; out vec4 color; void main(){color=vec4(vColor,1);}";
GLuint compile(const char* v, const char* f) {
    GLuint vs = glCreateShader(GL_VERTEX_SHADER), fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(vs, 1, &v, NULL); glCompileShader(vs);
    int ok; glGetShaderiv(vs, GL_COMPILE_STATUS, &ok);
    if (!ok) { char e[512]; glGetShaderInfoLog(vs, 512, 0, e); std::cout << "VS: " << e << "\n"; }
    glShaderSource(fs, 1, &f, NULL); glCompileShader(fs); glGetShaderiv(fs, GL_COMPILE_STATUS, &ok);
    if (!ok) { char e[512]; glGetShaderInfoLog(fs, 512, 0, e); std::cout << "FS: " << e << "\n"; }
    GLuint p = glCreateProgram(); glAttachShader(p, vs); glAttachShader(p, fs); glLinkProgram(p);
    glDeleteShader(vs); glDeleteShader(fs); return p;
}

int main() {
    srand(time(0));
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* win = glfwCreateWindow(SCRW, SCRH, "Plates Fast", 0, 0);
    glfwMakeContextCurrent(win);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    glfwSetInputMode(win, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPosCallback(win, mouse_cb);
    glEnable(GL_DEPTH_TEST);

    GLuint prog = compile(vsh, fsh);
    GLuint vao, vbo, instvbo;
    glGenVertexArrays(1, &vao); glBindVertexArray(vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, 0, 8 * sizeof(float), (void*)0); glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, 0, 8 * sizeof(float), (void*)(3 * sizeof(float))); glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, 0, 8 * sizeof(float), (void*)(6 * sizeof(float))); glEnableVertexAttribArray(2);

    glGenBuffers(1, &instvbo);
    glBindBuffer(GL_ARRAY_BUFFER, instvbo);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, 0, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_FLOAT, 0, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glVertexAttribDivisor(3, 1); glVertexAttribDivisor(4, 1);

    glm::vec3 fixed[3] = { {0,.7,.7},{.7,.7,0},{1,.28,0} };
    for (int i = 0; i < 12; i++) {
        glm::vec2 s(rand() % 300 - 150, rand() % 300 - 150);
        glm::vec3 c = fixed[i % 3];
        glm::vec2 v((rand() % 200 - 100) / 100.f, (rand() % 200 - 100) / 100.f);
        plates.emplace_back(s, c, v);
    }

    while (!glfwWindowShouldClose(win)) {
        float now = glfwGetTime(); dt = now - lastTime; lastTime = now;
        handle(win);
        if (!onGround)vel.y -= 9.8f * dt;
        camPos += vel * dt; collision();

        struct D { glm::vec3 pos, col; };
        std::vector<D> data; data.reserve(100000); // help avoid realloc
        int cx = int(camPos.x), cz = int(camPos.z);
        for (int x = cx - RADIUS; x <= cx + RADIUS; x++) {
            for (int z = cz - RADIUS; z <= cz + RADIUS; z++) {
                float stress; glm::vec3 col;
                int h = int(terrainHeightWithStress(x, z, stress, col));
                if (h > maxHeight) h = maxHeight; // limit height cap
                float s = glm::clamp(fabs(stress) / 2.f, 0.f, 1.f);
                glm::vec3 blend = glm::mix(col, glm::vec3(1, 0.3, 0), s);
                for (int y = 0; y <= h; y++) {
                    data.push_back({ glm::vec3(x,y,z), blend });
                }
            }
        }
        glBindBuffer(GL_ARRAY_BUFFER, instvbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(D) * data.size(), data.data(), GL_DYNAMIC_DRAW);

        glClearColor(0.5, 0.7, 0.9, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glm::vec3 front(glm::cos(glm::radians(yaw)) * glm::cos(glm::radians(pitch)),
            sin(glm::radians(pitch)),
            glm::sin(glm::radians(yaw)) * glm::cos(glm::radians(pitch)));
        glm::mat4 view = glm::lookAt(camPos, camPos + front, { 0,1,0 });
        glm::mat4 proj = glm::perspective(glm::radians(75.f), float(SCRW) / SCRH, .1f, 600.f);

        glUseProgram(prog);
        glUniformMatrix4fv(glGetUniformLocation(prog, "model"), 1, 0, glm::value_ptr(glm::mat4(1)));
        glUniformMatrix4fv(glGetUniformLocation(prog, "view"), 1, 0, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(prog, "proj"), 1, 0, glm::value_ptr(proj));

        glBindVertexArray(vao);
        glDrawArraysInstanced(GL_TRIANGLES, 0, 36, data.size());

        glfwSwapBuffers(win);
        glfwPollEvents();
    }

    glfwTerminate();
}

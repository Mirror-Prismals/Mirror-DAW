#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>
#include <cmath>

// Window
const unsigned int SCRW = 1280, SCRH = 720;
const int RADIUS = 30;

float hash(int x, int z) {
    int n = x * 73856093 ^ z * 19349663; n = (n << 13) ^ n;
    return 1.f - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.f;
}
float noise2d(float x, float z) {
    int xi = floor(x), zi = floor(z);
    float xf = x - xi, zf = z - zi;
    float v00 = hash(xi, zi), v10 = hash(xi + 1, zi), v01 = hash(xi, zi + 1), v11 = hash(xi + 1, zi + 1);
    float u = xf * xf * (3 - 2 * xf), v = zf * zf * (3 - 2 * zf);
    return (v00 * (1 - u) + v10 * u) * (1 - v) + (v01 * (1 - u) + v11 * u) * v;
}
float fbm(float x, float z) {
    float t = 0, a = 1, f = 1, max = 0;
    for (int i = 0; i < 5; i++) { t += noise2d(x * f, z * f) * a; max += a; a *= 0.5; f *= 2; }
    return t / max;
}

// Camera state
glm::vec3 camPos(0, 50, 0), vel(0);
float yaw = -90, pitch = 0, delta = 0, lastTime = 0;
bool onGround = false, firstMouse = true; float lastX = SCRW / 2, lastY = SCRH / 2;

void mouse(GLFWwindow*, double x, double y) {
    if (firstMouse) { lastX = x; lastY = y; firstMouse = false; }
    float dx = x - lastX, dy = lastY - y; lastX = x; lastY = y;
    yaw += dx * 0.1f; pitch += dy * 0.1f;
    if (pitch > 89)pitch = 89; if (pitch < -89)pitch = -89;
}
void keys(GLFWwindow* w) {
    glm::vec3 forward(glm::cos(glm::radians(yaw)) * glm::cos(glm::radians(pitch)), 0,
        glm::sin(glm::radians(yaw)) * glm::cos(glm::radians(pitch)));
    forward = glm::normalize(forward);
    glm::vec3 right = glm::normalize(glm::cross(forward, { 0,1,0 }));
    glm::vec3 m(0);
    if (glfwGetKey(w, GLFW_KEY_ESCAPE) == GLFW_PRESS)glfwSetWindowShouldClose(w, true);
    if (glfwGetKey(w, GLFW_KEY_W) == GLFW_PRESS)m += forward;
    if (glfwGetKey(w, GLFW_KEY_S) == GLFW_PRESS)m -= forward;
    if (glfwGetKey(w, GLFW_KEY_A) == GLFW_PRESS)m -= right;
    if (glfwGetKey(w, GLFW_KEY_D) == GLFW_PRESS)m += right;
    if (glm::length(m) > 0.01)m = glm::normalize(m);
    vel.x = m.x * 10; vel.z = m.z * 10;
    if (glfwGetKey(w, GLFW_KEY_SPACE) == GLFW_PRESS && onGround) { vel.y = 7; onGround = false; }
}
void collision() {
    int cx = int(camPos.x + 0.5f), cz = int(camPos.z + 0.5f);
    float n = fbm(cx * 0.05, cz * 0.05);
    int h = int(20 + n * 10);
    float feet = camPos.y - 1.6;
    if (feet < h) { camPos.y = h + 1.6; vel.y = 0; onGround = true; }
}

// FULL cube vertex data
float cube[] = {
-0.5f,-0.5f,-0.5f,0,0,-1,0,0,  0.5f,-0.5f,-0.5f,0,0,-1,1,0,  0.5f, 0.5f,-0.5f,0,0,-1,1,1,
0.5f,0.5f,-0.5f,0,0,-1,1,1,  -0.5f,0.5f,-0.5f,0,0,-1,0,1,  -0.5f,-0.5f,-0.5f,0,0,-1,0,0,
// Front
-0.5f,-0.5f,0.5f,0,0,1,0,0,  0.5f,-0.5f,0.5f,0,0,1,1,0,  0.5f,0.5f,0.5f,0,0,1,1,1,
0.5f,0.5f,0.5f,0,0,1,1,1,  -0.5f,0.5f,0.5f,0,0,1,0,1,  -0.5f,-0.5f,0.5f,0,0,1,0,0,
// Left
-0.5f,0.5f,0.5f,-1,0,0,1,0,  -0.5f,0.5f,-0.5f,-1,0,0,1,1,  -0.5f,-0.5f,-0.5f,-1,0,0,0,1,
-0.5f,-0.5f,-0.5f,-1,0,0,0,1,  -0.5f,-0.5f,0.5f,-1,0,0,0,0,  -0.5f,0.5f,0.5f,-1,0,0,1,0,
// Right
0.5f,0.5f,0.5f,1,0,0,1,0,   0.5f,0.5f,-0.5f,1,0,0,1,1,  0.5f,-0.5f,-0.5f,1,0,0,0,1,
0.5f,-0.5f,-0.5f,1,0,0,0,1,  0.5f,-0.5f,0.5f,1,0,0,0,0,  0.5f,0.5f,0.5f,1,0,0,1,0,
// Bottom
-0.5f,-0.5f,-0.5f,0,-1,0,0,1,  0.5f,-0.5f,-0.5f,0,-1,0,1,1,  0.5f,-0.5f,0.5f,0,-1,0,1,0,
0.5f,-0.5f,0.5f,0,-1,0,1,0,  -0.5f,-0.5f,0.5f,0,-1,0,0,0,  -0.5f,-0.5f,-0.5f,0,-1,0,0,1,
// Top
-0.5f,0.5f,-0.5f,0,1,0,0,1,  0.5f,0.5f,-0.5f,0,1,0,1,1,  0.5f,0.5f,0.5f,0,1,0,1,0,
0.5f,0.5f,0.5f,0,1,0,1,0,  -0.5f,0.5f,0.5f,0,1,0,0,0,  -0.5f,0.5f,-0.5f,0,1,0,0,1,
};

const char* vertsrc = R"(#version 330 core
layout(location=0)in vec3 pos;
layout(location=1)in vec3 norm;
layout(location=2)in vec2 uv;
layout(location=3)in vec3 offset;
out vec3 color;
uniform mat4 model,view,proj;
void main(){
color=vec3(0.3+0.7*uv,0.3);
gl_Position=proj*view*model*vec4(pos+offset,1);
})";

const char* fragsrc = R"(#version 330 core
in vec3 color;out vec4 FragColor;
void main(){FragColor=vec4(color,1);}
)";

GLuint compile(const char* v, const char* f) {
    GLuint vs = glCreateShader(GL_VERTEX_SHADER), fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(vs, 1, &v, NULL); glCompileShader(vs);
    int ok; glGetShaderiv(vs, GL_COMPILE_STATUS, &ok); if (!ok) { char e[512]; glGetShaderInfoLog(vs, 512, nullptr, e); puts(e); }
    glShaderSource(fs, 1, &f, NULL); glCompileShader(fs);
    glGetShaderiv(fs, GL_COMPILE_STATUS, &ok); if (!ok) { char e[512]; glGetShaderInfoLog(fs, 512, nullptr, e); puts(e); }
    GLuint p = glCreateProgram(); glAttachShader(p, vs); glAttachShader(p, fs); glLinkProgram(p);
    glDeleteShader(vs); glDeleteShader(fs); return p;
}

int main() {
    srand(time(0));
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* win = glfwCreateWindow(SCRW, SCRH, "Block terrain", nullptr, nullptr);
    glfwMakeContextCurrent(win);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    glfwSetInputMode(win, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPosCallback(win, mouse);
    glEnable(GL_DEPTH_TEST);

    GLuint prog = compile(vertsrc, fragsrc);
    GLuint vao, vbo, ibo, instanceVBO;
    glGenVertexArrays(1, &vao); glBindVertexArray(vao);
    glGenBuffers(1, &vbo); glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0); glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float))); glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float))); glEnableVertexAttribArray(2);

    glGenBuffers(1, &instanceVBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
    glEnableVertexAttribArray(3);
    glVertexAttribDivisor(3, 1);

    while (!glfwWindowShouldClose(win)) {
        float now = glfwGetTime(); delta = now - lastTime; lastTime = now;
        keys(win);
        if (!onGround)vel.y -= 9.81f * delta;
        camPos += vel * delta;
        collision();

        std::vector<glm::vec3> blockPos;
        int cx = int(camPos.x), cz = int(camPos.z);
        for (int x = cx - RADIUS; x <= cx + RADIUS; x++) {
            for (int z = cz - RADIUS; z <= cz + RADIUS; z++) {
                float n = fbm(x * 0.05, z * 0.05);
                int h = int(20 + n * 10);
                for (int y = 0; y <= h; y++) { blockPos.emplace_back(x, y, z); }
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * blockPos.size(), blockPos.data(), GL_DYNAMIC_DRAW);

        glClearColor(0.5, 0.7, 0.9, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::vec3 front(cos(glm::radians(yaw)) * cos(glm::radians(pitch)),
            sin(glm::radians(pitch)),
            sin(glm::radians(yaw)) * cos(glm::radians(pitch)));
        glm::mat4 V = glm::lookAt(camPos, camPos + front, { 0,1,0 });
        glm::mat4 P = glm::perspective(glm::radians(100.f), float(SCRW) / SCRH, 0.1f, 1000.f);

        glUseProgram(prog);
        GLuint mLoc = glGetUniformLocation(prog, "model");
        GLuint vLoc = glGetUniformLocation(prog, "view");
        GLuint pLoc = glGetUniformLocation(prog, "proj");
        glm::mat4 M(1);
        glUniformMatrix4fv(mLoc, 1, GL_FALSE, glm::value_ptr(M));
        glUniformMatrix4fv(vLoc, 1, GL_FALSE, glm::value_ptr(V));
        glUniformMatrix4fv(pLoc, 1, GL_FALSE, glm::value_ptr(P));

        glBindVertexArray(vao);
        glDrawArraysInstanced(GL_TRIANGLES, 0, 36, blockPos.size());

        glfwSwapBuffers(win); glfwPollEvents();
    }
    glfwTerminate();
}

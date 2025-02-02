#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <queue>
#include <map>
#include <functional>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstdint>

// Helper function to check if a string ends with a given suffix (C++17 compatible)
bool endsWith(const std::string& str, const std::string& suffix) {
    if (suffix.size() > str.size())
        return false;
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

const double PI = 3.14159265358979323846;

// Interface for anything that can be saved/restored
class Saveable {
public:
    virtual ~Saveable() = default;
    virtual void save(const std::string& path) const = 0;
    virtual bool restore(const std::string& path) = 0;
    virtual std::string getStateAsString() const = 0;
    virtual bool loadStateFromString(const std::string& state) = 0;
};

struct MidiEvent {
    enum class Type {
        NoteOn,
        NoteOff,
        ControlChange
    };

    Type type;
    uint8_t channel;
    uint8_t note;
    uint8_t velocity;
    double timestamp;
};

class AudioBuffer : public Saveable {
public:
    AudioBuffer(size_t capacity = 44100) : capacity_(capacity) {
        buffer_.resize(capacity);
    }

    void write(float sample) {
        if (writePos_ < capacity_) {
            buffer_[writePos_++] = sample;
        }
    }

    float read() {
        if (readPos_ < writePos_) {
            float sample = buffer_[readPos_++];
            return sample;
        }
        return 0.0f;
    }

    // Added isEmpty() method
    bool isEmpty() const {
        return readPos_ >= writePos_;
    }

    void reset() {
        readPos_ = 0;
    }

    void clear() {
        readPos_ = 0;
        writePos_ = 0;
        buffer_.assign(capacity_, 0.0f);
    }

    // Save interface implementation
    void save(const std::string& path) const override {
        if (endsWith(path, ".wav")) {
            saveAsWav(path);
        }
        else {
            saveBinary(path);
        }
    }

    bool restore(const std::string& path) override {
        if (endsWith(path, ".wav")) {
            return restoreFromWav(path);
        }
        return restoreFromBinary(path);
    }

    std::string getStateAsString() const override {
        std::stringstream ss;
        ss << writePos_ << " ";
        for (size_t i = 0; i < writePos_; ++i) {
            ss << buffer_[i] << " ";
        }
        return ss.str();
    }

    bool loadStateFromString(const std::string& state) override {
        std::stringstream ss(state);
        size_t size;
        ss >> size;
        if (size > capacity_) return false;

        clear();
        float sample;
        while (ss >> sample && writePos_ < size) {
            write(sample);
        }
        return true;
    }

private:
    void saveAsWav(const std::string& path) const {
        std::ofstream file(path, std::ios::binary);
        if (!file) return;

        std::vector<int16_t> intSamples;
        for (size_t i = 0; i < writePos_; ++i) {
            intSamples.push_back(static_cast<int16_t>(buffer_[i] * 32767.0f));
        }

        // WAV header generation
        uint32_t dataSize = intSamples.size() * sizeof(int16_t);
        uint32_t fileSize = 36 + dataSize;
        uint32_t sampleRate = 44100;
        uint16_t channels = 1;

        file.write("RIFF", 4);
        file.write(reinterpret_cast<const char*>(&fileSize), 4);
        file.write("WAVE", 4);
        file.write("fmt ", 4);
        uint32_t fmtSize = 16;
        file.write(reinterpret_cast<const char*>(&fmtSize), 4);
        uint16_t format = 1;
        file.write(reinterpret_cast<const char*>(&format), 2);
        file.write(reinterpret_cast<const char*>(&channels), 2);
        file.write(reinterpret_cast<const char*>(&sampleRate), 4);
        uint32_t byteRate = sampleRate * channels * sizeof(int16_t);
        file.write(reinterpret_cast<const char*>(&byteRate), 4);
        uint16_t blockAlign = channels * sizeof(int16_t);
        file.write(reinterpret_cast<const char*>(&blockAlign), 2);
        uint16_t bitsPerSample = 16;
        file.write(reinterpret_cast<const char*>(&bitsPerSample), 2);
        file.write("data", 4);
        file.write(reinterpret_cast<const char*>(&dataSize), 4);
        file.write(reinterpret_cast<const char*>(intSamples.data()), dataSize);
    }

    void saveBinary(const std::string& path) const {
        std::ofstream file(path, std::ios::binary);
        if (!file) return;
        file.write(reinterpret_cast<const char*>(&writePos_), sizeof(writePos_));
        file.write(reinterpret_cast<const char*>(buffer_.data()), writePos_ * sizeof(float));
    }

    bool restoreFromWav(const std::string& path) {
        std::ifstream file(path, std::ios::binary);
        if (!file) return false;

        char header[44];
        file.read(header, 44);

        uint32_t dataSize;
        file.seekg(40);
        file.read(reinterpret_cast<char*>(&dataSize), 4);

        std::vector<int16_t> intSamples(dataSize / sizeof(int16_t));
        file.read(reinterpret_cast<char*>(intSamples.data()), dataSize);

        clear();
        for (int16_t sample : intSamples) {
            write(sample / 32768.0f);
        }
        return true;
    }

    bool restoreFromBinary(const std::string& path) {
        std::ifstream file(path, std::ios::binary);
        if (!file) return false;

        size_t size;
        file.read(reinterpret_cast<char*>(&size), sizeof(size));
        if (size > capacity_) return false;

        clear();
        file.read(reinterpret_cast<char*>(buffer_.data()), size * sizeof(float));
        writePos_ = size;
        return true;
    }

    std::vector<float> buffer_;
    size_t capacity_;
    size_t readPos_ = 0;
    size_t writePos_ = 0;
};

class Plugin : public Saveable {
public:
    virtual ~Plugin() = default;
    virtual std::string getName() const = 0;
    virtual std::string getType() const = 0;
    virtual void process(std::vector<float>& buffer) = 0;
    virtual void setParameter(const std::string& name, float value) = 0;
    virtual float getParameter(const std::string& name) = 0;
    virtual void loadScript(const std::string& language, const std::string& code) {}
};

class Distortion : public Plugin {
public:
    std::string getName() const override { return "Distortion"; }
    std::string getType() const override { return "Effect"; }

    void process(std::vector<float>& buffer) override {
        for (auto& sample : buffer) {
            sample = std::tanh(drive_ * sample);
        }
    }

    void setParameter(const std::string& name, float value) override {
        if (name == "drive") {
            drive_ = value;
            std::cout << "Distortion drive set to: " << drive_ << std::endl;
        }
    }

    float getParameter(const std::string& name) override {
        if (name == "drive") return drive_;
        return 0.0f;
    }

    void save(const std::string& path) const override {
        std::ofstream file(path);
        file << getStateAsString();
    }

    bool restore(const std::string& path) override {
        std::ifstream file(path);
        if (!file) return false;
        std::string state;
        std::getline(file, state);
        return loadStateFromString(state);
    }

    std::string getStateAsString() const override {
        return std::to_string(drive_);
    }

    bool loadStateFromString(const std::string& state) override {
        try {
            drive_ = std::stof(state);
            return true;
        }
        catch (...) {
            return false;
        }
    }

private:
    float drive_ = 1.0f;
};

class ScriptPlugin : public Plugin {
public:
    ScriptPlugin(const std::string& name) : name_(name) {}

    std::string getName() const override { return name_; }
    std::string getType() const override { return "Script"; }

    void process(std::vector<float>& buffer) override {
        if (processFunc_) {
            processFunc_(buffer);
        }
    }

    void setParameter(const std::string& name, float value) override {}
    float getParameter(const std::string& name) override { return 0.0f; }

    void loadScript(const std::string& language, const std::string& code) override {
        script_ = code;
        if (language == "cpp") {
            processFunc_ = [](std::vector<float>& buf) {
                for (auto& sample : buf) {
                    sample *= 0.5f;
                }
                };
            std::cout << "Loaded script for " << name_ << std::endl;
        }
    }

    void save(const std::string& path) const override {
        std::ofstream file(path);
        file << getStateAsString();
    }

    bool restore(const std::string& path) override {
        std::ifstream file(path);
        if (!file) return false;
        std::string state;
        std::getline(file, state);
        return loadStateFromString(state);
    }

    std::string getStateAsString() const override {
        return name_ + "\n" + script_;
    }

    bool loadStateFromString(const std::string& state) override {
        std::stringstream ss(state);
        std::string name, script;
        std::getline(ss, name);
        std::getline(ss, script);
        name_ = name;
        loadScript("cpp", script);
        return true;
    }

private:
    std::string name_;
    std::string script_;
    std::function<void(std::vector<float>&)> processFunc_;
};

class AudioTrack : public Saveable {
public:
    AudioTrack(int id) : id_(id), volume_(1.0f), isPlaying_(false) {
        buffer_ = std::make_shared<AudioBuffer>();
        std::cout << "Created track " << id_ << std::endl;
    }

    void addPlugin(std::shared_ptr<Plugin> plugin) {
        plugins_.push_back(plugin);
        std::cout << "Added plugin " << plugin->getName() << " to track " << id_ << std::endl;
    }

    void removePlugin(const std::string& name) {
        plugins_.erase(
            std::remove_if(plugins_.begin(), plugins_.end(),
                [&name](const auto& p) { return p->getName() == name; }),
            plugins_.end()
        );
    }

    void generateTestTone(float frequency, float duration, float sampleRate = 44100) {
        float amplitude = 0.5f;
        size_t numSamples = static_cast<size_t>(duration * sampleRate);

        for (size_t i = 0; i < numSamples; i++) {
            float t = static_cast<float>(i) / sampleRate;
            float sample = amplitude * std::sin(2.0f * PI * frequency * t);
            buffer_->write(sample);
        }
        std::cout << "Generated " << duration << "s test tone at " << frequency << "Hz on track " << id_ << std::endl;
    }

    void addMidiEvent(const MidiEvent& event) {
        std::cout << "Track " << id_ << " received MIDI event: ";
        if (event.type == MidiEvent::Type::NoteOn) {
            float freq = 440.0f * std::pow(2.0f, (event.note - 69) / 12.0f);
            std::cout << "Note On, frequency: " << freq << "Hz" << std::endl;
            generateTestTone(freq, 1.0f);
        }
    }

    float getSample() {
        if (!isPlaying_ || buffer_->isEmpty()) {
            return 0.0f;
        }

        float sample = buffer_->read() * volume_;
        std::vector<float> buf = { sample };
        for (auto& plugin : plugins_) {
            plugin->process(buf);
        }
        return buf[0];
    }

    void play() {
        isPlaying_ = true;
        buffer_->reset();
        std::cout << "Track " << id_ << " started playing" << std::endl;
    }

    void stop() {
        isPlaying_ = false;
        std::cout << "Track " << id_ << " stopped playing" << std::endl;
    }

    void setVolume(float volume) {
        volume_ = volume;
        std::cout << "Track " << id_ << " volume set to " << volume_ << std::endl;
    }

    // Save interface implementation
    void save(const std::string& path) const override {
        if (endsWith(path, ".wav")) {
            buffer_->save(path);
            std::cout << "Saved track " << id_ << " audio to: " << path << std::endl;
        }
        else {
            std::ofstream file(path);
            file << getStateAsString();
            std::cout << "Saved track " << id_ << " state to: " << path << std::endl;
        }
    }

    bool restore(const std::string& path) override {
        if (endsWith(path, ".wav")) {
            bool success = buffer_->restore(path);
            if (success) {
                std::cout << "Restored track " << id_ << " audio from: " << path << std::endl;
            }
            return success;
        }
        else {
            std::ifstream file(path);
            if (!file) return false;
            std::string state;
            std::getline(file, state);
            bool success = loadStateFromString(state);
            if (success) {
                std::cout << "Restored track " << id_ << " state from: " << path << std::endl;
            }
            return success;
        }
    }

    std::string getStateAsString() const override {
        std::stringstream ss;
        ss << id_ << " " << volume_ << " " << isPlaying_ << "\n";
        ss << plugins_.size() << "\n";
        for (const auto& plugin : plugins_) {
            ss << plugin->getStateAsString() << "\n";
        }
        ss << buffer_->getStateAsString();
        return ss.str();
    }

    bool loadStateFromString(const std::string& state) override {
        std::stringstream ss(state);
        ss >> id_ >> volume_ >> isPlaying_;

        size_t pluginCount;
        ss >> pluginCount;
        plugins_.clear();

        // Load plugins based on a simple string match
        for (size_t i = 0; i < pluginCount; ++i) {
            std::string pluginState;
            std::getline(ss, pluginState);
            if (pluginState.empty())
                std::getline(ss, pluginState); // Skip empty lines if necessary

            if (pluginState.find("Distortion") != std::string::npos) {
                auto plugin = std::make_shared<Distortion>();
                plugin->loadStateFromString(pluginState);
                plugins_.push_back(plugin);
            }
            else if (pluginState.find("Script") != std::string::npos) {
                auto plugin = std::make_shared<ScriptPlugin>("restored_script");
                plugin->loadStateFromString(pluginState);
                plugins_.push_back(plugin);
            }
        }

        std::string bufferState;
        std::getline(ss, bufferState);
        return buffer_->loadStateFromString(bufferState);
    }

private:
    int id_;
    float volume_;
    bool isPlaying_;
    std::shared_ptr<AudioBuffer> buffer_;
    std::vector<std::shared_ptr<Plugin>> plugins_;
};

class DAW : public Saveable {
public:
    DAW() {
        std::cout << "DAW initialized" << std::endl;
    }

    std::shared_ptr<AudioTrack> createTrack() {
        auto track = std::make_shared<AudioTrack>(nextTrackId_++);
        tracks_.push_back(track);
        return track;
    }

    float processNextFrame() {
        float mixedSample = 0.0f;
        for (auto& track : tracks_) {
            mixedSample += track->getSample();
        }
        return std::max(-1.0f, std::min(1.0f, mixedSample));
    }

    void playAll() {
        std::cout << "Playing all tracks" << std::endl;
        for (auto& track : tracks_) {
            track->play();
        }
    }

    void stopAll() {
        std::cout << "Stopping all tracks" << std::endl;
        for (auto& track : tracks_) {
            track->stop();
        }
    }

    float getNextSample() {
        return processNextFrame();
    }

    void save(const std::string& path) const override {
        std::cout << "Saving DAW state to: " << path << std::endl;
        std::ofstream file(path);
        file << getStateAsString();
    }

    bool restore(const std::string& path) override {
        std::cout << "Restoring DAW state from: " << path << std::endl;
        std::ifstream file(path);
        if (!file) return false;
        std::string state;
        std::getline(file, state, '\0');
        return loadStateFromString(state);
    }

    std::string getStateAsString() const override {
        std::stringstream ss;
        ss << tracks_.size() << "\n" << nextTrackId_ << "\n";
        for (const auto& track : tracks_) {
            ss << track->getStateAsString() << "\n===TRACK_END===\n";
        }
        return ss.str();
    }

    bool loadStateFromString(const std::string& state) override {
        std::stringstream ss(state);
        size_t trackCount;
        ss >> trackCount >> nextTrackId_;

        tracks_.clear();
        for (size_t i = 0; i < trackCount; ++i) {
            auto track = std::make_shared<AudioTrack>(0); // Temporary ID
            std::string trackState;
            std::string line;
            while (std::getline(ss, line) && line != "===TRACK_END===") {
                trackState += line + "\n";
            }
            if (track->loadStateFromString(trackState)) {
                tracks_.push_back(track);
            }
        }
        return true;
    }

private:
    std::vector<std::shared_ptr<AudioTrack>> tracks_;
    int nextTrackId_ = 0;
};

int main() {
    std::cout << "Starting DAW application with full save/restore system\n" << std::endl;

    DAW daw;

    // Create and set up tracks
    auto track1 = daw.createTrack();
    auto track2 = daw.createTrack();

    // Add plugins
    auto distortion = std::make_shared<Distortion>();
    distortion->setParameter("drive", 2.0f);
    track1->addPlugin(distortion);

    auto scriptPlugin = std::make_shared<ScriptPlugin>("CustomGain");
    scriptPlugin->loadScript("cpp", "gain_reduction");
    track2->addPlugin(scriptPlugin);

    // Test MIDI event
    MidiEvent noteOn{
        MidiEvent::Type::NoteOn,
        0,   // channel
        60,  // middle C
        100, // velocity
        0.0  // timestamp
    };
    track1->addMidiEvent(noteOn);

    // Generate some audio on track2
    track2->generateTestTone(554.37f, 2.0f);

    // Set volumes
    track1->setVolume(0.8f);
    track2->setVolume(0.5f);

    // Demonstrate saving
    std::cout << "\nDemonstrating save system..." << std::endl;
    distortion->save("distortion_preset.txt");
    track2->save("track2.wav");
    track1->save("track1_state.txt");
    daw.save("session.daw");

    // Demonstrate restore functionality
    std::cout << "\nDemonstrating restore system..." << std::endl;
    DAW daw2;
    if (daw2.restore("session.daw")) {
        std::cout << "Successfully restored DAW session" << std::endl;
    }

    auto newDistortion = std::make_shared<Distortion>();
    if (newDistortion->restore("distortion_preset.txt")) {
        std::cout << "Successfully restored distortion preset" << std::endl;
    }

    // Play and process original DAW
    daw.playAll();

    std::cout << "\nProcessing first 10 samples:" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < 10; i++) {
        float sample = daw.getNextSample();
        std::cout << "Sample " << i << ": " << sample << std::endl;
    }

    std::cout << "\nPress Enter to stop playback..." << std::endl;
    std::cin.get();

    daw.stopAll();
    return 0;
}

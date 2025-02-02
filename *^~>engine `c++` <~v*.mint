0 #include <iostream>
1 #include <vector>
2 #include <memory>
3 #include <string>
4 #include <cmath>
5 #include <queue>
6 #include <map>
7 #include <functional>
8 #include <fstream>
9 #include <sstream>
10 #include <iomanip>
11 #include <algorithm>
12 #include <cstdint>
13 bool endsWith(const std::string& str, const std::string& suffix) {
14 if (suffix.size() > str.size())
15 return false;
16 return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
17 }
18 const double PI = 3.14159265358979323846;
19 class Saveable {
20 public:
21 virtual ~Saveable() = default;
22 virtual void save(const std::string& path) const = 0;
23 virtual bool restore(const std::string& path) = 0;
24 virtual std::string getStateAsString() const = 0;
25 virtual bool loadStateFromString(const std::string& state) = 0;
26 };
27 struct MidiEvent {
28 enum class Type {
29 NoteOn,
30 NoteOff,
31 ControlChange
32 };
33 Type type;
34 uint8_t channel;
35 uint8_t note;
36 uint8_t velocity;
37 double timestamp;
38 };
39 class AudioBuffer : public Saveable {
40 public:
41 AudioBuffer(size_t capacity = 44100) : capacity_(capacity) {
42 buffer_.resize(capacity);
43 }
44 void write(float sample) {
45 if (writePos_ < capacity_) {
46 buffer_[writePos_++] = sample;
47 }
48 }
49 float read() {
50 if (readPos_ < writePos_) {
51 float sample = buffer_[readPos_++];
52 return sample;
53 }
54 return 0.0f;
55 }
56 bool isEmpty() const {
57 return readPos_ >= writePos_;
58 }
59 void reset() {
60 readPos_ = 0;
61 }
62 void clear() {
63 readPos_ = 0;
64 writePos_ = 0;
65 buffer_.assign(capacity_, 0.0f);
66 }
67 void save(const std::string& path) const override {
68 if (endsWith(path, ".wav")) {
69 saveAsWav(path);
70 }
71 else {
72 saveBinary(path);
73 }
74 }
75 bool restore(const std::string& path) override {
76 if (endsWith(path, ".wav")) {
77 return restoreFromWav(path);
78 }
79 return restoreFromBinary(path);
80 }
81 std::string getStateAsString() const override {
82 std::stringstream ss;
83 ss << writePos_ << " ";
84 for (size_t i = 0; i < writePos_; ++i) {
85 ss << buffer_[i] << " ";
86 }
87 return ss.str();
88 }
89 bool loadStateFromString(const std::string& state) override {
90 std::stringstream ss(state);
91 size_t size;
92 ss >> size;
93 if (size > capacity_) return false;
94 clear();
95 float sample;
96 while (ss >> sample && writePos_ < size) {
97 write(sample);
98 }
99 return true;
100 }
101 private:
102 void saveAsWav(const std::string& path) const {
103 std::ofstream file(path, std::ios::binary);
104 if (!file) return;
105 
106 std::vector<int16_t> intSamples;
107 for (size_t i = 0; i < writePos_; ++i) {
108 intSamples.push_back(static_cast<int16_t>(buffer_[i] * 32767.0f));
109 }
110 uint32_t dataSize = intSamples.size() * sizeof(int16_t);
111 uint32_t fileSize = 36 + dataSize;
112 uint32_t sampleRate = 44100;
113 uint16_t channels = 1;
114 
115 file.write("RIFF", 4);
116 file.write(reinterpret_cast<const char*>(&fileSize), 4);
117 file.write("WAVE", 4);
118 file.write("fmt ", 4);
119 uint32_t fmtSize = 16;
120 file.write(reinterpret_cast<const char*>(&fmtSize), 4);
121 uint16_t format = 1;
122 file.write(reinterpret_cast<const char*>(&format), 2);
123 file.write(reinterpret_cast<const char*>(&channels), 2);
124 file.write(reinterpret_cast<const char*>(&sampleRate), 4);
125 uint32_t byteRate = sampleRate * channels * sizeof(int16_t);
126 file.write(reinterpret_cast<const char*>(&byteRate), 4);
127 uint16_t blockAlign = channels * sizeof(int16_t);
128 file.write(reinterpret_cast<const char*>(&blockAlign), 2);
129 uint16_t bitsPerSample = 16;
130 file.write(reinterpret_cast<const char*>(&bitsPerSample), 2);
131 file.write("data", 4);
132 file.write(reinterpret_cast<const char*>(&dataSize), 4);
133 file.write(reinterpret_cast<const char*>(intSamples.data()), dataSize);
134 }
135 void saveBinary(const std::string& path) const {
136 std::ofstream file(path, std::ios::binary);
137 if (!file) return;
138 file.write(reinterpret_cast<const char*>(&writePos_), sizeof(writePos_));
139 file.write(reinterpret_cast<const char*>(buffer_.data()), writePos_ * sizeof(float));
140 }
141 bool restoreFromWav(const std::string& path) {
142 std::ifstream file(path, std::ios::binary);
143 if (!file) return false;
144 char header[44];
145 file.read(header, 44);
146 uint32_t dataSize;
147 file.seekg(40);
148 file.read(reinterpret_cast<char*>(&dataSize), 4);
149 std::vector<int16_t> intSamples(dataSize / sizeof(int16_t));
150 file.read(reinterpret_cast<char*>(intSamples.data()), dataSize);
151 clear();
152 for (int16_t sample : intSamples) {
153 write(sample / 32768.0f);
154 }
155 return true;
156 }
157 bool restoreFromBinary(const std::string& path) {
158 std::ifstream file(path, std::ios::binary);
159 if (!file) return false;
160 size_t size;
161 file.read(reinterpret_cast<char*>(&size), sizeof(size));
162 if (size > capacity_) return false;
163 clear();
164 file.read(reinterpret_cast<char*>(buffer_.data()), size * sizeof(float));
165 writePos_ = size;
166 return true;
167 }
168 std::vector<float> buffer_;
169 size_t capacity_;
170 size_t readPos_ = 0;
171 size_t writePos_ = 0;
172 };
173 class Plugin : public Saveable {
174 public:
175 virtual ~Plugin() = default;
176 virtual std::string getName() const = 0;
177 virtual std::string getType() const = 0;
178 virtual void process(std::vector<float>& buffer) = 0;
179 virtual void setParameter(const std::string& name, float value) = 0;
180 virtual float getParameter(const std::string& name) = 0;
181 virtual void loadScript(const std::string& language, const std::string& code) {}
182 };
183 class Distortion : public Plugin {
184 public:
185 std::string getName() const override { return "Distortion"; }
186 std::string getType() const override { return "Effect"; }
187 void process(std::vector<float>& buffer) override {
188 for (auto& sample : buffer) {
189 sample = std::tanh(drive_ * sample);
190 }
191 }
192 void setParameter(const std::string& name, float value) override {
193 if (name == "drive") {
194 drive_ = value;
195 std::cout << "Distortion drive set to: " << drive_ << std::endl;
196 }
197 }
198 float getParameter(const std::string& name) override {
199 if (name == "drive") return drive_;
200 return 0.0f;
201 }
202 void save(const std::string& path) const override {
203 std::ofstream file(path);
204 file << getStateAsString();
205 }
206 bool restore(const std::string& path) override {
207 std::ifstream file(path);
208 if (!file) return false;
209 std::string state;
210 std::getline(file, state);
211 return loadStateFromString(state);
212 }
213 std::string getStateAsString() const override {
214 return std::to_string(drive_);
215 }
216 bool loadStateFromString(const std::string& state) override {
217 try {
218 drive_ = std::stof(state);
219 return true;
220 }
221 catch (...) {
222 return false;
223 }
224 }
225 private:
226 float drive_ = 1.0f;
227 };
228 class ScriptPlugin : public Plugin {
229 public:
230 ScriptPlugin(const std::string& name) : name_(name) {}
231 std::string getName() const override { return name_; }
232 std::string getType() const override { return "Script"; }
233 void process(std::vector<float>& buffer) override {
234 if (processFunc_) {
235 processFunc_(buffer);
236 }
237 }
238 void setParameter(const std::string& name, float value) override {}
239 float getParameter(const std::string& name) override { return 0.0f; }
240 void loadScript(const std::string& language, const std::string& code) override {
241 script_ = code;
242 if (language == "cpp") {
243 processFunc_ = [](std::vector<float>& buf) {
244 for (auto& sample : buf) {
245 sample *= 0.5f;
246 }
247 };
248 std::cout << "Loaded script for " << name_ << std::endl;
249 }
250 }
251 void save(const std::string& path) const override {
252 std::ofstream file(path);
253 file << getStateAsString();
254 }
255 bool restore(const std::string& path) override {
256 std::ifstream file(path);
257 if (!file) return false;
258 std::string state;
259 std::getline(file, state);
260 return loadStateFromString(state);
261 }
262 std::string getStateAsString() const override {
263 return name_ + "\n" + script_;
264 }
265 bool loadStateFromString(const std::string& state) override {
266 std::stringstream ss(state);
267 std::string name, script;
268 std::getline(ss, name);
269 std::getline(ss, script);
270 name_ = name;
271 loadScript("cpp", script);
272 return true;
273 }
274 private:
275 std::string name_;
276 std::string script_;
277 std::function<void(std::vector<float>&)> processFunc_;
278 };
279 class AudioTrack : public Saveable {
280 public:
281 AudioTrack(int id) : id_(id), volume_(1.0f), isPlaying_(false) {
282 buffer_ = std::make_shared<AudioBuffer>();
283 std::cout << "Created track " << id_ << std::endl;
284 }
285 void addPlugin(std::shared_ptr<Plugin> plugin) {
286 plugins_.push_back(plugin);
287 std::cout << "Added plugin " << plugin->getName() << " to track " << id_ << std::endl;
288 }
289 void removePlugin(const std::string& name) {
290 plugins_.erase(
291 std::remove_if(plugins_.begin(), plugins_.end(),
292 [&name](const auto& p) { return p->getName() == name; }),
293 plugins_.end()
294 );
295 }
296 void generateTestTone(float frequency, float duration, float sampleRate = 44100) {
297 float amplitude = 0.5f;
298 size_t numSamples = static_cast<size_t>(duration * sampleRate);
299 
300 for (size_t i = 0; i < numSamples; i++) {
301 float t = static_cast<float>(i) / sampleRate;
302 float sample = amplitude * std::sin(2.0f * PI * frequency * t);
303 buffer_->write(sample);
304 }
305 std::cout << "Generated " << duration << "s test tone at " << frequency << "Hz on track " << id_ << std::endl;
306 }
307 void addMidiEvent(const MidiEvent& event) {
308 std::cout << "Track " << id_ << " received MIDI event: ";
309 if (event.type == MidiEvent::Type::NoteOn) {
310 float freq = 440.0f * std::pow(2.0f, (event.note - 69) / 12.0f);
311 std::cout << "Note On, frequency: " << freq << "Hz" << std::endl;
312 generateTestTone(freq, 1.0f);
313 }
314 }
315 float getSample() {
316 if (!isPlaying_ || buffer_->isEmpty()) {
317 return 0.0f;
318 }
319 
320 float sample = buffer_->read() * volume_;
321 std::vector<float> buf = { sample };
322 for (auto& plugin : plugins_) {
323 plugin->process(buf);
324 }
325 return buf[0];
326 }
327 void play() {
328 isPlaying_ = true;
329 buffer_->reset();
330 std::cout << "Track " << id_ << " started playing" << std::endl;
331 }
332 void stop() {
333 isPlaying_ = false;
334 std::cout << "Track " << id_ << " stopped playing" << std::endl;
335 }
336 void setVolume(float volume) {
337 volume_ = volume;
338 std::cout << "Track " << id_ << " volume set to " << volume_ << std::endl;
339 }
340 // Save interface implementation
341 void save(const std::string& path) const override {
342 if (endsWith(path, ".wav")) {
343 buffer_->save(path);
344 std::cout << "Saved track " << id_ << " audio to: " << path << std::endl;
345 }
346 else {
347 std::ofstream file(path);
348 file << getStateAsString();
349 std::cout << "Saved track " << id_ << " state to: " << path << std::endl;
350 }
351 }
352 bool restore(const std::string& path) override {
353 if (endsWith(path, ".wav")) {
354 bool success = buffer_->restore(path);
355 if (success) {
356 std::cout << "Restored track " << id_ << " audio from: " << path << std::endl;
357 }
358 return success;
359 }
360 else {
361 std::ifstream file(path);
362 if (!file) return false;
363 std::string state;
364 std::getline(file, state);
365 bool success = loadStateFromString(state);
366 if (success) {
367 std::cout << "Restored track " << id_ << " state from: " << path << std::endl;
368 }
369 return success;
370 }
371 }
372 std::string getStateAsString() const override {
373 std::stringstream ss;
374 ss << id_ << " " << volume_ << " " << isPlaying_ << "\n";
375 ss << plugins_.size() << "\n";
376 for (const auto& plugin : plugins_) {
377 ss << plugin->getStateAsString() << "\n";
378 }
379 ss << buffer_->getStateAsString();
380 return ss.str();
381 }
382 bool loadStateFromString(const std::string& state) override {
383 std::stringstream ss(state);
384 ss >> id_ >> volume_ >> isPlaying_;
385 size_t pluginCount;
386 ss >> pluginCount;
387 plugins_.clear();
388 // Load plugins based on a simple string match
389 for (size_t i = 0; i < pluginCount; ++i) {
390 std::string pluginState;
391 std::getline(ss, pluginState);
392 if (pluginState.empty())
393 std::getline(ss, pluginState); // Skip empty lines if necessary
394 
395 if (pluginState.find("Distortion") != std::string::npos) {
396 auto plugin = std::make_shared<Distortion>();
397 plugin->loadStateFromString(pluginState);
398 plugins_.push_back(plugin);
399 }
400 else if (pluginState.find("Script") != std::string::npos) {
401 auto plugin = std::make_shared<ScriptPlugin>("restored_script");
402 plugin->loadStateFromString(pluginState);
403 plugins_.push_back(plugin);
404 }
405 }
406 std::string bufferState;
407 std::getline(ss, bufferState);
408 return buffer_->loadStateFromString(bufferState);
409 }
410 private:
411 int id_;
412 float volume_;
413 bool isPlaying_;
414 std::shared_ptr<AudioBuffer> buffer_;
415 std::vector<std::shared_ptr<Plugin>> plugins_;
416 };
417 class DAW : public Saveable {
418 public:
419 DAW() {
420 std::cout << "DAW initialized" << std::endl;
421 }
422 std::shared_ptr<AudioTrack> createTrack() {
423 auto track = std::make_shared<AudioTrack>(nextTrackId_++);
424 tracks_.push_back(track);
425 return track;
426 }
427 float processNextFrame() {
428 float mixedSample = 0.0f;
429 for (auto& track : tracks_) {
430 mixedSample += track->getSample();
431 }
432 return std::max(-1.0f, std::min(1.0f, mixedSample));
433 }
434 void playAll() {
435 std::cout << "Playing all tracks" << std::endl;
436 for (auto& track : tracks_) {
437 track->play();
438 }
439 }
440 void stopAll() {
441 std::cout << "Stopping all tracks" << std::endl;
442 for (auto& track : tracks_) {
443 track->stop();
444 }
445 }
446 float getNextSample() {
447 return processNextFrame();
448 }
449 void save(const std::string& path) const override {
450 std::cout << "Saving DAW state to: " << path << std::endl;
451 std::ofstream file(path);
452 file << getStateAsString();
453 }
454 bool restore(const std::string& path) override {
455 std::cout << "Restoring DAW state from: " << path << std::endl;
456 std::ifstream file(path);
457 if (!file) return false;
458 std::string state;
459 std::getline(file, state, '\0');
460 return loadStateFromString(state);
461 }
462 std::string getStateAsString() const override {
463 std::stringstream ss;
464 ss << tracks_.size() << "\n" << nextTrackId_ << "\n";
465 for (const auto& track : tracks_) {
466 ss << track->getStateAsString() << "\n===TRACK_END===\n";
467 }
468 return ss.str();
469 }
470 bool loadStateFromString(const std::string& state) override {
471 std::stringstream ss(state);
472 size_t trackCount;
473 ss >> trackCount >> nextTrackId_;
474 
475 tracks_.clear();
476 for (size_t i = 0; i < trackCount; ++i) {
477 auto track = std::make_shared<AudioTrack>(0); // Temporary ID
478 std::string trackState;
479 std::string line;
480 while (std::getline(ss, line) && line != "===TRACK_END===") {
481 trackState += line + "\n";
482 }
483 if (track->loadStateFromString(trackState)) {
484 tracks_.push_back(track);
485 }
486 }
487 return true;
488 }
489 private:
490 std::vector<std::shared_ptr<AudioTrack>> tracks_;
491 int nextTrackId_ = 0;
492 };
493 int main() {
494 std::cout << "Starting DAW application with full save/restore system\n" << std::endl;
495 DAW daw;
496 // Create and set up tracks
497 auto track1 = daw.createTrack();
498 auto track2 = daw.createTrack();
499 // Add plugins
500 auto distortion = std::make_shared<Distortion>();
501 distortion->setParameter("drive", 2.0f);
502 track1->addPlugin(distortion);
503 auto scriptPlugin = std::make_shared<ScriptPlugin>("CustomGain");
504 scriptPlugin->loadScript("cpp", "gain_reduction");
505 track2->addPlugin(scriptPlugin);
506 // Test MIDI event
507 MidiEvent noteOn{
508 MidiEvent::Type::NoteOn,
509 0,   // channel
510 60,  // middle C
511 100, // velocity
512 0.0  // timestamp
513 };
514 track1->addMidiEvent(noteOn);
515 // Generate some audio on track2
516 track2->generateTestTone(554.37f, 2.0f);
517 // Set volumes
518 track1->setVolume(0.8f);
519 track2->setVolume(0.5f);
520 // Demonstrate saving
521 std::cout << "\nDemonstrating save system..." << std::endl;
522 distortion->save("distortion_preset.txt");
523 track2->save("track2.wav");
524 track1->save("track1_state.txt");
525 daw.save("session.daw");
526 // Demonstrate restore functionality
527 std::cout << "\nDemonstrating restore system..." << std::endl;
528 DAW daw2;
529 if (daw2.restore("session.daw")) {
530 std::cout << "Successfully restored DAW session" << std::endl;
531 }
532 auto newDistortion = std::make_shared<Distortion>();
533 if (newDistortion->restore("distortion_preset.txt")) {
534 std::cout << "Successfully restored distortion preset" << std::endl;
535 }
536 // Play and process original DAW
537 daw.playAll();
538 std::cout << "\nProcessing first 10 samples:" << std::endl;
539 std::cout << std::fixed << std::setprecision(6);
540 for (int i = 0; i < 10; i++) {
541 float sample = daw.getNextSample();
542 std::cout << "Sample " << i << ": " << sample << std::endl;
543 }
544 std::cout << "\nPress Enter to stop playback..." << std::endl;
545 std::cin.get();
546 daw.stopAll();
547 return 0;
548 }
549
550
551
552
553
554
555
556
557
558
559
560
561
562
563
564
565
566
567
568
569
570
571
572
573
574
575
576
577
578
579
580
581
582
583
584
585
586
587
588
589
590
591
592
593
594
595
596
597
598
599
600
601
602
603
604
605
606
607
608
609
610
611
612
613
614
615
616
617
618
619
620
621
622
623
624
625
626
627
628
629
630
631
632
633
634
635
636
637
638
639
640
641
642
643
644
645
646
647
648
649
650
651
652
653
654
655
656
657
658
659
660
661
662
663
664
665
666

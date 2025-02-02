0 0 #!/usr/bin/env python3
1 0 1 
2 1 0 def generate_tree(base_start, base_end, layers_config):
3 2 3 || """
4 3 4 Generate a list of lines representing the tree.
5 4 5 ||
6 5 6 Each line starts with the base number. For each additional layer (from layers_config)
7 6 7 we check if the current number n falls into one of the defined segments.
8 7 8 If it does, we append the computed value (n - segment_start) to that line.
9 8 9 ||
10 9 10 |>
11 10 11 layers_config is a list where each element represents one layer.
12 11 12 Each layer is itself a list of segments.
13 12 13 A segment is a tuple: (segment_start, segment_end) (both inclusive).
14 13 14 """ |>
15 14 15 lines = []
16 15 16 for n in range(base_start, base_end + 1):
17 16 17 # Start with the base number (layer 0)
18 17 18 parts = [str(n)]
19 18 19 # Process each additional layer in order
20 19 20 for layer in layers_config:
21 20 21 appended = None
22 21 22 for seg in layer:
23 22 23 seg_start, seg_end = seg
24 23 24 if seg_start <= n <= seg_end:
25 24 25 appended = n - seg_start
26 25 26 break  # only one segment applies per layer
27 26 27 if appended is not None:
28 27 28 parts.append(str(appended))
29 28 29 lines.append(" ".join(parts))
30 29 30 return lines
31 30 31
32 31 32 def main():
33 32 33 # Set the base range (for example, 0 to 14)
34 33 34 base_start = 0
35 34 35 base_end = 14
36 35 36
37 36 37 # Define layer configurations.
38 37 38 # Each layer is a list of segments given as (segment_start, segment_end).
39 38 39 # For the example tree:
40 39 40 #   Layer 1 segments: [(1,1), (3,6), (8,13)]
41 40 41 #   Layer 2 segments: [(4,5), (9,12)]
42 41 41 #   Layer 3 segments: [(11,11)]
43 42 42 layers_config = [
44 43 44 [(1, 1), (3, 6), (8, 13)],  # Layer 1
45 44 45 [(4, 5), (9, 12)],          # Layer 2
46 45 46 [(11, 11)]                  # Layer 3
47 46 47 ]
48 47 48
49 48 49 # Generate the tree lines
50 49 50 tree_lines = generate_tree(base_start, base_end, layers_config)
51 50 51
52 51 52 # Print the result
53 52 53 for line in tree_lines:
54 53 54 print(line)
55 54 55
56 55 56 # Save the tree to a text file
57 56 57 output_filename = "tree.txt"
58 57 58 with open(output_filename, "w") as f:
59 58 59 for line in tree_lines:
60 59 60 f.write(line + "\n")
61 60 61 print(f"Tree saved to '{output_filename}'.")
62 61 62
63 62 63 if __name__ == "__main__":
64 63 63 main()
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80 0 0 #include <iostream>
81 1 1 #include <fstream>
82 2 2 #include <sstream>
83 3 3 #include <string>
84 4 4 #include <algorithm>
85 5 5 #include <cctype>
86 6 6 
87 7 7 // Trim leading whitespace from a string.
88 8 8 std::string trimLeading(const std::string& s) {
89 9 9 size_t start = s.find_first_not_of(" \t");
90 10 10 return (start != std::string::npos) ? s.substr(start) : "";
91 11 11 }
92 12 12 
93 13 13 // Mode "t": Convert any text file to midathonic format.
94 14 14 void toMidathonic(const std::string& inputFile, const std::string& outputFile) {
95 15 15 std::ifstream in(inputFile);
96 16 16 if (!in) {
97 17 17 std::cerr << "Error: Unable to open input file '" << inputFile << "'." << std::endl;
98 18 18 return;
99 19 19 }
100 20 20 
101 21 21 std::ofstream out(outputFile);
102 22 22 if (!out) {
103 23 23 std::cerr << "Error: Unable to open output file '" << outputFile << "' for writing." << std::endl;
104 24 24 return;
105 25 25 }
106 26 26 
107 27 27 // Write a header line; change this as desired.
108 28 28 out << "{|;;;;|}" << std::endl;
109 29 29 
110 30 30 std::string line;
111 31 31 int lineNumber = 0;
112 32 32 while (std::getline(in, line)) {
113 33 33 std::string trimmed = trimLeading(line);
114 34 34 out << lineNumber << " " << trimmed << std::endl;
115 35 35 ++lineNumber;
116 36 36 }
117 37 37 
118 38 38 std::cout << "Conversion to midathonic format complete. Output written to '"
119 39 39 << outputFile << "'." << std::endl;
120 40 40 }
121 41 41 
122 42 42 // Helper: Check if a string consists solely of digits.
123 43 43 bool isDigits(const std::string& s) {
124 44 44 return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
125 45 45 }
126 46 46 
127 47 47 // Mode "s": Strip the leading number and following space from each line.
128 48 48 void stripLineNumbers(const std::string& inputFile, const std::string& outputFile) {
129 49 49 std::ifstream in(inputFile);
130 50 50 if (!in) {
131 51 51 std::cerr << "Error: Unable to open input file '" << inputFile << "'." << std::endl;
132 52 52 return;
133 53 53 }
134 54 54 
135 55 55 std::ofstream out(outputFile);
136 56 56 if (!out) {
137 57 57 std::cerr << "Error: Unable to open output file '" << outputFile << "' for writing." << std::endl;
138 58 58 return;
139 59 59 }
140 60 60 
141 61 61 std::string line;
142 62 62 while (std::getline(in, line)) {
143 63 63 size_t pos = line.find(' ');
144 64 64 if (pos != std::string::npos) {
145 65 65 std::string possibleNumber = line.substr(0, pos);
146 66 66 if (isDigits(possibleNumber)) {
147 67 67 out << line.substr(pos + 1) << std::endl;
148 68 68 continue;
149 69 69 }
150 70 70 }
151 71 71 out << line << std::endl;
152 72 72 }
153 73 73 
154 74 74 std::cout << "Stripping line numbers complete. Output written to '"
155 75 75 << outputFile << "'." << std::endl;
156 76 76 }
157 77 77 
158 78 78 // Print usage instructions.
159 79 79 void printUsage(const char* progName) {
160 80 80 std::cout << "Usage:\n"
161 81 81 << "  " << progName << " mode input_file output_file\n\n"
162 82 82 << "Modes:\n"
163 83 83 << "  t or to_midathonic : Convert a text file to midathonic format.\n"
164 84 84 << "  s or strip         : Remove leading line numbers from a midathonic file.\n";
165 85 85 }
166 86 86 
167 87 87 int main(int argc, char* argv[]) {
168 88 88 if (argc != 4) {
169 89 89 printUsage(argv[0]);
170 90 90 return 1;
171 91 91 }
172 92 92 
173 93 93 std::string mode = argv[1];
174 94 94 std::string inputFile = argv[2];
175 95 95 std::string outputFile = argv[3];
176 96 96 
177 97 97 if (mode == "t" || mode == "to_midathonic") {
178 98 98 toMidathonic(inputFile, outputFile);
179 99 99 }
180 100 100 else if (mode == "s" || mode == "strip") {
181 101 101 stripLineNumbers(inputFile, outputFile);
182 102 102 }
183 103 103 else {
184 104 104 std::cerr << "Error: Unknown mode '" << mode << "'." << std::endl;
185 105 105 printUsage(argv[0]);
186 106 106 return 1;
187 107 107 }
188 108 108 
189 109 109 return 0;
190 110 110 }
191 111
192 112
193 113
194 114
195 115
196
197
198
199
200 0 0 #include <iostream>
201 1 1 #include <fstream>
202 2 2 #include <sstream>
203 3 3 #include <vector>
204 4 4 #include <string>
205 5 5 #include <cstdlib>
206 6 6 
207 7 7 // Prints usage instructions.
208 8 8 void printUsage(const char* progName) {
209 9 9 std::cout << "Usage: " << progName << " map_file code_file insertion_line output_file\n"
210 10 10 << "  map_file:      The premade map (tree) file\n"
211 11 11 << "  code_file:     The file containing the code snippet to insert\n"
212 12 12 << "  insertion_line:The line number (from the map file) where insertion begins\n"
213 13 13 << "  output_file:   The file to write the merged output\n";
214 14 14 }
215 15 15 
216 16 16 // Reads an entire file into a vector of strings (one per line).
217 17 17 std::vector<std::string> readFileLines(const std::string& filename) {
218 18 18 std::vector<std::string> lines;
219 19 19 std::ifstream in(filename);
220 20 20 if (!in) {
221 21 21 std::cerr << "Error: Unable to open file '" << filename << "' for reading.\n";
222 22 22 return lines;
223 23 23 }
224 24 24 std::string line;
225 25 25 while (std::getline(in, line)) {
226 26 26 lines.push_back(line);
227 27 27 }
228 28 28 return lines;
229 29 29 }
230 30 30 
231 31 31 // Writes a vector of strings (lines) into a file.
232 32 32 bool writeFileLines(const std::string& filename, const std::vector<std::string>& lines) {
233 33 33 std::ofstream out(filename);
234 34 34 if (!out) {
235 35 35 std::cerr << "Error: Unable to open file '" << filename << "' for writing.\n";
236 36 36 return false;
237 37 37 }
238 38 38 for (const auto& line : lines) {
239 39 39 out << line << "\n";
240 40 40 }
241 41 41 return true;
242 42 42 }
243 43 43 
244 44 44 int main(int argc, char* argv[]) {
245 45 45 if (argc != 5) {
246 46 46 printUsage(argv[0]);
247 47 47 return 1;
248 48 48 }
249 49 49 
250 50 50 // Get command-line arguments.
251 51 51 std::string mapFile = argv[1];
252 52 52 std::string codeFile = argv[2];
253 53 53 int insertionLine = std::atoi(argv[3]);
254 54 54 std::string outputFile = argv[4];
255 55 55 
256 56 56 // Read the map and code files.
257 57 57 std::vector<std::string> mapLines = readFileLines(mapFile);
258 58 58 if (mapLines.empty()) {
259 59 59 std::cerr << "Error: Map file is empty or could not be read.\n";
260 60 60 return 1;
261 61 61 }
262 62 62 
263 63 63 std::vector<std::string> codeLines = readFileLines(codeFile);
264 64 64 if (codeLines.empty()) {
265 65 65 std::cerr << "Error: Code file is empty or could not be read.\n";
266 66 66 return 1;
267 67 67 }
268 68 68 
269 69 69 // Check that insertionLine is non-negative.
270 70 70 if (insertionLine < 0) {
271 71 71 std::cerr << "Error: insertion_line must be non-negative.\n";
272 72 72 return 1;
273 73 73 }
274 74 74 
275 75 75 // (Optional) If insertionLine is beyond the current map size, fill the gap with new lines.
276 76 76 // This ensures we have a valid target even if inserting at the very end.
277 77 77 while (insertionLine > static_cast<int>(mapLines.size())) {
278 78 78 // Create a new line with its base number.
279 79 79 mapLines.push_back(std::to_string(mapLines.size()));
280 80 80 }
281 81 81 
282 82 82 // For each code snippet line, append it to the corresponding map line.
283 83 83 // If the target line exists, append " " + codeLine.
284 84 84 // Otherwise, create a new line with the base number and the code.
285 85 85 for (size_t i = 0; i < codeLines.size(); ++i) {
286 86 86 int targetIndex = insertionLine + static_cast<int>(i);
287 87 87 if (targetIndex < static_cast<int>(mapLines.size())) {
288 88 88 // Append a space and the code line.
289 89 89 mapLines[targetIndex] += " " + codeLines[i];
290 90 90 }
291 91 91 else {
292 92 92 // Create a new line. Its base number is the targetIndex.
293 93 93 mapLines.push_back(std::to_string(targetIndex) + " " + codeLines[i]);
294 94 94 }
295 95 95 }
296 96 96 
297 97 97 // Write the merged output to the output file.
298 98 98 if (!writeFileLines(outputFile, mapLines)) {
299 99 99 std::cerr << "Error writing output file.\n";
300 100 100 return 1;
301 101 101 }
302 102 102 
303 103 103 std::cout << "Insertion complete. Output written to '" << outputFile << "'.\n";
304 104 104 return 0;
305 105 105 }
306 106
307 107
308 108
309 109
310 110
311 111
312 112
313 113
314 114

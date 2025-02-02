{|;;;;|}
0 #include <iostream>
1 #include <fstream>
2 #include <sstream>
3 #include <vector>
4 #include <string>
5 #include <cstdlib>
6 
7 // Prints usage instructions.
8 void printUsage(const char* progName) {
9 std::cout << "Usage: " << progName << " map_file code_file insertion_line output_file\n"
10 << "  map_file:      The premade map (tree) file\n"
11 << "  code_file:     The file containing the code snippet to insert\n"
12 << "  insertion_line:The line number (from the map file) where insertion begins\n"
13 << "  output_file:   The file to write the merged output\n";
14 }
15 
16 // Reads an entire file into a vector of strings (one per line).
17 std::vector<std::string> readFileLines(const std::string& filename) {
18 std::vector<std::string> lines;
19 std::ifstream in(filename);
20 if (!in) {
21 std::cerr << "Error: Unable to open file '" << filename << "' for reading.\n";
22 return lines;
23 }
24 std::string line;
25 while (std::getline(in, line)) {
26 lines.push_back(line);
27 }
28 return lines;
29 }
30 
31 // Writes a vector of strings (lines) into a file.
32 bool writeFileLines(const std::string& filename, const std::vector<std::string>& lines) {
33 std::ofstream out(filename);
34 if (!out) {
35 std::cerr << "Error: Unable to open file '" << filename << "' for writing.\n";
36 return false;
37 }
38 for (const auto& line : lines) {
39 out << line << "\n";
40 }
41 return true;
42 }
43 
44 int main(int argc, char* argv[]) {
45 if (argc != 5) {
46 printUsage(argv[0]);
47 return 1;
48 }
49 
50 // Get command-line arguments.
51 std::string mapFile = argv[1];
52 std::string codeFile = argv[2];
53 int insertionLine = std::atoi(argv[3]);
54 std::string outputFile = argv[4];
55 
56 // Read the map and code files.
57 std::vector<std::string> mapLines = readFileLines(mapFile);
58 if (mapLines.empty()) {
59 std::cerr << "Error: Map file is empty or could not be read.\n";
60 return 1;
61 }
62 
63 std::vector<std::string> codeLines = readFileLines(codeFile);
64 if (codeLines.empty()) {
65 std::cerr << "Error: Code file is empty or could not be read.\n";
66 return 1;
67 }
68 
69 // Check that insertionLine is non-negative.
70 if (insertionLine < 0) {
71 std::cerr << "Error: insertion_line must be non-negative.\n";
72 return 1;
73 }
74 
75 // (Optional) If insertionLine is beyond the current map size, fill the gap with new lines.
76 // This ensures we have a valid target even if inserting at the very end.
77 while (insertionLine > static_cast<int>(mapLines.size())) {
78 // Create a new line with its base number.
79 mapLines.push_back(std::to_string(mapLines.size()));
80 }
81 
82 // For each code snippet line, append it to the corresponding map line.
83 // If the target line exists, append " " + codeLine.
84 // Otherwise, create a new line with the base number and the code.
85 for (size_t i = 0; i < codeLines.size(); ++i) {
86 int targetIndex = insertionLine + static_cast<int>(i);
87 if (targetIndex < static_cast<int>(mapLines.size())) {
88 // Append a space and the code line.
89 mapLines[targetIndex] += " " + codeLines[i];
90 }
91 else {
92 // Create a new line. Its base number is the targetIndex.
93 mapLines.push_back(std::to_string(targetIndex) + " " + codeLines[i]);
94 }
95 }
96 
97 // Write the merged output to the output file.
98 if (!writeFileLines(outputFile, mapLines)) {
99 std::cerr << "Error writing output file.\n";
100 return 1;
101 }
102 
103 std::cout << "Insertion complete. Output written to '" << outputFile << "'.\n";
104 return 0;
105 }

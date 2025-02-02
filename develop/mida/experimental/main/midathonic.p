{|;;;;|}
0 #include <iostream>
1 #include <fstream>
2 #include <sstream>
3 #include <string>
4 #include <algorithm>
5 #include <cctype>
6 
7 // Trim leading whitespace from a string.
8 std::string trimLeading(const std::string& s) {
9 size_t start = s.find_first_not_of(" \t");
10 return (start != std::string::npos) ? s.substr(start) : "";
11 }
12 
13 // Mode "t": Convert any text file to midathonic format.
14 void toMidathonic(const std::string& inputFile, const std::string& outputFile) {
15 std::ifstream in(inputFile);
16 if (!in) {
17 std::cerr << "Error: Unable to open input file '" << inputFile << "'." << std::endl;
18 return;
19 }
20 
21 std::ofstream out(outputFile);
22 if (!out) {
23 std::cerr << "Error: Unable to open output file '" << outputFile << "' for writing." << std::endl;
24 return;
25 }
26 
27 // Write a header line; change this as desired.
28 out << "{|;;;;|}" << std::endl;
29 
30 std::string line;
31 int lineNumber = 0;
32 while (std::getline(in, line)) {
33 std::string trimmed = trimLeading(line);
34 out << lineNumber << " " << trimmed << std::endl;
35 ++lineNumber;
36 }
37 
38 std::cout << "Conversion to midathonic format complete. Output written to '"
39 << outputFile << "'." << std::endl;
40 }
41 
42 // Helper: Check if a string consists solely of digits.
43 bool isDigits(const std::string& s) {
44 return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
45 }
46 
47 // Mode "s": Strip the leading number and following space from each line.
48 void stripLineNumbers(const std::string& inputFile, const std::string& outputFile) {
49 std::ifstream in(inputFile);
50 if (!in) {
51 std::cerr << "Error: Unable to open input file '" << inputFile << "'." << std::endl;
52 return;
53 }
54 
55 std::ofstream out(outputFile);
56 if (!out) {
57 std::cerr << "Error: Unable to open output file '" << outputFile << "' for writing." << std::endl;
58 return;
59 }
60 
61 std::string line;
62 while (std::getline(in, line)) {
63 size_t pos = line.find(' ');
64 if (pos != std::string::npos) {
65 std::string possibleNumber = line.substr(0, pos);
66 if (isDigits(possibleNumber)) {
67 out << line.substr(pos + 1) << std::endl;
68 continue;
69 }
70 }
71 out << line << std::endl;
72 }
73 
74 std::cout << "Stripping line numbers complete. Output written to '"
75 << outputFile << "'." << std::endl;
76 }
77 
78 // Print usage instructions.
79 void printUsage(const char* progName) {
80 std::cout << "Usage:\n"
81 << "  " << progName << " mode input_file output_file\n\n"
82 << "Modes:\n"
83 << "  t or to_midathonic : Convert a text file to midathonic format.\n"
84 << "  s or strip         : Remove leading line numbers from a midathonic file.\n";
85 }
86 
87 int main(int argc, char* argv[]) {
88 if (argc != 4) {
89 printUsage(argv[0]);
90 return 1;
91 }
92 
93 std::string mode = argv[1];
94 std::string inputFile = argv[2];
95 std::string outputFile = argv[3];
96 
97 if (mode == "t" || mode == "to_midathonic") {
98 toMidathonic(inputFile, outputFile);
99 }
100 else if (mode == "s" || mode == "strip") {
101 stripLineNumbers(inputFile, outputFile);
102 }
103 else {
104 std::cerr << "Error: Unknown mode '" << mode << "'." << std::endl;
105 printUsage(argv[0]);
106 return 1;
107 }
108 
109 return 0;
110 }

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>

// Prints usage instructions.
void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " map_file code_file insertion_line output_file\n"
        << "  map_file:      The premade map (tree) file\n"
        << "  code_file:     The file containing the code snippet to insert\n"
        << "  insertion_line:The line number (from the map file) where insertion begins\n"
        << "  output_file:   The file to write the merged output\n";
}

// Reads an entire file into a vector of strings (one per line).
std::vector<std::string> readFileLines(const std::string& filename) {
    std::vector<std::string> lines;
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Error: Unable to open file '" << filename << "' for reading.\n";
        return lines;
    }
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    return lines;
}

// Writes a vector of strings (lines) into a file.
bool writeFileLines(const std::string& filename, const std::vector<std::string>& lines) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Unable to open file '" << filename << "' for writing.\n";
        return false;
    }
    for (const auto& line : lines) {
        out << line << "\n";
    }
    return true;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printUsage(argv[0]);
        return 1;
    }

    // Get command-line arguments.
    std::string mapFile = argv[1];
    std::string codeFile = argv[2];
    int insertionLine = std::atoi(argv[3]);
    std::string outputFile = argv[4];

    // Read the map and code files.
    std::vector<std::string> mapLines = readFileLines(mapFile);
    if (mapLines.empty()) {
        std::cerr << "Error: Map file is empty or could not be read.\n";
        return 1;
    }

    std::vector<std::string> codeLines = readFileLines(codeFile);
    if (codeLines.empty()) {
        std::cerr << "Error: Code file is empty or could not be read.\n";
        return 1;
    }

    // Check that insertionLine is non-negative.
    if (insertionLine < 0) {
        std::cerr << "Error: insertion_line must be non-negative.\n";
        return 1;
    }

    // (Optional) If insertionLine is beyond the current map size, fill the gap with new lines.
    // This ensures we have a valid target even if inserting at the very end.
    while (insertionLine > static_cast<int>(mapLines.size())) {
        // Create a new line with its base number.
        mapLines.push_back(std::to_string(mapLines.size()));
    }

    // For each code snippet line, append it to the corresponding map line.
    // If the target line exists, append " " + codeLine.
    // Otherwise, create a new line with the base number and the code.
    for (size_t i = 0; i < codeLines.size(); ++i) {
        int targetIndex = insertionLine + static_cast<int>(i);
        if (targetIndex < static_cast<int>(mapLines.size())) {
            // Append a space and the code line.
            mapLines[targetIndex] += " " + codeLines[i];
        }
        else {
            // Create a new line. Its base number is the targetIndex.
            mapLines.push_back(std::to_string(targetIndex) + " " + codeLines[i]);
        }
    }

    // Write the merged output to the output file.
    if (!writeFileLines(outputFile, mapLines)) {
        std::cerr << "Error writing output file.\n";
        return 1;
    }

    std::cout << "Insertion complete. Output written to '" << outputFile << "'.\n";
    return 0;
}

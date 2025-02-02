#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// For convenience, define a type alias for a segment: a pair of (start, end)
using Segment = std::pair<int, int>;

// The generateTree function builds the map as a vector of strings.
// Each string represents one line: the base number followed by one value per layer 
// if the current base number falls within one of that layer's segments.
std::vector<std::string> generateTree(int baseStart, int baseEnd,
    const std::vector<std::vector<Segment>>& layersConfig) {
    std::vector<std::string> lines;

    // Loop over each base number
    for (int n = baseStart; n <= baseEnd; ++n) {
        std::ostringstream line;
        // Base layer (layer 0) is always printed
        line << n;

        // For each additional layer in the configuration:
        for (const auto& layer : layersConfig) {
            // Check each segment in this layer; if the number falls within one,
            // compute the extra value as (n - segment_start) and break out.
            bool segmentFound = false;
            for (const auto& seg : layer) {
                int segStart = seg.first;
                int segEnd = seg.second;
                if (n >= segStart && n <= segEnd) {
                    line << " " << (n - segStart);
                    segmentFound = true;
                    break; // Only one segment applies per layer.
                }
            }
            // If no segment applies in this layer, nothing is appended.
        }

        lines.push_back(line.str());
    }

    return lines;
}

int main() {
    // Customize the base layer range. For example, 0 to 14 or 0 to 299 for larger maps.
    int baseStart = 0;
    int baseEnd = 14; // Change this to 299 (or any value) for a longer base layer.

    // Configure additional layers.
    // Each layer is defined as a vector of segments (each segment is a pair: (start, end)).
    // For the given example:
    //   Layer 1: segments: (1,1), (3,6), (8,13)
    //   Layer 2: segments: (4,5), (9,12)
    //   Layer 3: segments: (11,11)
    std::vector<std::vector<Segment>> layersConfig;
    layersConfig.push_back({ {1, 1}, {3, 6}, {8, 13} });  // Layer 1
    layersConfig.push_back({ {4, 5}, {9, 12} });           // Layer 2
    layersConfig.push_back({ {11, 11} });                  // Layer 3

    // Generate the tree/map.
    std::vector<std::string> treeLines = generateTree(baseStart, baseEnd, layersConfig);

    // Output the map to the console.
    for (const auto& line : treeLines) {
        std::cout << line << std::endl;
    }

    // Save the map to a text file.
    std::ofstream outFile("tree.txt");
    if (outFile) {
        for (const auto& line : treeLines) {
            outFile << line << "\n";
        }
        outFile.close();
        std::cout << "Map saved to 'tree.txt'." << std::endl;
    }
    else {
        std::cerr << "Error: Could not open file for writing." << std::endl;
    }

    return 0;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>

// Trim leading whitespace from a string.
std::string trimLeading(const std::string& s) {
    size_t start = s.find_first_not_of(" \t");
    return (start != std::string::npos) ? s.substr(start) : "";
}

// Mode "t": Convert any text file to midathonic format.
void toMidathonic(const std::string& inputFile, const std::string& outputFile) {
    std::ifstream in(inputFile);
    if (!in) {
        std::cerr << "Error: Unable to open input file '" << inputFile << "'." << std::endl;
        return;
    }

    std::ofstream out(outputFile);
    if (!out) {
        std::cerr << "Error: Unable to open output file '" << outputFile << "' for writing." << std::endl;
        return;
    }

    // Write a header line; change this as desired.
    out << "{|;;;;|}" << std::endl;

    std::string line;
    int lineNumber = 0;
    while (std::getline(in, line)) {
        std::string trimmed = trimLeading(line);
        out << lineNumber << " " << trimmed << std::endl;
        ++lineNumber;
    }

    std::cout << "Conversion to midathonic format complete. Output written to '"
        << outputFile << "'." << std::endl;
}

// Helper: Check if a string consists solely of digits.
bool isDigits(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

// Mode "s": Strip the leading number and following space from each line.
void stripLineNumbers(const std::string& inputFile, const std::string& outputFile) {
    std::ifstream in(inputFile);
    if (!in) {
        std::cerr << "Error: Unable to open input file '" << inputFile << "'." << std::endl;
        return;
    }

    std::ofstream out(outputFile);
    if (!out) {
        std::cerr << "Error: Unable to open output file '" << outputFile << "' for writing." << std::endl;
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        size_t pos = line.find(' ');
        if (pos != std::string::npos) {
            std::string possibleNumber = line.substr(0, pos);
            if (isDigits(possibleNumber)) {
                out << line.substr(pos + 1) << std::endl;
                continue;
            }
        }
        out << line << std::endl;
    }

    std::cout << "Stripping line numbers complete. Output written to '"
        << outputFile << "'." << std::endl;
}

// Print usage instructions.
void printUsage(const char* progName) {
    std::cout << "Usage:\n"
        << "  " << progName << " mode input_file output_file\n\n"
        << "Modes:\n"
        << "  t or to_midathonic : Convert a text file to midathonic format.\n"
        << "  s or strip         : Remove leading line numbers from a midathonic file.\n";
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printUsage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    std::string inputFile = argv[2];
    std::string outputFile = argv[3];

    if (mode == "t" || mode == "to_midathonic") {
        toMidathonic(inputFile, outputFile);
    }
    else if (mode == "s" || mode == "strip") {
        stripLineNumbers(inputFile, outputFile);
    }
    else {
        std::cerr << "Error: Unknown mode '" << mode << "'." << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    return 0;
}

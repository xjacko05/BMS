/** 
 * BMS 2023
 * Project 1
 * Implementation and Frequency analysis of Affine cipher
 * by Martin Jacko <xjacko05@stud.fit.vutbr.cz>
 **/

#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <sstream>

#define ASCII_OFFSET_UPPERCASE 65
#define ASCII_OFFSET_LOWERCASE 97
#define ALPHABET_LENGTH 26

#define FROM_ASCII(c) ((c < ASCII_OFFSET_LOWERCASE) ? (c - ASCII_OFFSET_UPPERCASE) : (c - ASCII_OFFSET_LOWERCASE))
#define TO_ASCII(c) (((c < 0) ? c + ALPHABET_LENGTH : c) + ASCII_OFFSET_UPPERCASE)

#define BUFFER_SIZE 4096

std::vector<int> AsciiToVector(const std::string& input) {

    std::vector<int> binaryASCII;

    for (char ch : input) {
        // Ignore non-alphanumeric characters
        if (std::isalnum(ch)) {
            // Convert the character to its ASCII value
            int asciiValue = static_cast<int>(ch);

            int binary[8] = {0,0,0,0,0,0,0,0};
            int j = 0;

            while (asciiValue > 0) { 

                binary[j++] = asciiValue % 2;
                asciiValue = asciiValue / 2;
            }

            for (int k = 8 - 1; k >= 0; --k){
                binaryASCII.push_back(binary[k]);
            }
        }
    }

    return binaryASCII;
}

std::vector<int> BinaryToVector(const std::string& input) {

    std::vector<int> binaryVector;

    for (char ch : input) {
        if (ch == '0' || ch == '1') {
            binaryVector.push_back(ch - '0');  // Convert character '0' or '1' to integer 0 or 1
        }
    }

    return binaryVector;
}

std::string VectorToBinary(const std::vector<int>& inputVector) {

    std::stringstream ss;

    for (int value : inputVector) {
        ss << value;
    }

    return ss.str();
}

std::string VectorToAscii(const std::vector<int>& binaryVector) {
    // Ensure the binary vector size is a multiple of 8
    size_t size = binaryVector.size();

    // Convert binary vector to ASCII string
    std::string asciiString;
    for (size_t i = 0; i < size; i += 8) {
        int asciiValue = 0;
        for (size_t j = 0; j < 8; ++j) {
            asciiValue = (asciiValue << 1) | binaryVector[i + j];
        }
        asciiString.push_back(static_cast<char>(asciiValue));
    }

    return asciiString;
}

std::vector<std::vector<int>> csvToMatrix(const std::string& filename) {
    std::vector<std::vector<int>> result;

    std::ifstream file(filename);
    if (!file.is_open()) {
        // Handle the case where the file could not be opened
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return result;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::vector<int> row;

        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            // Split each cell using both comma and semicolon as delimiters
            if (cell.empty()){
                continue;
            }
            std::istringstream cellStream(cell);
            std::string value;
            while (std::getline(cellStream, value, ';')) {
                if (value.empty()){
                    continue;
                }
                try {
                    int intValue = std::stoi(value);
                    row.push_back(intValue);
                } catch (...) {
                    // Handle the case where conversion to integer fails
                    std::cerr << "Error: Invalid value in file " << filename << ": " << value << ", value treated as 0" << std::endl;
                    row.push_back(0);
                }
            }
        }

        result.push_back(row);
    }

    return result;
}

void MatrixToCsv(const std::vector<std::vector<int>>& data) {
    std::ofstream file('matica.csv');
    if (!file.is_open()) {
        // Handle the case where the file could not be opened
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            // Add a comma if not the last element in the row
            if (i < row.size() - 1) {
                file << ',';
            }
        }
        file << '\n';  // Add a newline character after each row
    }
}

enum Mode {
    ENCODE,
    DECODE,
    NONE
};

/* handles profram flow and argument parsing */
int main(int argc, char *argv[]){

    std::string matrixFileName;
    Mode mode;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-m" && i + 1 < argc) {
            matrixFileName = argv[i + 1];
            i++;
        } else if (arg == "-e") {
            mode = Mode::ENCODE;
        } else if (arg == "-d") {
            mode = Mode::DECODE;
        }
    }

    if (matrixFileName.empty() && mode == Mode::DECODE){
        std::cerr << "Decoding requires -m argument" << std::endl;
        return 1;
    }

    std::string input;

    getline(std::cin, input);

    std::vector<std::vector<int>> parityCheckMatrix;

    if (mode == Mode::ENCODE){

        auto binary = AsciiToVector(input);

        if (matrixFileName.empty()){
            parityCheckMatrix = parity_check_matrix(binary.size(), 0);
        }

        if (!binary.empty()){
            
            std::cout << VectorToBinary(binary) << std::endl;
        }
    }else if (mode == Mode::DECODE){

        auto binary = BinaryToVector(input);

        parityCheckMatrix = csvToMatrix(matrixFileName);
        
        if (!binary.empty()){
            
            std::cout << VectorToAscii(binary) << std::endl;
        }
    }

    return 0;
}
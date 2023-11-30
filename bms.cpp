/** 
 * BMS 2023
 * Project 1
 * Implementation of LDPC encoder and decoder
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
#include <chrono>
#include <random>

std::vector<int> extractNthElements(const std::vector<std::vector<int>>& input, int n) {
    std::vector<int> result;
    for (const auto& row : input) {
        result.push_back(row[n]);
    }
    return result;
}

std::vector<std::vector<int>> generateIdentityMatrix(int n) {
    // Initialize the identity matrix with zeros
    std::vector<std::vector<int>> identityMatrix(n, std::vector<int>(n, 0));

    // Set diagonal elements to 1
    for (int i = 0; i < n; ++i) {
        identityMatrix[i][i] = 1;
    }

    return identityMatrix;
}

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> gaussjordan(const std::vector<std::vector<int>>& X) {
    std::vector<std::vector<int>> A = X;
    int m = A.size();//rows
    int n = A[0].size();//cols

    std::vector<std::vector<int>> P = generateIdentityMatrix(m);

    int pivot_old = -1;
    for (int j = 0; j < n; ++j) {
        //std::vector<int> filtre_down(A.begin() + pivot_old + 1, A.end());
        std::vector<int> one_col = extractNthElements(A, j);
        std::vector<int> filtre_down(one_col.begin() + pivot_old + 1, one_col.end());
        auto max_it = std::max_element(filtre_down.begin(), filtre_down.end());
        int pivot = std::distance(filtre_down.begin(), max_it) + pivot_old + 1;

        if (A[pivot][j]) {
            pivot_old += 1;
            if (pivot_old != pivot) {
                std::swap(A[pivot], A[pivot_old]);
                std::swap(P[pivot], P[pivot_old]);
            }

            for (int i = 0; i < m; ++i) {
                if (i != pivot_old && A[i][j]) {
                    for (int k = 0; k < n; ++k) {
                        A[i][k] = abs(A[i][k] - A[pivot_old][k]);
                        //std::cout << i << "/" << m << "\t" << k << "/" << n << "\n";
                        if (static_cast<unsigned long>(k) < P[0].size()){
                            P[i][k] = abs(P[i][k] - P[pivot_old][k]);
                            //std::cout << i << "/" << P.size() << "\t" << k << "/" << P[0].size() << "\n";
                        }
                    }
                    if (n < m){
                        for (int k = n; k < m; ++k) {
                            P[i][k] = abs(P[i][k] - P[pivot_old][k]);
                        }
                    }
                }
            }
        }

        if (pivot_old == m - 1) {
            break;
        }
    }
    //std::cout << "HERE\n" << m << " " << n << std::endl;

    return std::make_pair(A, P);
}

std::vector<std::vector<int>> binaryproduct(const std::vector<std::vector<int>>& matrix1, const std::vector<std::vector<int>>& matrix2) {
    // Check if matrices are compatible for multiplication
    if (matrix1.empty() || matrix2.empty() || matrix1[0].size() != matrix2.size()) {
        // Matrices cannot be multiplied
        return {};
    }

    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int cols2 = matrix2[0].size();

    // Initialize the result matrix with zeros
    std::vector<std::vector<int>> result(rows1, std::vector<int>(cols2, 0));

    // Perform matrix multiplication
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            for (int k = 0; k < cols1; ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
            // Apply mod 2 to the result
            result[i][j] %= 2;
        }
    }

    return result;
}

std::vector<std::vector<int>> transposeMatrix(const std::vector<std::vector<int>>& matrix) {
    // Get the number of rows and columns in the original matrix
    int rows = matrix.size();
    int cols = matrix[0].size();

    // Initialize the transposed matrix with zeros
    std::vector<std::vector<int>> transposedMatrix(cols, std::vector<int>(rows, 0));

    // Perform the transpose operation
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposedMatrix[j][i] = matrix[i][j];
        }
    }

    return transposedMatrix;
}

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> coding_matrix_systematic(const std::vector<std::vector<int>>& H) {
    int n_equations = H.size();
    int n_code = H[0].size();

    // Create an identity matrix
    std::vector<std::vector<int>> P1 = generateIdentityMatrix(n_code);

    // Compute row-reduced form of H
    auto Hrowreduced = gaussjordan(H);

    // Compute the number of bits
    int n_bits = n_code - std::count_if(Hrowreduced.first.begin(), Hrowreduced.first.end(),
                                        [](const std::vector<int>& a) {
                                            return std::any_of(a.begin(), a.end(), [](int x) { return x != 0; });
                                        });

    // Permutation loop
    while (true) {
        std::vector<int> zeros;
        for (int i = 0; i < std::min(n_equations, n_code); ++i) {
            if (Hrowreduced.first[i][i] == 0) {
                zeros.push_back(i);
            }
        }

        if (!zeros.empty()) {
            int indice_colonne_a = *std::min_element(zeros.begin(), zeros.end());
            std::vector<int> list_ones;
            for (int j = indice_colonne_a + 1; j < n_code; ++j) {
                if (Hrowreduced.first[indice_colonne_a][j] != 0) {
                    list_ones.push_back(j);
                }
            }

            if (!list_ones.empty()) {
                int indice_colonne_b = *std::min_element(list_ones.begin(), list_ones.end());

                // Swap columns in Hrowreduced
                for (int i = 0; i < n_equations; ++i) {
                    std::swap(Hrowreduced.first[i][indice_colonne_a], Hrowreduced.first[i][indice_colonne_b]);
                }

                // Swap columns in P1
                for (int i = 0; i < n_code; ++i) {
                    std::swap(P1[i][indice_colonne_a], P1[i][indice_colonne_b]);
                }
            } else {
                break;
            }
        } else {
            break;
        }
    }

    // Transpose P1
    std::vector<std::vector<int>> P1_T = transposeMatrix(P1);

    // Permutation matrix P2
    std::vector<int> identity_vec(n_code, 0);
    for (int i = 0; i < n_code; ++i) {
        identity_vec[i] = i;
    }
    std::vector<int> sigma(identity_vec.begin() + (n_code - n_bits), identity_vec.end());
    sigma.insert(sigma.end(), identity_vec.begin(), identity_vec.begin() + (n_code - n_bits));

    std::vector<std::vector<int>> P2(n_code, std::vector<int>(n_code, 0));
    for (int i = 0; i < n_code; ++i) {
        P2[identity_vec[i]][sigma[i]] = 1;
    }

    // Compute permutation matrix P
    std::vector<std::vector<int>> P = binaryproduct(P2, P1_T);

    // Compute H_new
    std::vector<std::vector<int>> H_new = binaryproduct(H, transposeMatrix(P));

    // Compute G_systematic
    std::vector<std::vector<int>> G_systematic(n_bits, std::vector<int>(n_code, 0));
    for (int i = 0; i < n_bits; ++i) {
        G_systematic[i][i] = 1;
        for (int j = 0; j < (n_code - n_bits); ++j) {
            G_systematic[i][n_bits + j] = Hrowreduced.first[j][n_code - n_bits + i];
        }
    }

    return std::make_pair(H_new, transposeMatrix(G_systematic));
}

std::vector<std::vector<int>> coding_matrix(const std::vector<std::vector<int>>& H) {
    int n_code = H[0].size();

    // DOUBLE GAUSS-JORDAN:

    auto Href_colonnes_tQ = gaussjordan(transposeMatrix(H));

    auto Href_diag = gaussjordan(transposeMatrix(Href_colonnes_tQ.first));

    //auto Q = Href_colonnes_tQ.second;
    auto Q = transposeMatrix(Href_colonnes_tQ.second);

    //int n_bits = n_code - std::accumulate(Href_diag.begin(), Href_diag.end(), 0);
    int n_bits = n_code - static_cast<int>(std::accumulate(Href_diag.first.begin(), Href_diag.first.end(), 0, 
        [](int sum, const std::vector<int>& row) {
            return sum + std::accumulate(row.begin(), row.end(), 0);
        }));

    std::vector<std::vector<int>> Y(n_code, std::vector<int>(n_bits, 0));
    for (int i = n_code - n_bits; i < n_code; ++i) {
        Y[i][i - (n_code - n_bits)] = 1;
    }

    auto tG = binaryproduct(Q, Y);

    //return transposeMatrix(tG);
    return tG;
}

// Function to generate a regular Parity-Check Matrix H
std::vector<std::vector<int>> parity_check_matrix(int length, int seed) {
    //std::default_random_engine rng = check_random_state(seed);
    if (seed == 0){
        seed = std::chrono::steady_clock::now().time_since_epoch().count();
    }
    std::mt19937 rng(seed);

    int n_code = 2 * length;
    int d_v = length - 1;
    int d_c = length;

    int n_equations = (n_code * d_v) / d_c;

    std::vector<std::vector<int>> block(n_equations / d_v, std::vector<int>(n_code, 0));
    std::vector<std::vector<int>> H(n_equations, std::vector<int>(n_code));
    int block_size = n_equations / d_v;

    // Filling the first block with consecutive ones in each row of the block
    for (int i = 0; i < block_size; ++i) {
        for (int j = i * d_c; j < (i + 1) * d_c; ++j) {
            block[i][j] = 1;
        }
        H[i] = block[i];
    }

    // Create remaining blocks by permutations of the first block's columns
    for (int i = 1; i < d_v; ++i) {
        std::vector<std::vector<int>> permuted_block = transposeMatrix(block);
        std::shuffle(permuted_block.begin(), permuted_block.end(), rng);
        permuted_block = transposeMatrix(permuted_block);



        /*
        std::cout << "Block:\n";
        for (const auto& row : block) {
            for (int elem : row) {
                std::cout << elem << ' ';
            }
            std::cout << '\n';
        }
        std::cout << "Permuted block:\n";
        for (const auto& row : permuted_block) {
            for (int elem : row) {
                std::cout << elem << ' ';
            }
            std::cout << '\n';
        }
        */



        for (int j = 0; j < block_size; ++j) {
            H[i * block_size + j] = permuted_block[j];
        }
    }

    return H;
}


std::vector<int> encode(const std::vector<std::vector<int>>& tG, const std::vector<int>& v) {

    std::vector<std::vector<int>> v_2;
    v_2.push_back(v);

    std::vector<std::vector<int>> d = binaryproduct(v_2, transposeMatrix(tG));

    return d[0];
}

void checkEncode(const std::vector<int>& enc, const std::vector<std::vector<int>>& H){

    auto check = binaryproduct({enc}, transposeMatrix(H));
    bool ok = true;

    for (const auto& row : check) {
        for (int value : row) {
            if (value != 0) {
                ok = false;  // If any value is not zero, return false
            }
        }
    }

    if (ok){
        std::cout << "Check successful" << std::endl;
    }else{
        std::cout << "Check failed" << std::endl;
    }
}

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
    std::ofstream file("matica.csv");
    if (!file.is_open()) {
        // Handle the case where the file could not be opened
        std::cerr << "Error: Unable to open file matica.csv for writing." << std::endl;
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

void printMatrix(const std::vector<std::vector<int>>& H){
    for (const auto& row : H) {
        for (int elem : row) {
            std::cout << elem << ' ';
        }
        std::cout << '\n';
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
        }else{
            parityCheckMatrix = csvToMatrix(matrixFileName);
        }

        auto codingMatrix = coding_matrix(parityCheckMatrix);

        auto encodedBinary = encode(codingMatrix, binary);

        checkEncode(encodedBinary, parityCheckMatrix);

        if (matrixFileName.empty()){
            MatrixToCsv(parityCheckMatrix);
        }

        if (!encodedBinary.empty()){
            
            std::cout << VectorToBinary(encodedBinary) << std::endl;
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
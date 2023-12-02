/** 
 * BMS 2023
 * Project
 * Implementation of LDPC encoder and decoder
 * by Martin Jacko <xjacko05@stud.fit.vutbr.cz>
 **/

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <random>

enum Mode {
    ENCODE,
    DECODE,
    NONE
};

//print matrix (debug purposes)
void printMatrix(const std::vector<std::vector<int>>& H){

    for (const auto& row : H){
        for (int elem : row){
            std::cout << elem << ' ';
        }
        std::cout << std::endl;
    }
}

//extract n-th column of matrix
std::vector<int> extractNthColumn(const std::vector<std::vector<int>>& input, int n){

    std::vector<int> result;

    for (const auto& row : input){
        result.push_back(row[n]);
    }

    return result;
}

//generate identity matrix of given size
std::vector<std::vector<int>> generateIdentityMatrix(int n){

    //initialize the identity matrix with zeros
    std::vector<std::vector<int>> identityMatrix(n, std::vector<int>(n, 0));

    //set diagonal elements to 1
    for (int i = 0; i < n; ++i){
        identityMatrix[i][i] = 1;
    }

    return identityMatrix;
}

//perform Gauss-Jordan elimination, e.g. compute the binary row reduced echelon form of X
//returns 2 matrices, row reduced form of input and transformations applied to identity matrix
//based on pyldpc/utils.py/gaussjordan
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> gaussJordan(const std::vector<std::vector<int>>& X){

    std::vector<std::vector<int>> A = X;

    //get the number of rows and columns in the matrix
    int m = A.size();
    int n = A[0].size();

    //initialize matrix of transformations
    std::vector<std::vector<int>> P = generateIdentityMatrix(m);

    int pivot_old = -1;
    
    for (int j = 0; j < n; ++j){
        std::vector<int> one_col = extractNthColumn(A, j);

        //extract elements below the previous pivot row
        std::vector<int> filtre_down(one_col.begin() + pivot_old + 1, one_col.end());

        //find the index of the maximum element in the filtered column
        auto max_it = std::max_element(filtre_down.begin(), filtre_down.end());
        int pivot = std::distance(filtre_down.begin(), max_it) + pivot_old + 1;

        //check if the pivot element is non-zero
        if (A[pivot][j]){
            //update pivot_old and swap rows if necessary
            pivot_old += 1;
            if (pivot_old != pivot){
                std::swap(A[pivot], A[pivot_old]);
                std::swap(P[pivot], P[pivot_old]);
            }

            //eliminate non-zero elements in the current column
            for (int i = 0; i < m; ++i){
                if (i != pivot_old && A[i][j]){
                    for (int k = 0; k < n; ++k){
                        A[i][k] = abs(A[i][k] - A[pivot_old][k]);

                        //update transformation matrix
                        if (static_cast<unsigned long>(k) < P[0].size()){
                            P[i][k] = abs(P[i][k] - P[pivot_old][k]);
                        }
                    }
                    //if the matrix has more rows than columns, update remaining elements in P
                    if (n < m){
                        for (int k = n; k < m; ++k){
                            P[i][k] = abs(P[i][k] - P[pivot_old][k]);
                        }
                    }
                }
            }
        }

        //check if the whole matrix was processed
        if (pivot_old == m - 1){
            break;
        }
    }

    return std::make_pair(A, P);
}

//perform Gauss elimination on system of linear equations
//returns a pair of matrices, first matrix is the upper triangular form of input, second matrix is modified vector b
//based on pyldpc/utils.py/gausselimination
std::pair<std::vector<std::vector<int>>, std::vector<int>> gaussElimination(const std::vector<std::vector<int>>& A, const std::vector<int>& b){

    std::vector<std::vector<int>> A_copy = A;
    std::vector<int> b_copy = b;

    //get the number of rows and columns in the matrix
    int n = A_copy.size();
    int k = A_copy[0].size();

    for (int j = 0; j < std::min(k, n); ++j){
        //potential pivot rows
        std::vector<int> listedepivots;
        //find non-zero elements in the current column below the diagonal
        for (int i = j; i < n; ++i){
            if (A_copy[i][j]){
                listedepivots.push_back(i);
            }
        }

        //check if current column contains non-zero elements
        if (!listedepivots.empty()){
            int pivot = *std::min_element(listedepivots.begin(), listedepivots.end());

            //swap rows if the pivot row is not the current row
            if (pivot != j){
                std::swap(A_copy[j], A_copy[pivot]);
                std::swap(b_copy[j], b_copy[pivot]);
            }

            //eliminate non-zero elements below the pivot in the current column
            for (int i = j + 1; i < n; ++i){
                if (A_copy[i][j]){
                    //update A matrix and b vector
                    for (int col = 0; col < k; ++col){
                        A_copy[i][col] = abs(A_copy[i][col] - A_copy[j][col]);
                    }
                    b_copy[i] = abs(b_copy[i] - b_copy[j]);
                }
            }
        }
    }

    return std::make_pair(A_copy, b_copy);
}

//matrix multiplication with dimensions checking
std::vector<std::vector<int>> binaryproduct(const std::vector<std::vector<int>>& matrix1, const std::vector<std::vector<int>>& matrix2){

    //check if matrices are compatible for multiplication
    if (matrix1.empty() || matrix2.empty() || matrix1[0].size() != matrix2.size()){
        // Matrices cannot be multiplied
        std::cerr << "Error: Matrices cannot be multiplied. Matrix 1: " << matrix1.size() << "x" << matrix1[0].size() << " Matrix 2: " << matrix2.size() << "x" << matrix2[0].size() << std::endl;
        return {};
    }

    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int cols2 = matrix2[0].size();

    //initialize the result matrix with zeros
    std::vector<std::vector<int>> result(rows1, std::vector<int>(cols2, 0));

    //perform matrix multiplication
    for (int i = 0; i < rows1; ++i){
        for (int j = 0; j < cols2; ++j){
            for (int k = 0; k < cols1; ++k){
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
            //apply mod 2 to the result as we are working in mod 2 domain
            result[i][j] %= 2;
        }
    }

    return result;
}

//matrix transposition
std::vector<std::vector<int>> transposeMatrix(const std::vector<std::vector<int>>& matrix){

    // Get the number of rows and columns in the original matrix
    int rows = matrix.size();
    int cols = matrix[0].size();

    // Initialize the transposed matrix with zeros
    std::vector<std::vector<int>> transposedMatrix(cols, std::vector<int>(rows, 0));

    // Perform the transpose operation
    for (int i = 0; i < rows; ++i){
        for (int j = 0; j < cols; ++j){
            transposedMatrix[j][i] = matrix[i][j];
        }
    }

    return transposedMatrix;
}

//generate a pseudo-random regular parity-check matrix
//based on pyldpc/code.py/parity_check_matrix
std::vector<std::vector<int>> generateParityCheckMatrix(int length){

    //initialize pseudo-random generator with current time
    int seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);

    //values as pre assignment
    int n_code = 2 * length;
    int d_v = length - 1;
    int d_c = length;

    int n_equations = (n_code * d_v) / d_c;

    std::vector<std::vector<int>> block(n_equations / d_v, std::vector<int>(n_code, 0));
    std::vector<std::vector<int>> H(n_equations, std::vector<int>(n_code));
    int block_size = n_equations / d_v;

    //filling the first block with consecutive ones in each row of the block
    for (int i = 0; i < block_size; ++i){
        for (int j = i * d_c; j < (i + 1) * d_c; ++j){
            block[i][j] = 1;
        }
        H[i] = block[i];
    }

    //create remaining blocks by permutations of the first block's columns
    for (int i = 1; i < d_v; ++i){
        std::vector<std::vector<int>> permuted_block = transposeMatrix(block);
        std::shuffle(permuted_block.begin(), permuted_block.end(), rng);
        permuted_block = transposeMatrix(permuted_block);

        for (int j = 0; j < block_size; ++j){
            H[i * block_size + j] = permuted_block[j];
        }
    }

    return H;
}

//generate LDPC transposed coding matrix from parity-check matrix using Gauss elimination
//based on pyldpc/code.py/coding_matrix
std::vector<std::vector<int>> generateCodingMatrix(const std::vector<std::vector<int>>& H){

    //get the number of encoded bits
    int n_code = H[0].size();

    //double Gauss elimination
    auto [Href_colonnes, tQ] = gaussJordan(transposeMatrix(H));
    auto Href_diag = gaussJordan(transposeMatrix(Href_colonnes)).first;

    auto Q = transposeMatrix(tQ);

    //calculate the number of information bits in the code
    int n_bits = n_code - static_cast<int>(std::accumulate(Href_diag.begin(), Href_diag.end(), 0, 
        [](int sum, const std::vector<int>& row){
            return sum + std::accumulate(row.begin(), row.end(), 0);
        }));

    //initialize a matrix Y used in construction of coding matrix
    std::vector<std::vector<int>> Y(n_code, std::vector<int>(n_bits, 0));
    for (int i = n_code - n_bits; i < n_code; ++i){
        Y[i][i - (n_code - n_bits)] = 1;
    }

    //obtain transposed coding matrix
    auto tG = binaryproduct(Q, Y);

    return tG;
}


//encoding done by multiplication of input and coding matrix
std::vector<int> encode(const std::vector<std::vector<int>>& parityCheckMatrix, const std::vector<int>& binary){

    auto tG = generateCodingMatrix(parityCheckMatrix);

    std::vector<std::vector<int>> d = binaryproduct({binary}, transposeMatrix(tG));

    return d[0];
}


//compute original message from corrected encoded data
//based on pyldpc/encoder.py/get_message
std::vector<int> getMessage(const std::vector<std::vector<int>>& tG, const std::vector<int>& binary){

    //get the number of bits in the original message
    int n_bits = tG[0].size();

    auto [rtG, rx] = gaussElimination(tG, binary);

    std::vector<int> message(n_bits, 0);

    //perform back-substitution to compute the message
    message[n_bits - 1] = rx[n_bits - 1];
    for (int i = n_bits - 2; i >= 0; --i){
        message[i] = rx[i];
        for (int j = i + 1; j < n_bits; ++j){
            //update the current bit by subtracting the product of the upper triangular matrix and the message
            message[i] -= rtG[i][j] * message[j];
        }
        //enforce mod 2 domain
        message[i] %= 2;
    }

    //force all bits to be positive
    for (int& bit : message){
        bit = abs(bit);
    }

    return message;
}

//check if given input is codeword by checking if product of input and parity check matrix is equal to zero
bool isCodeWord(const std::vector<int>& enc, const std::vector<std::vector<int>>& H){

    auto product = binaryproduct({enc}, transposeMatrix(H));
    bool check = true;

    //check every element of product matrix
    for (const auto& row : product){
        for (int value : row){
            if (value != 0){
                check = false;
            }
        }
    }

    return check;
}

//converts filtered program ASCII input into binary vector form
std::vector<int> AsciiToVector(const std::string& input){

    std::vector<int> binaryASCII;

    for (char ch : input){
        //ignore non-alphanumeric characters
        if (std::isalnum(ch)){
            //convert the character to its ASCII value
            int asciiValue = static_cast<int>(ch);

            int binary[8] = {0,0,0,0,0,0,0,0};
            int j = 0;

            while (asciiValue > 0){ 

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

//converts filtered program binary input into binary vector form
std::vector<int> BinaryToVector(const std::string& input){

    std::vector<int> binaryVector;

    for (char ch : input){
        //ignore non-binary characters
        if (ch == '0' || ch == '1'){
            binaryVector.push_back(ch - '0');
        }
    }

    return binaryVector;
}

//convert binary vector to binary string
std::string VectorToBinary(const std::vector<int>& inputVector){

    std::stringstream ss;

    for (int value : inputVector){
        ss << value;
    }

    return ss.str();
}

//convert binary vector to ASCII string
std::string VectorToAscii(const std::vector<int>& binaryVector){

    size_t size = binaryVector.size();

    std::string asciiString;
    for (size_t i = 0; i < size; i += 8){
        int asciiValue = 0;
        for (size_t j = 0; j < 8; ++j){
            asciiValue = (asciiValue << 1) | binaryVector[i + j];
        }
        asciiString.push_back(static_cast<char>(asciiValue));
    }

    return asciiString;
}

//load matrix from a given file
std::vector<std::vector<int>> csvToMatrix(const std::string& filename){

    std::vector<std::vector<int>> matrix;

    std::ifstream file(filename);
    if (!file.is_open()){
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return matrix;
    }

    std::string line;
    while (std::getline(file, line)){
        std::istringstream lineStream(line);
        std::vector<int> row;

        //split each cell using both comma and semicolon as delimiters
        std::string cell;
        while (std::getline(lineStream, cell, ',')){
            if (cell.empty()){
                continue;
            }
            std::istringstream cellStream(cell);
            std::string value;
            while (std::getline(cellStream, value, ';')){
                if (value.empty()){
                    continue;
                }
                try {
                    int intValue = std::stoi(value);
                    row.push_back(intValue);
                } catch (...){
                    //cells with non-integer values are treated as zero
                    std::cerr << "Error: Invalid value in file " << filename << ": " << value << ", value treated as 0" << std::endl;
                    row.push_back(0);
                }
            }
        }

        matrix.push_back(row);
    }

    return matrix;
}

//save generated parity check matrix to file in csv format
void MatrixToCsv(const std::vector<std::vector<int>>& data){

    std::ofstream file("matica.csv");
    if (!file.is_open()){
        std::cerr << "Error: Unable to open file matica.csv for writing." << std::endl;
        return;
    }

    for (const auto& row : data){
        for (size_t i = 0; i < row.size(); ++i){
            file << row[i];
            if (i < row.size() - 1){
                file << ',';
            }
        }
        file << '\n';
    }
}

//compute and perform minimal adjustments to ASCII character so that it will become alphanumeric
std::vector<int> minimalBitFlipsChar(const std::vector<int>& binaryRepresentation){

    //convert the binary representation to ASCII character
    char originalChar = 0;
    for (int i = 7; i >= 0; --i){
        originalChar |= binaryRepresentation[i] << (7 - i);
    }

    //check if the character is already alphanumeric
    if (std::isalnum(originalChar)){
        return binaryRepresentation;
    }

    // Find the closest alphanumeric character with the least number of bit flips
    char closestAlphanumeric = originalChar;
    int minFlips = 8;


    for (char target = '0'; target <= '9'; ++target){
        int flips = 0;
        for (int i = 0; i < 8; ++i){
            if (((target >> i) & 1) != ((originalChar >> i) & 1)){
                flips++;
            }
        }

        if (flips < minFlips){
            minFlips = flips;
            closestAlphanumeric = target;
        }
    }

    for (char target = 'A'; target <= 'Z'; ++target){
        int flips = 0;
        for (int i = 0; i < 8; ++i){
            if (((target >> i) & 1) != ((originalChar >> i) & 1)){
                flips++;
            }
        }

        if (flips < minFlips){
            minFlips = flips;
            closestAlphanumeric = target;
        }
    }

    for (char target = 'a'; target <= 'z'; ++target){
        int flips = 0;
        for (int i = 0; i < 8; ++i){
            if (((target >> i) & 1) != ((originalChar >> i) & 1)){
                flips++;
            }
        }

        if (flips < minFlips){
            minFlips = flips;
            closestAlphanumeric = target;
        }
    }

    //save found character in binary form
    std::vector<int> closestBinaryRepresentation(8, 0);
    for (int i = 7; i >= 0; --i){
        closestBinaryRepresentation.push_back((closestAlphanumeric >> i) & 1);
    }

    return closestBinaryRepresentation;
}

//compute and perform minimal adjustments to ASCII string so that it will become alphanumeric
std::vector<int> minimalBitFlipsString(const std::vector<int>& asciiString){

    std::vector<int> result;

    //process consecutive 8-element subvectors
    for (size_t i = 0; i < asciiString.size(); i += 8){
        //get the subvector of size 8
        std::vector<int> subVector(asciiString.begin() + i, asciiString.begin() + i + 8);

        std::vector<int> subResult = minimalBitFlipsChar(subVector);

        result.insert(result.end(), subResult.begin(), subResult.end());
    }

    return result;
}

//check if given value represents alphanumeric ASCII character
bool checkAscii(const std::vector<int>& binaryVector){

    size_t size = binaryVector.size();

    for (size_t i = 0; i < size; i += 8){
        int asciiValue = 0;
        for (size_t j = 0; j < 8; ++j){
            asciiValue = (asciiValue << 1) | binaryVector[i + j];
        }
        if (!std::isalnum(static_cast<char>(asciiValue))){
            return false;
        }
    }

    return true;
}


//implementation of Hard-decision decoding using bit nodes and check nodes in Tanner graph 
//inspired by section 3.1 of https://www.bernh.net/media/download/papers/ldpc.pdf
std::vector<int> hardDecision(const std::vector<std::vector<int>>& parityCheckMatrix, const std::vector<int>& codeword){

    //get the size of check and bit nodes
    int checkNodeSize = parityCheckMatrix.size();
    int bitNodeSize = parityCheckMatrix[0].size();

    //coding matrix is required for decoding corrected data
    auto codingMatrix = generateCodingMatrix(parityCheckMatrix);

    //initialize the estimate with the received codeword
    std::vector<int> estimate = codeword;

    //variables for message passing
    std::vector<std::vector<int>> msgCheckToBit(checkNodeSize, std::vector<int>(bitNodeSize, 0));
    std::vector<int> flipped(bitNodeSize, 0);

    //iterative correcting
    for (int iter = 0; iter < 50; ++iter){

        //check if current estimate is a valid codeword
        if (isCodeWord(estimate, parityCheckMatrix)){
            //retrieve the original message from the estimate
            auto message = getMessage(codingMatrix, estimate);
            //check if the decoded message consists of valid alphanumeric ASCII characters
            if (checkAscii(message)){
                return message;
            }
        }

        //score table, keeps track of numver of conflicts
        std::vector<int> score(bitNodeSize, 0);
        
        //update messages from check nodes to bit nodes
        for (int i = 0; i < checkNodeSize; ++i){
            for (int j = 0; j < bitNodeSize; ++j){
                //kip if the bit has already been flipped 2 times (2 beacause, original value is preferred)
                if (flipped[j] == 2){
                    continue;
                }
                //check if the check node is connected to the bit node
                if (parityCheckMatrix[i][j]){
                    //compute sum of other bits connected to the check node
                    int sum = 0;
                    for (int k = 0; k < bitNodeSize; ++k){
                        if (parityCheckMatrix[i][k]){
                            if (k != j){
                                sum += estimate[k];
                            }
                        }
                    }
                    sum %= 2;
                    //increase conflict count if actual bit does not match expected one
                    if (estimate[j] != sum){
                        score[j]++;
                    }
                }
            }
        }

        //flip bit with highest score
        auto max_it = std::max_element(score.begin(), score.end());
        int index = std::distance(score.begin(), max_it);
        if (estimate[index] == 0){
            estimate[index] = 1;
        }else{
            estimate[index] = 0;
        }
        flipped[index]++;
        
    }

    //if the correcting process did not yield valid ASCII codeword, attempt to decode anyways
    auto message = getMessage(codingMatrix, codeword);
    //check if decoded message consists only of alphanumeric ASCII characters
    if (!checkAscii(message)){
        //correct no alphanumeric ASCII characters to closest viable alternative (least bits flipped required to do so)
        return minimalBitFlipsString(message);
    }

    return message;
}


//main function, handles primarily program flow and argument checking
int main(int argc, char *argv[]){

    std::string matrixFileName;
    Mode mode;

    //argument checking
    for (int i = 1; i < argc; ++i){
        std::string arg = argv[i];

        if (arg == "-m" && i + 1 < argc){
            matrixFileName = argv[i + 1];
            i++;
        } else if (arg == "-e"){
            mode = Mode::ENCODE;
        } else if (arg == "-d"){
            mode = Mode::DECODE;
        }
    }

    if (matrixFileName.empty() && mode == Mode::DECODE){
        std::cerr << "Error: Decoding requires -m argument" << std::endl;
        return 1;
    }

    //reading stdin
    std::string input;
    getline(std::cin, input);

    std::vector<std::vector<int>> parityCheckMatrix;

    if (mode == Mode::ENCODE){

        //transforming input to vector
        auto binary = AsciiToVector(input);

        //parity-check matrix generation
        if (matrixFileName.empty()){
            parityCheckMatrix = generateParityCheckMatrix(binary.size());
            MatrixToCsv(parityCheckMatrix);
        //parity-check matrix load
        }else{
            parityCheckMatrix = csvToMatrix(matrixFileName);
            //input and matrix size compatibility check
            if (parityCheckMatrix[0].size() / 2 != binary.size()){
                std::cerr << "Error: Incompatible matrix and input sizes" << std::endl;
                return 1;
            }
        }

        //encoding
        auto encodedBinary = encode(parityCheckMatrix, binary);

        //outputting vector as binary string
        std::cout << VectorToBinary(encodedBinary) << std::endl;
        
    }else if (mode == Mode::DECODE){

        //transforming input to vector
        auto binary = BinaryToVector(input);

        //parity-check matrix load
        parityCheckMatrix = csvToMatrix(matrixFileName);

        //decoding using hard decision bit flipping
        auto decodedBinary = hardDecision(parityCheckMatrix, binary);
        
        //outputting vector as ASCII string
        std::cout << VectorToAscii(decodedBinary) << std::endl;
    }

    return 0;
}
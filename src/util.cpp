#include "util.h"

void setDims(matrix& mat, const int numRow, const int numCol)
{
    mat.resize(numRow);
    for(int i = 0; i < numRow; i++)
    {
        mat[i].resize(numCol);
    }
}

void randMatrix(matrix& matrix, int numRow, int numCol, uint64_t modulus) 
{
    // Seed for the random number generator
    random_device rd;
    mt19937 gen(rd());

    // Define the distribution for random integers
    uniform_int_distribution<int> distribution(0, modulus - 1); // You can adjust the range as needed

    matrix.resize(numRow);
    for (int i = 0; i < numRow; ++i) 
    {
        matrix[i].resize(numCol);
        for (int j = 0; j < numCol; ++j) 
        {
            matrix[i][j] = distribution(gen);
        }
    }
}


void randMatrix(ringMatrix& matrix, int numRow, int numCol, int degree, uint64_t modulus) 
{
    // Seed for the random number generator
    random_device rd;
    mt19937 gen(rd());

    // Define the distribution for random integers
    uniform_int_distribution<uint64_t> distribution(0, modulus - 1); // Use uint64_t for distribution range

    matrix.resize(numRow);
    for (int i = 0; i < numRow; ++i) 
    {
        matrix[i].resize(numCol);
        for (int j = 0; j < numCol; ++j) 
        {
            poly p(degree);
            for (int k = 0; k < degree; ++k) 
            {
                p[k] = distribution(gen);
            }
            matrix[i][j] = p;
        }
    }
}


void randVector(vector<int64_t>& vector, int length, uint64_t modulus)
{
    // Seed for the random number generator
    random_device rd;
    mt19937 gen(rd());

    // Define the distribution for random integers
    uniform_int_distribution<uint64_t> distribution(0, modulus - 1); // You can adjust the range as needed

    vector.resize(length);
    for (int i = 0; i < length; ++i) 
    {
        vector[i] = distribution(gen);
    }
}


// Function to perform matrix addition over Z_q
void matrixAdd(const matrix& matrix1, const matrix& matrix2, int modulus, matrix& resultMatrix) 
{
    int numRows = matrix1.size();
    int numCols = matrix1[0].size();

    // Check if matrices have the same dimensions
    if (numRows != matrix2.size() || numCols != matrix2[0].size()) 
    {
        cerr << "Error: Matrices must have the same dimensions for addition." << endl;
        exit(EXIT_FAILURE);
    }

    // Resize the result matrix
    resultMatrix.resize(numRows);
    for (int i = 0; i < numRows; ++i) 
    {
        resultMatrix[i].resize(numCols);
    }

    // Perform matrix addition over Z_q
    for (int i = 0; i < numRows; ++i) 
    {
        for (int j = 0; j < numCols; ++j)
        {
            resultMatrix[i][j] = (matrix1[i][j] + matrix2[i][j]) % modulus;
        }
    }
}



void matrixMultiply(const matrix& matrix1, const matrix& matrix2, const uint64_t modulus, matrix& resultMatrix) 
{
    int numRows1 = matrix1.size();
    int numCols1 = matrix1[0].size();
    int numRows2 = matrix2.size();
    int numCols2 = matrix2[0].size();

    // Check if matrices are compatible for multiplication
    if (numCols1 != numRows2) {
        cerr << "Error: Incompatible matrix dimensions for multiplication." << endl;
        cerr << "column of matrix 1: " << numCols1 << endl;
        cerr << "row of matrix 2: " << numRows2 << endl;
        exit(EXIT_FAILURE);
    }

    // Resize the result matrix
    resultMatrix.resize(numRows1, vector<int64_t>(numCols2));
    // for (int i = 0; i < numRows1; ++i) 
    // { 
    //     resultMatrix[i].resize(numCols2);
    // }

    // Perform matrix multiplication over Z_q
    for (int i = 0; i < numRows1; ++i) 
    {
        for (int j = 0; j < numCols2; ++j) 
        {
            resultMatrix[i][j] = 0;
            for (int k = 0; k < numCols1; ++k) 
            {
                // resultMatrix[i][j] = (resultMatrix[i][j] + matrix1[i][k] * matrix2[k][j]) % modulus;
                resultMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
            resultMatrix[i][j] %= modulus;
        }
    }
}


void matrixMultiply(const matrix& matrix_input, const vector<int64_t>& vector_input, uint64_t modulus, vector<int64_t>& resultVector) {
    int numRows = matrix_input.size();
    int numCols = matrix_input[0].size();

    // Check if matrix and vector dimensions are compatible
    if (numCols != vector_input.size()) 
    {
        // cerr << "Error: Incompatible matrix and vector dimensions for multiplication." << endl;
        // exit(EXIT_FAILURE);
        cerr << "Error: Incompatible matrix and vector dimensions for multiplication." << endl;
        cerr << "column of matrix: " << numCols << endl;
        cerr << "row of vector: " << vector_input.size() << endl;
        exit(EXIT_FAILURE);
    }

    // Resize the result vector
    resultVector.resize(numRows, 0);

    // Perform matrix-vector multiplication over Z_q
    for (int i = 0; i < numRows; ++i) 
    {
        for (int j = 0; j < numCols; ++j) 
        {
            resultVector[i] = (resultVector[i] + matrix_input[i][j] * vector_input[j]) % modulus;
        }
    }
}


// Function to display a matrix
void print(const matrix& matrix) 
{
    for (const auto& row : matrix) 
    {
        std::cout << "[";
        for (uint64_t value : row) 
        {
            cout << value << " ";
        }
        std::cout << "]" << std::endl;
    }
}


void print(const ringMatrix& matrix) 
{
    for (const auto& row : matrix) 
    {
        std::cout << "[";
        for (const auto& poly_vector : row) 
        {
            std::cout << "[";
            for (const auto& value : poly_vector) 
            {
                std::cout << value << " ";
            }
            std::cout << "]";
        }
        std::cout << "]" << std::endl;
    }
}

void print(const poly vector)
{
    std::cout << "[";
    for (const auto& value : vector) {
        std::cout << value << " ";
    }
    std::cout << "]" << std::endl;
}

// Function to check if a number is probably prime using the Miller-Rabin primality test
bool isProbablePrime(uint64_t n, int k = 5) {
    if (n <= 1 || n == 4)
        return false;
    if (n <= 3)
        return true;

    // Find r such that n = 2^d * r + 1 for some r >= 1
    uint64_t d = n - 1;
    while (d % 2 == 0)
        d /= 2;

    // Witness loop
    for (int i = 0; i < k; ++i) {
        uint64_t a = 2 + rand() % (n - 3);
        uint64_t x = static_cast<uint64_t>(pow(a, d)) % n;

        if (x == 1 || x == n - 1)
            continue;

        while (d != n - 1) {
            x = (x * x) % n;
            d *= 2;

            if (x == 1)
                return false;
            if (x == n - 1)
                break;
        }

        if (x != n - 1)
            return false;
    }

    return true;
}

// Function to generate an n-bit prime number
void generatePrime(uint64_t& prime, int n) {
    srand(static_cast<unsigned int>(time(nullptr)));

    uint64_t lowerBound = 1LL << (n - 1);
    uint64_t upperBound = (1LL << n) - 1;

    do {
        prime = lowerBound + rand() % (upperBound - lowerBound + 1);
    } while (!isProbablePrime(prime));

}


int64_t modExp(int64_t base, int64_t exp, int64_t mod) 
{
    int64_t result = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1) // If exp is odd, multiply base with result
            result = (result * base) % mod;
        exp = exp >> 1; // exp = exp / 2
        base = (base * base) % mod; // base = base^2
    }
    return result;
}


void ntt(poly& a, const uint64_t modulus, const int64_t root, bool invert) 
{
    size_t n = a.size();
    assert((n & (n - 1)) == 0); // n must be a power of 2

    for (size_t i = 1, j = 0; i < n; ++i) 
    {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (int len = 2; len <= n; len <<= 1) 
    {
        int64_t wlen = invert ? modExp(root, (modulus - 1 - (modulus - 1) / len), modulus) : modExp(root, (modulus - 1) / len, modulus);
        for (int i = 0; i < n; i += len) 
        {
            int64_t w = 1;
            for (int j = 0; j < len / 2; ++j) 
            {
                int64_t u = a[i + j];
                int64_t v = (a[i + j + len / 2] * w) % modulus;
                a[i + j] = (u + v) % modulus;
                a[i + j + len / 2] = (u + modulus - v) % modulus;
                w = (w * wlen) % modulus;
            }
        }
    }

    if (invert) 
    {
        int64_t n_inv = modExp(n, modulus - 2, modulus);
        for (int64_t &x : a) {
            x = (x * n_inv) % modulus;
        }
    }
}


// void multiply_ntt(const poly& src1, const poly &src2, poly& dest, const uint64_t modulus, const int64_t root) 
// {
//     size_t n = 1;
//     while (n < src1.size() + src2.size())    
//         n <<= 1;
//     poly fa(src1.begin(), src1.end()), fb(src2.begin(), src2.end());
//     fa.resize(n);
//     fb.resize(n);

//     ntt(fa, modulus, root, false);
//     ntt(fb, modulus, root, false);
//     for (size_t i = 0; i < n; ++i) 
//     {
//         fa[i] = (fa[i] * fb[i]) % modulus;
//     }

//     ntt(fa, modulus, root, true);
//     dest.resize(n / 2);
//     for (size_t i = 0; i < n / 2; ++i) 
//     {
//         dest[i] = (fa[i] + modulus - (i + n / 2 < fa.size() ? fa[i + n / 2] : 0)) % modulus;
//     }
// }

void multiply_ntt(const poly& src1, const poly &src2, poly& dest, const uint64_t modulus, const uint64_t root, bool is_ntt_form) 
{
    size_t n;
    poly fa, fb;

    if(!is_ntt_form)
    {
        n = 1;
        while (n < src1.size() + src2.size())    
            n <<= 1;
        fa = poly(src1.begin(), src1.end());
        fb = poly(src2.begin(), src2.end());
        fa.resize(n);
        fb.resize(n);
        ntt(fa, modulus, root, false);
        ntt(fb, modulus, root, false);
    }
    else
    {
        n = src1.size();
        fa = src1;
        fb = src2;
    } 

    for (size_t i = 0; i < n; ++i) 
    {
        fa[i] = (fa[i] * fb[i]) % modulus;
    }

    ntt(fa, modulus, root, true);
    dest.resize(n / 2);
    for (size_t i = 0; i < n / 2; ++i) 
    {
        dest[i] = (fa[i] + modulus - (i + n / 2 < fa.size() ? fa[i + n / 2] : 0)) % modulus;
    }
}



void invert_ntt(const poly& src, poly& dest, const uint64_t modulus, const int64_t root)
{
    int n = src.size();
    poly tmp = src;
    dest.resize(n/2, 0);

    ntt(tmp, modulus, root, true);
    // We are in the ring Z_p[x] / <x^n + 1>, so reduce modulo x^n + 1
    for (size_t i = 0; i < dest.size(); ++i) 
    {
        dest[i] = (tmp[i] + modulus - (i + dest.size() < tmp.size() ? tmp[i + dest.size()] : 0)) % modulus;
    }
}
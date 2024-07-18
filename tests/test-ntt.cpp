#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>

#include "util.h"

#define __DEBUG 0

// Function to compute n-th root of unity
int64_t get_root_of_unity(int64_t n, int64_t mod, int64_t primitive_root) 
{
    return modExp(primitive_root, (mod - 1) / n, mod);
}

// Naive polynomial multiplication function
void multiply_naive(const poly& src1, const poly& src2, poly& dest, int64_t modulus) 
{
    size_t n = src1.size();
    poly result(2 * n - 1, 0);
    result.resize(2 * n - 1, 0);
    for (size_t i = 0; i < n; ++i) 
    {
        for (size_t j = 0; j < n; ++j) 
        {
            result[i + j] = (result[i + j] + src1[i] * src2[j]) % modulus;
        }
    }
    // Reduce modulo x^n + 1
    // poly final_result(n, 0);
    dest.resize(n, 0);
    for (size_t i = 0; i < n; ++i) 
    {
        dest[i] = (result[i] + modulus - (i + n < result.size() ? result[i + n] : 0)) % modulus;
    }
}

// int64_t mod = 998244353; // A prime number (2^23 * 7 * 17 + 1) 

int main() 
{
    uint64_t mod = 998244353;
    // uint64_t mod = 17; // p = k 2^m + 1
    // 2m-th root of unity
    int64_t primitive_root = 3; // A primitive root modulo mod
    int64_t degree = 32; // Degree of polynomials (n-1)
    // int64_t root = get_root_of_unity(degree * 2, mod, primitive_root);
    // int64_t root = 3;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> input_distribution(0, mod - 1); 

    poly input1, input2, ntt_input1, ntt_input2, result_ntt_1, result_ntt_2, result_naive;
    for(size_t i = 0; i < degree; i++) 
    {
        input1.push_back(input_distribution(gen));
        input2.push_back(input_distribution(gen));
    }

    // Print inputs
    std::cout << "\nParameter"  << std::endl;
    std::cout << " - modulus: " << mod << std::endl;
    std::cout << " - degree: " << degree << std::endl;
    std::cout << " - primitive root: " << primitive_root << std::endl;
    // std::cout << " - root: " << root << std::endl;

    ntt_input1 = input1; 
#if __DEBUG == 1
    std::cout << "Input1: "; print(input1);
#endif
    ntt(ntt_input1, mod, primitive_root);
    // ntt(ntt_input1, mod, root);
#if __DEBUG == 1
    std::cout << " -> NTT: "; print(ntt_input1);
#endif
    ntt(ntt_input1, mod, primitive_root, true);
    // ntt(ntt_input1, mod, root, true);
#if __DEBUG == 1
    std::cout << " -> iNTT: "; print(ntt_input1);   
#endif
    // Check if results match
    assert(input1 == ntt_input1);
    std::cout << "> Test passed: NTT and iNTT works well!" << std::endl;

    multiply_naive(input1, input2, result_naive, mod);
    multiply_ntt(input1, input2, result_ntt_1, mod, primitive_root);

    // transform ntt form
    size_t n = 1;
    while(n < 2 * degree)
        n <<= 1;
    ntt_input1 = poly(input1.begin(), input1.end());
    ntt_input2 = poly(input2.begin(), input2.end());
    ntt_input1.resize(n);
    ntt_input2.resize(n);

    ntt(ntt_input1, mod, primitive_root);
    ntt(ntt_input2, mod, primitive_root);

    multiply_ntt(ntt_input1, ntt_input2, result_ntt_2, mod, primitive_root, true);

    
#if __DEBUG == 2
    std::cout << "\nInput1: "; print(input1);
    std::cout << "Input2: "; print(input2);
    // Print results
    std::cout << " -> NTT Multiplication: "; print(result_ntt_1);
    std::cout << " -> Naive Multiplication: "; print(result_naive);
#endif
    
    // Check if results match
    assert(result_ntt_1 == result_naive);
    assert(result_ntt_1 == result_ntt_2);
    std::cout << "> Test passed: multiplication results match!" << std::endl;

    return 0;
}
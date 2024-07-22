#ifndef __UTIL
#define __UTIL
#include <cstdlib>
#include <iostream>
#include <vector>
#include <random>
#include <cstdint>
#include <cmath>
#include <cassert>

#include "omp.h"

using namespace std;

typedef vector<int64_t> poly;
typedef vector<vector<int64_t> > matrix;
typedef vector<vector<poly> > ringMatrix;

void setDims(matrix& matrix, const int numRow, const int numCol);
void randMatrix(matrix& matrix, int numRow, int numCol, uint64_t modulus);
void randMatrix(ringMatrix& matrix, int numRow, int numCol, int degree, uint64_t modulus);
void randVector(vector<int64_t>& vector, int length, uint64_t modulus);

void matrixAdd(const matrix& matrix1, const matrix& matrix2, int modulus, matrix& resultMatrix);
void matrixMultiply(const matrix& matrix1, const matrix& matrix2, const uint64_t modulus, matrix& resultMatrix);
void matrixMultiply(const matrix& matrix_input, const vector<int64_t>& vector_input, uint64_t modulus, vector<int64_t>& resultVector);

void print(const matrix& matrix);
void print(const ringMatrix& matrix);
void print(const poly vector);

bool isProbablePrime(uint64_t n, int k);
void generatePrime(uint64_t& prime, int n);
int generateDiscreteGaussian(int mean, double stddev);

int64_t modExp(int64_t base, int64_t exp, int64_t mod);
void ntt(poly& a, const uint64_t modulus, const int64_t root, bool invert = false);
void multiply_ntt(const poly& src1, const poly &src2, poly& dest, const uint64_t modulus, const uint64_t root, bool is_ntt_form = false);
void invert_ntt(const poly& src, poly& dest, const uint64_t modulus, const int64_t root);

#endif
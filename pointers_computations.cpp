//
// Created by edaravig on 26/1/2017.
//
#include <iostream>
#include <complex>

using namespace std;

/**
 * Prints a row vector of order equal to size
 * @param vector: the vector to be printed
 * @param size: the order of the vector
 */
void print_vector(double *vector, const size_t size) {
    cerr<<"[";
    for(int i=0; i<size; ++i) {
        cerr<<vector[i]<<endl;
    }
    cerr<<"]"<<endl;
}

/**
 * Prints a matrix of order equal to size
 * @param matrix: the matrix to be printed
 * @param size: the order of the matrix
 */
void print_matrix(double **matrix, const size_t size) {
    cerr<<"[";
    for(int i=0; i<size; ++i) {
        cerr<<"[";
        for(int j=0; j<size; ++j) {
            cerr<<matrix[i][j]<<" ";
        }
        cerr<<"]"<<endl;
    }
    cerr<<"]"<<endl;
}

/**
 * Addition of two row vectors of order equal to size
 * @param x: the first vector of the summation
 * @param y: the second vector of the summation
 * @param size: the order of the vectors
 * @return add: the summation of the two vectors (add = x+y)
 */
double* addition(const double *x, const double *y, const size_t size) {
    double *add = new double[size];
    for(int i=0; i<size; ++i) {
        add[i] = x[i] + y[i];
    }
    return add;
}

/**
 * Subtraction of two row vectors of order equal to size
 * @param x: the vector from which the elements are subtracted
 * @param y: the vector of which the elements values are subtracted
 * @param size: the order of the vectors
 * @return sub: the subtraction of the two vector (sub = x-y)
 */
double* subtraction(const double *x, const double *y, const size_t size) {
    double *sub = new double[size];
    for(int i=0; i<size; ++i) {
        sub[i] = x[i] - y[i];
    }
    return sub;
}

/**
 * Multiplication of the vector with a number
 * @param x: the vector which elements are multiplied
 * @param num: the number of the multiplication
 * @param size: the order of the vectors
 * @return mul_n: the multiplied with number vector (mul_n = num*x)
 */
double* multiply_with_number(const double* x, const double num, const size_t size) {
    double *mul_n = new double[size];
    for(int i=0; i<size; ++i) {
        mul_n[i] = num * x[i];
    }
    return mul_n;
}

/**
 * Transpose multiplication of two vectors: transposed with dimensions (1xn) and vector with dimensions (nx1).
 * This results always to a number, since the other-way multiplication -resulting in a matrix of (nxn)- is not supported
 * @param vector: the vector of dimensions (nx1)
 * @param transposed: the vector of dimensions (1xn)
 * @param size: the order of the vectors
 * @return value: the value of the row vector and the column vector (value = r'*r)
 */
double transpose_multiplication(const double *vector, const double *transposed, const size_t size) {
    double value = 0;
    for(int i=0; i<size; ++i) {
        value += vector[i] * transposed[i];
    }
    return value;
}

/**
 * Vector norm computes and returns the euclidean norm, or norm of size 2 of the input vector
 * @param x: the input vector
 * @param size: the order of the vector
 * @return norm: the euclidean norm
 */
double vector_norm_2(const double *x, const size_t size) {
    double norm = 0;
    for(int i=0; i<size; ++i) {
        norm += x[i]*x[i];
    }
    return sqrt(norm);
}
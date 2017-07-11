#include <iostream>
#include <functional>
#include <fstream>
#include <vector>
#include <iterator>
#include <random>
#include <cfloat>

using namespace std;

/*these functions may not need to be here, an include would be enough*/
void print_vector(double *vector, const size_t size);
void print_matrix(double **matrix, const size_t size);
double* addition(const double *x, const double *y, const size_t size);
double* subtraction(const double *x, const double *y, const size_t size);
double* multiply_with_number(const double *x, const double num, const size_t size);
double transpose_multiplication(const double *vector, const double *transposed, const size_t size);
double vector_norm_2(const double *x, const size_t size);

/**
 * Initializes (allocates memory) for matrix A
 */
double** initialize_A(size_t dimension) {
    double** A = new double*[dimension];
    for(int i=0; i<dimension; ++i) {
        A[i] = new double[dimension];
    }
    return A;
}

/**
 * Generates random numbers for the initialization of vector x
 * @param x: the vector with the random numbers
 */
void randomize_x(double *x, size_t dimension) {
    random_device rd; ///random device used for initializing the distribution
    double random_num; ///random number for seed
    double lower_bound = 1; ///lower bound for seed
    double upper_bound = DBL_MAX; ///upper bound for seed
    uniform_int_distribution<unsigned int> unif(lower_bound, upper_bound); ///declaration of uniform distribution

    x = new double[dimension];
    for(int i=0; i<dimension; ++i){
        random_num = unif(rd); ///get a random number
        x[i] = random_num;
    }
}

/**
 * Reads the data of the input file. Format of the file should be as the given one.
 * @credits: the initial function was found here: http://stackoverflow.com/questions/17392218/how-to-parse-string-file-txt-into-array-with-c
 * @param filename: the input filename to read data from
 * @return -1: if there was any issue with the input file, 0: if everything was done properly
 */
vector<tuple<double**, double*, size_t>> read_file (string filename /*,size_t dimension*/) {
    ifstream in(filename);///read from the file
    vector<tuple<double**, double*, size_t>> examples; ///used to store the instances
    double **A;

    if(!in.good()) {///check if the stream was read properly
        cerr<<"An error occured while trying to parse the file"<<endl;
        return examples;
    }

    vector<string> lines {istream_iterator<string>(in), istream_iterator<string>() }; ///initialize the vector from the values in the file

    int order = 0; ///represents the order of the input problem
    int c_line = 0; ///represents the counter of the input lines
    int number_of_example = 0; ///represents the read number of examples
    double *b; ///represents the input vector

    while(c_line < lines.size()) {

        if(lines[c_line]=="system:") {
            ++number_of_example;///increase number of read examples
            string dimensions = lines[c_line-1];///get the line with the dimensions
            int end = dimensions.find('x');///find the number corresponding to the dimension
            order = stoi(dimensions.substr(0, end));///type cast the data type
        }

        if(lines[c_line]=="A") {
            //dimension = order;
            A = initialize_A(order);///initilize A
            c_line+=3;///skip the next three lines
            ///iteratively fill the matrix
            for(int i=0; i<order; ++i) {
                for(int j=0; j<order; ++j) {
                    double num = stod(lines[c_line]);///type cast the data type
                    A[i][j] = num;///append it to the matrix
                    ++c_line;///move to the next line
                }
            }
        }

        if(lines[c_line]=="b") {
            b = new double[order];
            c_line+=3;///skip the next three lines
            ///iteratively fill the vector
            for(int i=0; i<order; ++i) {
                double num = stod(lines[c_line]);///type cast the data type
                b[i] = num;///append it to the matrix
                ++c_line;///move to the next line
            }
            examples.push_back(make_tuple(A, b, order));///append read example to the vector of examples
        }
        ++c_line;///move to the next line
    }
    in.close();///close stream
    return examples;
}

/**
 * "a function which calculates the matrix/vector product (y = A*x)"
 * @param x: the vector multiplied with the matrix
 * @param y: the result of the multiplication
 */
void my_cool_matrix_vector_product(const double *x, double *y, double **A, size_t dimension) {
    for(int i=0; i<dimension; ++i) {
        double value = 0;
        for(int j=0; j<dimension; ++j) {
            value += A[i][j]*x[j]; ///summing the (row x column) * row element multiplication -inner product
        }
        y[i] = value;
    }
}

/**
 * Calculates the result of the equation Ax = b, for an unknow vector x
 * @param epsilon: the margin of error
 * @param k_max: the maximum number of error
 * @param n: the order of the problem
 * @param b: the vector of the equation
 * @param x: the unknown vector, contains the solution in the end
 * @param matvec: the function for computing the inner product (y = A*x)
 * @return results: a pair showing whether the method converged and the total number of iterations taken place
 */
pair<bool, size_t> conjugate_gradient(double epsilon, size_t k_max, size_t n, const double *b, double *x, size_t dimension, double **A, function<void(const double*, double*,  double**, size_t)> matvec) {

    const size_t order = n; ///constant variable equal to the order of the problem
    size_t k = 0; ///represents the number of iterations taken place
    bool converged = false; ///variable representing whether the method converged or not

    double alpha; ///represents alpha_k

    double *p = new double[order]; ///represents p
    double *r = new double[order]; ///represents r
    double *w = new double[order]; ///represents w

    double *y; ///represents a temporary vector for computations

    matvec(x, r, A, n); ///(assisted computation) r = A*x_0
    r = subtraction(b, r, order); ///(pseudocode) r = b - A*x_0
    copy(r, r+order, p);///(pseudocode) p = r

    double rsold = transpose_multiplication(r, r, order); ///rsold initial value set to ||r||^2

    while(k<k_max) {///the first condition of possible end of method execution

        matvec(p, w, A, n); ///(pseudocode) w = A*p
        alpha = rsold / transpose_multiplication(p, w, order); ///(pseudocode) a_k = rho_{k-1}/p'*w

        y = multiply_with_number(p, alpha, order); ///(assisted computation) y = a_k*p
        x = addition(x, y, order); ///(pseudocode) x = x + a_k*p

        y = multiply_with_number(w, alpha, order); ///(assisted computation) y = a_k*w
        r = subtraction(r, y, order); ///(pseudocode) r = r - a_k*p

        double rsnew = transpose_multiplication(r, r, order); /// (pseudocode) rho_k = ||r||^2

        if(rsnew<epsilon) { ///the other condition of possible end of the method execution
            cerr<<"x = "; print_vector(x, order); ///print the vector of solution
            converged = true;
            return make_pair(converged, k);
        }

        y = multiply_with_number(p, rsnew/rsold, order); ///(assisted computation) y = beta_k*p (beta_k = rsnew/rsold)
        p = addition(r, y, order); ///(pseudocode) p = r + beta_k*p

        rsold = rsnew; ///updating rsold variable
        ++k; ///getting to the next execution
    }

    cerr<<"x = "; print_vector(x, order);///print the vector of solution
    return make_pair(converged, k);
}


int main() {

    vector<tuple<double**, double*, size_t>> read_examples = read_file("test_data.txt"); ///read data from file for testing
    if(read_examples.size() == 0) {
        return -1;
    }

    double epsilon = 0.000001; ///choosing an epsilon of error equal to 1.0e-6
    double *x; ///the pointer for the solution of the method
    size_t k_max = 420; ///choosing a specific number of maximum iterations

    vector<tuple<double**, double*, size_t>>::iterator it = read_examples.begin();///executing the method for each example
    while(it!= read_examples.end()) {

        double **A = get<0>(*it); ///set A to the matrix
        double *b = get<1>(*it);///set b to the vector
        size_t dim  = get<2>(*it);///set dimension to the order(n)
        size_t dimension = dim;
        randomize_x(x, dimension);///random initialization of vector x

        cerr<<"Executing example of dimension equal to: "<<dim<<"x"<<dim<<endl;
        cerr<<"A = "; print_matrix(A, dim); ///print the input matrix
        cerr<<"b = "; print_vector(b, dim); ///print the input vector

        pair<bool, size_t> results = make_pair(false, 0); ///initialization of results -used for storing the results
        results = conjugate_gradient(epsilon, k_max, dim, b, x, dimension, A, my_cool_matrix_vector_product); ///execute the method
        cerr<<"Method converged (1=yes, 0=no): "<<results.first<<" number of total iterations: "<<results.second<<endl;

        ++it;///move to the next example
    }
    return 0;
}
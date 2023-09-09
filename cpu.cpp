#include <iostream>
#include <cmath>

void matmul1(const float *a, const float *b, float *c, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
}

void matmul2(const float *a, const float *_b, float *c, int n) {
    float *b = new float[n * n];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            b[i * n + j] = _b[j * n + i];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[j * n + k]; // <- note the indices
}

// a vector of 256 / 32 = 8 floats
typedef float vec __attribute__ (( vector_size(32) ));

// a helper function that allocates n vectors and initializes them with zeros
vec* alloc(int n) {
    vec* ptr = (vec*) std::aligned_alloc(32, 32 * n);
    memset(ptr, 0, 32 * n);
    return ptr;
}

void matmul3(const float *_a, const float *_b, float *c, int n) {
    int nB = (n + 7) / 8; // number of 8-element vectors in a row (rounded up)

    vec *a = alloc(n * nB);
    vec *b = alloc(n * nB);

    // move both matrices to the aligned region
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i * nB + j / 8][j % 8] = _a[i * n + j];
            b[i * nB + j / 8][j % 8] = _b[j * n + i]; // <- b is still transposed
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vec s{}; // initialize the accumulator with zeros

            // vertical summation
            for (int k = 0; k < nB; k++)
                s += a[i * nB + k] * b[j * nB + k];

            // horizontal summation
            for (int k = 0; k < 8; k++)
                c[i * n + j] += s[k];
        }
    }

    std::free(a);
    std::free(b);
}


int main() {
    //int n = 1920;
    int n = 1024;
    float *a = (float*) malloc(n*n*sizeof(float));
    float *b = (float*) malloc(n*n*sizeof(float));
    float *c = (float*) malloc(n*n*sizeof(float));
    auto t1 = std::chrono::high_resolution_clock::now();
    matmul3(a, b, c, n);
    auto t2 = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0;
    double GHz = 1e9;
    double fma_clock = 0.0625;
    double freq = 3.2*GHz;
    double measured = t * freq / (std::pow((double)n, 3));
    double percent_peak = fma_clock / measured * 100;
    std::cout << "Size (n x n): n = " << n << std::endl;
    std::cout << "Size MB: " << 4.*n*n/(1024*1024) << std::endl;
    std::cout << "Time: " << t << "s" << std::endl;
    std::cout << "Clock cycles per element:" << std::endl;
    std::cout << "Theoretical performance peak: " << fma_clock << " cycles" << std::endl;
    std::cout << "Measured: " << measured << " cycles" << std::endl;
    std::cout << "Percent peak: " << percent_peak << " %" << std::endl;
    return 0;
}

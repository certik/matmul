#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>

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

// update 6x16 submatrix C[x:x+6][y:y+16]
// using A[x:x+6][l:r] and B[l:r][y:y+16]
void kernel(float *a, vec *b, vec *c, int x, int y, int l, int r, int n) {
    vec t[6][2]{}; // will be zero-filled and stored in ymm registers

    for (int k = l; k < r; k++) {
        for (int i = 0; i < 6; i++) {
            // broadcast a[x + i][k] into a register
            vec alpha = vec{} + a[(x + i) * n + k]; // converts to a broadcast
            // multiply b[k][y:y+16] by it and update t[i][0] and t[i][1]
            for (int j = 0; j < 2; j++)
                t[i][j] += alpha * b[(k * n + y) / 8 + j]; // converts to an fma
        }
    }

    // write the results back to C
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 2; j++)
            c[((x + i) * n + y) / 8 + j] += t[i][j];
}

// a helper function that allocates n vectors and initializes them with zeros
float* alloc4(int n) {
    float* ptr = (float*) std::aligned_alloc(32, 32 * n);
    memset(ptr, 0, 32 * n);
    return ptr;
}

void matmul4(const float *_a, const float *_b, float *_c, int n) {
    // to simplify the implementation, we pad the height and width
    // so that they are divisible by 6 and 16 respectively
    int nx = (n + 5) / 6 * 6;
    int ny = (n + 15) / 16 * 16;

    float *a = alloc4(nx * ny);
    float *b = alloc4(nx * ny);
    float *c = alloc4(nx * ny);

    for (int i = 0; i < n; i++) {
        memcpy(&a[i * ny], &_a[i * n], 4 * n);
        memcpy(&b[i * ny], &_b[i * n], 4 * n); // we don't need to transpose b this time
    }

    for (int x = 0; x < nx; x += 6)
        for (int y = 0; y < ny; y += 16)
            kernel(a, (vec*) b, (vec*) c, x, y, 0, n, ny);

    for (int i = 0; i < n; i++)
        memcpy(&_c[i * n], &c[i * ny], 4 * n);

    std::free(a);
    std::free(b);
    std::free(c);
}

void matmul5(const float *_a, const float *_b, float *_c, int n) {
    // to simplify the implementation, we pad the height and width
    // so that they are divisible by 6 and 16 respectively
    int nx = (n + 5) / 6 * 6;
    int ny = (n + 15) / 16 * 16;

    float *a = alloc4(nx * ny);
    float *b = alloc4(nx * ny);
    float *c = alloc4(nx * ny);

    for (int i = 0; i < n; i++) {
        memcpy(&a[i * ny], &_a[i * n], 4 * n);
        memcpy(&b[i * ny], &_b[i * n], 4 * n); // we don't need to transpose b this time
    }

    const int s3 = 64;  // how many columns of B to select
    const int s2 = 120; // how many rows of A to select
    const int s1 = 240; // how many rows of B to select

    for (int i3 = 0; i3 < ny; i3 += s3)
        // now we are working with b[:][i3:i3+s3]
        for (int i2 = 0; i2 < nx; i2 += s2)
            // now we are working with a[i2:i2+s2][:]
            for (int i1 = 0; i1 < ny; i1 += s1)
                // now we are working with b[i1:i1+s1][i3:i3+s3]
                // and we need to update c[i2:i2+s2][i3:i3+s3] with [l:r] = [i1:i1+s1]
                for (int x = i2; x < std::min(i2 + s2, nx); x += 6)
                    for (int y = i3; y < std::min(i3 + s3, ny); y += 16)
                        kernel(a, (vec*) b, (vec*) c, x, y, i1, std::min(i1 + s1, n), ny);

    for (int i = 0; i < n; i++)
        memcpy(&_c[i * n], &c[i * ny], 4 * n);

    std::free(a);
    std::free(b);
    std::free(c);
}



const int B = 8; // number of elements in a vector
typedef float vector __attribute__ (( vector_size(4 * B) ));

float* alloc6(int n) {
    float* ptr = (float*) std::aligned_alloc(64, 4 * n);
    memset(ptr, 0, 4 * n);
    return ptr;
}

// c: 6 x 16
// a: 6 x k
// b: k x 16
// c[x:x+6][y:y+16] += a[x:x+6][l:r] * b[l:r][y:y+16]

void kernel6(const float *a, const vector *b, vector *c, int x, int y, int l, int r, int n) {
    vector t00, t01,
           t10, t11,
           t20, t21,
           t30, t31,
           t40, t41,
           t50, t51;

    t00 = c[((x + 0) * n + y) / 8 + 0];
    t01 = c[((x + 0) * n + y) / 8 + 1];

    t10 = c[((x + 1) * n + y) / 8 + 0];
    t11 = c[((x + 1) * n + y) / 8 + 1];

    t20 = c[((x + 2) * n + y) / 8 + 0];
    t21 = c[((x + 2) * n + y) / 8 + 1];

    t30 = c[((x + 3) * n + y) / 8 + 0];
    t31 = c[((x + 3) * n + y) / 8 + 1];

    t40 = c[((x + 4) * n + y) / 8 + 0];
    t41 = c[((x + 4) * n + y) / 8 + 1];

    t50 = c[((x + 5) * n + y) / 8 + 0];
    t51 = c[((x + 5) * n + y) / 8 + 1];

    for (int k = l; k < r; k++) {
        vector a0 = vector{} + a[(x + 0) * n + k];
        t00 += a0 * b[(k * n + y) / 8];
        t01 += a0 * b[(k * n + y) / 8 + 1];

        vector a1 = vector{} + a[(x + 1) * n + k];
        t10 += a1 * b[(k * n + y) / 8];
        t11 += a1 * b[(k * n + y) / 8 + 1];

        vector a2 = vector{} + a[(x + 2) * n + k];
        t20 += a2 * b[(k * n + y) / 8];
        t21 += a2 * b[(k * n + y) / 8 + 1];

        vector a3 = vector{} + a[(x + 3) * n + k];
        t30 += a3 * b[(k * n + y) / 8];
        t31 += a3 * b[(k * n + y) / 8 + 1];

        vector a4 = vector{} + a[(x + 4) * n + k];
        t40 += a4 * b[(k * n + y) / 8];
        t41 += a4 * b[(k * n + y) / 8 + 1];

        vector a5 = vector{} + a[(x + 5) * n + k];
        t50 += a5 * b[(k * n + y) / 8];
        t51 += a5 * b[(k * n + y) / 8 + 1];
    }

    c[((x + 0) * n + y) / 8 + 0] = t00;
    c[((x + 0) * n + y) / 8 + 1] = t01;

    c[((x + 1) * n + y) / 8 + 0] = t10;
    c[((x + 1) * n + y) / 8 + 1] = t11;

    c[((x + 2) * n + y) / 8 + 0] = t20;
    c[((x + 2) * n + y) / 8 + 1] = t21;

    c[((x + 3) * n + y) / 8 + 0] = t30;
    c[((x + 3) * n + y) / 8 + 1] = t31;

    c[((x + 4) * n + y) / 8 + 0] = t40;
    c[((x + 4) * n + y) / 8 + 1] = t41;

    c[((x + 5) * n + y) / 8 + 0] = t50;
    c[((x + 5) * n + y) / 8 + 1] = t51;
}

extern "C" {
void matmul6(const float *_a, const float *_b, float *_c, int n) {
    int nx = (n + 5) / 6 * 6;
    int ny = (n + 15) / 16 * 16;

    constexpr int MAXN = 1920 * 1920 *4*4; // ~15MB each
    alignas(64) static float a[MAXN], b[MAXN], c[MAXN];
    ptrdiff_t n4 = n * ptrdiff_t(sizeof(float));
    for (ptrdiff_t i = 0; i < n; i++)
        memcpy(&a[i * nx], &_a[i * n], n4);
    for (ptrdiff_t i = 0; i < n; i++)
        memcpy(&b[i * ny], &_b[i * n], n4);
    // for (ptrdiff_t i = 0; i < n; i++) 
        // memcpy(&c[i * ny], &_c[i * n], n4);
    std::memset(c, 0, n4 * ny);
    const int s3 = 64;
    const int s2 = 120;
    const int s1 = 240;

    for (int i3 = 0; i3 < ny; i3 += s3)
        // now we are working with b[:][i3:i3+s3]
        for (int i2 = 0; i2 < nx; i2 += s2)
            // now we are working with a[i2:i2+s2][:]
            for (int i1 = 0; i1 < ny; i1 += s1)
                // now we are working with b[i1:i1+s1][i3:i3+s3]
                // this equates to updating c[i2:i2+s2][i3:i3+s3]
                // with [l:r] = [i1:i1+s1]
                for (int x = i2; x < i2 + s2; x += 6)
                    for (int y = i3; y < i3 + s3; y += 16)
                        kernel6(a, (vector*) b, (vector*) c, x, y, i1, i1 + s1, ny);

    for (ptrdiff_t i = 0; i < n; i++)
       memcpy(&_c[i * n], &c[i * ny], n4);
}
}


int main() {
    //int n = 1920;
    int n = 15*128*2*2;
    int iter = 1;
    float *a = (float*) malloc(n*n*sizeof(float));
    float *b = (float*) malloc(n*n*sizeof(float));
    float *c = (float *)malloc(n * n * sizeof(float));
    for (int i=0; i < n*n; i++) {
        a[i] = 0.125F;
        b[i] = 0.125F;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i=0; i < iter; i++) {
        matmul6(a, b, c, n);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0;
    t = t / iter;
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

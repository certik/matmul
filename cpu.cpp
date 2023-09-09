#include <iostream>

void matmul(const float *a, const float *b, float *c, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
}

int main() {
    int n = 1024;
    float *a = (float*) malloc(n*n*sizeof(float));
    float *b = (float*) malloc(n*n*sizeof(float));
    float *c = (float*) malloc(n*n*sizeof(float));
    matmul(a, b, c, n);
    std::cout << "Done." << std::endl;
    return 0;
}

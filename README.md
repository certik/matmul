# MatMul Benchmark

Run:
```
FC=gfortran cmake -DMATMUL_BLAS=OpenBLAS .
make
OMP_NUM_THREADS=1 ./matmul
```

# Benchmark results on Apple M1

The theoretical performance peak for matmul is just the cost of `fma`, which is
0.125 clock cycles per double precision matrix element (`fmla.2d v0, v0, v0`
takes 0.25 cycles), and 0.0625 per single precision element.

Single precison (f32) matmul

peak = 0.0625 clock cycles

    n    OpenBlas
    512  0.0768
    1024 0.0672
    2048 0.0640
    4096 0.0632
    8192 0.0631

To convert these clock cycles to seconds, multiply by n^3 and divide by 3.2GHz.
For example the n=8192 case gives 10.84s:

    >>> n = 8192; 0.0631*n**3 / 3.2e9
    10.840497455104

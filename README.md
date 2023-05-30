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

    n    OpenBlas  Percent Peak
    512  0.0768       81.4%
    1024 0.0672       93.0%
    2048 0.0640       97.7%
    4096 0.0632       98.9%
    8192 0.0631       99.0%

To convert these clock cycles to seconds, multiply by n^3 and divide by 3.2GHz.
For example the n=8192 case gives 10.84s:

    >>> n = 8192; 0.0631*n**3 / 3.2e9
    10.840497455104

## Rectangle

peak = 0.0625 cycles

 n1    n2   n3  OpenBlas  Percent Peak
32768  240  235  0.0768      81.4%
32768  480  470  0.0684      91.4%
32768  960  940  0.0658      95.0%
32768 2400 2350  0.0647      96.6%
32768 4800 4700  0.0636      98.3%
32768 9600 9400  0.0643      97.2%

To convert these clock cycles to seconds, multiply by `n1*n2*n3` and divide by
3.2GHz. For example the following case gives 3.74s:

    >>> n1 = 32768; n2 = 2400; n3 = 2350; 0.0647 * n1*n2*n3 / 3.2e9
    3.7366579199999994

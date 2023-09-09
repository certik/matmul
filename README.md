# f32 OpenBLAS (via NumPy) Benchmark

This benchmarks (n,n) matrix multiplication using NumPy (which uses OpenBLAS),
on single core:

   n    Time    Cycles    Percent Peak   Size in MB
  512   0.00s   0.07360      84.92%           1
 1024   0.02s   0.06942      90.04%           4
 2048   0.18s   0.06538      95.60%          16
 4096   1.37s   0.06397      97.70%          64
 8192  10.90s   0.06345      98.50%         256
16384  86.88s   0.06321      98.87%        1024

Peak            0.0625      100.00%

# f64 OpenBLAS (via NumPy) Benchmark

This benchmarks (n,n) matrix multiplication using NumPy (which uses OpenBLAS),
on single core:

   n     Time    Cycles    Percent Peak   Size in MB
  512    0.01s   0.15275      80.04%           2
 1024    0.05s   0.13771      90.77%           8
 2048    0.35s   0.13027      95.96%          32
 4096    2.75s   0.12845      97.31%         128
 8192   21.85s   0.12758      97.98%         512
16384  174.82s   0.12720      98.27%        2048

Peak             0.125       100.00%

# f32 C matmul1

   n    Time    Cycles    Percent Peak
  512   0.09s   2.12193       2.94%
 1024   0.88s   2.62856       2.38%
 2048  23.34s   8.69408       0.72%

Peak            0.0625      100.00%

# f32 C matmul2

   n    Time    Cycles    Percent Peak
  512   0.01s   0.238419     26.21%
 1024   0.08s   0.244379     25.58%
 2048   0.69s   0.258535     24.17%
 4096   5.35s   0.249129     25.09%
 8192  42.85s   0.249431     25.06%

Peak            0.0625      100.00%

# f32 C matmul3

   n    Time    Cycles    Percent Peak
  512   0.01s   0.309944     20.17%
 1024   0.12s   0.342727     18.24%
 2048   1.14s   0.426173     14.67%
 4096   9.98s   0.46487      13.44%

Peak            0.0625      100.00%

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

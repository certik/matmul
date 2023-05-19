# MatMul Benchmark

Run:
```
FC=gfortran cmake -DMATMUL_BLAS=OpenBLAS .
make
OMP_NUM_THREADS=1 ./matmul
```


#!/bin/bash

set -ex

FC=gfortran cmake -Bbuild -DMATMUL_BLAS=OpenBLAS
cmake --build build --parallel
OMP_NUM_THREADS=1 build/matmul

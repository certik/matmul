program matmul_benchmark

! Apple M1
! https://www.techpowerup.com/cpu-specs/core-i9-10980hk.c2276
! https://dougallj.github.io/applecpu/firestorm-simd.html
!
! Operation speeds Apple M1 (ARM64) per double
!
! R: 0.1665 - 0.2 - 0.25 (`ldr q0, [x1]` takes 0.333)
! W: 0.25 - 0.333 - 0.5   (Both `stp d0, d1, [x1]` and `str q1, [x1]` take 0.5 cycles; `stp q1, q2, [x1]` takes 1 cycle)
! *: 0.125  (`fmul.2d v0, v0, v0` takes 0.25)
! +: 0.125  (`fadd.2d v0, v0, v0` takes 0.25)
! fma: 0.125  (`fmla.2d v0, v0, v0` takes 0.25)
! min/max: 0.125 (`fmaxnm.2d` takes 0.25)
! abs: 0.125 (`fabs.2d` takes 0.25)
! int->double, double->int: 0.125 (`fcvtzs` and `scvtf` take 0.25 each)
!
! Note: There are two units for R, one unit for W and one unit can do both. The
! first number in 0.25 - 0.333 - 0.5 is only writing, then read/write sharing
! 50% of the common unit, and the last number is only the 1 unit for  W.
! Example: For array copy, we assume the unit gets used 50%, use the middle
! number and expect 0.333 for the peak.

use linalg, only: matmul_2d
use m, only: matmul1, matmul2, matmul3, matmul4
implicit none
integer, parameter :: dp = kind(0.d0)
integer, parameter :: sp = kind(0.0)
real(sp), allocatable :: A(:,:), B(:,:), C(:,:), C_ref(:,:)
real(dp) :: t1, t2, cycles, peak, GHz, freq
integer :: n
GHz = 1e9_dp
freq = 3.2_dp * GHz

n = 2048
allocate(A(n,n), B(n,n), C(n,n), C_ref(n,n))
call random_number(A)
call random_number(B)

peak = 0.125/2 ! cycles per array element, dominated by fma, single precision
print '("Benchmarking matmul; n = ",i0," peak =", f7.4, " cycles")', n, peak

call cpu_time(t1)
C = matmul(A, B)
call cpu_time(t2)
cycles = (t2-t1)*freq/real(n,dp)**3
print "('Fortran matmul: ',f7.3,' s =',f7.4,' cycles (',f4.1,'%)')", t2-t1, &
    cycles, 100*peak/cycles
C_ref = C
print *, "Error =", maxval(abs(C-C_ref))

call cpu_time(t1)
call matmul_2d(A, B, C)
call cpu_time(t2)
cycles = (t2-t1)*freq/real(n,dp)**3
print "('OpenBLAS matmul: ',f7.3,' s =',f7.4,' cycles (',f4.1,'%)')", t2-t1, &
    cycles, 100*peak/cycles
print *, "Error =", maxval(abs(C-C_ref))

call cpu_time(t1)
call matmul1(n, A, B, C)
call cpu_time(t2)
cycles = (t2-t1)*freq/real(n,dp)**3
print "('matmul 3 loops: ',f7.3,' s =',f7.4,' cycles (',f4.1,'%)')", t2-t1, &
    cycles, 100*peak/cycles
print *, "Error =", maxval(abs(C-C_ref))

call cpu_time(t1)
call matmul2(n, A, B, C)
call cpu_time(t2)
cycles = (t2-t1)*freq/real(n,dp)**3
print "('matmul blocks 1:',f7.3,' s =',f7.4,' cycles (',f4.1,'%)')", t2-t1, &
    cycles, 100*peak/cycles
print *, "Error =", maxval(abs(C-C_ref))

call cpu_time(t1)
call matmul3(n, A, B, C)
call cpu_time(t2)
cycles = (t2-t1)*freq/real(n,dp)**3
print "('matmul blocks 2:',f7.3,' s =',f7.4,' cycles (',f4.1,'%)')", t2-t1, &
    cycles, 100*peak/cycles
print *, "Error =", maxval(abs(C-C_ref))

call cpu_time(t1)
call matmul4(n, A, B, C)
call cpu_time(t2)
cycles = (t2-t1)*freq/real(n,dp)**3
print "('matmul blocks 3:',f7.3,' s =',f7.4,' cycles (',f4.1,'%)')", t2-t1, &
    cycles, 100*peak/cycles
print *, "Error =", maxval(abs(C-C_ref))

print *, "C_ref(:2,:2) =", C_ref(:2,:2)

end program

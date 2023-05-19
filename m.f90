module m

integer, parameter :: sp = kind(0.0)

contains

    subroutine matmul1(n, A, B, C)
    ! matmul using 3 loops
    real(sp), intent(in) :: A(n,n), B(n,n)
    real(sp), intent(out) :: C(n,n)
    integer :: i, j, k
    ! C_ij = A_i^k B_kj
    C = 0
    do j = 1, n
    do k = 1, n
    do i = 1, n
        C(i,j) = C(i,j) + A(i,k)*B(k,j)
    end do
    end do
    end do
    end subroutine

    subroutine matmul2(n, A, B, C)
    ! matmul using blocks
    real(sp), intent(in) :: A(n,n), B(n,n)
    real(sp), intent(out) :: C(n,n)
    integer, parameter :: M = 8 ! Block size
    integer :: i, j, k, ii, jj, kk
    ! C_ij = A_i^k B_kj
    C = 0
    do jj = 0, n-1, M
    do kk = 0, n-1, M
    do ii = 0, n-1, M
        do j = 1+jj, jj+M
        do k = 1+kk, kk+M
        do i = 1+ii, ii+M
            C(i,j) = C(i,j) + A(i,k)*B(k,j)
        end do
        end do
        end do
    end do
    end do
    end do
    end subroutine

    subroutine matmul3(n, A, B, C)
    ! matmul using blocks
    real(sp), intent(in) :: A(n,n), B(n,n)
    real(sp), intent(out) :: C(n,n)
    integer, parameter :: M = 64 ! Block size
    real(sp) :: BA(M,M), BB(M,M)
    integer :: i, j, k, ii, jj, kk
    ! C_ij = A_i^k B_kj
    C = 0
    do jj = 0, n-1, M
    do kk = 0, n-1, M
    do ii = 0, n-1, M
        BA(:,:) = A(ii+1:ii+M, kk+1:kk+M)
        BB(:,:) = B(kk+1:kk+M, jj+1:jj+M)
        do j = 1, M
        do k = 1, M
        do i = 1, M
            C(ii+i, jj+j) = C(ii+i, jj+j) + BA(i,k)*BB(k,j)
        end do
        end do
        end do
    end do
    end do
    end do
    end subroutine

end module

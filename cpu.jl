using LoopVectorization, LinearAlgebra

function kernel_turbo!(C,A,B,β)
  @turbo for n = indices((C,B),2), m = indices((C,A),1)
    Cmn = zero(eltype(C))
    for k = indices((A,B),(2,1))
      Cmn += A[m,k]*B[k,n]
    end
    C[m,n] = β===static(false) ? Cmn : Cmn + C[m,n]
  end
end
function matmul_turbo!(C,A,B)
  N = LinearAlgebra.checksquare(C)
  NA = LinearAlgebra.checksquare(A)
  NB = LinearAlgebra.checksquare(B)
  @assert N==NA==NB
  S3 = 64;
  S2 = 120;
  S1 = 240;
  @assert N%S3==0 && N % S1 == 0
  m = 0
  while m < N
    n = 0
    while n < N
      @views kernel_turbo!(
        C[1+m:m+S3,1+n:n+S2],
        A[1+m:m+S3,1:S1],
        B[1:S1,1+n:n+S2], static(false))
      k = S1
      while k < N
        @views kernel_turbo!(
          C[1+m:m+S3,1+n:n+S2],
          A[1+m:m+S3,1+k:k+S1],
          B[1+k:k+S1,1+n:n+S2], static(true))
        k += S1
      end      
      n += S2
    end    
    m += S3  
  end
end


run(`clang++ -std=c++17 -O3 -ffast-math -funroll-loops cpu.cpp -shared -fPIC -o libmatmul.so`)
const libcppmatmul = joinpath(pwd(), "libmatmul.so")
function matmul_cpp!(C,A,B)
  N = LinearAlgebra.checksquare(C)
  NA = LinearAlgebra.checksquare(A)
  NB = LinearAlgebra.checksquare(B)
  @assert N==NA==NB
  @ccall libcppmatmul.matmul6(B::Ptr{Float32}, A::Ptr{Float32}, C::Ptr{Float32}, (N%Int32)::Int32)::Cvoid
end

N=15*128*2*2;
T=Float32;M=K=N; A=rand(T,M,K);B=rand(T,K,N);C=Matrix{T}(undef,M,N);
@time matmul_turbo!(C,A,B);
C2 = similar(C);
@time matmul_cpp!(C2,A,B);
@assert C ≈ C2


using Kronecker
using CUDA


A = CUDA.randn(100, 100);
B = CUDA.rand(50, 50);

v = CUDA.rand(5000);

K = A ⊗ B
r = K * v


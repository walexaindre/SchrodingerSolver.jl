using SchrodingerSolver
using BenchmarkTools
using SparseArrays
using GLMakie
using Krylov
using LinearAlgebra
using CUDA
CUDA.allowscalar(false)

sz=8

PA = PeriodicAbstractMesh(Int, (sz,sz))
AI,AJ,AV = get_A_format_IJV(Float64, PA, (:ord4,:ord4))
opA = sparse(AI, AJ, AV|>Array{ComplexF64})

sparsity_pattern = sparsitypattern(AI, AJ)

b = rand(ComplexF64, sz^2)

@btime gmres($opA, $b)
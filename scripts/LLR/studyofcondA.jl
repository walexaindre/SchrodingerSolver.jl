using SchrodingerSolver
using BenchmarkTools
using SparseArrays
using GLMakie
using Krylov
using LinearAlgebra
using CUDA
CUDA.allowscalar(false)

sz=100

PA = PeriodicAbstractMesh(Int, (sz,sz,sz))
Grid = PeriodicGrid(Int,Float64,1.0,(sz,sz,sz),((0.0,1.0),(0.0,1.0),(0.0,1.0)))

AI,AJ,AV = get_A_format_IJV(Float64, PA, (:ord4,:ord4,:ord4))
opA = sparse(AI, AJ, AV|>Array{ComplexF64})

#sparsity_pattern = sparsitypattern(AI, AJ,AV)
DI,DJ,DV = get_D_format_IJV(Float64, Grid, (:ord4,:ord4,:ord4))
b = rand(ComplexF64, sz^3)
opD = sparse(DI, DJ, DV|>Array{ComplexF64})
I,J,V = drop(opB,PA,0.001,1:3)
rf = sparse(I,J,V)

rs,_ = cg(opD, b,M=rf)
rs,_ = gmres(opD, b,M=rf)
rs,_ = gmres(opB, b)
rs,_ = cg(opB, b)

opB = 4im*opA+opD

@btime gmres($opA, $b)
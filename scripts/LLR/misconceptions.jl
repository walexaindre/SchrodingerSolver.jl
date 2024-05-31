using SchrodingerSolver
using BenchmarkTools
using SparseArrays
using GLMakie
using LinearAlgebra
using CUDA
CUDA.allowscalar(false)

Grid = PeriodicGrid(Int,Float64,1.0,(1:0.5:4,1:0.5:4))
Grid1 = PeriodicGrid(Int,Float64,1.0,(1:0.5:4,))
PA = PeriodicAbstractMesh(Int, size(Grid))
PA1 = PeriodicAbstractMesh(Int, size(Grid1))

AI,AJ,AV = get_A_format_COO(Float64, PA, (:ord4,:ord4))

DI,DJ,DV = get_D_format_COO(Float64, Grid, (:ord4,:ord4))

A1I,A1J,A1V = get_A_format_IJV(Float64, PA1, (:ord4))
D1I,D1J,D1V = get_D_format_IJV(Float64, Grid1, (:ord4,))

opA= sparse(AI, AJ, AV)
opD = sparse(DI, DJ, DV)
op1A = sparse(A1I, A1J, A1V)
op1D = sparse(D1I, D1J, D1V)
opProd1A1D = op1A * op1D

B= 4im*opA+opD

BI,BJ,BV = findnz(B)

Ix = I(size(op1A,1))

kronsum = kron(op1A, op1D) + kron(op1D, op1A)

anotherkronsum = opA * opD

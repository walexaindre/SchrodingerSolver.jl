module SchrodingerSolver
    using CUDA
    using CUDA.CUSPARSE

    using Distributed
    using SparseMatricesCSR
    using StaticArrays

    using Printf
    using MKLSparse
    using SparseArrays

    using LinearAlgebra
    using GLMakie
    using Krylov
    using IncompleteLU


    using Base.Threads 
    using Base.Iterators
    using ProgressMeter
    using SparseMatricesCSR
    using Dictionaries
    using PrettyTables

    include("./core/index.jl")

end
using SchrodingerSolver
using Base.Iterators: product
using Test
using Aqua

@testset "SchrodingerSolver.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        #Aqua.test_all(SchrodingerSolver)
    end
    # Write your tests here.

    include("index.jl")
end

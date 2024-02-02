using BenchmarkTools

using SchrodingerSolver

Mesh2 = AbstractMesh(Int64,10,4,2)
function llvms()
    Mesh = AbstractMesh(Int64,10,4,2)

    @inbounds Mesh[:,99,:]
end

llvms()
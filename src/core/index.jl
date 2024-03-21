#Threaded CSR
include("ThreadedSparseCSR/ThreadedSparseCSR.jl")
using .ThreadedSparseCSR

#Types
include("AbstractMesh/AbstractMeshTypes.jl")
include("Grid/GridTypes.jl")
#include("Stencil/StencilTypes.jl")
#include("ComputeBackend/ComputeBackend.jl")
#include("FiniteDifferences/FiniteDifferenceTypes.jl")

#Methods and declarations
include("AbstractMesh/AbstractMesh.jl")
include("Grid/Grid.jl")
#include("FiniteDifferences/FiniteDifferences.jl")
#include("Stencil/Stencil.jl")
#include("SpaceTimeGrid/SpaceTimeGrid.jl")
#include("SchrodingerPDE/SchrodingerPDE.jl")
#include("SpaceDiscretization/SpaceDiscretization.jl")
#include("TimeDiscretization/TimeDiscretization.jl")

#Solver end
#include("Solver/Solver.jl")
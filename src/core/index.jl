#Threaded CSR
include("ThreadedSparseCSR/ThreadedSparseCSR.jl")
using .ThreadedSparseCSR

#Types
include("AbstractMesh/AbstractMeshTypes.jl")
include("Grid/GridTypes.jl")
include("SpaceDiscretization/SpaceDiscretizationTypes.jl")
include("TimeComposition/TimeCompositionTypes.jl")
include("SchrodingerPDE/SchrodingerPDETypes.jl")
include("PDEStats/PDEStatsTypes.jl")
include("Offset/OffsetTypes.jl")
include("Image/ImageTypes.jl")
include("Solver/SolverTypes.jl")
include("Solver/MethodTypes.jl")
include("Solver/Backend/BackendTypes.jl")
include("Solver/Procedure/indextypes.jl")
#include("Stencil/StencilTypes.jl")
#include("ComputeBackend/ComputeBackend.jl")
#include("FiniteDifferences/FiniteDifferenceTypes.jl")

#Methods and declarations
include("TupleSort/TupleSort.jl")
include("AbstractMesh/AbstractMesh.jl")
include("Grid/Grid.jl")
include("SpaceDiscretization/SpaceDiscretization.jl")
include("TimeComposition/TimeComposition.jl")
include("SchrodingerPDE/SchrodingerPDE.jl")
include("PDEStats/PDEStats.jl")
include("Preconditioner/Preconditioner.jl")
include("Offset/Offset.jl")
include("Image/Image.jl")
include("Solver/Solver.jl")
include("Solver/Method.jl")
include("Solver/Backend/Backend.jl")
include("Solver/Procedure/index.jl")




#Solver end
#include("Solver/Solver.jl")
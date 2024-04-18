using SchrodingerSolver
using BenchmarkTools

A = AssemblySymmetricOffset(AllZeroOffset,(1:2,1:3,1:2))

PA = PeriodicAbstractMesh(Int,(5,10,5))

I = (1,1,1)

offsets = apply_offsets(PA, I, A)

SchrodingerSolver.core_circulant_matrix_format_IJV(1:length(A), offsets, PA)
using BenchmarkTools

using SchrodingerSolver

import SchrodingerSolver:Stencil,StencilPattern, AbstractMesh,to_matrix,LinearOffsetStencil

mesh2 = AbstractMesh(Int,20,15)

pat = StencilPattern(2,1,0,1)

stencil = Stencil(pat,mesh2)

mat = to_matrix(Matrix{Int},stencil)


rr =  LinearOffsetStencil(Int,(3,4),(0,1))

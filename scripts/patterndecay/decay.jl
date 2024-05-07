using SchrodingerSolver
using GLMakie
using SparseArrays
using LinearAlgebra

savepath = "./scripts/patterndecay/plots/"
basefilename = "patterndecay2d_"

Mesh = PeriodicAbstractMesh(Int64,(6,6))

AI,AJ,AV = get_A_format_COO(Float64,Mesh,(:ord4,:ord4))

opAdense = sparse(AI,AJ,AV)|>Array
opAinv=inv(opAdense)

opAinvI,opAinvJ,opAinvV = findnz(opAinv|>sparse)
f=Figure(size=(800,800))

ax = Axis(f[1,1])
hidedecorations!(ax)
hidespines!(ax)

St  = GenerateOffset(OffsetUniqueZero,1,(1:1,1:1))

CI,CJ,CV =core_circulant_matrix_format_COO([1,2,2,3,3],St,Mesh)

sparsitypattern!(ax,CI,CJ,CV)

sparsitypattern!(ax,opAinvI,opAinvJ,opAinvV)
sparsitypattern!(ax,AI,AJ,AV,colormap=:grayC100,colorscale=identity)

approximatesparsitypattern!(ax,opAinvI,opAinvJ,opAinvV,colormap=:binary,colorscale=identity)
approximatesparsitypattern!(ax,AI,AJ,AV,colormap=:grayC100,colorscale=sqrt)

save(joinpath(savepath,basefilename*"circulantsparsitypattern.png"),f)

f
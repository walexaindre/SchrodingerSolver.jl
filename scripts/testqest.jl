using SchrodingerSolver
using BenchmarkTools
using SparseArrays
using GLMakie
using CUDA
CUDA.allowscalar(false)

A = GenerateOffset(OffsetUniqueZero,1, (1:2,1:3))

PA = PeriodicAbstractMesh(Int, (8,8))
Grid = PeriodicGrid(Int,Float64,1.0,(1:0.5:4,1:0.5:4))

SI,SJ,SV = core_circulant_matrix_format_COO(collect(1:length(A)),A,PA)
sparse(SI,SJ,SV)
sparsitypattern(SI,SJ,SV)


I = (1, 1)

offsets = apply_offsets(PA, I, A)

function gensym(sz)
    r = div(sz - 1, 2)
    vl = repeat(collect(1:r); inner = 2)
    return vcat(-1, vl)
end
M = SchrodingerSolver.core_circulant_matrix_format_IJV(gensym(length(A)), A, PA)
#@code_warntype SchrodingerSolver.core_circulant_matrix_format_IJV(gensym(length(A)), A, PA)
sparsitypattern(M[1], M[2], M[3])

for typev in Iterators.drop(subtypes(CUDA.CUSPARSE.AbstractCuSparseMatrix), 2)
    I = [1, 1]
    J = [1, 2]
    V = [1.0, 2.0]
    sarr = sparse(I, J, V)
    println(SchrodingerSolver.matrix_to_vector(typeof(typev(sarr))))
end

function rsz()
    n = 500
    c = 55444
    vecs = rand(1:(c - 2), n, 2)
    println("$n $c")
    for row in eachrow(vecs)
        println(join(row, " "))
    end
end
open("myfile.txt", "w") do io
    redirect_stdout(io) do
        rsz()
    end

end

sortperm
for i in 0:10
    c=11
    r = abs(i)
   println( min(r,c-r))
end



AI,AJ,AV = get_A_format_COO(Float64, PA, (:ord4,:ord4))

DI,DJ,DV = get_D_format_IJV(Float64, Grid, (:ord4,:ord4))

sparsitypattern(DI,DJ,DV)


Z=sparse(AI,AJ,AV)
BI,BJ,BV = drop(Z,PA,0.000001)
sparsitypattern(AI,AJ,AV)
f=Figure(size=(800,800))
MA = sparse(AI,AJ,AV)
kron(MA,MA)
ax = Axis(f[1,1])
hidedecorations!(ax)
hidespines!(ax)
sparsitypattern(BI,BJ,BV)
sparsitypattern!(ax,AI,AJ,AV,colormap=:magma)

save("sparsitypattern.png",f)

rows = Vector{Int}(undef, 0)
cols = Vector{Int}(undef, 0)
vals = Vector{Float64}(undef, 0)
idx = 1
sizehint!(rows, length(PA)*length(A))
sizehint!(cols, length(PA)*length(A))
sizehint!(vals, length(PA)*length(A))
for i in  CartesianIndices(PA)
    res = apply_offsets(PA,i,A)
    append!(rows,res)
    append!(cols,fill(idx,length(res)))
    append!(vals,fill(1.0,length(res)))
    idx+=1
end
f
sparsitypattern(rows,cols,vals)
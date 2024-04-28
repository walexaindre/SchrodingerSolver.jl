using SchrodingerSolver
using BenchmarkTools
using SparseArrays
using CUDA
CUDA.allowscalar(false)

A = AssemblySymmetricOffset(UniqueZeroOffset, (1:3,))

PA = PeriodicAbstractMesh(Int, (36,))

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



AI,AJ,AV = get_A_format_IJV(Float64, PA, (:ord8))
Z=sparse(AI,AJ,AV)
BI,BJ,BV = drop(Z,PA,0.000001)

sparsitypattern(BI,BJ,BV)
sparsitypattern(AI,AJ,AV)
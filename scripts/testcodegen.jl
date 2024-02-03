using BenchmarkTools

using SchrodingerSolver

Mesh2 = AbstractMesh(Int64,10,4,2)
function llvms(idx)
    Mesh = AbstractMesh(Int64,10,4,2)

    @inbounds Mesh[idx]
end

a=1
#@btime llvms($a)


function rev(i,n)
    mod(i - 1, n) + 1
end

function rews(i,n)
    if 1<=i<=n
        return i
    else
        return mod(i - 1, n) + 1
    end
end

n = rand(50000:3000000);
a = rand(-n:2*n,10000000);

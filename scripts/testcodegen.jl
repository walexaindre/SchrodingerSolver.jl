using BenchmarkTools

using SchrodingerSolver

import SchrodingerSolver:  PeriodicAbstractMesh, getindex2, getindex3


#stencil = Stencil(pat,mesh2)

#it = StencilIterator(Int, 3, (1,2), (1, 2, 3))
#it = StencilIterator(Int, 3, (1,2), (1, 3, 2))
#it = StencilIterator(Int, 3, (1,2), (2, 1, 3))
#it = StencilIterator(Int, 3, (1,2), (2, 3, 1)) 
#it = StencilIterator(Int, 3, (1,2), (3, 1, 2))
#it = StencilIterator(Int, 3, (1,2), (3, 2, 1))


#rr = LinearOffsetStencil(Int, (3, 4), (1, 0))
#=



av = rand(1:180,500)

function sample(a)
    mesh3 = AbstractMesh(Int, 3, 3, 3)
    stencil = Stencil(Int, mesh3, (1,))
    return stencil[a]
end


# Iteration no. 1, manual iteration, assuming that filter is 3x3
function conv_1(input, filter)
    input_r, input_c = size(input)
    filter_r, filter_c = size(filter)

    if filter_r != filter_c
        throw(DomainError(filter, "Filter row and column must be equals"))
    end

    # Because the filter is 3x3, then the conv matrix will be reduced by 2.
    result = zeros(input_r-2, input_c-2)
    result_r, result_c = size(result)

    for i in 1:result_r
        for j in 1:result_c
            result[i,j] += input[i,j]*filter[1,1] + input[i,j+1]*filter[1,2] + input[i,j+2]*filter[1,3]
            result[i,j] += input[i+1,j]*filter[2,1] + input[i+1,j+1]*filter[2,2] + input[i+1,j+2]*filter[2,3]
            result[i,j] += input[i+2,j]*filter[3,1] + input[i+2,j+1]*filter[3,2] + input[i+2,j+2]*filter[3,3]
        end
    end

    return result
end

input = [[3,2,4,2,7,6] [8,0,2,1,7,8] [2,2,10,4,1,9] [1,5,4,6,5,0] [5,4,1,7,5,6] [5,0,2,7,6,8]]
filterv = [[0,1,0] [1,2,1] [0,1,0]]
conv_1(input, filterv)

=#

#bmark

mesh3d = PeriodicAbstractMesh(Int, 30, 30, 30)

randvals = rand(-20*20*20:30*30*34,500,3)

function numleak(randvals)
    rx = 0;
    for i in 1:500
        rx += getindex(mesh3d,randvals[i,1],randvals[i,2],randvals[i,3])
    end
    rx
end

function rv(A,col)
    if col == 1
        return repeat(1:A.dims[1]; outer=size(A, 2) * size(A, 3))
    elseif col == 2
        return repeat(1:A.dims[2]; inner=size(A, 1), outer=size(A, 3))
    else
        return repeat(1:A.dims[3]; inner=size(A, 1) * size(A, 2))
    end
end

randvals = [rv(mesh3d,1) rv(mesh3d,2) rv(mesh3d,3)]

@btime getindex(mesh3d,$randvals[1,1],$randvals[1,2],$randvals[1,3])
@btime numleak($randvals)
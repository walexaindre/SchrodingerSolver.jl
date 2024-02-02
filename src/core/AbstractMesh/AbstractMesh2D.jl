#2D

@inline function linear_indexing(x::Int, y::Int, metadata::MetaMesh2D)
    (y - 1) * metadata.M + x
end

@inline function linear_indexing(index::Int, metadata::MetaMesh2D)
    col, row = divrem(index - 1, metadata.M)
    row + 1, col + 1
end




@inline function Base.getindex(Mesh::MetaMesh2D, index::Int)
    @boundscheck begin
        if !(1 <= index <= Mesh.MN)
            throw(BoundsError(Mesh, index))
        end
    end
    linear_indexing(index, Mesh)
end


@inline function Base.getindex(Mesh::MetaMesh2D, i::Int, j::Int)
    @boundscheck begin
        if !((1 <= i <= Mesh.M) && (1 <= j <= Mesh.N))
            throw(BoundsError(Mesh, (i, j)))
        end
    end
    linear_indexing(i, j, Mesh)
end






@inline function Base.lastindex(Mesh::MetaMesh2D, d::Int)
    if d == 1
        return Mesh.M
    elseif d == 2
        return Mesh.N
    else
        throw(BoundsError(Mesh, d))
    end
end

    
function Base.iterate(AbstractMesh::MetaMesh2D, state::Int = 1)
    state > length(AbstractMesh) ? nothing : (linear_indexing(state,AbstractMesh),state+1)
end

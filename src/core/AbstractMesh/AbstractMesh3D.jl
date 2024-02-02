#3D

@inline function linear_indexing(index::Int, metadata::MetaMesh3D)
    z, remaining_position = divrem(index - 1, metadata.MN)
    y, x = divrem(remaining_position, metadata.M)
    x + 1, y + 1, z + 1
end

@inline function linear_indexing(x::Int, y::Int, z::Int, metadata::MetaMesh3D)
    (z - 1) * metadata.MN + (y - 1) * metadata.M + x
end






@inline function Base.getindex(Mesh::MetaMesh3D, index::Int)
    @boundscheck begin
        if !(1 <= index <= linear_size(Mesh))
            throw(BoundsError(Mesh, index))
        end
    end
    linear_indexing(index, Mesh)
end




@inline function Base.getindex(Mesh::MetaMesh3D, i::Int, j::Int, k::Int)
    @boundscheck begin
        if !((1 <= i <= Mesh.M) && (1 <= j <= Mesh.N) && (1 <= k <= Mesh.L))
            throw(BoundsError(Mesh, (i, j, k)))
        end
    end
    linear_indexing(i, j, k, Mesh)
end




@inline function Base.lastindex(Mesh::MetaMesh3D, d::Int)
    if d == 1
        return Mesh.M
    elseif d == 2
        return Mesh.N
    elseif d == 3
        return Mesh.L
    else
        throw(BoundsError(Mesh, d))
    end
end

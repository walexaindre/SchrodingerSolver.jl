#Collection of generic functions for ND meshes
####################################################################


@inline function Base.firstindex(Mesh::T) where {T<:MetaMesh}
    1
end

@inline function Base.lastindex(Mesh::T) where {T<:MetaMesh}
    length(Mesh)
end


# N dims

@inline function Base.ndims(Mesh::MetaMesh1D)
    1
end

@inline function Base.ndims(Mesh::MetaMesh2D)
    2
end

@inline function Base.ndims(Mesh::MetaMesh3D)
    3
end


# Size

@inline function Base.size(Mesh::MetaMesh1D)
    Mesh.M
end

@inline function Base.size(Mesh::MetaMesh2D)
    (Mesh.M, Mesh.N)
end

@inline function Base.size(Mesh::MetaMesh2D, dim::Int)
    @boundscheck begin
        if !(1 <= dim <= 2)
            throw(BoundsError(Mesh, dim))
        end
    end

    if dim == 1
        return Mesh.M
    end
    Mesh.N
end

@inline function Base.size(Mesh::MetaMesh3D)
    (Mesh.M, Mesh.N, Mesh.L)
end

@inline function Base.size(Mesh::MetaMesh3D, dim::Int)
    @boundscheck begin
        if !(1 <= dim <= 3)
            throw(BoundsError(Mesh, dim))
        end
    end

    if dim == 1
        return Mesh.M
    elseif dim == 2
        return Mesh.N
    end
    Mesh.L
end

# length

@inline Base.length(AbstractMesh::MetaMesh1D) = AbstractMesh.M
@inline Base.length(AbstractMesh::MetaMesh2D) = AbstractMesh.MN
@inline Base.length(AbstractMesh::MetaMesh3D) = AbstractMesh.MNL


# eltype

@inline Base.eltype(::Type{MetaMesh1D}) = Int64
@inline Base.eltype(::Type{MetaMesh2D}) = (Int64,Int64)
@inline Base.eltype(::Type{MetaMesh3D}) = (Int64,Int64,Int64)


# checkbounds_indices

@inline function checkbounds_indices(Mesh::MetaMesh1D,i::Int)
    Base.checkindex(Bool,1:Mesh.M,i)
end


@inline function checkbounds_indices(Mesh::MetaMesh2D,i::Int, j::Int)
    Base.checkbounds_indices(Bool,(1:Mesh.M,1:Mesh.N),(i,j))
end

@inline function checkbounds_indices(Mesh::MetaMesh3D,i::Int, j::Int, k::Int)
    Base.checkbounds_indices(Bool,(1:Mesh.M,1:Mesh.N,1:Mesh.L),(i,j,k))
end

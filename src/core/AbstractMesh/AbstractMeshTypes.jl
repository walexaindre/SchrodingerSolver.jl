# An implementation of an abstract N dimensional mesh.
# Here we can map from the N dimensional index to the 1D index and vice versa.

struct AbstractMesh{T<:Int,N} <: AbstractArray{T,N}
    ghost_cells_depth::Int
    dims::NTuple{N,Int}
    multiplied_dims::NTuple{N,Int}
end

@inline AbstractMesh(::Type{T}, dims::NTuple{N,Int}) where {T<:Int,N} = AbstractMesh{T,N}(dims, cumprod(dims))
@inline AbstractMesh(::Type{T}, dims::Int...) where {T<:Int} = AbstractMesh(T, dims)

@inline Base.similar(A::AbstractMesh, ::Type{T}, dims::Dims) where {T<:Int} = AbstractMesh(T, dims)
@inline Base.copy(A::AbstractMesh{T}) where {T<:Int} = AbstractMesh(T,A.dims)
@inline Base.size(A::AbstractMesh) = A.dims
@inline Base.length(A::AbstractMesh) = A.multiplied_dims[end]


@inline function Base.getindex(A::AbstractMesh{T,N}, I::Vararg{Int,N}) where {T<:Int,N}
    I = @. mod1(I-1,A.dims)
    index = I[1]
    @simd for i in 2:N
        @inbounds index += (I[i]-1)*A.multiplied_dims[i-1]
    end
    index
end

@inline function getindex2(A::AbstractMesh{T,N}, I::Vararg{Int,N}) where {T<:Int,N}
    index = mod1(I[1]-1,A.dims[1])
    @simd for i in 2:N
        @inbounds index += (mod(I[i],A.dims[i]))*A.multiplied_dims[i-1]
    end
    index
end
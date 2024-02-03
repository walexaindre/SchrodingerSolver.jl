struct AbstractMesh{T,N} <: AbstractArray{Int,N}
    dims::NTuple{N,Int}
    multiplied_dims::NTuple{N,Int}
end

@inline AbstractMesh(::Type{T}, dims::Int...) where {T<:Int} = AbstractMesh(T, dims)
@inline AbstractMesh(::Type{T}, dims::NTuple{N,Int}) where {T<:Int,N} = AbstractMesh{T,N}(dims, cumprod(dims))


@inline Base.similar(A::AbstractMesh, ::Type{T}, dims::Dims) where {T<:Int} = AbstractMesh(T, dims)
@inline Base.copy(A::AbstractMesh{T}) where {T<:Int} = AbstractMesh(T,A.dims)
@inline Base.size(A::AbstractMesh) = A.dims
@inline Base.length(A::AbstractMesh) = A.multiplied_dims[end]

@inline function Base.getindex(A::AbstractMesh{T,N}, I::Vararg{Int,N}) where {T<:Int,N}
    I = @. mod(I-1,A.dims) + 1
    index = I[1]
    for i in 2:N
        index += (I[i]-1)*A.multiplied_dims[i-1]
    end
    index
end

@inline function Base.getindex(A::AbstractMesh{T, N}, index::V) where {T<:Int, V<:Int, N}
    @boundscheck begin
        if !(1 <= index <= length(A))
            throw(BoundsError(A, index))
        end
    end

    indices = MVector{N,Int}(undef)
    remainder = index - 1

    for idx in N:-1:2
        indices[idx], remainder = divrem(remainder, A.multiplied_dims[idx-1])
    end

    indices[1] = remainder

    Tuple(indices.+1)
end

@inline function Base.getindex(A::AbstractMesh{T,1}, ::Colon) where {T<:Int}
    collect(1:A.dims[1])
end

@inline function Base.getindex(A::AbstractMesh{T,2}, ::Colon,col::V) where {T<:Int,V<:Int}
    collect(getindex(A,1,col):getindex(A,A.dims[1],col))
end

@inline function Base.getindex(A::AbstractMesh{T,2},row::V,::Colon) where {T<:Int,V<:Int}
    collect(getindex(A,row,1):getindex(A,row,A.dims[2]))
end

@inline function Base.getindex(A::AbstractMesh{T,3},::Colon,col::V,::Colon) where {T<:Int,V<:Int}

    @boundscheck begin
        if !(1 <= col <= A.dims[2])
            throw(BoundsError(A, (:,col,:)))
        end
    end

    out  = Vector{T}(undef,A.dims[1]*A.dims[3])

    for zidx = 1:A.dims[3]
        display((zidx-1)*A.dims[1]+1:zidx*A.dims[1])
        out[(zidx-1)*A.dims[1]+1:zidx*A.dims[1]] = collect(getindex(A,1,col,zidx):getindex(A,A.dims[1],col,zidx))
    end
    out
end








################################################################
#               Implement Stencil Like Operations
#
################################################################
struct Stencil{T<:Int,N} <: AbstractArray{Int,N}
    dims::NTuple{N,Int}
    depth::T
end

@inline Base.IndexStyle(::Type{<:Stencil}) = IndexLinear()

@inline function Base.getindex(A::Stencil{T,N}, I::Vararg{Int,N}) where {T<:Int,N} 

end


struct LinearStencil{T<:Int,N}
    up::NTuple{N,T}
    down::NTuple{N,T}
    left::NTuple{N,T}
    right::NTuple{N,T}
    center::T
end

struct NDStencil{T<:Int,N,M}
    up::NTuple{M,NTuple{N,T}}
    down::NTuple{M,NTuple{N,T}}
    left::NTuple{M,NTuple{N,T}}
    right::NTuple{M,NTuple{N,T}}
    center::NTuple{N,T}
end
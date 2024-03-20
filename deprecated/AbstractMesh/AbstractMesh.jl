@inline AbstractMesh(::Type{T}, dims::Int...) where {T<:Int} = AbstractMesh(T, dims)
@inline AbstractMesh(::Type{T}, dims::NTuple{N,Int}) where {T<:Int,N} = AbstractMesh{T,N}(dims, cumprod(dims))

@inline Base.similar(A::AbstractMesh, ::Type{T}, dims::Dims) where {T<:Int} = AbstractMesh(T, dims)
@inline Base.copy(A::AbstractMesh{T}) where {T<:Int} = AbstractMesh(T,A.dims)
@inline Base.size(A::AbstractMesh) = A.dims
@inline Base.length(A::AbstractMesh) = A.multiplied_dims[end]

@inline function Base.getindex(A::AbstractMesh{T,N}, I::Vararg{Int,N}) where {T<:Int,N}
    I = @. mod1(I-1,A.dims)
    index = I[1]
    for i in 2:N
        index += (I[i]-1)*A.multiplied_dims[i-1]
    end
    index
end

@inline function Base.getindex(A::AbstractMesh{T, N}, index::V) where {T<:Int, V<:Int, N}
    
    index = mod1(index-1,length(A))
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

@inline function Base.getindex(A::AbstractMesh{T,2},::Colon) where {T<:Int}
    [repeat(1:A.dims[1]; outer=size(A, 2)) repeat(1:A.dims[2]; inner=size(A, 1))]
end

@inline function Base.getindex(A::AbstractMesh{T,3},::Colon) where {T<:Int}
    [repeat(1:A.dims[1]; outer=size(A, 2) * size(A, 3)) repeat(1:A.dims[2]; inner=size(A, 1), outer=size(A, 3)) repeat(1:A.dims[3]; inner=size(A, 1) * size(A, 2))]
end

@inline function Base.getindex(A::AbstractMesh{T,N},::Colon) where {T<:Int,N}
    error("Unimplemented function")
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
        out[(zidx-1)*A.dims[1]+1:zidx*A.dims[1]] = collect(getindex(A,1,col,zidx):getindex(A,A.dims[1],col,zidx))
    end
    out
end
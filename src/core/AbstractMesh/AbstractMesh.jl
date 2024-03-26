################################################################################################
##                                                                                         
##     Generic Methods for AbstractMesh                             
##
################################################################################################

@inline Base.size(A::T) where {T<:AbstractMesh} = A.dims
@inline Base.length(A::T) where {T<:AbstractMesh} = A.multiplied_dims[end]

################################################################################################
##
##     Periodic Abstract Mesh Methods
##
################################################################################################

@inline PeriodicAbstractMesh(::Type{T}, dims::T) where {T<:Integer} = PeriodicAbstractMesh(T,
                                                                                           dims)

@inline PeriodicAbstractMesh(::Type{T}, dims::NTuple{N,T}) where {T<:Integer,N} = PeriodicAbstractMesh{T,
                                                                                                       N}(dims,
                                                                                                          cumprod(dims))

@inline Base.similar(A::PeriodicAbstractMesh, ::Type{T}, dims::Dims) where {T<:Integer} = PeriodicAbstractMesh(T,
                                                                                                               dims)
@inline Base.copy(A::PeriodicAbstractMesh{T}) where {T<:Integer} = PeriodicAbstractMesh(T,
                                                                                        A.dims)

@inline function Base.getindex(A::PeriodicAbstractMesh{T,N},
                               I::Vararg{Int,N}) where {T<:Integer,N}
    return @inbounds LinearIndices(A.dims)[mod1.(I, A.dims)...]
end

@inline function Base.getindex(A::PeriodicAbstractMesh{T,1}, ::Colon) where {T<:Integer}
    return collect(1:A.dims[1])
end

@inline function Base.getindex(A::PeriodicAbstractMesh{T,2}, ::Colon,
                               col::V) where {T<:Integer,V<:Integer}
    @boundscheck begin
        if !(1 <= col <= 2)
            throw(BoundsError(A, (:, col)))
        end
    end

    return collect(getindex(A, 1, col):getindex(A, A.dims[1], col))
end

@inline function Base.getindex(A::PeriodicAbstractMesh{T,N}, ::Colon,
                               col::V) where {T<:Integer,V<:Integer,N}
    @boundscheck begin
        if !(1 <= col <= A.dims[2])
            throw(BoundsError(A, (:, col)))
        end
    end

    output_size = div(prod(size(A)), size(A, 2))
    out = Vector{T}(undef, output_size)

    for (index, valid_indexes) in enumerate(product(axes(A)[3:N]...))
        @inbounds out[((index - 1) * A.dims[1] + 1):(index * A.dims[1])] = collect(getindex(A,
                                                                                            1,
                                                                                            col,
                                                                                            valid_indexes...):getindex(A,
                                                                                                                       A.dims[1],
                                                                                                                       col,
                                                                                                                       valid_indexes...))
    end
    return out
end

@inline function getindex2(A::PeriodicAbstractMesh{T,N},
                           I::Vararg{T,N}) where {T<:Integer,N}
    index = mod1(I[1], A.dims[1])
    @simd for i in 2:N
        @inbounds index += (mod1(I[i], A.dims[i]) - 1) * A.multiplied_dims[i - 1]
    end
    return index
end

@inline function getindex3(A::PeriodicAbstractMesh{T,N},
                           I::Vararg{T,N}) where {T<:Integer,N}
    I = @. mod1(I, A.dims)
    index = I[1]
    @simd for i in 2:N
        @inbounds index += (I[i] - 1) * A.multiplied_dims[i - 1]
    end
    return index
end
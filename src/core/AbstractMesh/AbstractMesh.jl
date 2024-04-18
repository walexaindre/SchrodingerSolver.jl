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

@inline PeriodicAbstractMesh(::Type{T}, dim::T) where {T<:Integer} = PeriodicAbstractMesh{T,
                                                                                          1}((dim,),
                                                                                             (dim,))

@inline PeriodicAbstractMesh(::Type{T}, dims::NTuple{N,T}) where {T<:Integer,N} = PeriodicAbstractMesh{T,
                                                                                                       N}(dims,
                                                                                                          cumprod(dims))

@inline Base.similar(A::PeriodicAbstractMesh, ::Type{T}, dims::Dims) where {T<:Integer} = PeriodicAbstractMesh(T,
                                                                                                               dims)
@inline Base.copy(A::PeriodicAbstractMesh{T}) where {T<:Integer} = PeriodicAbstractMesh(T,
                                                                                        A.dims)

#@inline Base.getindex(A::PeriodicAbstractMesh{T,1},
#                      I::T) where {T<:Integer} = @inbounds LinearIndices(A.dims)[I]

@inline function Base.getindex(A::PeriodicAbstractMesh{T,N},
                               I::Vararg{Int,N}) where {T<:Integer,N}
    return @inbounds LinearIndices(A.dims)[mod1.(I, A.dims)...]
end

@inline function Base.getindex(A::PeriodicAbstractMesh{T,1},
                               ::Colon) where {T<:Integer}
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

@inline function extract_every_dimension(A::PeriodicAbstractMesh{T,N}) where {T<:Integer,
                                                                              N}
    return (PeriodicAbstractMesh(T, A.dims[idx]) for idx in 1:N)
end


@inline function core_circulant_matrix_format_IJV(col::Vec, offsets::VecO,
                                                  AMesh::AM) where {Vt<:Integer,
                                                                    VecO<:AbstractVector{Vt},
                                                                    Vec<:AbstractVector,
                                                                    AM<:AbstractMesh{Vt}}
    if (length(col) != length(offsets))
        throw(DimensionMismatch("The length of the column vector must be equal to the length of the offsets vector"))
    end

    sz = length(col) * length(AMesh)
    sz_offset = length(offsets)

    I = Vector{V}(undef, sz)
    J = similar(I)
    V = Vector{eltype(col)}(undef, sz)

    @threads for idx in 1:length(Amesh)
        b = sz_offset * idx
        a = b - sz_offset + 1

        pos_to_modify = a:b

        I[pos_to_modify] .= apply_offsets(AMesh, idx, offsets)
        J[pos_to_modify] .= idx
        V[pos_to_modify] .= col
    end

    I, J, V
end

@inline function assembly_circulant_matrix_format_IJV(col::Vec, offsets_vector::VecO,
                                                      AMesh::AM) where {Vt<:Integer,
                                                                        VecO<:Vector{Vt},
                                                                        Vec<:AbstractVector,
                                                                        AM<:AbstractMesh{Vt}}
    offsets = offset_generator(offsets_vector)

    core_circulant_matrix_format_IJV(col, offsets, AMesh)
end

@inline function assembly_circulant_matrix_format_IJV(col::Vec, offsets_range::R,
                                                      AMesh::AM) where {Vt<:Integer,
                                                                        R<:AbstractRange{Vt},
                                                                        Vec<:AbstractVector,
                                                                        AM<:AbstractMesh{Vt}}
    offsets = offset_generator(Vt, offsets_range)

    core_circulant_matrix_format_IJV(col, offsets, AMesh)
end

export apply_offsets, offset_generator, extract_every_dimension,assembly_circulant_matrix_format_IJV,core_circulant_matrix_format_IJV

#Conversions

@inline PeriodicAbstractMesh(P::PeriodicGrid{V,T,R,N}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N} = PeriodicAbstractMesh(V,
                                                                                                                               P.dims)

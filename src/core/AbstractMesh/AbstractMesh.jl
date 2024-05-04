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

@inline function extract_every_dimension(A::PeriodicAbstractMesh{T,N}) where {T<:Integer,
                                                                              N}
    return (PeriodicAbstractMesh(T, A.dims[idx]) for idx in 1:N)
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

export apply_offsets, offset_generator, extract_every_dimension,
       assembly_circulant_matrix_format_IJV, core_circulant_matrix_format_IJV

#Conversions

@inline PeriodicAbstractMesh(P::PeriodicGrid{V,T,R,N}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N} = PeriodicAbstractMesh(V,
                                                                                                                               P.dims)
################################################################################################

"""
    S = matrix_to_vector(M)

Return the dense vector storage type `S` related to the dense matrix storage type `M`.
"""
function matrix_to_vector(::Type{M}) where {M<:DenseMatrix}
    T = hasproperty(M, :body) ? M.body : M
    par = T.parameters
    npar = length(par)
    (2 ≤ npar ≤ 3) || error("Type $M is not supported.")
    if npar == 2
        S = T.name.wrapper{par[1],1}
    else
        S = T.name.wrapper{par[1],1,par[3]}
    end
    return S
end

matrix_to_vector(::Type{M}) where {M<:AbstractCuSparseMatrix} = M.types[3]

matrix_to_vector(::Type{M}) where M<:SparseMatrixCSC = M.types[5]



"""
    v = vundef(S, n)

Create an uninitialized vector of storage type `S` of length `n`.
"""
vundef(S, n) = S(undef, n)


"""
    v = vzeros(S, n)

Create a vector of storage type `S` of length `n` only composed of zero.
"""
vzeros(S, n) = fill!(S(undef, n), zero(eltype(S)))

"""
    v = vones(S, n)

Create a vector of storage type `S` of length `n` only composed of one.
"""
vones(S, n) = fill!(S(undef, n), one(eltype(S)))

"""
    v = vcanonical(S, n, idx)

Create a vector of storage type `S` of length `n` where the position `idx` is one and any other position is zero.
"""
function vcanonical(S, n, idx = 1)
    v = vzeros(S, n)
    v[idx:idx] .= 1
    return v
end

"""
    v = vseq(S, n)

Create a vector of storage type `S` of length `n` where the position `i` is `i`.
"""
function vseq(S,n)
    cumsum(vones(S,n))
end

"""
    v = test_if_zero(V,idx)

Return `true` if the element at position `idx` of the vector `V` is zero.
"""
function test_if_zero end


function test_if_zero(Vec::GPUV,idx) where {GPUV<:CuVector}
    CUDA.@allowscalar if Vec[idx] == zero(eltype(Vec))
        return true
    end
    return false
end

function test_if_zero(Vec::GenV,idx) where {GenV<:AbstractVector}
    if Vec[idx] == zero(eltype(Vec))
        return true
    end
    return false
end

"""
    v = vfirst(V)

Return the first element of the vector `V`.
"""
function vfirst end

function vfirst(Vec::GenV) where {GenV<:AbstractVector}
    first(Vec)
end


function vfirst(Vec::GPUV) where {GPUV<:CuVector}
    CUDA.@allowscalar first(Vec)
end

"""
    v = vlast(V)

Return the last element of the vector `V`.
"""
function vlast end

function vlast(Vec::GenV) where {GenV<:AbstractVector}
    last(Vec)
end

function vlast(Vec::GPUV) where {GPUV<:CuVector}
    CUDA.@allowscalar last(Vec)
end


#=
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
=#
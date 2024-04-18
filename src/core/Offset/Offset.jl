Base.length(o::BaseOffset{V,S}) where {V,S} = length(o.offsets)
Base.size(o::BaseOffset{V,S}) where {V,S} = size(o.offsets)
Base.eltype(o::BaseOffset{V,S}) where {V,S} = V
Base.getindex(o::BaseOffset{V,S}, idx::V) where {V<:Integer,S} = o.offsets[idx]

@inline function check_offset_min(voff)
    if minimum(voff) < 1
        throw(ArgumentError("Start index must be non-negative nor zero"))
    end
end

function AssemblyBaseOffset(Rank::R) where {V<:Integer,R<:AbstractRange{V}}
    check_offset_min(Rank)
    BaseOffset(vcat(0, collect(flatten(zip(Rank, -Rank)))))
end

function AssemblyBaseOffset(::Type{NonZeroOffset},
                            Rank::R) where {R<:AbstractRange{V}} where {V<:Integer}
    check_offset_min(Rank)
    BaseOffset(collect(flatten(zip(Rank, -Rank))))
end

function AssemblyBaseOffset(vec::Vec) where {V<:Integer,Vec<:AbstractVector{V}}
    check_offset_min(vec)
    BaseOffset(vcat(0, collect(flatten(zip(vec, -vec)))))
end

function AssemblyBaseOffset(::Type{NonZeroOffset},
                            vec::Vec) where {V<:Integer,Vec<:AbstractVector{V}}
    check_offset_min(vec)
    BaseOffset(collect(flatten(zip(vec, -vec))))
end

##############################################################################################

Base.ndims(::SymmetricOffset{V,N,OffsetTuple}) where {V<:Integer,N,OffsetTuple} = N
Base.length(o::SymmetricOffset{V,N,OffsetTuple}) where {V<:Integer,N,OffsetTuple} = sum(o.dims)
Base.size(o::SymmetricOffset{V,N,OffsetTuple}) where {V<:Integer,N,OffsetTuple} = o.dims
Base.getindex(o::SymmetricOffset{V,N,OffsetTuple}, odim::V, elemidx::V) where {V<:Integer,N,OffsetTuple} = o.offsets[odim][elemidx]
Base.iterate(S::SymmetricOffset{V,N,OffsetTuple}) where {V<:Integer,N,OffsetTuple} = flatten(S.offsets)

@inline cdims(offsets) =
    ntuple(length(offsets)) do i
        length(offsets[i])
    end

@inline function iterative_only_first_zero(ranks)
    ntuple(length(ranks)) do i
        if i != 1
            return AssemblyBaseOffset(NonZeroOffset, ranks[i])
        else
            return AssemblyBaseOffset(ranks[i])
        end
    end
end

@inline function iterative_dim_sum(dims)
    prev = 1

    ntuple(length(dims)) do i
        if (i == 1)
            return 1
        else
            prev += dims[i - 1]
            return prev
        end
    end
end

@inline function AssemblySymmetricOffset(::Type{NonZeroOffset},
                                         Ranks::Tup) where {Tup}
    offsets = Tuple(AssemblyBaseOffset(NonZeroOffset, rank) for rank in Ranks)
    dims = cdims(offsets)
    dsum = iterative_dim_sum(dims)
    SymmetricOffset(dims, dsum, offsets)
end

@inline function AssemblySymmetricOffset(Ranks::Tup) where {Tup}
    offsets = Tuple(AssemblyBaseOffset(rank) for rank in Ranks)
    dims = cdims(offsets)
    dsum = iterative_dim_sum(dims)
    SymmetricOffset(dims, dsum, offsets)
end

@inline function AssemblySymmetricOffset(::Type{AllZeroOffset},
                                         Ranks::Tup) where {Tup}
    AssemblySymmetricOffset(Ranks)
end

@inline function AssemblySymmetricOffset(::Type{UniqueZeroOffset},
                                         Ranks::Tup) where {Tup}
    offsets = iterative_only_first_zero(Ranks)
    dims = cdims(offsets)
    dsum = iterative_dim_sum(dims)
    SymmetricOffset(dims, dsum, offsets)
end

@inline function apply_offset_by_dim!(out, idx, I, A, offset, dim)
    midx = idx
    for oidx in offset
        tmp = I[dim] + oidx
        J = Base.setindex(I, tmp, dim)
        out[midx] = getindex(A, J...)
        midx += 1
    end
end

@inline function apply_offsets!(out::Vec, start_idx::V, A::PeriodicAbstractMesh{V,N},
                        I::NTuple{N,V},
                        offsets::SymmetricOffset{V,N,OTup}) where {V<:Integer,
                                                                   N,
                                                                   OTup,
                                                                   Vec<:AbstractVector}

    # 0 based offset start
    sidx = start_idx - 1

    for dim in 1:N
        apply_offset_by_dim!(out, sidx + offsets.dsum[dim], I, A,
                             offsets.offsets[dim], dim)
    end
    out
end

@inline function apply_offsets(A::PeriodicAbstractMesh{V,N}, I::NTuple{N,V},
                       offsets::SymmetricOffset{V,N,OTup}) where {V<:Integer,
                                                                  N,
                                                                  OTup}
    offlen = length(offsets)
    out = Vector{V}(undef, offlen)

    apply_offsets!(out, 1, A, I, offsets)
end

##############################################################################################

@inline function core_circulant_matrix_format_IJV(col::Vec,
                                                  SOff::SymmetricOffset{V,N,OTup},
                                                  AMesh::PeriodicAbstractMesh{V,N}) where {V<:Integer,
                                                                                           N,
                                                                                           Vec,
                                                                                           OTup}
    if (length(col) != length(SOff))
        throw(DimensionMismatch("The length of the column vector must be equal to the length of the offsets vector"))
    end

    sz = length(col) * length(AMesh)
    sz_offset = length(SOff)

    _I = Vector{V}(undef, sz)
    _J = similar(I)
    _V = Vector{eltype(col)}(undef, sz)
    
    @threads for idx in CartesianIndices(AMesh)
        b = sz_offset * idx
        a = b - sz_offset + 1

        pos_to_modify = a:b

        _I[pos_to_modify] .= apply_offsets(AMesh, idx, SOff)
        _J[pos_to_modify] .= idx
        _V[pos_to_modify] .= col
    end

    _I, _J, _V
end

export BaseOffset, SymmetricOffset, AssemblyBaseOffset, AssemblySymmetricOffset,
       apply_offsets
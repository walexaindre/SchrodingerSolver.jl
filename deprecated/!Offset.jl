Base.length(o::BaseOffset{V,S}) where {V,S} = length(o.offsets)
Base.size(o::BaseOffset{V,S}) where {V,S} = size(o.offsets)
Base.eltype(o::BaseOffset{V,S}) where {V,S} = V
Base.getindex(o::BaseOffset{V,S}, idx::V) where {V<:Integer,S} = o.offsets[idx]

@inline function check_offset_min(voff)
    if minimum(voff; init = 1) < 1
        throw(ArgumentError("Start index must be non-negative nor zero"))
    end
end

@inline function check_value_range_compatibility(::Type{NonZeroOffset}, rankvalues,
                                                 values)
    if length(rankvalues) * 2 != length(values)
        throw(DimensionMismatch("The length of the offsets vector must be equal to the length of the values vector"))
    end
end

@inline function check_value_range_compatibility(rankvalues, values)
    if length(rankvalues) * 2 + 1 != length(values)
        throw(DimensionMismatch("The length of the offsets vector must be equal to the length of the values vector"))
    end
end

BaseOffset(S::Vec) where {V<:Integer,Vec<:AbstractVector{V}} = BaseOffset(S,
                                                                          vzeros(Vec,
                                                                                 length(S)))

##############################################################################################

function AssemblyBaseOffset(::Type{VecStorage}, ::Type{NonZeroOffset},
                            Data::Container,
                            Values::Vec) where {Tv,V<:Integer,VecStorage,
                                                Container<:VectorOrRank{V},
                                                Vec<:VectorOrRank{Tv}}
    check_offset_min(Data)
    check_value_range_compatibility(NonZeroOffset, Data, Values)
    SType = VecStorage{V}
    SVType = VecStorage{Tv}
    offset = SType(collect(flatten(zip(Data, -Data))))
    BaseOffset(offset, SVType(Values))
end

function AssemblyBaseOffset(::Type{VecStorage}, ::Type{NonZeroOffset},
                            Data::Container,
                            f::Fun) where {V<:Integer,VecStorage,
                                           Container<:VectorOrRank{V},Fun}
    check_offset_min(Data)
    SType = VecStorage{V}
    SVType = VecStorage{eltype(f(vfirst(Data)))}
    offset = SType(collect(flatten(zip(Data, -Data))))
    Values = SVType(f.(offset))
    BaseOffset(offset, Values)
end

function AssemblyBaseOffset(::Type{VecStorage}, Data::Container,
                            Values::Vec) where {Tv,V<:Integer,VecStorage,
                                                Container<:VectorOrRank{V},
                                                Vec<:VectorOrRank{Tv}}
    check_offset_min(Data)
    check_value_range_compatibility(Data, Values)
    SType = VecStorage{V}
    SVType = VecStorage{Tv}
    offset = SType(collect(flatten(zip(Data, -Data))))
    BaseOffset(SType(offset), SVType(Values))
end

function AssemblyBaseOffset(::Type{VecStorage}, Data::Container,
                            f::Fun) where {VecStorage,Fun,V<:Integer,
                                           Container<:VectorOrRank{V}}
    check_offset_min(Data)
    SType = VecStorage{V}
    SVType = VecStorage{eltype(f(0))}
    offset = SType(vcat(0, collect(flatten(zip(Data, -Data)))))
    Values = SVType(f.(offset))
    BaseOffset(offset, Values)
end

##############################################################################################

function AssemblyBaseOffset(::Type{NonZeroOffset},
                            Data::Container) where {V<:Integer,
                                                    Container<:VectorOrRank{V}}
    AssemblyBaseOffset(Vector, NonZeroOffset, Data, Returns(0))
end

function AssemblyBaseOffset(::Type{NonZeroOffset}, Data::Container,
                            Value) where {V<:Integer,Container<:VectorOrRank{V}}
    AssemblyBaseOffset(Vector, NonZeroOffset, Data, Value)
end

function AssemblyBaseOffset(Data::Container) where {V<:Integer,
                                                    Container<:VectorOrRank{V}}
    AssemblyBaseOffset(Vector, Data, Returns(0))
end

function AssemblyBaseOffset(Data::Container,
                            Value) where {V<:Integer,Container<:VectorOrRank{V}}
    AssemblyBaseOffset(Vector, Data, Value)
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

@inline function SymmetricOffset(offsets)
    dims = cdims(offsets)
    dsum = iterative_dim_sum(dims)
    SymmetricOffset(dims, dsum, offsets)
end

@inline function AssemblySymmetricOffset(::Type{NonZeroOffset},
                                         Ranks::Tup) where {Tup}
    offsets = Tuple(AssemblyBaseOffset(NonZeroOffset, rank) for rank in Ranks)
    SymmetricOffset(offsets)
end

@inline function AssemblySymmetricOffset(Ranks::Tup) where {Tup}
    offsets = Tuple(AssemblyBaseOffset(rank) for rank in Ranks)
    SymmetricOffset(offsets)
end

@inline function AssemblySymmetricOffset(::Type{AllZeroOffset},
                                         Ranks::Tup) where {Tup}
    AssemblySymmetricOffset(Ranks)
end

@inline function AssemblySymmetricOffset(::Type{UniqueZeroOffset},
                                         Ranks::Tup) where {Tup}
    offsets = iterative_only_first_zero(Ranks)
    SymmetricOffset(offsets)
end


@inline function extract_all_symmetric_offsets(offsets::SymmetricOffset{V,N,OTup}) where {V<:Integer,
                                                                                       N,
                                                                                       OTup}
    ntuple(N) do i
        SymmetricOffset((offsets.offsets[i],))
    end
end


@inline function apply_offset_by_dim!(out, idx, I, A, offset, dim)
    midx = idx
    for oidx in offset
        tmp = I[dim] + oidx
        J = Base.setindex(I, tmp, dim)
        out[midx] = getindex(A, CartesianIndex(J))
        midx += 1
    end
end

@inline function get_offsets_vector(offsets::SymmetricOffset{V,N,OTup}) where {V<:Integer,
                                                                               N,
                                                                               OTup}
    foldr(vcat, offsets.offsets)
end

@inline function apply_offsets!(out::Vec, start_idx::V, A::PeriodicAbstractMesh{V,N},
                                I::Ind,
                                offsets::SymmetricOffset{V,N,OTup}) where {V<:Integer,
                                                                           N,
                                                                           OTup,
                                                                           Ind<:TupleOrCartesianIndex{N,
                                                                                                      V},
                                                                           Vec<:AbstractVector}

    # 0 based offset start
    sidx = start_idx - 1

    for dim in 1:N
        apply_offset_by_dim!(out, sidx + offsets.dsum[dim], I, A,
                             offsets.offsets[dim], dim)
    end
    out
end

@inline function apply_offsets(A::PeriodicAbstractMesh{V,N}, I::Ind,
                               offsets::SymmetricOffset{V,N,OTup}) where {V<:Integer,
                                                                          N,
                                                                          Ind<:TupleOrCartesianIndex{N,
                                                                                                     V},
                                                                          OTup}
    offlen = length(offsets)
    out = Vector{V}(undef, offlen)

    apply_offsets!(out, 1, A, I, offsets)
end

@inline function offset_generator(::Type{ZOffsetRule}, AMesh::AM,
                                  R::TRank) where {N,V<:Integer,
                                                   ZOffsetRule<:ZeroOffsetRule,
                                                   TRank<:AbstractRange{V},
                                                   AM<:AbstractMesh{V,N}}
    AssemblySymmetricOffset(ZOffsetRule, ntuple(Returns(R), N))
end

@inline function offset_generator(::Type{ZOffsetRule}, AMesh::AM,
                                  RankTuple::TupTRank) where {N,V<:Integer,
                                                              ZOffsetRule<:ZeroOffsetRule,
                                                              TupTRank<:NTuple{N},
                                                              AM<:AbstractMesh{V,N}}
    AssemblySymmetricOffset(ZOffsetRule, RankTuple)
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
    _J = similar(_I)
    _V = Vector{eltype(col)}(undef, sz)
    LinInd = LinearIndices(AMesh)

    @threads for idx in CartesianIndices(AMesh)
        lidx = LinInd[idx]
        b = sz_offset * lidx
        a = b - sz_offset + 1

        pos_to_modify = a:b

        _I[pos_to_modify] .= apply_offsets(AMesh, idx, SOff)
        _J[pos_to_modify] .= lidx
        _V[pos_to_modify] .= col
    end

    _I, _J, _V
end

@inline function span_by_dim(dim::V,
                             SOff::SymmetricOffset{V,N,OTup}) where {V<:Integer,
                                                                     N,
                                                                     OTup}
    (SOff.dsum[dim], SOff.dsum[dim] + SOff.dims[dim] - 1)
end

@inline function offset_to_vector(offset::BaseOffset{V,Tv,S,SV}) where {V,Tv,S,SV}
    haszero = test_if_zero(offset.offsets, 1)

    if haszero
        len = div(length(offset) - 1, 2) + 1
        out = vzeros(S, len)
        out[2:end] .= offset[2:2:end]
    else
        len = div(length(offset), 2)
        out = offset[1:2:end]
    end

    return out
end

"""


"""
@inline function infer_minimal_offsets(oindices::Vec,
                                       SOff::SymmetricOffset{V,N,OTup}) where {V<:Integer,
                                                                               N,
                                                                               OTup,
                                                                               Vec<:AbstractVector{V}}
    offset_vec = get_offsets_vector(SOff)

    offsets = ntuple(N) do dim
        start, stop = span_by_dim(dim, SOff)
        candidate_indices = oindices .>= start .&& oindices .<= stop
        candidates = oindices[candidate_indices]
        normalized_indices = offset_vec[candidates]

        uniqueoffset = Vector(offset_to_vector(SOff.offsets[dim]))

        out = vzeros(Vec, 0)

        for uidx in uniqueoffset
            if uidx == 0
                continue
            end

            leftv = any(i -> i == uidx, normalized_indices)
            rightv = any(i -> i == -uidx, normalized_indices)

            if leftv && rightv
                push!(out, uidx)
            elseif leftv
                @warn "Asymetrical offset detected: $uidx at dimension: $dim minifying step... Dropping offset"
            elseif rightv
                @warn "Asymetrical offset detected: $(-uidx) at dimension: $dim at minifying step... Dropping offset"
            end
        end

        if test_if_zero(uniqueoffset, 1)
            return AssemblyBaseOffset(out)
        else
            return AssemblyBaseOffset(NonZeroOffset, out)
        end
    end

    SymmetricOffset(offsets)
end

export BaseOffset, SymmetricOffset, AssemblyBaseOffset, AssemblySymmetricOffset,
       apply_offsets
@inline Base.ndims(o::BaseOffset{Ti,Tv,M,S,SV}) where {Ti<:Integer,Tv,M,S<:NTuple{M,Ti},SV<:TupleOrVector{Tv}} = 1
@inline Base.length(o::BaseOffset{Ti,Tv,M,S,SV}) where {Ti<:Integer,Tv,M,S<:NTuple{M,Ti},SV<:TupleOrVector{Tv}} = M
@inline Base.size(o::BaseOffset{Ti,Tv,M,S,SV}) where {Ti<:Integer,Tv,M,S<:NTuple{M,Ti},SV<:TupleOrVector{Tv}} = (M,)
@inline Base.getindex(o::BaseOffset{Ti,Tv,M,S,SV}, idx::Ti) where {Ti<:Integer,Tv,M,S<:NTuple{M,Ti},SV<:TupleOrVector{Tv}} = o.offsets[idx]
@inline Base.getindex(o::BaseOffset{Ti,Tv,M,S,SV}, idx::CartesianIndex{1}) where {Ti<:Integer,Tv,M,S<:NTuple{M,Ti},SV<:TupleOrVector{Tv}} = o.offsets[idx]
@inline Base.eltype(o::BaseOffset{Ti,Tv,M,S,SV}) where {Ti<:Integer,Tv,M,S<:NTuple{M,Ti},SV<:TupleOrVector{Tv}} = Ti

function check_symmetry(offset_values)
    fidx = firstindex(offset_values)
    lidx = lastindex(offset_values)

    if length(offset_values) > 1
        #Checking for existence of zero element
        if offset_values[1] == 0
            if iseven(length(offset_values))
                throw(ArgumentError("Symmetric offsets with zero element must have an odd number of elements"))
            end
            fidx += 1
        end
        #Checking for symmetry
        for i in fidx:2:lidx
            if offset_values[i] != -offset_values[i + 1]
                throw(ArgumentError("Symmetric offsets must be symmetrically paired. Status: [$i  $(i + 1)]"))
            end
        end
        #Checking for uniqueness
        if length(unique(offset_values)) != length(offset_values)
            throw(ArgumentError("Symmetric offsets must be unique"))
        end
    end
end

function check_validty(values, offsets)
    if length(values) != length(offsets)
        throw(ArgumentError("Values and offsets must have the same length"))
    end
end

@inline Base.ndims(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {N,Ti<:Integer,Tv,OffsetTuple} = N
@inline Base.length(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {N,Ti<:Integer,Tv,OffsetTuple} = sum(s.dims)
@inline Base.size(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {N,Ti<:Integer,Tv,OffsetTuple} = s.dims
@inline Base.iterate(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {N,Ti<:Integer,Tv,OffsetTuple} = flatten(s.offsetbydim)
@inline Base.eltype(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {N,Ti<:Integer,Tv,OffsetTuple} = Ti

@inline Base.getindex(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}, dim::Ti, index::Ti) where {N,Ti<:Integer,Tv,OffsetTuple} = s.offsetbydim[dim][index]
@inline Base.getindex(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}, dimindex::NTuple{2,Ti}) where {N,Ti<:Integer,Tv,OffsetTuple} = s.offsetbydim[dimindex[1]][dimindex[2]]
@inline Base.getindex(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}, dimindex::CartesianIndex{2}) where {N,Ti<:Integer,Tv,OffsetTuple} = s.offsetbydim[dimindex[1]][dimindex[2]]
@inline Base.getindex(s::SymmetricOffset{N,Ti,Tv,OffsetTuple}, dim::Ti) where {N,Ti<:Integer,Tv,OffsetTuple} = s.offsetbydim[dim]
function Base.show(io::IO,
                   s::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {N,Ti<:Integer,Tv,
                                                                   OffsetTuple}
    header = (["SymmetricOffset"], ["{Dimensions:$N,Ti:$Ti,Tv:$Tv}"], ["Offsets"])
    pretty_table(io, collect(s.offsetbydim); header = header)
end

function (s::SymmetricOffset{N,Ti,Tv,OffsetTuple})(dim::Ti,
                                                   index::Ti) where {N,Ti<:Integer,
                                                                     Tv,
                                                                     OffsetTuple}
    CartesianIndex(ntuple(d -> d != dim ? zero(Ti) : s.offsetbydim[dim][index],
                          Val{N}))
end

function (s::SymmetricOffset{N,Ti,Tv,OffsetTuple})(dimindex::CartesianIndex{2}) where {N,
                                                                                       Ti<:Integer,
                                                                                       Tv,
                                                                                       OffsetTuple}
    CartesianIndex(ntuple(d -> d != dimindex[1] ? zero(Ti) :
                               s.offsetbydim[dimindex[1]][dimindex[2]], Val{N}))
end

@inline cdims(offsets) =
    ntuple(length(offsets)) do i
        length(offsets[i])
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

@inline function GenerateRank(R::RVec, ndims::Ti) where {Ti,RVec<:VectorOrRank{Ti}}
    ntuple(Returns(R), ndims)
end

@inline function GenerateRank(Rank::RTup, ndims::Ti) where {Ti,RTup<:Tuple}
    if ndims != length(Rank)
        throw(ArgumentError("The rank must have the same length as the number of dimensions"))
    end

    Rank
end

@inline function BaseOffset(Values::T,Offsets::T) where {T<:Tuple{}}
    BaseOffset{Int,Int,0,Tuple{},Tuple{}}((),())
end

@inline function SymmetricOffset(values::Tup1,
                                 offsets::Tup2) where {Tup1<:Tuple,Tup2<:Tuple}
    dims = cdims(offsets)
    dsum = iterative_dim_sum(dims)
    offsets = swapsort.(offsets, by = normalize)
    check_symmetry.(offsets)
    check_validty.(values, offsets)

    constructed_base = ntuple(length(dims)) do i
        BaseOffset(offsets[i], values[i])
    end

    SymmetricOffset(dims, dsum, constructed_base)
end

@inline function SymmetricOffset(offsets::Tup1) where {Tup1<:Tuple}
    dims = cdims(offsets)
    dsum = iterative_dim_sum(dims)
    offsets = swapsort.(offsets, by = normalize)
    check_symmetry.(offsets)

    constructed_base = ntuple(length(dims)) do i
        BaseOffset(offsets[i], offsets[i])
    end

    SymmetricOffset(dims, dsum, constructed_base)
end

function GenerateOffset(::Type{OffsetAllZero}, Ranges::Tup1,
                        ValueTup::Tup2) where {Tup1<:Tuple,Tup2<:Tuple}
    offsets = ntuple(length(Ranges)) do i
        vals = collect(Ranges[i])
        return Tuple(unique(vcat(0, vals, -vals)))
    end

    SymmetricOffset(offsets, ValueTup)
end

function GenerateOffset(::Type{OffsetAllZero}, Ranges::Tup1) where {Tup1<:Tuple}
    offsets = ntuple(length(Ranges)) do i
        vals = collect(Ranges[i])
        return Tuple(unique(vcat(0, vals, -vals)))
    end

    SymmetricOffset(offsets)
end

function GenerateOffset(::Type{OffsetNonZero}, Ranges::Tup1) where {Tup1<:Tuple}
    offsets = ntuple(length(Ranges)) do i
        vals = collect(Ranges[i])
        filter!(x -> x != 0, vals)
        return Tuple(unique(vcat(vals, -vals)))
    end
    SymmetricOffset(offsets)
end

function GenerateOffset(::Type{OffsetNonZero}, Ranges::Tup1,
                        ValueTup::Tup2) where {Tup1<:Tuple,Tup2<:Tuple}
    offsets = ntuple(length(Ranges)) do i
        vals = collect(Ranges[i])
        filter!(x -> x != 0, vals)
        return Tuple(unique(vcat(vals, -vals)))
    end

    SymmetricOffset(offsets, ValueTup)
end

function GenerateOffset(::Type{OffsetUniqueZero}, dim,
                        Ranges::Tup1) where {Tup1<:Tuple}
    offsets = ntuple(length(Ranges)) do i
        vals = collect(Ranges[i])

        if i == dim
            return Tuple(unique(vcat(0, vals, -vals)))
        else
            filter!(x -> x != 0, vals)
            return Tuple(unique(vcat(vals, -vals)))
        end
    end

    SymmetricOffset(offsets)
end

function GenerateOffset(::Type{OffsetUniqueZero}, dim,
                        Ranges::Tup1,
                        ValuesTup::Tup2) where {Tup1<:Tuple,Tup2<:Tuple}
    offsets = ntuple(length(Ranges)) do i
        vals = collect(Ranges[i])

        if i == dim
            return Tuple(unique(vcat(0, vals, -vals)))
        else
            filter!(x -> x != 0, vals)
            return Tuple(unique(vcat(vals, -vals)))
        end
    end

    SymmetricOffset(offsets, ValuesTup)
end

@inline function extract_all_symmetric_offsets(offsets::SymmetricOffset{N,Ti,Tv,
                                                                        OTup}) where {N,
                                                                                      Ti<:Integer,
                                                                                      Tv,
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

@inline function apply_offsets!(out::Vec, start_idx::Ti,
                                A::PeriodicAbstractMesh{Ti,N},
                                I::Ind,
                                offsets::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {Ti<:Integer,
                                                                                      Tv,
                                                                                      N,
                                                                                      OffsetTuple,
                                                                                      Ind<:TupleOrCartesianIndex{N,
                                                                                                                 Ti},
                                                                                      Vec<:AbstractVector}

    # 0 based offset start
    sidx = start_idx - 1

    for dim in 1:N
        apply_offset_by_dim!(out, sidx + offsets.dsum[dim], I, A,
                             offsets.offsetbydim[dim], dim)
    end

    out
end

@inline function apply_offsets(A::PeriodicAbstractMesh{Ti,N},
                               I::Ind,
                               offsets::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {Ti<:Integer,
                                                                                     Tv,
                                                                                     N,
                                                                                     OffsetTuple,
                                                                                     Ind<:TupleOrCartesianIndex{N,
                                                                                                                Ti}}
    offlen = length(offsets)
    out = Vector{Ti}(undef, offlen)

    apply_offsets!(out, 1, A, I, offsets)
end

@inline function offset_to_vector(dim::Ti,
                                  offsets::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {Ti<:Integer,
                                                                                        Tv,
                                                                                        N,
                                                                                        OffsetTuple}
    collect(offsets[dim])
end

@inline function span_by_dim(dim::Ti,
                             SOff::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {Ti<:Integer,
                                                                                Tv,
                                                                                N,
                                                                                OffsetTuple}
    (SOff.dsum[dim], SOff.dsum[dim] + SOff.dims[dim] - 1)
end

@inline function infer_minimal_offsets(oindices::Vec,
                                       offsets::SymmetricOffset{N,Ti,Tv,OffsetTuple}) where {Ti<:Integer,
                                                                                             Tv,
                                                                                             N,
                                                                                             OffsetTuple,
                                                                                             Vec<:AbstractVector{Ti}}
    minimal_offsets = ntuple(N) do dim
        vodim = offset_to_vector(dim, offsets)
        start, stop = span_by_dim(dim, offsets)
        candidate_indices = oindices .>= start .&& oindices .<= stop
        candidates = oindices[candidate_indices]
        candidates = candidates .- (start-1)
        normalized_indices = vodim[candidates]
        return Tuple(normalized_indices) #BaseOffset(normalized_indices, candidates)
    end
    SymmetricOffset(minimal_offsets)
end

@inline function core_circulant_matrix_format_COO(col::Vec,
                                                  SOff::SymmetricOffset{N,Ti,Tv,
                                                                        OTup},
                                                  AMesh::PeriodicAbstractMesh{Tv,N}) where {Ti<:Integer,
                                                                                            Tv,
                                                                                            N,
                                                                                            Vec,
                                                                                            OTup}
    if (length(col) != length(SOff))
        throw(DimensionMismatch("The length of the column vector must be equal to the length of the offsets vector"))
    end

    sz = length(col) * length(AMesh)
    sz_offset = length(SOff)

    _I = Vector{Ti}(undef, sz)
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

export BaseOffset, SymmetricOffset, OffsetNonZero, OffsetAllZero, OffsetUniqueZero,
       GenerateOffset, apply_offsets!, infer_minimal_offsets,
       extract_all_symmetric_offsets, core_circulant_matrix_format_COO
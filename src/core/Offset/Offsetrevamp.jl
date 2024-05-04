@inline Base.ndims(o::BaseOffset{M,Ti,Tv,S,SV}) where {M,Ti<:Integer,Tv,S<:NTuple{M,V},SV<:AbstractVector{Tv}} = 1
@inline Base.length(o::BaseOffset{M,Ti,Tv,S,SV}) where {M,Ti<:Integer,Tv,S<:NTuple{M,V},SV<:AbstractVector{Tv}} = M
@inline Base.size(o::BaseOffset{M,Ti,Tv,S,SV}) where {M,Ti<:Integer,Tv,S<:NTuple{M,V},SV<:AbstractVector{Tv}} = (M,)
@inline Base.getindex(o::BaseOffset{M,Ti,Tv,S,SV}, idx::Ti) where {M,Ti<:Integer,Tv,S<:NTuple{M,V},SV<:AbstractVector{Tv}} = o.offsets[idx]
@inline Base.getindex(o::BaseOffset{M,Ti,Tv,S,SV}, idx::CartesianIndex{1}) where {M,Ti<:Integer,Tv,S<:NTuple{M,V},SV<:AbstractVector{Tv}} = o.offsets[idx]
@inline Base.eltype(o::BaseOffset{M,Ti,Tv,S,SV}) where {M,Ti<:Integer,Tv,S<:NTuple{M,V},SV<:AbstractVector{Tv}} = Ti

function check_symmetry(offset_values)
end

function check_offset(offset_values::S) where {M,V,S<:NTuple{M,V}}
end

@inline Base.ndims(s::SymmetricOffset{N,Ti,OffsetTuple}) where {N,Ti<:Integer,OffsetTuple} = N
@inline Base.length(s::SymmetricOffset{N,Ti,OffsetTuple}) where {N,Ti<:Integer,OffsetTuple} = sum(s.dims)
@inline Base.size(s::SymmetricOffset{N,Ti,OffsetTuple}) where {N,Ti<:Integer,OffsetTuple} = s.dims
@inline Base.iterate(s::SymmetricOffset{N,Ti,OffsetTuple}) where {N,Ti<:Integer,OffsetTuple} = flatten(s.offsetbydim)
@inline Base.eltype(s::SymmetricOffset{N,Ti,OffsetTuple}) where {N,Ti<:Integer,OffsetTuple} = Ti

@inline Base.getindex(s::SymmetricOffset{N,Ti,OffsetTuple}, dim::Ti, index::Ti) where {N,Ti<:Integer,OffsetTuple} = s.offsetbydim[dim][index]
@inline Base.getindex(s::SymmetricOffset{N,Ti,OffsetTuple}, dimindex::NTuple{2,Ti}) where {N,Ti<:Integer,OffsetTuple} = s.offsetbydim[dimindex[1]][dimindex[2]]
@inline Base.getindex(s::SymmetricOffset{N,Ti,OffsetTuple}, dimindex::CartesianIndex{2}) where {N,Ti<:Integer,OffsetTuple} = s.offsetbydim[dimindex[1]][dimindex[2]]

function (s::SymmetricOffset{N,Ti,OffsetTuple})(dim::Ti,
                                                index::Ti) where {N,Ti<:Integer,
                                                                  OffsetTuple}
    CartesianIndex(ntuple(d -> d != dim ? zero(Ti) : s.offsetbydim[dim][index],
                          Val{N}))
end

function (s::SymmetricOffset{N,Ti,OffsetTuple})(dimindex::CartesianIndex{2}) where {N,
                                                                                    Ti<:Integer,
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

@inline function SymmetricOffset(offsets)
    dims = cdims(offsets)
    dsum = iterative_dim_sum(dims)
    SymmetricOffset(dims, dsum, offsets)
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
                                offsets::SymmetricOffset{N,Ti,OffsetTuple}) where {V<:Integer,
                                                                                   N,
                                                                                   OffsetTuple,
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
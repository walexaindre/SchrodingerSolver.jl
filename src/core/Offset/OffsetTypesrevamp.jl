TupleOrCartesianIndex{N,V} = Union{NTuple{N,V},
                                   CartesianIndex{N}} where {N,V<:Integer}
TupleOrRange{N,V} = Union{Tuple,AbstractRange{V}} where {V<:Integer}

VectorOrRank = Union{AbstractVector{V},AbstractRange{V}} where {V<:Integer}

abstract type Offset{V} <: AbstractVector{V} end

#Base struct for 1D offsets
struct BaseOffset{M,Ti<:Integer,Tv,S<:NTuple{M,Ti},SV<:AbstractVector{Tv}} <: Offset{Ti}
    offsets::S
    values::SV
end

struct SymmetricOffset{N,Ti<:Integer,OffsetTuple}
    dims::NTuple{N,Ti}
    dsum::NTuple{N,Ti}
    offsetbydim::OffsetTuple
end


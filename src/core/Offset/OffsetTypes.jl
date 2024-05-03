TupleOrCartesianIndex{N,V} = Union{NTuple{N,V},
                                   CartesianIndex{N}} where {N,V<:Integer}
TupleOrRange{N,V} = Union{Tuple,AbstractRange{V}} where {V<:Integer}

VectorOrRank = Union{AbstractVector{V},AbstractRange{V}} where {V<:Integer}

abstract type Offset{V} <: AbstractArray{V,1} end

abstract type ZeroOffsetRule end

"Every dimension has the zero offset"
struct AllZeroOffset <: ZeroOffsetRule end
"No dimension has the zero offset"
struct NonZeroOffset <: ZeroOffsetRule end
"Only the first dimension has the zero offset"
struct UniqueZeroOffset <: ZeroOffsetRule end

#Base struct for 1D offsets
struct BaseOffset{V<:Integer,Tv,S<:AbstractVector{V},SV<:AbstractVector{Tv}} <: Offset{V}
    offsets::S
    values::SV
end


struct SymmetricOffset{V<:Integer,N,OffsetTuple}
    dims::NTuple{N,V}
    dsum::NTuple{N,V}
    offsets::OffsetTuple
end

export Offset, BaseOffset, SymmetricOffset, AllZeroOffset, UniqueZeroOffset

include("OffsetTypesrevamp.jl")
include("Offsetrevamp.jl")
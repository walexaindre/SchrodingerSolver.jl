TupleOrCartesianIndex{N,V} = Union{NTuple{N,V},
                                   CartesianIndex{N}} where {N,V<:Integer}
TupleOrRange{N,V} = Union{Tuple,AbstractRange{V}} where {V<:Integer}

TupleOrVector{Tv} = Union{Tuple{Vararg{Tv}},AbstractVector{Tv}}
 
VectorOrRank = Union{AbstractVector{V},AbstractRange{V}} where {V<:Integer}

abstract type Offset{V} <: AbstractVector{V} end

#Base struct for 1D offsets
struct BaseOffset{Ti<:Integer,Tv,M,S<:NTuple{M,Ti},SV<:TupleOrVector{Tv}} <: Offset{Ti}
    offsets::S
    values::SV
end


struct SymmetricOffset{N,Ti<:Integer,Tv,OffsetTuple<:NTuple{N,BaseOffset{Ti,Tv}}}
    dims::NTuple{N,Ti}
    dsum::NTuple{N,Ti}
    offsetbydim::OffsetTuple
end

abstract type OffsetRule end

"Every dimension has the zero offset"
struct OffsetNonZero <: OffsetRule end
"No dimension has the zero offset"
struct OffsetAllZero <: OffsetRule end
"Only the first dimension has the zero offset"
struct OffsetUniqueZero <: OffsetRule end
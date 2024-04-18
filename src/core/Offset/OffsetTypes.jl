abstract type Offset{V} <: AbstractArray{V,1} end

"Every dimension has the zero offset"
struct AllZeroOffset end
"No dimension has the zero offset"
struct NonZeroOffset end
"Only the first dimension has the zero offset"
struct UniqueZeroOffset end

#Base struct for 1D offsets
struct BaseOffset{V<:Integer,S<:AbstractVector{V}} <: Offset{V} 
    offsets::S
end

struct SymmetricOffset{V<:Integer,N,OffsetTuple}
    dims::NTuple{N, V}
    dsum::NTuple{N, V}
    offsets::OffsetTuple
end

export Offset, BaseOffset, SymmetricOffset, AllZeroOffset, UniqueZeroOffset

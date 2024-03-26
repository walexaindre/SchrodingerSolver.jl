abstract type AbstractGrid{T,N} <: AbstractArray{T,N} end

struct PeriodicGrid{V<:Integer,T<:Real,Range<:AbstractRange{T},N} <:
       AbstractGrid{Tuple{V,V},N}
    ranges::NTuple{N,Range}
    dims::NTuple{N,V}
    h::NTuple{N,T}
    τ::T    
end

struct GhostPeriodicGrid{V<:Integer,T<:Real,Range<:AbstractRange{T},N} <: AbstractGrid{V,N}
    ranges::NTuple{N,Range}
    dims::NTuple{N,V}
    h::NTuple{N,T}
    τ::T
    depth::T
end

export AbstractGrid, PeriodicGrid, GhostPeriodicGrid
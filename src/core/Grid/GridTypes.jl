abstract type AbstractGrid{T<:Integer,N} <: AbstractArray{T,N} end

struct PeriodicGrid{T<:Integer,Range<:AbstractRange{T},N} <: AbstractGrid{T,N}
    ranges::NTuple{N,Range}
    dims::NTuple{N,T}
    h::NTuple{N,T}
    τ::T
end

struct GhostPeriodicGrid{T<:Integer,Range<:AbstractRange{T},N} <: AbstractGrid{T,N}
    ranges::NTuple{N,Range}
    dims::NTuple{N,T}
    h::NTuple{N,T}
    τ::T
    depth::T
end
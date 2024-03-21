abstract type AbstractMesh{T<:Integer,N} <: AbstractArray{T,N} end

# An implementation of an abstract N dimensional mesh.
# Here we can map from the N dimensional index to the 1D index and vice versa.

struct PeriodicAbstractMesh{T<:Integer,N} <: AbstractMesh{T,N}
    dims::NTuple{N,T}
    multiplied_dims::NTuple{N,T}
end

struct GhostPeriodicAbstractMesh{T<:Integer,N} <: AbstractMesh{T,N}
    dims::NTuple{N,T}
    multiplied_dims::NTuple{N,T}
    depth::T
end
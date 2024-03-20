# An implementation of an abstract N dimensional mesh.
# Here we can map from the N dimensional index to the 1D index and vice versa.
# 

struct AbstractMesh{T<:Int,N} <: AbstractArray{T,N}
    dims::NTuple{N,Int}
    multiplied_dims::NTuple{N,Int}
end
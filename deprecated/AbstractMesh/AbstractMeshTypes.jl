struct AbstractMesh{T<:Int,N} <: AbstractArray{T,N}
    dims::NTuple{N,Int}
    multiplied_dims::NTuple{N,Int}
end
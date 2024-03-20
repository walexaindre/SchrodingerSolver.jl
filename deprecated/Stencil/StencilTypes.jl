struct LinearOffsetStencil{T<:Int,N} <: AbstractArray{T,N}
    dims::NTuple{N,T}
    multiplied_dims::NTuple{N,T}
    offset_dims::NTuple{N,T}
end

const Forward = true
const Backward = false

struct StencilIterator{T<:Int,M,R} <: AbstractArray{T,3}
    dims::NTuple{3,T}

    linear_index::NTuple{M,T}
    iteration_sequence::NTuple{3,T} # N,M,R
    iteration_order::NTuple{3,Bool} # Forward/Backward
    linear_index_multipliers::NTuple{R,T} # Integer multipliers for linear index

    dim_count::T # Numbers of dimensions
    r::T # Elems per linear index
end


struct Stencil{T<:Int,N,M,R} <: AbstractArray{T,2}
    dims::NTuple{2,T} #Stencil dims in some given permutation of N,M,R (By default 3,2,1)
    mesh::AbstractMesh{T,N} #Abstract Mesh to calculate offsets with given dims
    stencil_iterator::StencilIterator{T,M,R} #Iterator to get offsets in any dimension based on linear indexes.
end

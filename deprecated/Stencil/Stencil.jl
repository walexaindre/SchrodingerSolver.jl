#LinearOffsetStencil with broadcasting

Base.IndexStyle(::Type{LinearOffsetStencil}) = IndexLinear()
@inline LinearOffsetStencil(::Type{T}, dims::NTuple{N,T}, offset_dims::NTuple{N,T}) where {T<:Int,N} = LinearOffsetStencil(dims,
                                                                                                                           cumprod(dims),
                                                                                                                           offset_dims)
@inline Base.copy(A::LinearOffsetStencil{T}) where {T<:Int} = LinearOffsetStencil(T, A.dims,
                                                                                  A.offset_dims)
@inline Base.size(A::LinearOffsetStencil) = A.dims
@inline Base.length(A::LinearOffsetStencil) = A.multiplied_dims[end]

@inline function Base.getindex(A::LinearOffsetStencil{T,N},
                               I::Vararg{T,N}) where {T<:Int,N}
    Rank = mod1.(I .+ A.offset_dims .- 1, A.dims)
    index = Rank[1]

    for i in 2:N
        index += (Rank[i] - 1) * A.multiplied_dims[i - 1]
    end

    return index
end

#Generic offset kernel to make unitary calculations with given offset
@inline function offsetkernel(mesh::AbstractMesh{T,N}, offset::NTuple{N,T},
                              I::NTuple{N,T}) where {T<:Int,N}
    Rank = mod1.(I .+ offset .- 1, mesh.dims)

    index = Rank[1]

    @inbounds for i in 2:N
        index += (Rank[i] - 1) * mesh.multiplied_dims[i - 1]
    end

    return index
end

#Generic offset kernel to make unitary calculations with given offset
@inline function offsetkernel2(mesh::AbstractMesh{T,N}, offset::NTuple{N,T},
                               I::NTuple{N,T}) where {T<:Int,N}
    index = mod1(I[1] + offset[1] - 1, mesh.dims[1])

    @inbounds for i in 2:N
        index += mod(I[i] + offset[i] - 1, mesh.dims[i]) * mesh.multiplied_dims[i - 1]
    end

    return index
end

#StencilIterator with inmutable struct

@inline function calculatedims(dim_count::T, linear_index_multipliers::NTuple{R,T},
                               iteration_sequence::NTuple{3,T},
                               linear_index::NTuple{M,T}) where {T<:Int,M,R}
    result = MVector{3,T}(undef)
    #Defining N
    result[iteration_sequence[1]] = dim_count
    #Defining M
    result[iteration_sequence[2]] = M
    #Defining R
    result[iteration_sequence[3]] = R
    return Tuple(result)
end

@inline StencilIterator(::Type{T}, dim_count::T, linear_index::NTuple{M,T}, iteration_sequence::NTuple{3,T}=ntuple(i -> i, 3), iteration_order::NTuple{3,Bool}=ntuple(i -> Forward, 3), linear_index_multipliers::NTuple{R,T}=tuple(-1, 1)) where {T<:Int,M,R} = StencilIterator(calculatedims(dim_count,
                                                                                                                                                                                                                                                                                               linear_index_multipliers,
                                                                                                                                                                                                                                                                                               iteration_sequence,
                                                                                                                                                                                                                                                                                               linear_index),
                                                                                                                                                                                                                                                                                 linear_index,
                                                                                                                                                                                                                                                                                 iteration_sequence,
                                                                                                                                                                                                                                                                                 iteration_order,
                                                                                                                                                                                                                                                                                 linear_index_multipliers,
                                                                                                                                                                                                                                                                                 dim_count,
                                                                                                                                                                                                                                                                                 R)

@inline Base.size(A::StencilIterator) = A.dims

@inline Base.IndexStyle(::Type{StencilIterator}) = IndexLinear()

@inline function set_value(cidx::T, Nv::T, Mv::T, Rv::T) where {T<:Int}
    if cidx != Nv
        return zero(T)
    else
        return Mv * Rv
    end
end

@inline function Base.getindex(A::StencilIterator{T,M,R}, I::Vararg{T,3}) where {T<:Int,M,R}
    ON = A.iteration_order[1]
    OM = A.iteration_order[2]
    OR = A.iteration_order[3]

    Nidx = A.iteration_sequence[1]
    Midx = A.iteration_sequence[2]
    Ridx = A.iteration_sequence[3]

    Nmax = A.dims[Nidx]
    Mmax = A.dims[Midx]
    Rmax = A.dims[Ridx]

    Nout = ON == Forward ? I[Nidx] : Nmax + 1 - I[Nidx]
    Mvalue = OM == Forward ? I[Midx] : Mmax + 1 - I[Midx]
    Rvalue = OR == Backward ? I[Ridx] : Rmax + 1 - I[Ridx]

    Mout = A.linear_index[Mvalue]
    Rout = A.linear_index_multipliers[Rvalue]

    return ntuple(index -> set_value(index, Nout, Mout, Rout), A.dim_count)
end

#Stencil with inmutable struct

@inline function validate_sequence(iteration_sequence::NTuple{3,T}) where {T<:Int}
    if iteration_sequence[1] == iteration_sequence[2] ||
       iteration_sequence[1] == iteration_sequence[3] ||
       iteration_sequence[2] == iteration_sequence[3]
        throw(ArgumentError("Iteration sequence must have unique values"))
    end

    if all(i -> 1 <= i <= 3, iteration_sequence) == false
        throw(ArgumentError("Iteration sequence must have values between 1 and 3"))
    end

    return iteration_sequence
end

@inline Stencil(::Type{T}, mesh::AbstractMesh{T,N}, linear_index::NTuple{M,T}, iteration_sequence::NTuple{3,T}=(3, 2, 1), iteration_order::NTuple{3,Bool}=(Forward, Forward, Forward), linear_index_multipliers::NTuple{R,T}=(-1, 1)) where {T<:Int,N,M,R} = Stencil((length(mesh),
                                                                                                                                                                                                                                                                      R *
                                                                                                                                                                                                                                                                      M *
                                                                                                                                                                                                                                                                      N +
                                                                                                                                                                                                                                                                      1),
                                                                                                                                                                                                                                                                     mesh,
                                                                                                                                                                                                                                                                     StencilIterator(T,
                                                                                                                                                                                                                                                                                     N,
                                                                                                                                                                                                                                                                                     linear_index,
                                                                                                                                                                                                                                                                                     validate_sequence(iteration_sequence),
                                                                                                                                                                                                                                                                                     iteration_order,
                                                                                                                                                                                                                                                                                     linear_index_multipliers))

@inline Base.size(A::Stencil) = A.dims
@inline Base.IndexStyle(::Type{Stencil}) = IndexLinear()
@inline Base.length(A::Stencil) = length(A.mesh) * (length(A.stencil_iterator) + 1)

@inline function Base.getindex(A::Stencil{T,N,M}, I::Vararg{T,2}) where {T<:Int,N,M}
    if I[2] == 1
        return I[1]
    else
        return offsetkernel(A.mesh, A.stencil_iterator[I[2] - 1],
                            getindex(A.mesh, I[1]))
    end
end

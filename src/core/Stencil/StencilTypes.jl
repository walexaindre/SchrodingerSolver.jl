struct LinearOffsetStencil{T<:Int,R<:AbstractRange{T},N} <: AbstractArray{T,N}
    dims::NTuple{N,T}
    ranges::NTuple{N,R}
    multiplied_dims::NTuple{N,T}
    offset_dims::NTuple{N,T}
end

function getrange(vend::T) where {T<:Int}
    return 1:vend
end

function getranges(dims::NTuple{N,T}) where {T<:Int,N}
    return getrange.(dims)
end

Base.IndexStyle(::Type{LinearOffsetStencil}) = IndexLinear()
@inline LinearOffsetStencil(::Type{T}, dims::NTuple{N,T}, offset_dims::NTuple{N,T}) where {T<:Int,N} = LinearOffsetStencil(dims,
                                                                                                                           getranges(dims),
                                                                                                                           cumprod(dims),
                                                                                                                           offset_dims)
@inline Base.copy(A::LinearOffsetStencil{T}) where {T<:Int} = LinearOffsetStencil(T, A.dims,
                                                                                  A.ranges,
                                                                                  A.offset_dims)
@inline Base.size(A::LinearOffsetStencil) = A.dims
@inline Base.length(A::LinearOffsetStencil) = A.multiplied_dims[end]

@inline function Base.getindex(A::LinearOffsetStencil{T,R,N},
                               I::Vararg{T,N}) where {T<:Int,R<:AbstractRange{T},N}
    Rank = getindex.(A.ranges, mod.(I .+ A.offset_dims .- 1, A.dims) .+ 1)

    index = Rank[1]
    
    for i in 2:N
        index += (Rank[i] - 1) * A.multiplied_dims[i - 1]
    end

    return index
end










struct StencilPattern{V<:Int,N}
    dims::V
    pattern_values::NTuple{N,V}
    with_center::Bool
    with_center_per_layer::Bool

    @inline StencilPattern(dims::V, one_direction_pattern::NTuple{N,V}, with_center::Bool=true, with_center_per_layer::Bool=false) where {V<:Int,N} = new{V,
                                                                                                                                                          N}(dims,
                                                                                                                                                             one_direction_pattern,
                                                                                                                                                             with_center,
                                                                                                                                                             with_center_per_layer)
end

@inline StencilPattern(dims::V, one_direction_pattern::V...) where {V<:Int} = StencilPattern(dims,
                                                                                             one_direction_pattern)

@inline function calculate()
end

@inline function getpattern(pattern::StencilPattern)
    return permutedims(findall(!iszero, pattern.pattern_values))
end

@inline function left(pattern::StencilPattern)
    dims = pattern.dims
    NTuple{dims,Int}(0)
end

@inline function getoffsets(pattern::StencilPattern)
    stencil_pattern = getpattern(pattern)
    padding = zeros(Int, length(stencil_pattern))
    #Offset return order is up, down, left and right in increasing order...



    product(stencil_pattern, padding)
end

@inline function Base.length(S::StencilPattern{V,N}) where {V<:Int,N}
    if S.with_center_per_layer
        return (S.dims + 1) * N
    elseif S.with_center
        return S.dims * N + 1
    end
    return S.dims * N
end

struct Stencil{T<:Int,N}
    pattern::StencilPattern{T}
    mesh::AbstractMesh{T,N}
end

function to_matrix(::Type{M}, stencil::Stencil) where {M<:AbstractMatrix}
    pattern = stencil.pattern
    mesh = stencil.mesh
    ncols = length(pattern)

    stencil_array = to_op = mesh[:]

    out = M(undef, length(mesh), ncols)

    if pattern.with_center_per_layer
        error("Missing Implementation")
    elseif pattern.with_center
        out[:, 1] .= collect(1:length(mesh))
    end

end
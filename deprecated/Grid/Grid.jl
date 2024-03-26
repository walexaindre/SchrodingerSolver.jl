struct Grid{T<:AbstractFloat,Range<:AbstractRange{T},N} <: AbstractArray{Int,N}
    ranges::NTuple{N,Range}
    dims::NTuple{N,Int}
    h::NTuple{N,T}
    τ::T
end

@inline range_start_step_stop(start, step, stop) = start:step:stop
@inline range_start_stop_length(start, stop, length) = range(start, stop; length=length)

@inline Grid(::Type{T}, τ::T, h::NTuple{N,T}, ranges::NTuple{N,R}, dims::NTuple{N,Int}) where {T,R<:AbstractRange{T},N} = Grid{T,
                                                                                                                               R,
                                                                                                                               N}(ranges,
                                                                                                                                  dims,
                                                                                                                                  h,
                                                                                                                                  τ)
@inline Grid(::Type{T}, ranges::NTuple{N,R}) where {T,R<:AbstractRange{T},N} = Grid(T,
                                                                                    zero(T),
                                                                                    Base.step.(ranges),
                                                                                    ranges,
                                                                                    length.(ranges))
@inline Grid(::Type{T}, ranges::R...) where {T,R<:AbstractRange{T}} = Grid(T, ranges)

@inline Grid(::Type{T}, τ::T, h::NTuple{N,T}, ranges::NTuple{N,R}) where {T,R<:AbstractRange{T},N} = Grid(T,
                                                                                                          τ,
                                                                                                          h,
                                                                                                          ranges,
                                                                                                          length.(ranges))
@inline Grid(::Type{T}, τ::T, h::NTuple{N,T}, ranges::R...) where {T,R<:AbstractRange{T},N} = Grid(T,
                                                                                                   τ,
                                                                                                   h,
                                                                                                   ranges)
@inline Grid(::Type{T}, h::NTuple{N,T}, ranges::R...) where {T,R<:AbstractRange{T},N} = Grid(T,
                                                                                             h,
                                                                                             ranges)
@inline Grid(::Type{T}, h::NTuple{N,T}, ranges::NTuple{N,R}) where {T,R<:AbstractRange{T},N} = Grid(T,
                                                                                                    zero(T),
                                                                                                    h,
                                                                                                    ranges)

@inline Grid(::Type{T}, τ::T, ranges::NTuple{N,R}) where {T,R<:AbstractRange{T},N} = Grid(T,
                                                                                          τ,
                                                                                          Base.step.(ranges),
                                                                                          ranges)
@inline Grid(::Type{T}, τ::T, ranges::R...) where {T,R<:AbstractRange{T}} = Grid(T, τ,
                                                                                 ranges)

@inline Grid(::Type{T}, boundary_lower::NTuple{N,T}, boundary_upper::NTuple{N,T}, h::NTuple{N,T}) where {T,N} = Grid(T,
                                                                                                                     h,
                                                                                                                     range_start_step_stop.(boundary_lower,
                                                                                                                                            h,
                                                                                                                                            boundary_upper))
@inline Grid(::Type{T}, boundary_lower::NTuple{N,T}, boundary_upper::NTuple{N,T}, steps::NTuple{N,Int}) where {T,N} = Grid(T,
                                                                                                                           range_start_stop_length.(boundary_lower,
                                                                                                                                                    boundary_upper,
                                                                                                                                                    steps))

@inline Base.size(A::Grid) = A.dims
@inline Base.copy(A::Grid{T}) where {T<:Int} = Grid(T, A.ranges, A.dims)
@inline Base.convert(::Type{AbstractMesh}, A::Grid) = AbstractMesh(Int, A.dims)

#Evaluate every range at the given index I
@inline function Base.getindex(A::Grid{T,R,N},
                               I::Vararg{Int,N}) where {T,R<:AbstractRange{T},N}
    return getindex.(A.ranges, mod1.(I .- 1, A.dims))
end

@inline function Base.getindex(A::Grid{T,R,1}, ::Colon) where {T,R<:AbstractRange{T}}
    return collect(A.ranges[1])
end

@inline function Base.getindex(A::Grid{T,R,2}, ::Colon,
                               col::V) where {T,R<:AbstractRange{T},V<:Integer}
    @boundscheck begin
        if !(1 <= col <= 2)
            throw(BoundsError(A, (:, col)))
        end
    end

    if col == 1
        return repeat(A.ranges[1]; outer=size(A, 2))
    else
        return repeat(A.ranges[2]; inner=size(A, 1))
    end
end

mod

@inline function Base.getindex(A::Grid{T,R,3}, ::Colon,
                               col::V) where {T,R<:AbstractRange{T},V<:Integer}
    @boundscheck begin
        if !(1 <= col <= 3)
            throw(BoundsError(A, (:, col)))
        end
    end

    if col == 1
        return repeat(A.ranges[1]; outer=size(A, 2) * size(A, 3))
    elseif col == 2
        return repeat(A.ranges[2]; inner=size(A, 1), outer=size(A, 3))
    else
        return repeat(A.ranges[3]; inner=size(A, 1) * size(A, 2))
    end
end

@inline function Base.getindex(A::Grid{T,R,N}, ::Colon,
                               col::V) where {T,R<:AbstractRange{T},N,V<:Integer}
    @boundscheck begin
        if !(1 <= col <= N)
            throw(BoundsError(A, (:, col)))
        end
    end

    return repeat(A.ranges[col]; inner=prod(size(A)[begin:(col - 1)]),
                  outer=prod(size(A)[(col + 1):end]))
end

export Grid
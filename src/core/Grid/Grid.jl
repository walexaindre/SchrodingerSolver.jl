################################################################################################
##                                                                                         
##     Generic Methods for AbstractGrid                           
##
################################################################################################

@inline Base.size(A::T) where {T<:AbstractGrid} = A.dims
@inline Base.length(A::T) where {T<:AbstractGrid} = A.multiplied_dims[end]

@inline range_start_step_stop(start, step, stop) = start:step:stop
@inline range_start_stop_length(start, stop, length) = range(start, stop; length=length)
@inline range_length(rank) = Base.step(rank)

################################################################################################
##                                                                                         
##     Periodic Grid Methods                           
##
################################################################################################

@inline function Base.getindex(A::PeriodicGrid{T,R,N},
                               I::Vararg{Int,N}) where {T<:Integer,R<:AbstractRange{T},N}
    return @inbounds getindex.(A.ranges, mod1.(I .- 1, A.dims))
end

@inline function Base.getindex(A::Grid{T,R,2}, ::Colon,
                               col::V) where {T,R<:AbstractRange{T},V<:Int}
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

@inline function Base.getindex(A::Grid{T,R,3}, ::Colon, col::V,
                               ::Colon) where {T,R<:AbstractRange{T},V<:Int}
    @boundscheck begin
        if !(1 <= col <= 3)
            throw(BoundsError(A, (:, col, :)))
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
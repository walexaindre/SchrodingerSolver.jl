################################################################################################
##                                                                                         
##     Generic Methods for AbstractGrid                           
##
################################################################################################

@inline Base.size(A::T) where {T<:AbstractGrid} = A.dims
@inline Base.length(A::T) where {T<:AbstractGrid} = prod(A.dims)

@inline range_start_step_stop(start, step, stop) = start:step:stop
@inline range_start_stop_step(start, stop, step) = start:step:stop
@inline range_start_stop_step(start_stop, step) = start_stop[1]:step:start_stop[2]
@inline range_start_stop_length(start, stop, length) = range(start, stop; length=length)
@inline range_length(rank) = Base.step(rank)

################################################################################################
##                                                                                         
##     Periodic Grid Methods                           
##
################################################################################################

#TODO Perform integrity checks...

@inline PeriodicGrid(::Type{V}, ::Type{T}, τ::T, ranges::NTuple{N,R}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N} = PeriodicGrid(ranges,
                                                                                                                                      size.(ranges,
                                                                                                                                            1),
                                                                                                                                      range_length.(ranges),
                                                                                                                                      τ)
@inline PeriodicGrid(::Type{V}, ::Type{T}, τ::T, ranges::R...) where {V<:Integer,T<:Real,R<:AbstractRange{T}} = PeriodicGrid(V,
                                                                                                                             T,
                                                                                                                             τ,
                                                                                                                             ranges)

@inline function PeriodicGrid(::Type{V}, ::Type{T}, τ::T, h::NTuple{N,T},
                              boundaries::NTuple{N,Tuple{T,T}}) where {V<:Integer,T<:Real,N}
    ranges = Tuple(range_start_stop_step(boundaries[idx], h[idx]) for idx in 1:N)
    dims = Tuple(size(ranges[idx], 1) for idx in 1:N)

    return PeriodicGrid(ranges,
                        dims,
                        h,
                        τ)
end

@inline PeriodicGrid(::Type{V}, ::Type{T}, τ::T, h::NTuple{N,T}, boundaries::Tuple{T,T}...) where {V<:Integer,T<:Real,N} = PeriodicGrid(V,
                                                                                                                                        T,
                                                                                                                                        τ,
                                                                                                                                        h,
                                                                                                                                        boundaries)


                                                                                                                                        
@inline function Base.getindex(A::PeriodicGrid{V,T,R,N},
                               I::Vararg{V,N}) where {V<:Integer,T<:Real,
                                                      R<:AbstractRange{T},N}
    #println("IS: ", I)
    @inbounds getindex.(A.ranges, mod1.(I, A.dims))
end


@inline function Base.getindex(A::PeriodicGrid{V,T,R,N}, ::Colon,
                               col::V) where {V<:Integer,T<:Real,R<:AbstractRange{T},N}
    @boundscheck begin
        if !(1 <= col <= N)
            throw(BoundsError(A, (:, col)))
        end
    end

    #total_product = prod(size(A))
    left_slice = 1:(col - 1)
    right_slice = (col + 1):N

    #This product is handy because prod([]) = 1 (product of an empty collection is 1) 
    #and a range has dimension 1
    inner_elems = prod(size(A)[left_slice])
    outer_elems = prod(size(A)[right_slice])

    return repeat(A.ranges[col]; inner=inner_elems, outer=outer_elems)
end

@inline Base.convert(::Type{PeriodicAbstractMesh{V,N}}, A::PeriodicGrid{V,T,R,N}) where {V<:Integer,T<:Real,R<:AbstractRange{T},N} = PeriodicAbstractMesh(V,A.dims)

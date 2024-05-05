
@inline normalize(x::T) where {T} = x<=zero(T) ? one(T)-x-x : x+x
@inline min_max(x,y,lt=isless,by=identity) = lt(by(x),by(y)) ? (x,y) : (y,x)

const LIMIT_SIZE = 16


swapsort(x::T;lt=isless,by=identity) where {T} = x
swapsort(x::Tuple{1,T};lt=isless,by=identity) where {T} = (x,)

@inline function swapsort(vec::Vec;lt=isless,by=identity) where {Vec<:AbstractVector}
   n = length(vec)
   n>LIMIT_SIZE && throw(ArgumentError("Vector size $n exceeds the limit $LIMIT_SIZE"))
   swapsort(Val{n},vec;lt=lt,by=by)
end

include("swapsort.jl")

@inline normalize(x::T) where {T} = x<=zero(T) ? one(T)-x-x : x+x
@inline min_max(x,y,lt=isless,by=identity) = lt(by(x),by(y)) ? (x,y) : (y,x)

@inline function normalize(x::T) where {T}
   y=x+x
   if x<=0
    return 1-y 
   end
   return y
end

swapsort(x::T) where {T} = x
swapsort(x::Tuple{1,T}) where {T} = x

include("swapsort.jl")
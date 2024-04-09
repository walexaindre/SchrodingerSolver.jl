abstract type SpaceDiscretization end

AbstractFloatOrRational{V} = Union{AbstractFloat,Rational{V}} where {V<:Integer}

struct SecondDerivativeCoefficients{V <: Integer,T <: AbstractFloatOrRational{V}}
    a::T
    b::T
    c::T
    α::T
    β::T
    order::V
end

export SecondDerivativeCoefficients, SpaceDiscretization
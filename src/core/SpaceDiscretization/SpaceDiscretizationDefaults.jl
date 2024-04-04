#=
    if order == 2
        α = 0 // 1
        β = 0 // 1

        a = 1 // 1
        b = 0 // 1
        c = 0 // 1

    elseif order == 4
        α = 1 // 10
        β = 0 // 1

        a = (4 // 3) * (1 - α)
        b = (1 // 3) * (-1 + 10 * α)
        c = 0 // 1
    elseif order == 6
        α = 2 // 11
        β = 0 // 1

        a = 12 // 11
        b = 3 // 11
        c = 0 // 1
    elseif order == 8
        α = 344 // 1179
        β = (38 * α - 9) // 214

        a = (696 - 1191 * α) // 428
        b = (2454 * α - 294) // 535
        c = (1179 * α - 344) // 2140

    elseif order == 10
        α = 334 // 899
        β = 43 // 1798

        a = 1065 // 1798
        b = 1038 // 899
        c = 79 // 1798
    end
=#

const SpaceDiscretizationDefaults::Dict{Symbol,SecondDerivativeCoefficients{Int,Rational{Int}}} = Dict(
    :ord2 => SecondDerivativeCoefficients(1, 0, 0, 0, 0, 2),
    :ord4 => SecondDerivativeCoefficients(4//3, 1//3, 0, 1//10, 0, 4),
    :ord6 => SecondDerivativeCoefficients(12//11, 3//11, 0, 2//11, 0, 6),
    :ord8 => SecondDerivativeCoefficients(696//428, 2454//535, 1179//2140, 344//1179, 38//214, 8),
    :ord10 => SecondDerivativeCoefficients(1065//1798, 1038//899, 79//1798, 334//899, 43//1798, 10)
)

function register(sym::Symbol, a::T, b::T, c::T, α::T, β::T,
                  order::V) where {T<:Real,V<:Integer}
    Type = Rational{Int}

    SpaceDiscretizationDefaults[sym] = SecondDerivativeCoefficients(Type(a), Type(b),
                                                                    Type(c), Type(α),
                                                                    Type(β),
                                                                    Int(order))
end

unregister(sym::Symbol) = delete!(SpaceDiscretizationDefaults, sym)

get_available_defaults() = keys(SpaceDiscretizationDefaults)

register()

export register, unregister, get_available_defaults

function default_derivative_coefficient_order(order)
    if order == 2
        α = 0 // 1
        β = 0 // 1

        a = 1 // 1
        b = 0 // 1
        c = 0 // 1

        return SpaceDiscretization(a, b, c, α, β, order)

    elseif order == 4
        α = 1 // 10
        β = 0 // 1

        a = (4 // 3) * (1 - α)
        b = (1 // 3) * (-1 + 10 * α)
        c = 0 // 1

        return SpaceDiscretization(a, b, c, α, β, order)

    elseif order == 6
        α = 2 // 11
        β = 0 // 1

        a = 12 // 11
        b = 3 // 11
        c = 0 // 1

        return SpaceDiscretization(a, b, c, α, β, order)

    elseif order == 8
        α = 344 // 1179
        β = (38 * α - 9) // 214

        a = (696 - 1191 * α) // 428
        b = (2454 * α - 294) // 535
        c = (1179 * α - 344) // 2140

        return SpaceDiscretization(a, b, c, α, β, order)

    elseif order == 10
        α = 334 // 899
        β = 43 // 1798

        a = 1065 // 1798
        b = 1038 // 899
        c = 79 // 1798

        return SpaceDiscretization(a, b, c, α, β, order)
    end
end

const SpaceDiscretizationDefaults::Dict{Symbol,SecondDerivativeCoefficients{Int,Rational{Int}}} = Dict(:ord2 => default_derivative_coefficient_order(2),
                                                                                                       :ord4 => default_derivative_coefficient_order(4),
                                                                                                       :ord6 => default_derivative_coefficient_order(6),
                                                                                                       :ord8 => default_derivative_coefficient_order(8),
                                                                                                       :ord10 => default_derivative_coefficient_order(10))

function register(sym::Symbol, a::T, b::T, c::T, α::T, β::T,
                  order::V) where {T<:AbstractFloatOrRational,V<:Integer}
    Type = Rational{Int}

    SpaceDiscretizationDefaults[sym] = SpaceDiscretization(Type(a), Type(b),
                                                           Type(c), Type(α),
                                                           Type(β),
                                                           Int(order))
end

@inline unregister(::Type{SpaceDiscretization}, sym::Symbol) = delete!(SpaceDiscretizationDefaults,
                                                                       sym)

@inline get_available(::Type{SpaceDiscretization}) = keys(SpaceDiscretizationDefaults)

@inline get_space_discretization(sym::Symbol) = get(SpaceDiscretizationDefaults, sym, nothing)

@inline get_coefficients(::Type{SpaceDiscretization}, sym::Symbol) = get_space_discretization(sym)

export register, unregister, get_available, get_space_discretization, get_coefficients
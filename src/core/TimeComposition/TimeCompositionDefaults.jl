
function time_discretization_ord_2(::Type{T},
                                   substeps::Int64,
                                   index::Int64) where {T<:AbstractFloatOrRational}
    return ConstructSymmetricTimeCompositionMethod(2, 1, convert(Array{T}, [1.0]))
end

function time_discretization_ord_4(::Type{T},
                                   substeps::Int64,
                                   index::Int64) where {T<:AbstractFloatOrRational}
    if substeps == 5 && index == 1
        γ₁ = (3.0 + sqrt(3.0)) / 6
        γ₂ = (3.0 - sqrt(3.0)) / 6
        γ₃ = -1.0
        return ConstructSymmetricTimeCompositionMethod(4, 5,
                                                       convert(Array{T},
                                                               [γ₁, γ₂, γ₃]))
    elseif substeps == 5 && index == 2
        γ₁ = 0.28
        γ₂ = 0.62546642846767004501
        γ₃ = 1 − 2(γ₁ + γ₂)
        return ConstructSymmetricTimeCompositionMethod(4, 5,
                                                       convert(Array{T},
                                                               [γ₁, γ₂, γ₃]))
    elseif substeps == 7 && index == 1
        γ₁ = 1 / (6 - cbrt(6))
        γ₂ = γ₁
        γ₃ = γ₁
        γ₄ = 1 - 6 * γ₁
        return ConstructSymmetricTimeCompositionMethod(4, 7,
                                                       convert(Array{T},
                                                               [γ₁, γ₂, γ₃, γ₄]))
    else
        throw(BoundsError([(5, 1), (5, 2), (7, 1)], (substeps, index)))
    end
end

function time_discretization_ord_6(::Type{T}, substeps::Int64,
                                   index::Int64) where {T<:AbstractFloatOrRational}
    if substeps == 7 && index == 1
        γ₁ = 0.7845136104775572638
        γ₂ = 0.2355732133593581336
        γ₃ = -1.1776799841788710069
        γ₄ = 1 - 2 * (γ₁ + γ₂ + γ₃)
        return ConstructSymmetricTimeCompositionMethod(6, 7,
                                                       convert(Array{T},
                                                               [γ₁, γ₂, γ₃, γ₄]))
    elseif substeps == 9 && index == 1
        γ₁ = 0.1867
        γ₂ = 0.55549702371247839916
        γ₃ =  0.12946694891347535806
        γ₄ = -0.84326562338773460855
        γ₅ = 1 - 2 * (γ₁ + γ₂ + γ₃ + γ₄)
        return ConstructSymmetricTimeCompositionMethod(6, 9,
                                                       convert(Array{T},
                                                               [γ₁, γ₂, γ₃, γ₄, γ₅]))
    elseif substeps == 9 && index == 2
        γ₁ = 0.3921614440073141392
        γ₂ = 0.3325991367893594381
        γ₃ = -0.7062461725576393598
        γ₄ = 0.0822135962935508002
        γ₅ = 0.7985439909348299631
        return ConstructSymmetricTimeCompositionMethod(6, 9,
                                                       convert(Array{T},
                                                               [γ₁, γ₂, γ₃, γ₄, γ₅]))
    else
        throw(BoundsError([(7, 1), (9, 1), (9, 2)], (substeps, index)))
    end
end

function time_discretization_ord_8(::Type{T}, substeps::Int64,
                                   index::Int64) where {T<:AbstractFloatOrRational}
    if substeps == 15 && index == 1
        γ₁ = 0.74167036435061295345
        γ₂ = -0.40910082580003159400
        γ₃ = 0.19075471029623837995
        γ₄ = -0.57386247111608226666
        γ₅ = 0.29906418130365592384
        γ₆ = 0.33462491824529818378
        γ₇ = 0.31529309239676659663
        γ₈ = -0.7968879393529163540
        return ConstructSymmetricTimeCompositionMethod(8,
                                                       15,
                                                       convert(Array{T},
                                                               [γ₁, γ₂, γ₃, γ₄, γ₅,
                                                                γ₆, γ₇, γ₈]))
    else
        throw(BoundsError([(15, 1)], (substeps, index)))
    end
end

function default_time_discretization_coefficients(::Type{T}, order::V, substeps::V,
                                                  index::V) where {V<:Integer,
                                                                   T<:AbstractFloatOrRational}
    if order == 2
        return time_discretization_ord_2(T, substeps, index)
    elseif order == 4
        return time_discretization_ord_4(T, substeps, index)
    elseif order == 6
        return time_discretization_ord_6(T, substeps, index)
    elseif order == 8
        return time_discretization_ord_8(T, substeps, index)
    else
        throw(BoundsError([2, 4, 6, 8], order))
    end
end

const TimeCompositionMethodDefaults::Dict{Symbol,SymmetricTimeCompositionMethod{Int,Rational{Int},Vector{Rational{Int}}}} = Dict(:tord2_1_1 => default_time_discretization_coefficients(Rational{Int},
                                                                                                                                                                                        2,
                                                                                                                                                                                        1,
                                                                                                                                                                                        1),
                                                                                                                                 :tord4_5_1 => default_time_discretization_coefficients(Rational{Int},
                                                                                                                                                                                        4,
                                                                                                                                                                                        5,
                                                                                                                                                                                        1),
                                                                                                                                 :tord4_5_2 => default_time_discretization_coefficients(Rational{Int},
                                                                                                                                                                                        4,
                                                                                                                                                                                        5,
                                                                                                                                                                                        2),
                                                                                                                                 :tord4_7_1 => default_time_discretization_coefficients(Rational{Int},
                                                                                                                                                                                        4,
                                                                                                                                                                                        7,
                                                                                                                                                                                        1),
                                                                                                                                 :tord6_7_1 => default_time_discretization_coefficients(Rational{Int},
                                                                                                                                                                                        6,
                                                                                                                                                                                        7,
                                                                                                                                                                                        1),
                                                                                                                                 :tord6_9_1 => default_time_discretization_coefficients(Rational{Int},
                                                                                                                                                                                        6,
                                                                                                                                                                                        9,
                                                                                                                                                                                        1),
                                                                                                                                 :tord6_9_2 => default_time_discretization_coefficients(Rational{Int},
                                                                                                                                                                                        6,
                                                                                                                                                                                        9,
                                                                                                                                                                                        2),
                                                                                                                                 :tord8_15_1 => default_time_discretization_coefficients(Rational{Int},
                                                                                                                                                                                         8,
                                                                                                                                                                                         15,
                                                                                                                                                                                         1))

function register(sym::Symbol, order::V, substeps::V,
                  coefficients::Array, atol::T = T(0),
                  rtol::T = atol == T(0) ? atol : T(√eps(T))) where {V<:Integer,
                                                                     T<:AbstractFloatOrRational,
                                                                     Array<:AbstractArray{T,
                                                                                          1}}
    Vtype = Int
    Rtype = Rational{Int}

    TimeCompositionMethodDefaults[sym] = ConstructSymmetricTimeCompositionMethod(Vtype(order),
                                                                                 Vtype(substeps),
                                                                                 Vector{Rtype}(coefficients))
end

@inline unregister(::Type{TimeComposition}, sym::Symbol) = delete!(TimeCompositionMethodDefaults,
                                                                   sym)

@inline get_time_composition(sym::Symbol) = TimeCompositionMethodDefaults[sym]

@inline get_available(::Type{TimeComposition}) = keys(TimeCompositionMethodDefaults)

@inline get_coefficients(::Type{TimeComposition}, sym::Symbol) = get_time_composition(sym)

export register, unregister, get_available, get_coefficients, get_time_composition
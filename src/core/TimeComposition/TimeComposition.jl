@inline function check_coefficients(CompositionMethod::SymmetricTimeCompositionMethod{V,
                                                                                      T,
                                                                                      C},
                                    atol::T = T(0),
                                    rtol::T = atol == T(0) ? atol : T(√eps(T))) where {V<:Integer,
                                                                                       T<:AbstractFloatOrRational{V},

                                                                                       C<:AbstractArray{T,
                                                                                                        1}}
    sum_coefficients = 2 * sum(CompositionMethod.coefficients) -
                       CompositionMethod.coefficients[end]
    if !isapprox(sum_coefficients, 1; atol = atol, rtol = rtol)
        throw(DomainError(sum_coefficients,
                          "Must be nearest to one the sum of substeps if you want a valid Symmetric Time Composition Method"))
    end
end

function ConstructSymmetricTimeCompositionMethod(order::V, substeps::V,
                                                 coefficients::Array, atol::T = T(0),
                                                 rtol::T = atol == T(0) ? atol :
                                                           T(√eps(T))) where {V<:Integer,
                                                                              T<:AbstractFloatOrRational{V},
                                                                              Array<:AbstractArray{T,
                                                                                                   1}}
    R = SymmetricTimeCompositionMethod(order, substeps, coefficients)
    check_coefficients(R, atol, rtol)
    return R
end

@inline Base.size(CompositionMethod::SymmetricTimeCompositionMethod{V,T,
C}) where {V<:Integer,T<:AbstractFloatOrRational{V},C<:AbstractArray{T,1}} = (CompositionMethod.substeps,)

@inline function Base.getindex(CompositionMethod::SymmetricTimeCompositionMethod{V,T,
                                                                                 C},
                               index::V) where {V<:Integer,
                                                T<:AbstractFloatOrRational{V},
                                                C<:AbstractArray{T,1}}
    @boundscheck begin
        if !(1 <= index <= length(CompositionMethod))
            throw(BoundsError(1:length(CompositionMethod), index))
        end
    end

    return CompositionMethod.coefficients[mod1(index,
                                               length(CompositionMethod.coefficients))]
end

include("TimeCompositionDefaults.jl")
@inline get_coefficient(CompositionMethod::SymmetricTimeCompositionMethod{T,V,C}, index::V) where {T<:AbstractFloat,V<:Integer,C<:AbstractArray{T,1}} = CompositionMethod.coefficients[mod(index - 1, length(CompositionMethod.coefficients)) + 1]

@inline function check_coefficients(CompositionMethod::SymmetricTimeCompositionMethod{T,V,
                                                                                      C}) where {T<:AbstractFloat,
                                                                                                 V<:Integer,
                                                                                                 C<:AbstractArray{T,
                                                                                                                  1}}
    sum_coefficients = 2 * sum(CompositionMethod.coefficients) -
                       CompositionMethod.coefficients[end]
    if !isapprox(sum_coefficients, 1)
        throw(DomainError(sum_coefficients,
                          "Must be nearest to one the sum of substeps if you want a valid Symmetric Time Composition Method"))
    end
end

function ConstructSymmetricTimeCompositionMethod(order::V, substeps::V,
                                                 coefficients::Array) where {T<:AbstractFloat,
                                                                             V<:Integer,
                                                                             Array<:AbstractArray{T,
                                                                                                  1}}
    R = SymmetricTimeCompositionMethod(order, substeps, coefficients)
    check_coefficients(R)
    return R
end

@inline Base.size(CompositionMethod::SymmetricTimeCompositionMethod{T,V,
C}) where {T<:AbstractFloat,V<:Integer,C<:AbstractArray{T,1}} = (CompositionMethod.substeps,)

@inline function Base.getindex(CompositionMethod::SymmetricTimeCompositionMethod{T,V,
                                                                                 C},
                               index::V) where {T<:AbstractFloat,V<:Integer,
                                                C<:AbstractArray{T,1}}
    @boundscheck begin
        if !(1 <= index <= length(CompositionMethod))
            throw(BoundsError(1:length(CompositionMethod), index))
        end
    end
    
    return get_coefficient(CompositionMethod, index)
end

include("TimeDiscretizationDefaults.jl")
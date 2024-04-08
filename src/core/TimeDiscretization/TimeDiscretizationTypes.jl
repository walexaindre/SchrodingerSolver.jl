struct SymmetricTimeCompositionMethod{V<:Integer,T<:AbstractFloatOrRational{V},C<:AbstractArray{T,1}} <: AbstractArray{T,1}
    order::V
    substeps::V
    coefficients::C
end
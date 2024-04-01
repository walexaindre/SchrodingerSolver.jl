struct SymmetricTimeCompositionMethod{T<:AbstractFloat,V<:Integer,C<:AbstractArray{T,1}} <: AbstractArray{T,1}
    order::V
    substeps::V
    coefficients::C
end
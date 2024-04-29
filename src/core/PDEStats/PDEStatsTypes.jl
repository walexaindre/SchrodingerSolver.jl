abstract type AbstractStats{IntType,FloatType,ArrayType} end 

mutable struct RuntimeStats{IntType,FloatType,ArrayType<:AbstractArray{FloatType}} <: AbstractStats{IntType,FloatType,ArrayType} 
    const system_energy::ArrayType
    const system_power::ArrayType
    const step_time::ArrayType
    const first::IntType
    const last::IntType
    const step::IntType

    current_iter::IntType
end
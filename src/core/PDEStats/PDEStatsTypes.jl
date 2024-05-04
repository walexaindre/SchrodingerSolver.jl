abstract type AbstractStats{IntType,FloatType,ArrayType} end 

struct ComponentPower{FloatType,ArrayType<:AbstractArray{FloatType}} <: AbstractVector{FloatType}
    power::ArrayType
end

mutable struct RuntimeStats{IntType,FloatType,ArrayType<:AbstractArray{FloatType},Power<:ComponentPower{FloatType,ArrayType}} <: AbstractStats{IntType,FloatType,ArrayType} 
    const system_energy::ArrayType
    const system_power::Power
    const step_time::ArrayType
    const first::IntType
    const last::IntType
    const step::IntType

    current_iter::IntType
end
abstract type AbstractStats{IntType,FloatType,ArrayType} end

mutable struct RuntimeStats{IntType,FloatType,ArrayType} <: AbstractStats{ArrayType}
    const system_energy::ArrayType{FloatType,1}
    const system_power::ArrayType{FloatType,1}
    const step_time::ArrayType{FloatType,1}
    current_iter::IntType
end
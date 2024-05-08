abstract type AbstractStats{IntType,FloatType,ArrayType} end 

struct ComponentPower{FloatType,ArrayType<:AbstractArray{FloatType}} <: AbstractVector{FloatType}
    power::ArrayType
end

mutable struct RuntimeStats{IntType,FloatType,ArrayType<:AbstractArray{FloatType},Power<:Tuple{Vararg{ComponentPower{FloatType,ArrayType}}}} <: AbstractStats{IntType,FloatType,ArrayType} 
    const system_energy::ArrayType
    const system_power::Power
    const step_time::ArrayType
    const solver_time::ArrayType
    const solver_iterations::ArrayType
    const log_frequency::IntType
    const τ::FloatType
    const locked::Bool
    log_data::Bool
    current_iteration::IntType
    store_index::IntType
end


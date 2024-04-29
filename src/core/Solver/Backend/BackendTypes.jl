abstract type AbstractBackend end

struct CPUBackend <: AbstractBackend end

struct CPUParallelBackend <: AbstractBackend end

abstract type GPUBackend <: AbstractBackend end

struct CUDABackend <: GPUBackend end

struct AMDBackend <: GPUBackend end

abstract type AbstractBackendMemory end

struct BackendMemory{Tv,Cv<:Complex{Tv},VectorType<:AbstractVector{Cv},ArrayType<:AbstractArray{Cv,2},SolverStorage,LinOp1,LinOp2,LinOp3} <: AbstractBackendMemory
    current_state::ArrayType
    component_temp::VectorType

    opA::LinOp1
    preA::LinOp2
    opD::LinOp3

    b_temp::VectorType
    b0_temp::VectorType

    solver_memory::SolverStorage
end


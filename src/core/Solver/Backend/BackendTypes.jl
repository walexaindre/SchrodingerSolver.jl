abstract type AbstractBackend end

struct CPUBackend <: AbstractBackend end

struct CPUParallelBackend <: AbstractBackend end

struct GPUBackend <: AbstractBackend end

struct CUDABackend <: GPUBackend end

struct AMDBackend <: GPUBackend end

abstract type AbstractBackendMemory end

struct BackendMemory{Tv,ArrayType,SolverStorage} <: AbstractBackendMemory
    current_state::ArrayType{Complex{Tv},2}
    component_temp::ArrayType{Complex{Tv},1}

    b_temp::ArrayType{Complex{Tv},1}
    b0_temp::ArrayType{Complex{Tv},1}

    solver_memory::SolverStorage
end


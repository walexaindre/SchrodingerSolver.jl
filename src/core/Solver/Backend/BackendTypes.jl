abstract type AbstractBackend{IntType,FloatType,ComplexType,VectorType,
                              VectorComplexType,MatrixType,MatrixComplexType} end

abstract type AbstractCPUBackend{IntType,FloatType,ComplexType,VectorType,
                                 VectorComplexType,MatrixType,MatrixComplexType} <:
              AbstractBackend{IntType,FloatType,ComplexType,VectorType,
                              VectorComplexType,
                              MatrixType,MatrixComplexType} end

abstract type AbstractGPUBackend{IntType,FloatType,ComplexType,VectorType,
                                 VectorComplexType,
                                 MatrixType,MatrixComplexType} <:
              AbstractBackend{IntType,FloatType,ComplexType,VectorType,
                              VectorComplexType,MatrixType,MatrixComplexType} end

struct CPUBackend{IntType,FloatType,ComplexType,VectorType,VectorComplexType,
                  MatrixType,MatrixComplexType} <:
       AbstractCPUBackend{IntType,FloatType,ComplexType,VectorType,VectorComplexType,
                          MatrixType,MatrixComplexType} end

struct CPUParallelBackend{IntType,FloatType,ComplexType,VectorType,VectorComplexType,
                          MatrixType,MatrixComplexType} <:
       AbstractCPUBackend{IntType,FloatType,ComplexType,VectorType,VectorComplexType,
                          MatrixType,MatrixComplexType} end

struct CUDABackend{IntType,FloatType,ComplexType,VectorType,VectorComplexType,
                   MatrixType,MatrixComplexType} <:
       AbstractGPUBackend{IntType,FloatType,ComplexType,VectorType,VectorComplexType,
                          MatrixType,MatrixComplexType}
end

struct AMDBackend{IntType,FloatType,ComplexType,VectorType,VectorComplexType,
                  MatrixType,MatrixComplexType} <:
       AbstractGPUBackend{IntType,FloatType,ComplexType,VectorType,VectorComplexType,
                          MatrixType,MatrixComplexType} end

################################################################################

const CPUBackendF32 = CPUBackend{Int32,Float32,ComplexF32,Vector{Float32},
                           Vector{ComplexF32},Array{Float32,2},Array{ComplexF32,2}}

const CPUBackendF64 = CPUBackend{Int64,Float64,ComplexF64,Vector{Float64},
                           Vector{ComplexF64},Array{Float64,2},Array{ComplexF64,2}}

const CPUParallelBackendF32 = CPUParallelBackend{Int32,Float32,ComplexF32,Vector{Float32},
                                           Vector{ComplexF32},Array{Float32,2},
                                           Array{ComplexF32,2}}

const CPUParallelBackendF64 = CPUParallelBackend{Int64,Float64,ComplexF64,Vector{Float64},
                                           Vector{ComplexF64},Array{Float64,2},
                                           Array{ComplexF64,2}}

const CUDABackendF32 = CUDABackend{Int32,Float32,ComplexF32,CuVector{Float32},
                             CuVector{ComplexF32},CuArray{Float32,2},
                             CuArray{ComplexF32,2}}

const CUDABackendF64 = CUDABackend{Int64,Float64,ComplexF64,CuVector{Float64},
                             CuVector{ComplexF64},CuArray{Float64,2},
                             CuArray{ComplexF64,2}}

################################################################################














abstract type AbstractBackendMemory end

struct BackendMemory{Tv,Cv<:Complex{Tv},VectorType<:AbstractVector{Cv},
                     ArrayType<:AbstractArray{Cv,2},SolverStorage,LinOp1,LinOp2,
                     LinOp3} <: AbstractBackendMemory
    current_state::ArrayType
    component_temp::VectorType

    opA::LinOp1
    preA::LinOp2
    opD::LinOp3

    b_temp::VectorType
    b0_temp::VectorType

    solver_memory::SolverStorage
end



export CPUBackend, CPUParallelBackend, CUDABackend, AMDBackend, CPUBackendF32,
       CPUBackendF64, CPUParallelBackendF32, CPUParallelBackendF64, CUDABackendF32,
       CUDABackendF64, AbstractBackend, AbstractCPUBackend, AbstractGPUBackend
is_cpu_backend(::Type{Backend}) where {Backend<:AbstractCPUBackend} = true
is_cpu_backend(::Type{Backend}) where {Backend<:AbstractGPUBackend} = false
is_gpu_backend(::Type{Backend}) where {Backend<:AbstractGPUBackend} = true
is_gpu_backend(::Type{Backend}) where {Backend<:AbstractCPUBackend} = false

function backend_real_vector_type(::Type{Backend}) where {IntType,FloatType,ComplexType,
                                                     VectorType,VectorComplexType,
                                                     MatrixType,MatrixComplexType,
                                                     Backend<:AbstractBackend{IntType,
                                                                              FloatType,
                                                                              ComplexType,
                                                                              VectorType,
                                                                              VectorComplexType,
                                                                              MatrixType,
                                                                              MatrixComplexType}}
    VectorType
end

function backend_complex_vector_type(::Type{Backend}) where {IntType,FloatType,
                                                        ComplexType,
                                                        VectorType,VectorComplexType,
                                                        MatrixType,MatrixComplexType,
                                                        Backend<:AbstractBackend{IntType,
                                                                                 FloatType,
                                                                                 ComplexType,
                                                                                 VectorType,
                                                                                 VectorComplexType,
                                                                                 MatrixType,
                                                                                 MatrixComplexType}}
    VectorComplexType
end

function backend_int_type(::Type{Backend}) where {IntType,FloatType,ComplexType,
                                             VectorType,VectorComplexType,
                                             MatrixType,MatrixComplexType,
                                             Backend<:AbstractBackend{IntType,
                                                                      FloatType,
                                                                      ComplexType,
                                                                      VectorType,
                                                                      VectorComplexType,
                                                                      MatrixType,
                                                                      MatrixComplexType}}
    IntType
end

function backend_float_type(::Type{Backend}) where {IntType,FloatType,ComplexType,
                                               VectorType,VectorComplexType,
                                               MatrixType,MatrixComplexType,
                                               Backend<:AbstractBackend{IntType,
                                                                        FloatType,
                                                                        ComplexType,
                                                                        VectorType,
                                                                        VectorComplexType,
                                                                        MatrixType,
                                                                        MatrixComplexType}}
    FloatType
end

function backend_complex_type(::Type{Backend}) where {IntType,FloatType,ComplexType,
                                                 VectorType,VectorComplexType,
                                                 MatrixType,MatrixComplexType,
                                                 Backend<:AbstractBackend{IntType,
                                                                          FloatType,
                                                                          ComplexType,
                                                                          VectorType,
                                                                          VectorComplexType,
                                                                          MatrixType,
                                                                          MatrixComplexType}}
    ComplexType
end

function backend_matrix_type(::Type{Backend}) where {IntType,FloatType,ComplexType,
                                                VectorType,VectorComplexType,
                                                MatrixType,MatrixComplexType,
                                                Backend<:AbstractBackend{IntType,
                                                                         FloatType,
                                                                         ComplexType,
                                                                         VectorType,
                                                                         VectorComplexType,
                                                                         MatrixType,
                                                                         MatrixComplexType}}
    MatrixType
end

function backend_complex_matrix_type(::Type{Backend}) where {IntType,FloatType,
                                                        ComplexType,
                                                        VectorType,VectorComplexType,
                                                        MatrixType,MatrixComplexType,
                                                        Backend<:AbstractBackend{IntType,
                                                                                 FloatType,
                                                                                 ComplexType,
                                                                                 VectorType,
                                                                                 VectorComplexType,
                                                                                 MatrixType,
                                                                                 MatrixComplexType}}
    MatrixComplexType
end

function initialize_krylov_memory(::Type{Tv}, ::Type{ArrayType},
                                  ::Type{SolverStorage},
                                  opA,
                                  preA,
                                  opD,
                                  ncomponents,
                                  element_count,
                                  memory_size) where {Tv,ArrayType,SolverStorage}
    state = vzeros(ArrayType{Complex{Tv},2}, (element_count, ncomponents))
    mem = vzeros(ArrayType{Complex{Tv},1}, element_count)

    BackendMemory(state, mem, opA, preA, opD, copy(mem), copy(mem),
                  SolverStorage(element_count, element_count, memory_size,
                                typeof(mem)))
end

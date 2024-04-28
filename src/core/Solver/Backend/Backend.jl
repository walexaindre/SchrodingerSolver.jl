function initialize_krylov_memory(::Type{Tv}, ::Type{ArrayType},
                                  ::Type{SolverStorage},
                                  ncomponents,
                                  element_count,
                                  memory_size) where {Tv,ArrayType,SolverStorage}
    state = vzeros(ArrayType{Complex{Tv},2}, (element_count, ncomponents))
    mem = vzeros(ArrayType{Complex{Tv},1}, element_count)

    BackendMemory(state, mem, copy(mem), copy(mem),
                  SolverStorage(element_count, element_count, memory_size,
                                typeof(mem)))
end


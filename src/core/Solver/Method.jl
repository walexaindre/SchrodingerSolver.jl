IterativeLinearSolver(::Type{ComputeBackend}) where {IntType,FloatType,ComplexType,VectorType,
VectorComplexType,MatrixType,
MatrixComplexType,
ComputeBackend<:AbstractBackend{IntType,
FloatType,
ComplexType,
VectorType,
VectorComplexType,
MatrixType,
MatrixComplexType}} = IterativeLinearSolver(700 * eps(FloatType),
                                            700 * eps(FloatType), IntType(10000))

NormBased(::Type{ComputeBackend}) where {IntType,FloatType,ComplexType,VectorType,
VectorComplexType,MatrixType,
MatrixComplexType,
ComputeBackend<:AbstractBackend{IntType,
FloatType,
ComplexType,
VectorType,
VectorComplexType,
MatrixType,
MatrixComplexType}} = NormBased(700 * eps(FloatType), 700 * eps(FloatType),
                                IntType(10000))

function PaulMethod1(::Type{ComputeBackend}, PDE::PDEq, Grid, Conf,
                     stopping_criteria = NormBased(ComputeBackend),
                     preconditioner_drop_tol = FloatType(1e-4),
                     preconditioner_solve_params = IterativeLinearSolver(ComputeBackend),
                     linear_solve_params = IterativeLinearSolver(ComputeBackend)) where {N,
                                                                                         IntType,
                                                                                         FloatType,
                                                                                         ComplexType,
                                                                                         VectorType,
                                                                                         VectorComplexType,
                                                                                         MatrixType,
                                                                                         MatrixComplexType,
                                                                                         ComputeBackend<:CPUBackend{IntType,
                                                                                                                    FloatType,
                                                                                                                    ComplexType,
                                                                                                                    VectorType,
                                                                                                                    VectorComplexType,
                                                                                                                    MatrixType,
                                                                                                                    MatrixComplexType},
                                                                                         PDEq<:SchrodingerPDE{N}}
    #Initialize basic structures
    Mesh = PeriodicAbstractMesh(Grid)

    AI, AJ, AV = get_A_format_COO(FloatType, Mesh, Conf.space_order)
    DI, DJ, DV = get_D_format_COO(FloatType, Grid, Conf.space_order)

    opA = sparse(AI, AJ, AV|>VectorComplexType)

    opD = sparse(DI, DJ, DV|>VectorComplexType)

    preAI, preAJ, preAV = drop(opA, Mesh, preconditioner_drop_tol)

    preA = sparse(preAI, preAJ, preAV)

    #Initialize time related structures
    time_comp = get_time_composition(Conf.time_order)

    σset = Set(get_σ(PDE))
    TimeMultipliers = Grid.τ * coefficients(time_comp)

    time_substeps = Tuple(Grid.τ * collect(time_comp))

    dictionary_keys = Array{Tuple{FloatType,FloatType},1}(undef, 0)
    
    sizehint!(dictionary_keys, length(σset) * length(TimeMultipliers))

    dictionary_values = Array{Kernel{SparseMatrixCSC{ComplexType,IntType},
                                     SparseMatrixCSC{ComplexType,IntType},
                                     SparseMatrixCSC{ComplexType,IntType}},1}(undef,
                                                                              0)
    sizehint!(dictionary_values, length(σset) * length(TimeMultipliers))

    for (σ, βτ) in product(σset,
                           TimeMultipliers)
        opB = (4im * opA + βτ * σ * opD)
        opC = (4im * opA - βτ * σ * opD)
        preBI, preBJ, preBV = drop(opB, Mesh, preconditioner_drop_tol)
        preB = sparse(preBI, preBJ, preBV)

        Ker = Kernel(opB, preB, opC)

        push!(dictionary_keys, (σ, βτ))
        push!(dictionary_values, Ker)
    end

    KernelDict = Dictionary(dictionary_keys, dictionary_values)
    #Initialize memory and Method
    Memory = initialize_krylov_memory(ComputeBackend, PDE, Mesh, opA, preA, opD)
    Method = PaulMethod(Grid, KernelDict, linear_solve_params, stopping_criteria,
                        preconditioner_drop_tol,
                        preconditioner_solve_params,time_substeps)

    Method, Memory
end

function Base.show(io::IO,
                   s::PaulMethod{FloatType,Grid,TKernel,ItSolver,StoppingCriteria,
                                 ItSolver2}) where {FloatType,Grid,TKernel,ItSolver,
                                                    StoppingCriteria,ItSolver2}
    println(io,
            "$(typeof(s))")
end

export PaulMethod1, IterativeLinearSolver, NormBased
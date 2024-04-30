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

    AI, AJ, AV = get_A_format_IJV(FloatType, Mesh, Conf.space_order)
    DI, DJ, DV = get_D_format_IJV(FloatType, Grid, Conf.space_order)

    opA = sparse(AI, AJ, AV)

    opD = sparse(DI, DJ, DV)

    preAI, preAJ, preAV = drop(opA, Mesh, preconditioner_drop_tol)

    preA = sparse(preAI, preAJ, preAV)

    #Initialize time related structures
    time_comp = get_time_composition(Conf.time_order)

    σset = Set(get_σ(PDE))
    TimeMultipliers = Grid.τ * coefficients(time_comp)

    dictionary_keys = Array{Tuple{FloatType,FloatType},1}(undef, 0)
    sizehint!(dictionary_keys, length(σset) * length(TimeMultipliers))

    dictionary_values = Array{Kernel{SparseMatrixCSC{ComplexType,IntType},
                                     SparseMatrixCSC{ComplexType,IntType},
                                     SparseMatrixCSC{ComplexType,IntType}},1}(undef,
                                                                              0)
    sizehint!(dictionary_values, length(σset) * length(TimeMultipliers))

    for (σ, βτ) in product(σset,
                           TimeMultipliers)
        opB = (4im * opA + βτ * σ * D)
        opC = (4im * opA - βτ * σ * D)
        preB = drop(opB, Mesh, preconditioner_drop_tol)

        Ker = Kernel(opB, preB, opC)

        push!(dictionary_keys, (σ, βτ))
        push!(dictionary_values, Ker)
    end

    KernelDict = Dictionary(dictionary_keys, dictionary_values)
    #Initialize memory and Method
    Memory = initialize_krylov_memory(ComputeBackend, PDE, Mesh, opA, preA, opD)
    Method = PaulMethod(Mesh, KernelDict, linear_solve_params, stopping_criteria,
                        preconditioner_drop_tol,
                        preconditioner_solve_params)

    Method, Memory
end

export PaulMethod1
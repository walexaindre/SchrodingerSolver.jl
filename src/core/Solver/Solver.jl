
function DefaultSolver(::Type{ComputeBackend},
                       dim::IntType) where {IntType,FloatType,ComplexType,VectorType,
                                            VectorComplexType,MatrixType,
                                            MatrixComplexType,
                                            ComputeBackend<:AbstractBackend{IntType,
                                                                            FloatType,
                                                                            ComplexType,
                                                                            VectorType,
                                                                            VectorComplexType,
                                                                            MatrixType,
                                                                            MatrixComplexType}}
    SolverConfig(:tord2_1_1,
                 ntuple(Returns(:ord2), dim),
                 ComputeBackend,
                 zero(IntType),
                 false,
                 "/data",
                 (true, IntType(5)),
                 (false, zero(IntType), zero(IntType)))
end

"""

    With this function you can advance one step in time.

"""
function step! end

function update_component!(SP::SchrodingerProblem{PDEq,SolverConf,MethodT,Storage,
                                                  Statistics},
                           component_index::IntType,
                           τ::FloatType,
                           σ::FloatType) where {IntType,FloatType,PDEq,SolverConf,
                                                MethodT<:PaulMethod,
                                                Storage,Statistics}

    #Memory temporary arrays
    Mem = SP.Memory
    Grid = SP.Method.Mesh
    current_state = Mem.current_state
    ψ = view(current_state, :, component_index)
    zₗ = Mem.component_temp
    #dst src
    copy!(zₗ, ψ)

    current_state_abs2 = Mem.current_state_abs2
    temporary_abs2 = Mem.temp_state_abs2

    stage1 = Mem.stage1
    stage2 = Mem.stage2

    b0_temp = Mem.b0_temp
    b_temp = Mem.b_temp

    SolverMem = Mem.solver_memory
    #End of Memory temporary arrays

    PDE = SP.PDE
    sqrt_measure = sqrt(measure(Grid))

    #Method operators
    Method = SP.Method
    stopping_criteria = Method.stopping_criteria #[TODO]
    solver_params = Method.linear_solve_params
    Kernel = Method.Kernel[(σ, τ)]

    opC = Kernel.opC
    opB = Kernel.opB
    preB = Kernel.preB
    opA = Mem.opA
    #End of Method operators

    #N optimized in PDE for easy calculations
    N = get_optimized(PDE)
    #End of N optimized

    b0_temp = opC * ψ

    ndiff = one(FloatType)
    pdiffn = zero(FloatType)

    for l in 1:1000
        @. current_state_abs2 = abs2(current_state)
        @. temporary_abs2 = abs2(zₗ)
        @. stage1 = zₗ + ψ
        stage2 .= N(current_state_abs2, temporary_abs2, component_index)
        @. b_temp = stage1 * stage2
        stage1 .= opA * b_temp
        @. b_temp = -τ * stage1 + b0_temp
        gmres!(SolverMem, opB, b_temp; M = preB,atol=solver_params.atol*1e-1, rtol=solver_params.rtol*1e-1)
        copy!(stage2,zₗ)
        copy!(zₗ,SolverMem.x )
        @. stage1 = stage2 - zₗ

        pdiffn = ndiff
        ndiff = sqrt_measure * norm(stage1)
        #println("Norm: ", norm(stage1)) 
        if ndiff < 1e-10 #|| abs(ndiff - pdiffn) < 1e-10
            #println("Converged at iteration: ", l)
            break
        end
    end
    copy!(ψ, zₗ)
    nothing
end

function step!(SP::SchrodingerProblem{PDEq,SolverConf,Method,Storage,Statistics}) where {PDEq,
                                                                                         SolverConf,
                                                                                         Method,
                                                                                         Storage,
                                                                                         Statistics}
    

    σ_forward = get_σ(SP.PDE)
    σ_backward = reverse(σ_forward)

    time_substeps = SP.Method.time_collection
    
    start_timer = time()
    for τ in time_substeps
        #Forward
        for (component_index, σ) in enumerate(σ_forward)
            update_component!(SP, component_index, τ, σ)
        end
        #Backward

        for (component_index, σ) in zip(length(σ_backward):-1:1, σ_backward)
            update_component!(SP, component_index, τ, σ)
        end
    end

    work_timer = time() - start_timer
    update_stats!(SP.Stats, work_timer, SP.PDE, SP.Method.Mesh, SP.Memory,
                  IterativeLinearSolver(SP.Config.backend_type))
    work_timer
end

export DefaultSolver, step!, update_component!

# module SchrodingerSolver
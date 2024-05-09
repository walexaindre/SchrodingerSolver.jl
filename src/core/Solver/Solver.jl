
function DefaultSolver(::Type{ComputeBackend},
                       dim::IntType, spaceorder = ntuple(Returns(:ord2), dim),
                       timeorder = :tord2_1_1) where {IntType,FloatType,ComplexType,
                                                      VectorType,
                                                      VectorComplexType,MatrixType,
                                                      MatrixComplexType,
                                                      ComputeBackend<:AbstractBackend{IntType,
                                                                                      FloatType,
                                                                                      ComplexType,
                                                                                      VectorType,
                                                                                      VectorComplexType,
                                                                                      MatrixType,
                                                                                      MatrixComplexType}}
    SolverConfig(timeorder,
                 spaceorder,
                 ComputeBackend,
                 zero(IntType),
                 false,
                 "/data",
                 NormBased(ComputeBackend),
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
        gmres!(SolverMem, opB, b_temp; M = preB, atol = solver_params.atol * 1e-1,
               rtol = solver_params.rtol * 1e-1)
        copy!(stage2, zₗ)
        copy!(zₗ, SolverMem.x)
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

function update_component!(PDE, Method, Mem, Stats,
                           style::NormBased{IntType,FloatType},
                           component_index::IntType, τ::FloatType,
                           σ::FloatType) where {IntType,FloatType}
    grid_measure = sqrt(measure(Method.Mesh))
    solved = false

    #Memory temporary arrays
    current_state = Mem.current_state
    ψ = view(current_state, :, component_index)
    zₗ = Mem.component_temp
    #dst src
    copy!(zₗ, ψ)

    current_state_abs2 = Mem.current_state_abs2
    temporary_abs2 = Mem.temp_state_abs2

    stage1 = Mem.stage1 #Is assumed that the norm of stage1 is the norm of the difference
    stage2 = Mem.stage2

    b0_temp = Mem.b0_temp
    b_temp = Mem.b_temp

    SolverMem = Mem.solver_memory
    #End of Memory temporary arrays

    #Method operators
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

    mul!(b0_temp, opC, ψ)
    for l in 1:(style.max_steps)
        @. current_state_abs2 = abs2(current_state)
        @. temporary_abs2 = abs2(zₗ)
        @. stage1 = zₗ + ψ
        stage2 .= N(current_state_abs2, temporary_abs2, component_index)
        @. b_temp = stage1 * stage2
        mul!(stage1, opA, b_temp)

        @. b_temp = -τ * stage1 + b0_temp

        gmres!(SolverMem, opB, b_temp; atol = get_atol(solver_params),
               rtol = get_rtol(solver_params),
               itmax = get_max_iterations(solver_params) )#M = preB, ldiv = false

        update_solver_info!(Stats, SolverMem.stats.timer, SolverMem.stats.niter)

        copy!(stage2, zₗ)
        copy!(zₗ, SolverMem.x)
        @. stage1 = stage2 - zₗ

        znorm = grid_measure * norm(stage1)
        solved = znorm <= style.atol + style.rtol * znorm
        if solved
            break
        end
    end
    copy!(ψ, zₗ)

    if !solved
        @warn "Convergence not reached in $(style.max_steps) iterations..."
    end
    nothing
end

function update_component!(PDE, Method, Mem, Stats, style::FixedSteps{IntType},
                           component_index::IntType, τ::FloatType,
                           σ::FloatType) where {IntType,FloatType}

    #Memory temporary arrays
    current_state = Mem.current_state
    ψ = view(current_state, :, component_index)
    zₗ = Mem.component_temp
    #dst src
    copy!(zₗ, ψ)

    current_state_abs2 = Mem.current_state_abs2
    temporary_abs2 = Mem.temp_state_abs2

    stage1 = Mem.stage1 #Is assumed that the norm of stage1 is the norm of the difference
    stage2 = Mem.stage2

    b0_temp = Mem.b0_temp
    b_temp = Mem.b_temp

    SolverMem = Mem.solver_memory
    #End of Memory temporary arrays

    #Method operators
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

    mul!(b0_temp, opC, ψ)
    for _ in 1:(style.nsteps)
        @. current_state_abs2 = abs2(current_state)
        @. temporary_abs2 = abs2(zₗ)
        @. stage1 = zₗ + ψ
        stage2 .= N(current_state_abs2, temporary_abs2, component_index)
        @. b_temp = stage1 * stage2
        mul!(stage1, opA, b_temp)
        @. b_temp = -τ * stage1 + b0_temp
        gmres!(SolverMem, opB, b_temp; N = preB, atol = get_atol(solver_params),
               rtol = get_rtol(solver_params),
               itmax = get_max_iterations(solver_params))
        copy!(zₗ, SolverMem.x)
        update_solver_info!(Stats, SolverMem.timer, SolverMem.niter)
    end
    copy!(ψ, zₗ)
end

function step!(SP::SchrodingerProblem{PDEq,SolverConf,Meth,Storage,Statistics}) where {PDEq,
                                                                                       SolverConf,
                                                                                       Meth,
                                                                                       Storage,
                                                                                       Statistics}
    PDE = SP.PDE
    Method = SP.Method
    Memory = SP.Memory
    Stats = SP.Stats
    Conf = SP.Config

    Style = Conf.stopping_criteria

    σ_forward = get_σ(SP.PDE)
    σ_backward = reverse(σ_forward)

    time_substeps = SP.Method.time_collection

    start_timer = time()
    for τ in time_substeps
        #Forward
        for (component_index, σ) in enumerate(σ_forward)
            update_component!(PDE, Method, Memory, Stats, Style, component_index, τ,
                              σ)
        end
        #Backward

        for (component_index, σ) in zip(length(σ_backward):-1:1, σ_backward)
            update_component!(PDE, Method, Memory, Stats, Style, component_index, τ,
                              σ)
        end
    end

    work_timer = time() - start_timer
    update_stats!(SP.Stats, work_timer, SP.PDE, SP.Method.Mesh, SP.Memory,
                  IterativeLinearSolver(SP.Config.backend_type))
    work_timer
end

export DefaultSolver, step!, update_component!

# module SchrodingerSolver
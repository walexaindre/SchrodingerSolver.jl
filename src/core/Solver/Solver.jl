function base_name(S::SolverConfig)
end

function DefaultSolver(dim::IntType) where {IntType<:Integer}
    SolverConfig(:tord2_1_1,
                 ntuple(Returns(:ord2), dim),
                 CPUBackend,
                 0,
                 1e-8,
                 1e-8,
                 1e-4,
                 false,
                 "/data",
                 (true, 5),
                 (false, 0, 0))
end

"""

    With this function you can advance one step in time.

"""
function step! end

function update_component!()
end

function step!(SP::SchrodingerProblem{PDEq,SolverConf,Storage,Statistics}) where {PDEq,SolverConf,Storage,Statistics}

    σ_forward = get_σ(SP.PDE)
    σ_backward = reverse(σ_forward)
    start_time = time()
    for τ in 1:3
        #Forward
        for (component_index, σ) in enumerate(σ_forward)
            update_component!()
        end
        #Backward

        for (component_index,σ) in zip(length(σ_backward):-1:1,σ_backward)
            update_component!()
        end
    end
    

    workt_time = time() - start_time
    
end
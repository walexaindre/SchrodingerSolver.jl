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
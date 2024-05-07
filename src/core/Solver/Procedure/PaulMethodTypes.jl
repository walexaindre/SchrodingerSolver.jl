struct PaulMethod{FloatType,Grid,TKernel,ItSolver,StoppingCriteria,ItSolver2,TTime} <: AbstractSolverMethod{FloatType}
    Mesh::Grid
    Kernel::TKernel
    linear_solve_params::ItSolver
    preconditioner_drop_tol::FloatType
    preconditioner_solve_params::ItSolver2
    time_collection::TTime
end

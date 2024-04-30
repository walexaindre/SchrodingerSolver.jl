abstract type AbstractStoppingCriterion end

struct FixedSteps{IntType}<:AbstractStoppingCriterion
    nsteps::IntType
end

struct NormBased{IntType,FloatType}<:AbstractStoppingCriterion
    atol::FloatType
    rtol::FloatType
    max_steps::IntType
end

abstract type AbstractLinearSolver end

struct IterativeLinearSolver{IntType,FloatType}<:AbstractLinearSolver
    rtol::FloatType
    atol::FloatType
    max_iteration::IntType
end

struct DirectLinearSolver
end

abstract type AbstractSolverMethod end

struct Kernel{LinOp1,LinOp2,LinOp3}
    opB::LinOp1
    preB::LinOp2
    opC::LinOp3
end

struct PaulMethod{FloatType,Grid,TKernel,ItSolver,StoppingCriteria,ItSolver2}
    Mesh::Grid
    Kernel::TKernel
    linear_solve_params::ItSolver
    stopping_criteria::StoppingCriteria
    preconditioner_drop_tol::FloatType
    preconditioner_solve_params::ItSolver2
end


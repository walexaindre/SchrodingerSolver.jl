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



struct SolverConfig{N,IntType,FloatType,Backend,Method,LinearSolver,StoppingCriteria}
    time_order::Symbol
    space_order::NTuple{N,Symbol}
    backend_type::Type{Backend}

    method::Method

    verbose::IntType

    linear_solver::LinearSolver

    stop_criteria::StoppingCriteria

    droptol::FloatType

    show_progress::Bool
    
    output_folder::String

    OpA
    PreA

    stats::Tuple{Bool,IntType}
    debug::Tuple{Bool,IntType,IntType}
end

struct SchrodingerProblem{PDEq,SolverConf,Storage,Statistics}
    PDE::PDEq
    Config::SolverConf
    Memory::Storage
    Stats::Statistics
end
struct SolverConfig{N,IntType,Backend,StoppingCriteria}
    time_order::Symbol
    space_order::NTuple{N,Symbol}
    backend_type::Type{Backend}
    verbose::IntType
    show_progress::Bool
    output_folder::String
    stopping_criteria::StoppingCriteria
    stats::Tuple{Bool,IntType}
    debug::Tuple{Bool,IntType,IntType}
end


struct SchrodingerProblem{PDEq,SolverConf,Method,Storage,Statistics}
    PDE::PDEq
    Config::SolverConf
    Method::Method
    Memory::Storage
    Stats::Statistics
end


export SchrodingerProblem, SolverConfig

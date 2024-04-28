struct SolverConfig{N,IntType,FloatType,Backend}
    time_order::Symbol
    space_order::NTuple{N,Symbol}
    backend_type::Type{Backend}
    
    verbose::IntType

    rtol::FloatType
    atol::FloatType
    droptol::FloatType

    show_progress::Bool
    
    output_folder::String

    stats::Tuple{Bool,IntType}
    debug::Tuple{Bool,IntType,IntType}
end
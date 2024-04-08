const TimeDiscretizationDefaults::Dict{Symbol,SymmetricTimeCompositionMethod{Int,Rational{Int},Vector{Rational{Int}}}} = Dict(:tord2 => default_derivative_coefficient_order(2),
                                                                                                       :tord4 => default_derivative_coefficient_order(4),
                                                                                                       :tord6 => default_derivative_coefficient_order(6),
                                                                                                       :tord8 => default_derivative_coefficient_order(8),
                                                                                                       :tord10 => default_derivative_coefficient_order(10))
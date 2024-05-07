IterativeLinearSolver(::Type{ComputeBackend}) where {IntType,FloatType,ComplexType,VectorType,
VectorComplexType,MatrixType,
MatrixComplexType,
ComputeBackend<:AbstractBackend{IntType,
FloatType,
ComplexType,
VectorType,
VectorComplexType,
MatrixType,
MatrixComplexType}} = IterativeLinearSolver(700 * eps(FloatType),
                                            700 * eps(FloatType), IntType(10000))

NormBased(::Type{ComputeBackend}) where {IntType,FloatType,ComplexType,VectorType,
VectorComplexType,MatrixType,
MatrixComplexType,
ComputeBackend<:AbstractBackend{IntType,
FloatType,
ComplexType,
VectorType,
VectorComplexType,
MatrixType,
MatrixComplexType}} = NormBased(700 * eps(FloatType), 700 * eps(FloatType),
                                IntType(10000))


@inline get_atol(param::NormBased{IntType,FloatType}) where {IntType,FloatType}= param.atol
@inline get_rtol(param::NormBased{IntType,FloatType}) where {IntType,FloatType}= param.rtol
@inline get_max_iterations(param::NormBased{IntType,FloatType}) where {IntType,FloatType}= param.max_steps

@inline get_atol(param::IterativeLinearSolver{IntType,FloatType}) where {IntType,FloatType}= param.atol
@inline get_rtol(param::IterativeLinearSolver{IntType,FloatType}) where {IntType,FloatType}= param.rtol
@inline get_max_iterations(param::IterativeLinearSolver{IntType,FloatType}) where {IntType,FloatType}= param.max_iterations

Base.show(io::IO, param::IterativeLinearSolver) = print(io, "IterativeLinearSolver\nrtol: $(param.rtol)\natol: $(param.atol)\nmaxiter: $(param.max_iteration)")
Base.show(io::IO, param::NormBased) = print(io, "StoppingCriteria\nNormBased\natol: $(param.atol)\nrtol: $(param.rtol)\nmaxiter: $(param.max_steps)")
Base.show(io::IO, param::FixedSteps) = print(io, "StoppingCriteria\nFixedSteps\nnsteps: $(param.nsteps)")

export IterativeLinearSolver, NormBased, FixedSteps, DirectLinearSolver
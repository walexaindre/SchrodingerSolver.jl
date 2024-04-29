abstract type AbstractSolverMethod end

struct PaulMethod{ArrayType,PreB}
    B::ArrayType
    preB::PreB
    D::ArrayType
end
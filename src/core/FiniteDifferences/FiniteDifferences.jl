function default_discretization_coefficients(::Type{T},
    order::Int)::FiniteDifferencesSecondDerivativeCoefficients{T} where {T <: AbstractFloat}
    if order == 2
        α = 0 // 1
        β = 0 // 1

        a = 1 // 1
        b = 0 // 1
        c = 0 // 1

    elseif order == 4
        α = 1 // 10
        β = 0 // 1

        a = (4 // 3) * (1 - α)
        b = (1 // 3) * (-1 + 10 * α)
        c = 0 // 1
    elseif order == 6
        α = 2 // 11
        β = 0 // 1

        a = 12 // 11
        b = 3 // 11
        c = 0 // 1
    elseif order == 8
        α = 344 // 1179
        β = (38 * α - 9) // 214

        a = (696 - 1191 * α) // 428
        b = (2454 * α - 294) // 535
        c = (1179 * α - 344) // 2140

    elseif order == 10
        α = 334 // 899
        β = 43 // 1798

        a = 1065 // 1798
        b = 1038 // 899
        c = 79 // 1798
    else
        throw(BoundsError([2, 4, 6, 8, 10], order))
    end

    FiniteDifferencesSecondDerivativeCoefficients{T}(T(a), T(b), T(c), T(α), T(β))
end

function get_A_format_IJV(::Type{T},
    Mesh::AbstractMesh{T,1},
    space_discretization_coefficients::FiniteDifferencesSecondDerivativeCoefficients{T}) where {T <: AbstractFloat}

    α = space_discretization.α
    β = space_discretization.β
    
    count = 0
    value = [one(T)]
    start_depth = 0
    end_depth = 0

    if (α != 0)
        count += 2
        start_depth = 1
        end_depth = 1
        append!(value,fill(α, 2))
    end

    if (β != 0)
        count += 2
        end_depth = 2
        append!(value,fill(β, 2))
    end

    if (α == 0 && β != 0)
        start_depth = 2
    end

    count += 1

    space_usage = count * length(Mesh)

    I = zeros(Int64, space_usage) #row idx
    J = zeros(Int64, space_usage) #column idx
    V = repeat(value, length(Mesh))

    @threads for idx in 1:length(Mesh)
        I[(count * (idx - 1) + 1):(count * idx)] .= idx
        if (start_depth != 0)
            J[(count * (idx - 1) + 1):(count * idx)] .= get_linear_stencil(idx,
                start_depth,
                end_depth,
                Mesh)
        else
            J[(count * (idx - 1) + 1):(count * idx)] .= idx
        end

        #V[(count * (idx - 1) + 1):(count * idx)] .= value
    end
    I, J, V
end
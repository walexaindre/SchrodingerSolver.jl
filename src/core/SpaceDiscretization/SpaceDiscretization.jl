
function estimate_order(a::T, b::T, c::T, α::T, β::T, atol::T = 0,
                        rtol::T = atol == 0 ? atol : √eps(T)) where {T<:AbstractFloatOrRational{Int}}
    if !isapprox(a + b + c, 1 + 2 * (α + β); atol = atol, rtol = rtol)
        return 0
    end

    for idx in 2:2:36
        pow2 = 2^idx
        pow3 = 3^idx

        fact = (idx + 1) * (idx + 2)

        left_side = a + pow2 * b + pow3 * c
        right_side = fact * (α + pow2 * β)

        if !isapprox(left_side, right_side; atol = atol, rtol = rtol)
            return idx
        end
    end

    throw(ArgumentError("Order higher than 36... "))
end

function validate_constraints(a::T, b::T, c::T, α::T, β::T, order::V, atol::T = 0,
                              rtol::T = atol == 0 ? atol : √eps(T)) where {V<:Integer,
                                                                           T<:AbstractFloatOrRational{V}}
    if order < 2
        throw(ArgumentError("Order must be greater than 1 and even at space discretization => [ order = $order ]"))
    end

    if mod(order, 2) == 1
        throw(ArgumentError("Order must be even at space discretization => [ order = $order ]"))
    end

    if !isapprox(a + b + c, 1 + 2 * (α + β); atol = atol, rtol = rtol)
        throw(ArgumentError("Constraint not satisfied => [ a + b + c - 1 - 2 * ( α + β ) = $(a + b + c - 1 - 2 * α - 2 * β) ] != 0"))
    end

    for idx in 2:2:(order - 2)
        pow2 = 2^idx
        pow3 = 3^idx

        fact = (idx + 1) * (idx + 2)

        left_side = a + pow2 * b + pow3 * c
        right_side = fact * (α + pow2 * β)

        if !isapprox(left_side, right_side; atol = atol, rtol = rtol)
            throw(ArgumentError("Constraint not satisfied => [ a + 2 ^ $idx * b + 3 ^ $idx * c -  ($(idx+2)!/$idx!)( α - 2 ^ $idx * β ) = $(abs(left_side-right_side)) ] !≈ 0"))
        end
    end
end

function validate_constraints(SpaceDiscretization::SecondDerivativeCoefficients{V,T},
                              atol::T = 0,
                              rtol::T = atol == 0 ? atol : √eps(T)) where {V<:Integer,
                                                                           T<:AbstractFloatOrRational{V}}
    return validate_constraints(SpaceDiscretization.a, SpaceDiscretization.b,
                                SpaceDiscretization.c, SpaceDiscretization.α,
                                SpaceDiscretization.β, SpaceDiscretization.order,
                                atol,
                                rtol)
end

function validate_positive_definite(α::T,
                                    β::T) where {T<:AbstractFloatOrRational{Int}}

    #By Gershgorin's theorem the generated matrix A is PSD
    #If |λ-1|≤ 2 ( α + β ) => -2 ( α + β ) ≤ λ - 1 ≤ 2 ( α + β ) => 1 - 2 ( α - β ) ≤ λ ≤ 1 + 2 ( α + β )
    # Then 1 - 2 ( α + β ) > 0 is required...
    if 1 - 2 * (α + β) <= 0
        throw(ArgumentError("Generated A matrix will not be positive definite..."))
    end
end

function validate_positive_definite(SpaceDiscretization::SecondDerivativeCoefficients{V,
                                                                                      T}) where {V<:Integer,
                                                                                                 T<:AbstractFloatOrRational{V}}
    validate_positive_definite(SpaceDiscretization.α, SpaceDiscretization.β)
end

function check_validity(SpaceDiscretization::SecondDerivativeCoefficients{V,T},
                        atol::T = 0,
                        rtol::T = atol == 0 ? atol : √eps(T)) where {T<:AbstractFloat,
                                                                     V<:Integer}
    validate_constraints(SpaceDiscretization, atol, rtol)
    validate_positive_definite(SpaceDiscretization)
end

function check_validity(a::T, b::T, c::T, α::T, β::T, order::V, atol::T = 0,
                        rtol::T = atol == 0 ? atol : √eps(T)) where {V<:Integer,
                                                                     T<:AbstractFloatOrRational{V}}
    validate_constraints(a, b, c, α, β, order, atol, rtol)
    validate_positive_definite(α, β)
end

function SpaceDiscretization(a::T, b::T, c::T, α::T, β::T, order::V, atol::T = T(0),
                             rtol::T = atol == T(0) ? atol : T(√eps(T))) where {V<:Integer,
                                                                                T<:AbstractFloatOrRational{V}}
    check_validity(a, b, c, α, β, order, atol, rtol)
    return SecondDerivativeCoefficients{V,T}(a, b, c, α, β, order)
end

include("SpaceDiscretizationDefaults.jl")

function get_A_format_IJV(::Type{T}, Mesh::AM,
                          Order::Symbol) where {V<:Integer,AM<:AbstractMesh{V,1},
                                                T<:AbstractFloatOrRational{V}}

    #Fetch discretization coefficients
    SD = get(Order, SpaceDiscretizationDefaults, nothing)

    if isnothing(SD)
        throw(ArgumentError("Order not found in SpaceDiscretizationDefaults..."))
    end

    offsets = [V(0)]

    value = [one(T)]

    sizehint!(offsets, 5)
    sizehint!(value, 5)

    if !isapprox(SD.α, 0)
        push!(offsets, V(1), V(-1))
        push!(value, SD.α, -SD.α)
    end

    if !isapprox(SD.β, 0)
        push!(offsets, V(2), V(-2))
        push!(value, SD.β, -SD.β)
    end

    count = size(offsets, 1)
    space_usage = count * length(Mesh)

    _I = zeros(V, space_usage) #row idx
    _J = zeros(V, space_usage) #column idx
    _V = repeat(value, length(Mesh))

    #Get the IJV format of the matrix A
    @threads for idx in 1:length(Mesh)
        _I[(count * (idx - 1) + 1):(count * idx)] .= idx
        _J[(count * (idx - 1) + 1):(count * idx)] .= apply_offsets(Mesh, idx,
                                                                   offsets)

        #V[(count * (idx - 1) + 1):(count * idx)] .= value
    end
    _I, _J, _V
end

function get_A_format_IJV(::Type{T}, Mesh::AM,
                          Order::NTuple{N,Symbol}) where {V<:Integer,
                                                          AM<:AbstractMesh{V,1},
                                                          T<:AbstractFloatOrRational{V},
                                                          N}

    #Get individual dimensions
    submesh = extract_every_dimension(Mesh)

    #Get matrix for every sub mesh
    _A = [sparse(get_A_format_IJV(T, submesh[idx], Order[idx])...,
                 length(submesh[idx])) for idx in N:-1:1]
    A = kron(_A...)
    #Get the IJV format of the matrix A
    findnz(A)
end

function get_Δ_format_IJV(::Type{T}, Grid::AG,
                          Order::NTuple{N,Symbol}) where {N,
                                                          T<:AbstractFloatOrRational,
                                                          AG<:AbstractGrid{NTuple{N,
                                                                                  T},
                                                                           N}}
    Mesh = PeriodicAbstractMesh(Grid)
    V = typeof(Mesh.dims[1])
    submeshes = extract_every_dimension(Mesh)

    temp = []
    for (h, ord, submesh) in zip(Grid.h, Order, submeshes)
        SD = get(SpaceDiscretizationDefaults, ord, nothing)

        if isnothing(SD)
            throw(ArgumentError("Order not found in SpaceDiscretizationDefaults..."))
        end

        values = [zero(T)]
        offsets = [V(0)]

        if !isapprox(SD.a, 0)
            values[1] += (-2 * SD.a) / (h * h)
            push!(values, SD.a / (h * h))
            offsets = [V(1), V(-1)]
        end

        if !isapprox(SD.b, 0)
            values[1] += (-2 * SD.b) / (h * h)
            push!(values, SD.b / (4 * h * h))
            offsets = [V(2), V(-2)]
        end

        if !isapprox(SD.c, 0)
            values[1] += (-2 * SD.c) / (h * h)
            push!(values, SD.c / (9 * h * h))
            offsets = [V(3), V(-3)]
        end

        count = size(offsets, 1)
        space_usage = count * length(submesh)

        _I = zeros(Int64, space_usage) #row idx
        _J = zeros(Int64, space_usage) #column idx
        _V = repeat(values, length(submesh))

        #@threads 
        for idx in 1:length(submesh)
            section = (count * (idx - 1) + 1):(count * idx)
            _I[section] .= idx
            _J[section] .= apply_offsets(submesh, idx, offsets)
        end
        push!(temp, sparse(_I, _J, _V))
    end

    reverse!(temp)

    sizes = [size(temp[idx], 1) for idx in 1:N]

    #1D Dx
    #2D => kron(Iy, Dx) + kron(Dy, Ix)
    #3D => kron(Iz, Iy, Dx) + kron(Iz, Dy, Ix) + kron(Dz, Iy, Ix)

    if N == 1
        return findnz(temp[1])
    elseif N == 2
        Ix = sparse(T(1)I, sizes[2], sizes[2])
        Iy = sparse(T(1)I, sizes[1], sizes[1])
        return findnz(kron(Iy, temp[2]) + kron(temp[1], Ix))

    elseif N == 3
        Ix = sparse(T(1)I, sizes[3], sizes[3])
        Iy = sparse(T(1)I, sizes[2], sizes[2])
        Iz = sparse(T(1)I, sizes[1], sizes[1])
        return findnz(kron(Iz, Iy, temp[3]) + kron(Iz, temp[2], Ix) +
                      kron(temp[1], Iy, Ix))
    end
end

function get_D_format_IJV()
end

export SpaceDiscretization, get_A_format_IJV, get_Δ_format_IJV, check_validity,
       validate_positive_definite, validate_constraints, estimate_order

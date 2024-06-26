
function estimate_order(a::T, b::T, c::T, α::T,
                        β::T) where {T<:AbstractFloatOrRational{Int}}
    if !isapprox(a + b + c, 1 + 2 * (α + β))
        return 0
    end

    for idx in 2:2:36
        pow2 = 2^idx
        pow3 = 3^idx

        fact = (idx + 1) * (idx + 2)

        left_side = a + pow2 * b + pow3 * c
        right_side = fact * (α + pow2 * β)

        if !isapprox(left_side, right_side)
            return idx
        end
    end

    throw(ArgumentError("Order higher than 36... "))
end

function validate_constraints(a::T, b::T, c::T, α::T, β::T,
                              order::V) where {V<:Integer,
                                               T<:AbstractFloatOrRational{V}}
    if order < 2
        throw(ArgumentError("Order must be greater than 1 and even at space discretization => [ order = $order ]"))
    end

    if mod(order, 2) == 1
        throw(ArgumentError("Order must be even at space discretization => [ order = $order ]"))
    end

    if !isapprox(a + b + c, 1 + 2 * (α + β))
        throw(ArgumentError("Constraint not satisfied => [ a + b + c - 1 - 2 * ( α + β ) = $(a + b + c - 1 - 2 * α - 2 * β) ] != 0"))
    end

    for idx in 2:2:(order - 2)
        pow2 = 2^idx
        pow3 = 3^idx

        fact = (idx + 1) * (idx + 2)

        left_side = a + pow2 * b + pow3 * c
        right_side = fact * (α + pow2 * β)

        if !isapprox(left_side, right_side)
            throw(ArgumentError("Constraint not satisfied => [ a + 2 ^ $idx * b + 3 ^ $idx * c -  ($(idx+2)!/$idx!)( α - 2 ^ $idx * β ) = $(abs(left_side-right_side)) ] !≈ 0"))
        end
    end
end

function validate_constraints(SpaceDiscretization::SecondDerivativeCoefficients{V,T}) where {V<:Integer,
                                                                                             T<:AbstractFloatOrRational{V}}
    return validate_constraints(SpaceDiscretization.a, SpaceDiscretization.b,
                                SpaceDiscretization.c, SpaceDiscretization.α,
                                SpaceDiscretization.β, SpaceDiscretization.order)
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

function check_validity(SpaceDiscretization::SecondDerivativeCoefficients{V,T}) where {T<:AbstractFloat,
                                                                                       V<:Integer}
    validate_constraints(SpaceDiscretization)
    validate_positive_definite(SpaceDiscretization)
end

function check_validity(a::T, b::T, c::T, α::T, β::T,
                        order::V) where {V<:Integer,
                                         T<:AbstractFloatOrRational{V}}
    validate_constraints(a, b, c, α, β, order)
    validate_positive_definite(α, β)
end

function SpaceDiscretization(a::T, b::T, c::T, α::T, β::T,
                             order::V) where {V<:Integer,
                                              T<:AbstractFloatOrRational{V}}
    check_validity(a, b, c, α, β, order)
    return SecondDerivativeCoefficients{V,T}(a, b, c, α, β, order)
end

include("SpaceDiscretizationDefaults.jl")

function get_A_format_COO(::Type{T}, Mesh::AM,
                          Order::Symbol) where {V<:Integer,AM<:AbstractMesh{V,1},
                                                T<:AbstractFloatOrRational{V}}

    #Fetch discretization coefficients
    SD = get_coefficients(SpaceDiscretization, Order)

    if isnothing(SD)
        throw(ArgumentError("Order not found in SpaceDiscretizationDefaults..."))
    end

    offsets = zeros(V, 0)

    value = [one(T)]

    sizehint!(offsets, 5)
    sizehint!(value, 5)

    if !isapprox(SD.α, 0)
        push!(offsets, V(1))
        push!(value, SD.α, SD.α)
    end

    if !isapprox(SD.β, 0)
        push!(offsets, V(2))
        push!(value, SD.β, SD.β)
    end

    count = size(value, 1)
    space_usage = count * length(Mesh)

    _I = zeros(V, space_usage) #row idx
    _J = zeros(V, space_usage) #column idx
    _V = repeat(value, length(Mesh))
    
    offsets = GenerateOffset(OffsetUniqueZero, 1, (offsets,))

    #Get the IJV format of the matrix A
    @threads for idx in 1:length(Mesh)

        _I[(count * (idx - 1) + 1):(count * idx)] .= idx
        _J[(count * (idx - 1) + 1):(count * idx)] .= apply_offsets(Mesh, CartesianIndex(idx),
                                                                   offsets)

        #V[(count * (idx - 1) + 1):(count * idx)] .= value
    end
    _I, _J, _V
end

function get_A_format_COO(::Type{T}, Mesh::AM,
                          Order::NTuple{N,Symbol}) where {V<:Integer,
                                                          T<:AbstractFloatOrRational{V},
                                                          N,AM<:AbstractMesh{V,N}}

    #Get individual dimensions
    submesh = extract_every_dimension(Mesh)
    #Get matrix for every sub mesh

    _A = ntuple(ndims(Mesh)) do idx #This can have bugs easily... [TODO]

        pos = N+1-idx

        cmesh,_ = iterate(submesh,pos-1)
        sz = length(cmesh)

        _I,_J,_V = get_A_format_COO(T, cmesh, Order[pos])

        sparse(_I,_J,_V,sz,sz)
    end

    #Kronecker product of all matrices in reverse order...
    A = kron(_A...)
    #Get the IJV format of the matrix A
    findnz(A)
end

function get_D_format_COO(::Type{Tv}, Grid::AG,
                          Order::NTuple{N,Symbol}) where {N,
                                                          Tv<:AbstractFloatOrRational,
                                                          AG<:AbstractGrid{NTuple{N,
                                                                                  Tv},
                                                                           N}}
    if N > 3
        throw(ArgumentError("Only 1D, 2D and 3D grids are supported..."))
    end

    Mesh = PeriodicAbstractMesh(Grid)
    V = typeof(Mesh.dims[1])
    submeshes = extract_every_dimension(Mesh)

    temp = Vector{SparseMatrixCSC{Tv,V}}(undef, 0)
    temp_A = Vector{SparseMatrixCSC{Tv,V}}(undef, 0)

    for (h, ord, submesh) in zip(Grid.h, Order, submeshes)
        SD = get_coefficients(SpaceDiscretization, ord)

        if isnothing(SD)
            throw(ArgumentError("Order not found in SpaceDiscretizationDefaults... Please use register function to add new orders."))
        end

        values = [zero(Tv)]
        offsets = zeros(V, 0)
        hsqr = h^2

        if !isapprox(SD.a, 0)
            nondiaga = SD.a / hsqr
            diaga = -2 * nondiaga

            values[1] += diaga
            push!(values, nondiaga, nondiaga)
            push!(offsets, V(1))
        end

        if !isapprox(SD.b, 0)
            nondiagb = SD.b / (4 * hsqr)
            diagb = -2 * nondiagb

            values[1] += diagb
            push!(values, nondiagb, nondiagb)
            push!(offsets, V(2))
        end

        if !isapprox(SD.c, 0)
            nondiagc = SD.c / (9 * hsqr)
            diagc = -2 * nondiagc

            values[1] += diagc
            push!(values, nondiagc, nondiagc)
            push!(offsets, V(3), V(-3))
        end

        count = size(values, 1)
        space_usage = count * length(submesh)

        sym_offsets = GenerateOffset(OffsetUniqueZero, 1, (offsets,))

        _I = zeros(Int64, space_usage) #row idx
        _J = zeros(Int64, space_usage) #column idx
        _V = repeat(values, length(submesh))

        #@threads 
        for idx in 1:length(submesh)
            section = (count * (idx - 1) + 1):(count * idx)
            _I[section] .= idx
            _J[section] .= apply_offsets(submesh, CartesianIndex(idx), sym_offsets)
        end
        push!(temp, sparse(_I, _J, _V))
        AI,AJ,AV = get_A_format_COO(Tv, submesh, ord)
        push!(temp_A,sparse(AI,AJ,AV))
    end

    reverse!(temp)

    #1D Dx
    #2D => kron(Ay, Dx) + kron(Dy, Ax)
    #3D => kron(Az, Ay, Dx) + kron(Az, Dy, Ax) + kron(Dz, Ay, Ax)

    if N == 1
        return findnz(temp[1])
    elseif N == 2
        Ax = temp_A[2]
        Ay = temp_A[1]
        return findnz(kron(Ay, temp[2]) + kron(temp[1], Ax))

    elseif N == 3
        Ax = temp_A[3]
        Ay = temp_A[2]
        Az = temp_A[1]
        return findnz(kron(Az, Ay, temp[3]) + kron(Az, temp[2], Ax) +
                      kron(temp[1], Ay, Ax))
    end
end

export SpaceDiscretization, get_A_format_COO, get_D_format_COO, check_validity,
       validate_positive_definite, validate_constraints, estimate_order

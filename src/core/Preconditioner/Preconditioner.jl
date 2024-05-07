function drop(M::Matrix, AMesh::AM, τ::Tv, offset_ranges::RTup = 1:2,
              rtol::TolT = 700 * eps(real(eltype(M))),
              atol::TolT = 700 * eps(real(eltype(M)))) where {N,V<:Integer,
                                                              Tv<:AbstractFloat,
                                                              TolT<:AbstractFloat,
                                                              RTup<:TupleOrRange{V},
                                                              Matrix,
                                                              AM<:AbstractMesh{V,N}}
    if size(M, 1) != size(M, 2)
        throw(DimensionMismatch("Matrix must be square"))
    end

    if size(M, 1) != length(AMesh)
        throw(DimensionMismatch("Matrix size must be equal to the mesh size"))
    end

    if ndims(M) != 2
        throw(DimensionMismatch("Matrix must be 2D"))
    end

    #Vector storage type
    VecS = matrix_to_vector(Matrix)

    #Generate first element of the mesh.
    Ind = ntuple(Returns(1), N)

    #Generate a list of offsets for given ranges.
    offsets = GenerateOffset(OffsetUniqueZero,1, GenerateRank(offset_ranges, N))

    #Apply offset list to first element (1,1,...,1) of the mesh.
    one_idx_offset = apply_offsets(AMesh, Ind, offsets)

    #Get the number of rows/cols of the matrix. (We assume here that M is square)
    rows = size(M, 1)

    #Generate a vector with a 1 at the first position and zeros elsewhere.
    e₁ = vcanonical(VecS, rows)

    #Solve the system Mx = e₁
    col, _ = gmres(M, e₁; atol = atol, rtol = rtol)

    #Prune the inverse with the desired drop tolerance τ.
    bit_index = abs2.(col) .> τ

    #Get indexes of non pruned elements.
    nonpruned_index = findall(bit_index)

    #Undropped offset positions.
    undropped_offset_positions = findall(in(nonpruned_index), one_idx_offset)

    #Infer offset values at used offset positions.
    minimal_offset = infer_minimal_offsets(undropped_offset_positions, offsets)
    minimal_index = apply_offsets(AMesh, Ind, minimal_offset)

    core_circulant_matrix_format_COO(col[minimal_index], minimal_offset, AMesh)
end

function drop_kron(M::Matrix, AMesh::AM, τ::Tv, offset_ranges::RTup = 1:2,
     rtol::TolT = 700 * eps(real(eltype(M))),
     atol::TolT = 700 * eps(real(eltype(M)))) where {N,Ti<:Integer,
                                                     Tv<:AbstractFloat,
                                                     TolT<:AbstractFloat,
                                                     RTup<:TupleOrRange{Ti},
                                                     Matrix,
                                                     AM<:AbstractMesh{Ti,N}}

    if size(M, 1) != size(M, 2)
        throw(DimensionMismatch("Matrix must be square"))
    end

    if size(M, 1) != length(AMesh)
        throw(DimensionMismatch("Matrix size must be equal to the mesh size"))
    end

    if ndims(M) != 2
        throw(DimensionMismatch("Matrix must be 2D"))
    end

    #Vector storage type
    VecS = matrix_to_vector(Matrix)

    #Generate a list of offsets for given ranges.
    offsets = offset_generator(UniqueZeroOffset, AMesh, offset_ranges)

    #TODO

    error("Not implemented yet")
end

export drop
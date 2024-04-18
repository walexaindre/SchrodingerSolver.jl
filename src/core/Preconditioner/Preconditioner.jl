function drop(M::Matrix, AMesh::AM, τ::Tv, offset_range::R = 0:2,
              rtol::AbstractFloat = 700 * eps(real(eltype(M))),
              atol::AbstractFloat = 700 * eps(real(eltype(M)))) where {V<:Integer,
                                                                        Tv<:AbstractFloat,
                                                                        R<:AbstractRange{V},
                                                                        Matrix<:AbstractArray{Tv},
                                                                        AM<:AbstractMesh{V}}
    offsets = offset_generator(V, offset_range)
    full_locations = apply_offsets(AMesh, 1, offsets)

    rows = size(M, 1)

    e₁ = zeros(eltype(M), rows)
    e₁[1] .= 1

    col, _ = gmres(M, e₁; atol = atol, rtol = rtol)

    #Prune the inverse with the desired drop tolerance τ.
    bit_index = abs2.(col) .> τ

    #Get indexes of non pruned elements.
    int_index = findall(bit_index) 

    #Undropped stencil positions.
    undropped_stencil_positions = findall(in(int_index), full_locations)

    #Here the not eliminated entries are chosen.
    used_stencil_idx = full_locations[undropped_stencil_positions]

    nnz_stencil_positions = length(undropped_stencil_positions)

    #Checking a_ij==a_ji this check is flawed. You need other interpretation for this... [TODO] rework this kind of validation
    if 1:nnz_stencil_positions != undropped_stencil_positions
        @warn "Problem with matrix entries (a_ij!=a_ji). Use a different drop tolerance..."
    end

    #e₁ values at used stencil positions.
    ê₁ = col[used_stencil_idx]
    println(used_stencil_idx)
    
    #core_circulant_matrix_format_IJV(ê₁, undropped_stencil_positions, AMesh)
    
end

export drop
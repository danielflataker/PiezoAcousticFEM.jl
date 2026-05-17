"""
    assemble_k_form_dense(grid, material, formulation, ip, qr)

Assemble dense global K-form blocks for a piezoelectric structure in vacuum.
This is a small reference path for tests and algebraic verification.

- The Ferrite part (`DofHandler`, `CellValues`, `CellIterator`) handles mesh,
  interpolation, DOF numbering, and Gauss integration.
- The physics part (`piezo_element_matrices`) builds Kocbach's local element
  blocks.
- The adapter part (`KFormAssembly`) preserves the mapping from Ferrite DOFs to
  compact Kocbach blocks.
"""
function assemble_k_form_dense(grid, material::AxisymmetricRZPiezoMaterial, formulation::AxisymmetricRZ, ip, qr)
    dh = piezo_dofhandler(grid, ip)
    cellvalues_u = CellValues(qr, ip)
    cellvalues_ϕ = CellValues(qr, ip)

    return assemble_k_form_dense(dh, cellvalues_u, cellvalues_ϕ, material, formulation)
end


function assemble_k_form_dense(
    dh::DofHandler,
    cellvalues_u,
    cellvalues_ϕ,
    material::AxisymmetricRZPiezoMaterial,
    formulation::AxisymmetricRZ,
)
    grid = Ferrite.get_grid(dh)
    dofmap = PiezoFieldDofMap(dh)
    u_dofs = dofmap.u_dofs
    ϕ_dofs = dofmap.ϕ_dofs

    TK, TE, Tε, Tρ = k_form_block_types_from_grid(material, grid)
    Kuu = zeros(TK, length(u_dofs), length(u_dofs))
    Kuϕ = zeros(TE, length(u_dofs), length(ϕ_dofs))
    Kϕu = zeros(TE, length(ϕ_dofs), length(u_dofs))
    Kϕϕ = zeros(Tε, length(ϕ_dofs), length(ϕ_dofs))
    Muu = zeros(Tρ, length(u_dofs), length(u_dofs))

    u_range = dof_range(dh, :u)
    ϕ_range = dof_range(dh, :ϕ)

    for (cellid, cell) in enumerate(CellIterator(grid))
        reinit!(cellvalues_u, cell)
        reinit!(cellvalues_ϕ, cell)

        elem = piezo_element_matrices(
            cellvalues_u,
            cellvalues_ϕ,
            getcoordinates(cell),
            material,
            formulation,
        )

        cell_dofs = celldofs(dh, cellid)
        u_local = [compact_displacement_dof(dofmap, dof) for dof in cell_dofs[u_range]]
        ϕ_local = [compact_potential_dof(dofmap, dof) for dof in cell_dofs[ϕ_range]]

        add_block!(Kuu, u_local, u_local, elem.Kuu)
        add_block!(Kuϕ, u_local, ϕ_local, elem.Kuϕ)
        add_block!(Kϕu, ϕ_local, u_local, elem.Kϕu)
        add_block!(Kϕϕ, ϕ_local, ϕ_local, elem.Kϕϕ)
        add_block!(Muu, u_local, u_local, elem.Muu)
    end

    system = KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    return KFormAssembly(system, dh, dofmap)
end


"""
    assemble_k_form_sparse(grid, material, formulation, ip, qr)

Assemble sparse global K-form blocks with triplet collection. Repeated element
contributions are coalesced by `sparse(I, J, V, m, n)`.
"""
function assemble_k_form_sparse(grid, material::AxisymmetricRZPiezoMaterial, formulation::AxisymmetricRZ, ip, qr)
    dh = piezo_dofhandler(grid, ip)
    cellvalues_u = CellValues(qr, ip)
    cellvalues_ϕ = CellValues(qr, ip)

    return assemble_k_form_sparse(dh, cellvalues_u, cellvalues_ϕ, material, formulation)
end


function assemble_k_form_sparse(
    dh::DofHandler,
    cellvalues_u,
    cellvalues_ϕ,
    material::AxisymmetricRZPiezoMaterial,
    formulation::AxisymmetricRZ,
)
    grid = Ferrite.get_grid(dh)
    dofmap = PiezoFieldDofMap(dh)
    u_dofs = dofmap.u_dofs
    ϕ_dofs = dofmap.ϕ_dofs

    TK, TE, Tε, Tρ = k_form_block_types_from_grid(material, grid)
    Kuu_triplets = matrix_triplets(TK)
    Kuϕ_triplets = matrix_triplets(TE)
    Kϕu_triplets = matrix_triplets(TE)
    Kϕϕ_triplets = matrix_triplets(Tε)
    Muu_triplets = matrix_triplets(Tρ)

    u_range = dof_range(dh, :u)
    ϕ_range = dof_range(dh, :ϕ)
    ncells = getncells(grid)
    nᵤ_cell = length(u_range)
    nϕ_cell = length(ϕ_range)
    sizehint_triplets!(Kuu_triplets, ncells * nᵤ_cell * nᵤ_cell)
    sizehint_triplets!(Kuϕ_triplets, ncells * nᵤ_cell * nϕ_cell)
    sizehint_triplets!(Kϕu_triplets, ncells * nϕ_cell * nᵤ_cell)
    sizehint_triplets!(Kϕϕ_triplets, ncells * nϕ_cell * nϕ_cell)
    sizehint_triplets!(Muu_triplets, ncells * nᵤ_cell * nᵤ_cell)

    u_local = Vector{Int}(undef, nᵤ_cell)
    ϕ_local = Vector{Int}(undef, nϕ_cell)

    for (cellid, cell) in enumerate(CellIterator(grid))
        reinit!(cellvalues_u, cell)
        reinit!(cellvalues_ϕ, cell)

        elem = piezo_element_matrices(
            cellvalues_u,
            cellvalues_ϕ,
            getcoordinates(cell),
            material,
            formulation,
        )

        cell_dofs = celldofs(dh, cellid)
        fill_compact_dofs!(u_local, dofmap, cell_dofs[u_range], :u)
        fill_compact_dofs!(ϕ_local, dofmap, cell_dofs[ϕ_range], :ϕ)

        append_block_triplets!(Kuu_triplets, u_local, u_local, elem.Kuu)
        append_block_triplets!(Kuϕ_triplets, u_local, ϕ_local, elem.Kuϕ)
        append_block_triplets!(Kϕu_triplets, ϕ_local, u_local, elem.Kϕu)
        append_block_triplets!(Kϕϕ_triplets, ϕ_local, ϕ_local, elem.Kϕϕ)
        append_block_triplets!(Muu_triplets, u_local, u_local, elem.Muu)
    end

    nᵤ = length(u_dofs)
    nᵩ = length(ϕ_dofs)
    Kuu = sparse_matrix(Kuu_triplets, nᵤ, nᵤ)
    Kuϕ = sparse_matrix(Kuϕ_triplets, nᵤ, nᵩ)
    Kϕu = sparse_matrix(Kϕu_triplets, nᵩ, nᵤ)
    Kϕϕ = sparse_matrix(Kϕϕ_triplets, nᵩ, nᵩ)
    Muu = sparse_matrix(Muu_triplets, nᵤ, nᵤ)

    system = KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    return KFormAssembly(system, dh, dofmap)
end


function add_block!(A, rows, cols, block)
    @boundscheck size(block) == (length(rows), length(cols)) ||
        throw(DimensionMismatch("block size does not match assembly indices"))

    @inbounds for (j_local, j_global) in pairs(cols)
        for (i_local, i_global) in pairs(rows)
            A[i_global, j_global] += block[i_local, j_local]
        end
    end

    return A
end


matrix_triplets(::Type{T}) where {T} = (Int[], Int[], T[])


function sizehint_triplets!(triplets, n::Integer)
    I, J, V = triplets
    sizehint!(I, n)
    sizehint!(J, n)
    sizehint!(V, n)

    return triplets
end


function fill_compact_dofs!(dest, dofmap::PiezoFieldDofMap, ferrite_dofs, field::Symbol)
    @boundscheck length(dest) == length(ferrite_dofs) ||
        throw(DimensionMismatch("compact DOF buffer length does not match cell DOFs"))

    if field === :u
        @inbounds for i in eachindex(dest, ferrite_dofs)
            dest[i] = compact_displacement_dof(dofmap, ferrite_dofs[i])
        end
    elseif field === :ϕ
        @inbounds for i in eachindex(dest, ferrite_dofs)
            dest[i] = compact_potential_dof(dofmap, ferrite_dofs[i])
        end
    else
        throw(ArgumentError("unknown piezo field $field"))
    end

    return dest
end


function append_block_triplets!(triplets, rows, cols, block)
    @boundscheck size(block) == (length(rows), length(cols)) ||
        throw(DimensionMismatch("block size does not match assembly indices"))

    I, J, V = triplets
    @inbounds for (j_local, j_global) in pairs(cols)
        for (i_local, i_global) in pairs(rows)
            push!(I, i_global)
            push!(J, j_global)
            push!(V, block[i_local, j_local])
        end
    end

    return triplets
end


function sparse_matrix(triplets, nrows, ncols)
    I, J, V = triplets
    return sparse(I, J, V, nrows, ncols)
end


function k_form_block_types_from_grid(material::AxisymmetricRZPiezoMaterial, grid)
    Tx = eltype(first(getnodes(grid)).x)

    return (
        promote_type(eltype(material.cᴱ), Tx),
        promote_type(eltype(material.e), Tx),
        promote_type(eltype(material.εˢ), Tx),
        promote_type(typeof(material.ρ), Tx),
    )
end

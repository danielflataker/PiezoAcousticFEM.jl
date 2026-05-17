"""
    reconstruct_fields(assembly, solution)

Reconstruct nodal fields from compact K-form solution vectors.
"""
function reconstruct_fields(assembly::KFormAssembly, solution::DirectVoltageSolution)
    grid = Ferrite.get_grid(assembly.dofhandler)

    return (
        displacement=nodal_displacement(assembly, solution.displacement),
        potential=nodal_potential(assembly, solution.potential),
        radius=[Ferrite.get_node_coordinate(grid, nodeid)[1] for nodeid in 1:getnnodes(grid)],
    )
end


reconstruct_fields(result::HarmonicVoltageResult) =
    reconstruct_fields(result.assembled.assembly, result.solution)


"""
    reconstruct_fields(result::ShortCircuitModalResult, mode_index)

Reconstruct nodal displacement and potential fields for one short-circuit mode.
The compact displacement mode is already stored on the result; the internal
electric potentials are recovered from the condensed short-circuit K-form.
"""
function reconstruct_fields(result::ShortCircuitModalResult, mode_index::Integer)
    1 <= mode_index <= size(result.modes, 2) ||
        throw(ArgumentError("mode_index must be in 1:$(size(result.modes, 2))"))

    assembly = result.assembled.assembly
    mode = result.modes[:, mode_index]
    potential = reconstruct_short_circuit_modal_potential(result, mode)
    grid = Ferrite.get_grid(assembly.dofhandler)

    return (
        displacement=nodal_displacement(assembly, mode),
        potential=nodal_potential(assembly, potential),
        radius=[Ferrite.get_node_coordinate(grid, nodeid)[1] for nodeid in 1:getnnodes(grid)],
        eigenvalue=result.eigenvalues[mode_index],
        angular_frequency=result.angular_frequencies[mode_index],
        frequency=result.frequencies[mode_index],
        normalization=result.normalization,
    )
end


function reconstruct_short_circuit_modal_potential(result::ShortCircuitModalResult, mode)
    assembly = result.assembled.assembly
    partition = result.reduction.partition
    system = assembly.system
    potential = zeros(eltype(mode), partition.nϕ)

    if !isempty(partition.internal)
        Kii = Matrix(system.Kϕϕ[partition.internal, partition.internal])
        Kiu = Matrix(system.Kϕu[partition.internal, :])
        potential[partition.internal] = -(Kii \ (Kiu * mode))
    end

    return potential
end


function nodal_displacement(assembly::KFormAssembly, displacement)
    dh = assembly.dofhandler
    grid = Ferrite.get_grid(dh)
    values = [zero(Vec{2, eltype(displacement)}) for _ in 1:getnnodes(grid)]
    dofmap = assembly.dofmap

    for nodeid in 1:getnnodes(grid)
        uᵣ = displacement[compact_displacement_dof(dofmap, nodeid, 1)]
        uz = displacement[compact_displacement_dof(dofmap, nodeid, 2)]
        values[nodeid] = Vec{2}((uᵣ, uz))
    end

    return values
end


function nodal_potential(assembly::KFormAssembly, potential)
    dh = assembly.dofhandler
    grid = Ferrite.get_grid(dh)
    values = zeros(eltype(potential), getnnodes(grid))
    dofmap = assembly.dofmap

    for nodeid in 1:getnnodes(grid)
        values[nodeid] = potential[compact_potential_node_dof(dofmap, nodeid)]
    end

    return values
end

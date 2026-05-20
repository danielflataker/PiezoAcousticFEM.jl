"""
    CompactFieldDofs(displacement, potential; metadata)

Field values in this package's compact displacement and potential indexing,
not grid-node indexing.
"""
struct CompactFieldDofs{D,P,M}
    displacement::D
    potential::P
    metadata::M
end

CompactFieldDofs(displacement, potential; metadata=(output_kind=:compact_field_dofs,)) =
    CompactFieldDofs(displacement, potential, metadata)

"""
    NodalFieldOutput(displacement, potential, radius; metadata)

Nodal fields reconstructed on the grid vertices for visualization and simple
postprocessing.
"""
struct NodalFieldOutput{D,P,R,M}
    displacement::D
    potential::P
    radius::R
    metadata::M
end

NodalFieldOutput(displacement, potential, radius; metadata=(output_kind=:nodal_fields,)) =
    NodalFieldOutput(displacement, potential, radius, metadata)

"""
    reconstruct_field_dofs(assembly, solution)

Return compact field-DOF values from a K-form solution.
"""
function reconstruct_field_dofs(::KFormAssembly, solution::DirectVoltageSolution)
    return CompactFieldDofs(solution.displacement, solution.potential)
end


reconstruct_field_dofs(result::HarmonicVoltageResult) =
    reconstruct_field_dofs(result.assembled.assembly, result.solution)


function reconstruct_field_dofs(result::ShortCircuitModalResult, mode_index::Integer)
    1 <= mode_index <= size(result.modes, 2) ||
        throw(ArgumentError("mode_index must be in 1:$(size(result.modes, 2))"))

    mode = result.modes[:, mode_index]

    return CompactFieldDofs(
        mode,
        reconstruct_short_circuit_modal_potential(result, mode);
        metadata=(
            output_kind=:compact_field_dofs,
            analysis=:short_circuit_modal,
            mode_index=mode_index,
            normalization=result.normalization,
        ),
    )
end


"""
    reconstruct_fields(assembly, solution)

Reconstruct grid-vertex fields from compact K-form solution vectors. This is a
simple visualization path; compact field DOFs remain the solver truth.
"""
function reconstruct_fields(assembly::KFormAssembly, solution::DirectVoltageSolution; metadata=(output_kind=:nodal_fields,))
    grid = ferrite_grid(assembly)

    return NodalFieldOutput(
        nodal_displacement(assembly, solution.displacement),
        nodal_potential(assembly, solution.potential),
        [Ferrite.get_node_coordinate(grid, nodeid)[1] for nodeid in 1:getnnodes(grid)];
        metadata,
    )
end


reconstruct_fields(result::HarmonicVoltageResult) =
    reconstruct_fields(
        result.assembled.assembly,
        result.solution;
        metadata=(
            output_kind=:nodal_fields,
            analysis=:harmonic_voltage,
            harmonic_convention=result.analysis.convention,
            quantity_interpretation="complex amplitude",
        ),
    )


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
    grid = ferrite_grid(assembly)

    return NodalFieldOutput(
        nodal_displacement(assembly, mode),
        nodal_potential(assembly, potential),
        [Ferrite.get_node_coordinate(grid, nodeid)[1] for nodeid in 1:getnnodes(grid)];
        metadata=(
            output_kind=:nodal_fields,
            analysis=:short_circuit_modal,
            mode_index=mode_index,
            eigenvalue=result.eigenvalues[mode_index],
            angular_frequency=result.angular_frequencies[mode_index],
            frequency=result.frequencies[mode_index],
            normalization=result.normalization,
            quantity_interpretation="mode shape",
        ),
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
    grid = ferrite_grid(dh)
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
    grid = ferrite_grid(dh)
    values = zeros(eltype(potential), getnnodes(grid))
    dofmap = assembly.dofmap

    for nodeid in 1:getnnodes(grid)
        values[nodeid] = potential[compact_potential_node_dof(dofmap, nodeid)]
    end

    return values
end

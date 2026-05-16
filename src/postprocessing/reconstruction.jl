"""
    reconstruct_fields(assembly, solution)

Reconstruct nodal fields from compact K-form solution vectors.
"""
function reconstruct_fields(assembly::KFormAssembly, solution::DirectVoltageSolution)
    grid = Ferrite.get_grid(assembly.dofhandler)

    return (
        displacement=nodal_displacement(assembly, solution.displacement),
        potential=nodal_potential(assembly, solution.potential),
        radius=[node.x[1] for node in getnodes(grid)],
    )
end


reconstruct_fields(result::HarmonicVoltageResult) =
    reconstruct_fields(result.assembled.assembly, result.solution)


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

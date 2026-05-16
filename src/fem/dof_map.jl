"""
    PiezoFieldDofMap

Én kompakt mapping mellom Ferrite sine globale DOF-er og K-form-blokkenes
feltvise DOF-er.
"""
struct PiezoFieldDofMap
    u_dofs::Vector{Int}
    ϕ_dofs::Vector{Int}
    u_index::Dict{Int,Int}
    ϕ_index::Dict{Int,Int}
    node_to_u::Dict{Tuple{Int,Int},Int}
    node_to_phi::Dict{Int,Int}

    function PiezoFieldDofMap(
        u_dofs,
        ϕ_dofs;
        node_to_u=Dict{Tuple{Int,Int},Int}(),
        node_to_phi=Dict{Int,Int}(),
    )
        u = collect(Int, u_dofs)
        ϕ = collect(Int, ϕ_dofs)

        return new(
            u,
            ϕ,
            Dict(dof => i for (i, dof) in pairs(u)),
            Dict(dof => i for (i, dof) in pairs(ϕ)),
            node_to_u,
            node_to_phi,
        )
    end
end


function PiezoFieldDofMap(dh::DofHandler)
    u_dofs = _field_dofs(dh, :u)
    ϕ_dofs = _field_dofs(dh, :ϕ)
    u_index = Dict(dof => i for (i, dof) in pairs(u_dofs))
    ϕ_index = Dict(dof => i for (i, dof) in pairs(ϕ_dofs))
    node_to_u = Dict{Tuple{Int,Int},Int}()
    node_to_phi = Dict{Int,Int}()

    u_range = dof_range(dh, :u)
    ϕ_range = dof_range(dh, :ϕ)
    grid = Ferrite.get_grid(dh)

    for cellid in 1:getncells(grid)
        cell = getcells(grid, cellid)
        cell_dofs = celldofs(dh, cellid)
        cell_u_dofs = cell_dofs[u_range]
        cell_ϕ_dofs = cell_dofs[ϕ_range]

        for (a, nodeid) in pairs(cell.nodes)
            node_to_u[(nodeid, 1)] = u_index[cell_u_dofs[2a - 1]]
            node_to_u[(nodeid, 2)] = u_index[cell_u_dofs[2a]]
            node_to_phi[nodeid] = ϕ_index[cell_ϕ_dofs[a]]
        end
    end

    return PiezoFieldDofMap(u_dofs, ϕ_dofs; node_to_u, node_to_phi)
end

compact_displacement_dof(dofmap::PiezoFieldDofMap, ferrite_dof::Integer) =
    dofmap.u_index[Int(ferrite_dof)]

compact_potential_dof(dofmap::PiezoFieldDofMap, ferrite_dof::Integer) =
    dofmap.ϕ_index[Int(ferrite_dof)]

compact_displacement_dof(dofmap::PiezoFieldDofMap, nodeid::Integer, component::Integer) =
    dofmap.node_to_u[(Int(nodeid), Int(component))]

compact_displacement_dof(dofmap::PiezoFieldDofMap, nodeid::Integer, component::Symbol) =
    compact_displacement_dof(dofmap, nodeid, _displacement_component_index(component))

compact_potential_node_dof(dofmap::PiezoFieldDofMap, nodeid::Integer) =
    dofmap.node_to_phi[Int(nodeid)]


"""
    _field_dofs(dh, field)

Ferrite adapter: get all global Ferrite DOFs for one field.
"""
function _field_dofs(dh::DofHandler, field::Symbol)
    range = dof_range(dh, field)
    dofs = Int[]

    for cellid in 1:getncells(Ferrite.get_grid(dh))
        append!(dofs, celldofs(dh, cellid)[range])
    end

    return sort!(unique!(dofs))
end

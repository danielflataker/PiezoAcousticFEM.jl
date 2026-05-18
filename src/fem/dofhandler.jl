"""
    potential_dofs_on_facets(assembly, facetset)

Ferrite adapter: find compact `ϕ` indices for potential DOFs on a Ferrite
facetset. These indices can be used directly in analysis-specific electrical
DOF partitions.
"""
function potential_dofs_on_facets(assembly::KFormAssembly, facetset)
    dh = assembly.dofhandler
    dofmap = assembly.dofmap
    ferrite_dofs = ferrite_potential_dofs_on_facets(dh, facetset)
    indices = [compact_potential_dof(dofmap, dof) for dof in ferrite_dofs]

    return sort!(unique!(indices))
end


"""
    displacement_component_dofs_on_facets(assembly, facetset, component)

Find compact `u` indices for a displacement component on a Ferrite facetset.
`component` may be `:r`/`:ur` or `:z`/`:uz`.
"""
function displacement_component_dofs_on_facets(assembly::KFormAssembly, facetset, component::Symbol)
    dh = assembly.dofhandler
    dofmap = assembly.dofmap
    ferrite_dofs = ferrite_displacement_component_dofs_on_facets(dh, facetset, component)
    indices = [compact_displacement_dof(dofmap, dof) for dof in ferrite_dofs]

    return sort!(unique!(indices))
end


"""
    axis_radial_displacement_constraint(assembly; axis, value=0)

Build a homogeneous Dirichlet condition for `u_r` on the symmetry axis. `axis`
must be an explicit `AxisBoundary`.
"""
function axis_radial_displacement_constraint(assembly::KFormAssembly; axis::AxisBoundary, value=0)
    axis_facets = resolve_facetset(ferrite_grid(assembly), axis)
    indices = displacement_component_dofs_on_facets(assembly, axis_facets, :r)

    return DirichletDofs(indices, fill(value, length(indices)))
end


"""
    potential_partition(assembly; driven, grounded)

Build Kocbach's potential partition from explicit electrode objects.
"""
function potential_partition(assembly::KFormAssembly; driven::FacetElectrode, grounded::FacetElectrode)
    grid = ferrite_grid(assembly)
    driven_facets = resolve_facetset(grid, driven)
    grounded_facets = resolve_facetset(grid, grounded)

    p = potential_dofs_on_facets(assembly, driven_facets)
    g = potential_dofs_on_facets(assembly, grounded_facets)
    nϕ = length(assembly.dofmap.ϕ_dofs)
    i = setdiff(collect(1:nϕ), union(p, g))

    return HarmonicVoltageDofPartition(nϕ, i, p, g)
end


"""
    short_circuit_partition(assembly, electrodes)

Build the electric DOF partition for short-circuit analyses. Both terminals are
homogeneous electric Dirichlet sets.
"""
function short_circuit_partition(assembly::KFormAssembly, electrodes::TwoTerminalElectrodes)
    grid = ferrite_grid(assembly)
    signal_facets = resolve_facetset(grid, electrodes.signal)
    reference_facets = resolve_facetset(grid, electrodes.reference)

    signal = potential_dofs_on_facets(assembly, signal_facets)
    reference = potential_dofs_on_facets(assembly, reference_facets)
    grounded = sort!(unique!(vcat(signal, reference)))
    nϕ = length(assembly.dofmap.ϕ_dofs)
    internal = setdiff(collect(1:nϕ), grounded)

    return ShortCircuitDofPartition(nϕ, internal, grounded)
end


"""
    mechanical_dirichlet(assembly, boundary_conditions)

Build homogeneous mechanical Dirichlet DOFs from the problem boundary-condition
container. Returns `nothing` when no mechanical constraints are present.
"""
function mechanical_dirichlet(assembly::KFormAssembly, bcs::AxisymmetricBoundaryConditions)
    constraints = [
        _mechanical_dirichlet(assembly, bc)
        for bc in bcs.mechanical
    ]

    return combine_dirichlet_dofs(constraints)
end

_mechanical_dirichlet(assembly::KFormAssembly, axis::AxisBoundary) =
    axis_radial_displacement_constraint(assembly; axis)

function combine_dirichlet_dofs(constraints)
    isempty(constraints) && return nothing

    active = [constraint for constraint in constraints if constraint !== nothing]
    isempty(active) && return nothing

    indices = reduce(vcat, (constraint.indices for constraint in active); init=Int[])
    values = reduce(vcat, (constraint.values for constraint in active))
    isempty(indices) && return nothing

    order = sortperm(indices)
    sorted_indices = indices[order]
    sorted_values = values[order]
    length(unique(sorted_indices)) == length(sorted_indices) ||
        throw(ArgumentError("mechanical constraints contain duplicate displacement DOFs"))

    return DirichletDofs(sorted_indices, sorted_values)
end

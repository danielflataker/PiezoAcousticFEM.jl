# Public adapter API ---------------------------------------------------------

"""
    KFormAssembly

Adapterresultat mellom Ferrite og Kocbach-blokkene. Ferrite eier mesh,
celler, interpolasjon og global DOF-nummerering. `system` er vår
fysikkrepresentasjon på K-form, med kompakte blokker.
"""
struct KFormAssembly{S,DH,M}
    system::S
    dofhandler::DH
    dofmap::M
end

"""
    piezo_dofhandler(grid, ip)

Ferrite-adapter: lag en `DofHandler` for feltene `:u` og `:ϕ`.
Interpolasjonen `ip` er skalar; forskyvningsfeltet bruker `ip^2`.
"""
function piezo_dofhandler(grid, ip)
    dh = DofHandler(grid)
    add!(dh, :u, ip^2)
    add!(dh, :ϕ, ip)
    close!(dh)
    return dh
end


"""
    potential_dofs_on_facets(assembly, facetset)

Ferrite-adapter: finn kompakte `ϕ`-indekser for potensial-DOF-er på et
Ferrite facetset. Disse indeksene kan brukes direkte i analyse-spesifikke
elektriske DOF-partisjoner.
"""
function potential_dofs_on_facets(assembly::KFormAssembly, facetset)
    dh = assembly.dofhandler
    dofmap = assembly.dofmap
    ferrite_dofs = _potential_ferrite_dofs_on_facets(dh, facetset)
    indices = [compact_potential_dof(dofmap, dof) for dof in ferrite_dofs]

    return sort!(unique!(indices))
end


"""
    displacement_component_dofs_on_facets(assembly, facetset, component)

Finn kompakte `u`-indekser for en displacement-komponent på et Ferrite
facetset. `component` kan være `:r`/`:ur` eller `:z`/`:uz`.
"""
function displacement_component_dofs_on_facets(assembly::KFormAssembly, facetset, component::Symbol)
    dh = assembly.dofhandler
    dofmap = assembly.dofmap
    ferrite_dofs = _displacement_component_ferrite_dofs_on_facets(dh, facetset, component)
    indices = [compact_displacement_dof(dofmap, dof) for dof in ferrite_dofs]

    return sort!(unique!(indices))
end


"""
    axis_radial_displacement_constraint(assembly; axis, value=0)

Lag en homogen Dirichlet-betingelse for `u_r` på symmetriaksen. `axis` kan
være en eksplisitt `AxisBoundary`.
"""
function axis_radial_displacement_constraint(assembly::KFormAssembly; axis::AxisBoundary, value=0)
    grid = Ferrite.get_grid(assembly.dofhandler)
    axis_facets = _resolve_facetset(grid, axis)
    indices = displacement_component_dofs_on_facets(assembly, axis_facets, :r)

    return DirichletDofs(indices, value)
end


"""
    potential_partition(assembly; driven, grounded)

Bygg Kocbachs potensialpartisjon fra eksplisitte elektrodeobjekter.
"""
function potential_partition(assembly::KFormAssembly; driven::FacetElectrode, grounded::FacetElectrode)
    grid = Ferrite.get_grid(assembly.dofhandler)
    driven_facets = _resolve_facetset(grid, driven)
    grounded_facets = _resolve_facetset(grid, grounded)

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
    grid = Ferrite.get_grid(assembly.dofhandler)
    signal_facets = _resolve_facetset(grid, electrodes.signal)
    reference_facets = _resolve_facetset(grid, electrodes.reference)

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
        mechanical_dirichlet(assembly, bc)
        for bc in bcs.mechanical
    ]

    return combine_dirichlet_dofs(constraints)
end

mechanical_dirichlet(assembly::KFormAssembly, axis::AxisBoundary) =
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

# Ferrite DOF helpers --------------------------------------------------------

function _potential_ferrite_dofs_on_facets(dh::DofHandler, facetset)
    interpolation = Ferrite.getfieldinterpolation(dh, Ferrite.find_field(dh, :ϕ))
    facet_dofs = Ferrite.dirichlet_facetdof_indices(interpolation)
    field_range = dof_range(dh, :ϕ)
    field_offset = first(field_range) - 1
    dofs = Int[]

    for facet in facetset
        cellid, facetid = facet.idx
        cell_dofs = celldofs(dh, cellid)
        local_dofs = field_offset .+ collect(facet_dofs[facetid])
        append!(dofs, cell_dofs[local_dofs])
    end

    return sort!(unique!(dofs))
end


function _displacement_component_ferrite_dofs_on_facets(dh::DofHandler, facetset, component::Symbol)
    component_index = _displacement_component_index(component)
    scalar_interpolation = _scalar_base_interpolation(dh, :u)
    facet_dofs = Ferrite.dirichlet_facetdof_indices(scalar_interpolation)
    field_range = dof_range(dh, :u)
    field_offset = first(field_range) - 1
    dofs = Int[]

    for facet in facetset
        cellid, facetid = facet.idx
        cell_dofs = celldofs(dh, cellid)
        local_dofs = [
            field_offset + 2 * (scalar_dof - 1) + component_index
            for scalar_dof in facet_dofs[facetid]
        ]
        append!(dofs, cell_dofs[local_dofs])
    end

    return sort!(unique!(dofs))
end


# Private helpers ------------------------------------------------------------

_resolve_facetset(grid, electrode::FacetElectrode{<:AbstractString}) =
    getfacetset(grid, String(electrode.facetset))
_resolve_facetset(_, electrode::FacetElectrode) = electrode.facetset
_resolve_facetset(grid, boundary::FacetBoundary{<:AbstractString}) =
    getfacetset(grid, String(boundary.facetset))
_resolve_facetset(_, boundary::FacetBoundary) = boundary.facetset
_resolve_facetset(grid, boundary::AxisBoundary) = _resolve_facetset(grid, boundary.boundary)

_displacement_component_index(component::Symbol) =
    component in (:r, :ur, :u_r) ? 1 :
    component in (:z, :uz, :u_z) ? 2 :
    throw(ArgumentError("unknown displacement component $component; expected :r/:ur or :z/:uz"))

function _scalar_base_interpolation(dh::DofHandler, field::Symbol)
    interpolation = Ferrite.getfieldinterpolation(dh, Ferrite.find_field(dh, field))

    return _scalar_base_interpolation(interpolation)
end

_scalar_base_interpolation(interpolation) = interpolation.ip

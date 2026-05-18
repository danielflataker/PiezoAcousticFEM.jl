# Ferrite backend boundary ---------------------------------------------------

"""
    KFormAssembly

Adapter result between Ferrite and the Kocbach blocks. Ferrite owns the mesh,
cells, interpolation, and global DOF numbering. `system` is this package's
physics representation in K-form, with compact blocks.
"""
struct KFormAssembly{S,DH,M}
    system::S
    dofhandler::DH
    dofmap::M
end

"""
    piezo_dofhandler(grid, ip)

Build the Ferrite `DofHandler` for the current piezoelectric fields. The
interpolation `ip` is scalar; the displacement field uses `ip^2`.
"""
function piezo_dofhandler(grid, ip)
    dh = DofHandler(grid)
    add!(dh, :u, ip^2)
    add!(dh, :ϕ, ip)
    close!(dh)
    return dh
end

ferrite_grid(assembly::KFormAssembly) = ferrite_grid(assembly.dofhandler)
ferrite_grid(dh::DofHandler) = Ferrite.get_grid(dh)

piezo_field_dof_range(dh::DofHandler, field::Symbol) = dof_range(dh, field)

function ferrite_field_dofs(dh::DofHandler, field::Symbol)
    range = piezo_field_dof_range(dh, field)
    dofs = Int[]

    for cell in CellIterator(dh)
        append!(dofs, celldofs(cell)[range])
    end

    return sort!(unique!(dofs))
end

function ferrite_potential_dofs_on_facets(dh::DofHandler, facetset)
    interpolation = Ferrite.getfieldinterpolation(dh, Ferrite.find_field(dh, :ϕ))
    facet_dofs = Ferrite.dirichlet_facetdof_indices(interpolation)
    field_range = piezo_field_dof_range(dh, :ϕ)
    field_offset = first(field_range) - 1
    dofs = Int[]
    facet_cache = FacetCache(dh)

    for facet in facetset
        facetid = facet[2]
        reinit!(facet_cache, facet)
        cell_dofs = celldofs(facet_cache)
        local_dofs = field_offset .+ collect(facet_dofs[facetid])
        append!(dofs, cell_dofs[local_dofs])
    end

    return sort!(unique!(dofs))
end

function ferrite_displacement_component_dofs_on_facets(dh::DofHandler, facetset, component::Symbol)
    component_index = displacement_component_index(component)
    scalar_interpolation = scalar_base_interpolation(dh, :u)
    facet_dofs = Ferrite.dirichlet_facetdof_indices(scalar_interpolation)
    field_range = piezo_field_dof_range(dh, :u)
    field_offset = first(field_range) - 1
    dofs = Int[]
    facet_cache = FacetCache(dh)

    for facet in facetset
        facetid = facet[2]
        reinit!(facet_cache, facet)
        cell_dofs = celldofs(facet_cache)
        local_dofs = [
            field_offset + 2 * (scalar_dof - 1) + component_index
            for scalar_dof in facet_dofs[facetid]
        ]
        append!(dofs, cell_dofs[local_dofs])
    end

    return sort!(unique!(dofs))
end

displacement_component_index(component::Symbol) =
    component in (:r, :ur, :u_r) ? 1 :
    component in (:z, :uz, :u_z) ? 2 :
    throw(ArgumentError("unknown displacement component $component; expected :r/:ur or :z/:uz"))

function scalar_base_interpolation(dh::DofHandler, field::Symbol)
    interpolation = Ferrite.getfieldinterpolation(dh, Ferrite.find_field(dh, field))

    return scalar_base_interpolation(interpolation)
end

scalar_base_interpolation(interpolation) = interpolation.ip

function resolve_facetset(grid, boundary_or_electrode; label="facet set")
    facetset = try
        _resolve_facetset(grid, boundary_or_electrode)
    catch err
        err isa KeyError || rethrow()
        throw(ArgumentError("$label $(err.key) does not exist"))
    end

    isempty(facetset) &&
        throw(ArgumentError("$label must contain at least one facet"))

    return facetset
end

_resolve_facetset(grid, electrode::FacetElectrode{<:AbstractString}) =
    getfacetset(grid, String(electrode.facetset))
_resolve_facetset(_, electrode::FacetElectrode) = electrode.facetset
_resolve_facetset(grid, boundary::FacetBoundary{<:AbstractString}) =
    getfacetset(grid, String(boundary.facetset))
_resolve_facetset(_, boundary::FacetBoundary) = boundary.facetset
_resolve_facetset(grid, boundary::AxisBoundary) = _resolve_facetset(grid, boundary.boundary)

function ferrite_facet_measure(grid, interpolation, facetset)
    T = eltype(Ferrite.get_node_coordinate(grid, 1))
    measure = zero(T)
    facetvalues = FacetValues(FacetQuadratureRule{Ferrite.getrefshape(interpolation)}(1), interpolation)

    for facet in FacetIterator(grid, facetset)
        reinit!(facetvalues, facet)
        for q_point in 1:getnquadpoints(facetvalues)
            measure += getdetJdV(facetvalues, q_point)
        end
    end

    return measure
end

function ferrite_facet_radius_range(grid, interpolation, facetset)
    rmin = Inf
    rmax = -Inf
    facetvalues = FacetValues(FacetQuadratureRule{Ferrite.getrefshape(interpolation)}(2), interpolation)

    for facet in FacetIterator(grid, facetset)
        coordinates = getcoordinates(facet)
        reinit!(facetvalues, facet)

        for q_point in 1:getnquadpoints(facetvalues)
            x = spatial_coordinate(facetvalues, q_point, coordinates)
            r = x[1]
            rmin = min(rmin, r)
            rmax = max(rmax, r)
        end
    end

    return rmin, rmax
end

function ferrite_geometry_tolerance(grid)
    scale = zero(eltype(Ferrite.get_node_coordinate(grid, 1)))

    for nodeid in 1:getnnodes(grid)
        scale = max(scale, maximum(abs, Ferrite.get_node_coordinate(grid, nodeid)))
    end

    return max(scale, one(scale)) * sqrt(eps(float(scale)))
end

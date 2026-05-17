"""
    validate(problem)

Preflight a piezoelectric problem before matrix assembly. The current supported
discretization is nodal equal-order Lagrange interpolation for displacement and
electric potential on a 2D axisymmetric `r-z` grid.
"""
function validate(problem::PiezoProblem)
    validate_supported_discretization(problem)
    validate_boundary_semantics(problem)

    return nothing
end


function validate_boundary_semantics(problem::PiezoProblem)
    grid = problem.grid
    interpolation = problem.interpolation
    tol = _geometry_tolerance(grid)

    signal = _validate_electrode(grid, interpolation, problem.electrodes.signal, "signal electrode", tol)
    reference = _validate_electrode(grid, interpolation, problem.electrodes.reference, "reference electrode", tol)
    isempty(intersect(signal, reference)) ||
        throw(ArgumentError("signal and reference electrodes must use disjoint facets"))

    for (i, boundary) in pairs(problem.boundary_conditions.mechanical)
        _validate_mechanical_boundary(grid, interpolation, boundary, "mechanical boundary $i", tol)
    end

    return nothing
end


function validate_supported_discretization(problem::PiezoProblem)
    grid = problem.grid
    interpolation = problem.interpolation
    quadrature = problem.quadrature

    getncells(grid) > 0 ||
        throw(ArgumentError("unsupported discretization: grid must contain at least one cell"))

    Ferrite.getspatialdim(grid) == 2 ||
        throw(ArgumentError("unsupported discretization: axisymmetric r-z models require 2D node coordinates"))

    _is_lagrange_interpolation(interpolation) ||
        throw(ArgumentError("unsupported discretization: interpolation must be a nodal Lagrange interpolation"))

    ip_shape = Ferrite.getrefshape(interpolation)
    qr_shape = Ferrite.getrefshape(quadrature)
    ip_shape == qr_shape ||
        throw(ArgumentError(
            "unsupported discretization: quadrature reference shape $qr_shape must match interpolation reference shape $ip_shape",
        ))

    nbase = Ferrite.getnbasefunctions(interpolation)

    for cellid in 1:getncells(grid)
        cell = getcells(grid, cellid)
        cell_shape = Ferrite.getrefshape(typeof(cell))
        nnodes = length(cell.nodes)

        cell_shape == ip_shape ||
            throw(ArgumentError(
                "unsupported discretization: cell $cellid has reference shape $cell_shape but interpolation uses $ip_shape",
            ))
        nbase == nnodes ||
            throw(ArgumentError(
                "unsupported discretization: cell $cellid has $nnodes mesh nodes but interpolation has $nbase scalar shape functions; higher-order or interior DOFs are not supported yet",
            ))
    end

    return nothing
end


_is_lagrange_interpolation(::Lagrange) = true
_is_lagrange_interpolation(_) = false


function _validate_electrode(grid, interpolation, electrode::FacetElectrode, label::AbstractString, tol)
    facetset = _validated_facetset(grid, electrode, label)
    measure = _facetset_measure(grid, interpolation, facetset)
    measure > tol ||
        throw(ArgumentError("$label must have nonzero boundary measure; observed measure=$measure with tolerance=$tol"))

    return facetset
end


function _validate_mechanical_boundary(grid, interpolation, boundary::AxisBoundary, label::AbstractString, tol)
    facetset = _validated_facetset(grid, boundary, label)
    rmin, rmax = _facetset_radius_range(grid, interpolation, facetset)
    max(abs(rmin), abs(rmax)) <= tol ||
        throw(ArgumentError(
            "$label must lie on the axis r=0; observed r range=[$rmin, $rmax] with tolerance=$tol",
        ))

    return nothing
end


function _validated_facetset(grid, boundary_or_electrode, label::AbstractString)
    facetset = try
        _resolve_facetset(grid, boundary_or_electrode)
    catch err
        err isa KeyError || rethrow()
        throw(ArgumentError("$label facet set $(err.key) does not exist"))
    end

    isempty(facetset) &&
        throw(ArgumentError("$label facet set must contain at least one facet"))

    return facetset
end


function _facetset_measure(grid, interpolation, facetset)
    measure = zero(eltype(first(getnodes(grid)).x))
    facet_dofs = Ferrite.dirichlet_facetdof_indices(interpolation)

    for facet in facetset
        cellid, facetid = facet.idx
        cell = getcells(grid, cellid)
        facet_nodes = cell.nodes[collect(facet_dofs[facetid])]
        x1 = Ferrite.get_node_coordinate(grid, first(facet_nodes))
        x2 = Ferrite.get_node_coordinate(grid, last(facet_nodes))
        measure += norm(x2 - x1)
    end

    return measure
end


function _facetset_radius_range(grid, interpolation, facetset)
    rmin = Inf
    rmax = -Inf
    facet_dofs = Ferrite.dirichlet_facetdof_indices(interpolation)

    for facet in facetset
        cellid, facetid = facet.idx
        cell = getcells(grid, cellid)

        for nodeid in cell.nodes[collect(facet_dofs[facetid])]
            r = Ferrite.get_node_coordinate(grid, nodeid)[1]
            rmin = min(rmin, r)
            rmax = max(rmax, r)
        end
    end

    return rmin, rmax
end


function _geometry_tolerance(grid)
    scale = zero(eltype(first(getnodes(grid)).x))

    for node in getnodes(grid)
        scale = max(scale, maximum(abs, node.x))
    end

    return max(scale, one(scale)) * sqrt(eps(float(scale)))
end

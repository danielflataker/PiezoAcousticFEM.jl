"""
    validate(problem)

Preflight a piezoelectric problem before matrix assembly. The current supported
discretization is equal-order scalar `Lagrange` or `Serendipity` interpolation
for displacement and electric potential on a 2D axisymmetric `r-z` grid.
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

    _is_supported_scalar_interpolation(interpolation) ||
        throw(ArgumentError("unsupported discretization: interpolation must be Lagrange or Serendipity"))

    ip_shape = Ferrite.getrefshape(interpolation)
    qr_shape = Ferrite.getrefshape(quadrature)
    ip_shape == qr_shape ||
        throw(ArgumentError(
            "unsupported discretization: quadrature reference shape $qr_shape must match interpolation reference shape $ip_shape",
        ))

    for cellid in 1:getncells(grid)
        cell = getcells(grid, cellid)
        cell_shape = Ferrite.getrefshape(typeof(cell))

        cell_shape == ip_shape ||
            throw(ArgumentError(
                "unsupported discretization: cell $cellid has reference shape $cell_shape but interpolation uses $ip_shape",
            ))
    end

    return nothing
end


_is_supported_scalar_interpolation(::Union{Lagrange,Serendipity}) = true
_is_supported_scalar_interpolation(_) = false


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


function _facetset_radius_range(grid, interpolation, facetset)
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


function _geometry_tolerance(grid)
    scale = zero(eltype(Ferrite.get_node_coordinate(grid, 1)))

    for nodeid in 1:getnnodes(grid)
        scale = max(scale, maximum(abs, Ferrite.get_node_coordinate(grid, nodeid)))
    end

    return max(scale, one(scale)) * sqrt(eps(float(scale)))
end

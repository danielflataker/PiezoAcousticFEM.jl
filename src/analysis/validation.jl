"""
    validate(problem)

Preflight a piezoelectric problem before matrix assembly. The current supported
discretization is nodal equal-order Lagrange interpolation for displacement and
electric potential on a 2D axisymmetric `r-z` grid.
"""
function validate(problem::PiezoProblem)
    validate_supported_discretization(problem)

    return nothing
end


function validate_supported_discretization(problem::PiezoProblem)
    grid = problem.grid
    interpolation = problem.interpolation
    quadrature = problem.quadrature

    getncells(grid) > 0 ||
        throw(ArgumentError("unsupported discretization: grid must contain at least one cell"))

    _coordinate_dimension(grid) == 2 ||
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


function _coordinate_dimension(grid)
    nodes = getnodes(grid)
    isempty(nodes) &&
        throw(ArgumentError("unsupported discretization: grid must contain at least one node"))

    return length(first(nodes).x)
end


_is_lagrange_interpolation(::Lagrange) = true
_is_lagrange_interpolation(_) = false

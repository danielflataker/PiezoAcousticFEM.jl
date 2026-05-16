"""
    FacetBoundary(facetset)

Geometric boundary described by a Ferrite facetset or facetset name.
"""
abstract type AbstractBoundary end

struct FacetBoundary{F} <: AbstractBoundary
    facetset::F
end


"""
    AxisBoundary(facetset)

Symmetry-axis boundary for axisymmetric r-z models.
"""
struct AxisBoundary{B<:FacetBoundary} <: AbstractBoundary
    boundary::B
end


"""
    AxisymmetricBoundaryConditions(mechanical)

Boundary-condition container for the current one-domain axisymmetric piezo
problem. `mechanical` is a tuple of mechanical boundary-condition objects,
such as `AxisBoundary(FacetBoundary("left"))`.
"""
struct AxisymmetricBoundaryConditions{M<:Tuple}
    mechanical::M
end

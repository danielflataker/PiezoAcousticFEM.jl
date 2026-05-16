# formulations.jl

"""
    AbstractFormulation

Supertype for field formulations, such as the axisymmetric r-z formulation or a
future full 3D/acoustic formulation.
"""
abstract type AbstractFormulation end

"""
    AxisymmetricRZ

Axisymmetric formulation in the r-z plane. A 2D mesh in coordinates
`x = [r, z]` represents a 3D body rotated around the z-axis.
"""
struct AxisymmetricRZ <: AbstractFormulation end

coordinate_measure(::AxisymmetricRZ, x) = x[1]
symmetry_measure(::AxisymmetricRZ) = 2π

"""
    integration_weight(::AxisymmetricRZ, x, dΩ)

Return the axisymmetric integration weight

    2π*r*dΩ

where `r = x[1]`. Here `dΩ` is Ferrite's area element in the r-z plane.
"""
function integration_weight(kin::AxisymmetricRZ, x, dΩ)
    r = coordinate_measure(kin, x)
    r > zero(r) || throw(ArgumentError("axisymmetric quadrature radius must be positive, got r=$r"))
    dΩ > zero(dΩ) || throw(ArgumentError("axisymmetric quadrature area measure must be positive, got dΩ=$dΩ"))

    return symmetry_measure(kin) * r * dΩ
end

"""
    strain(::AxisymmetricRZ, u, ∇u, x)

Return the Voigt strain for the axisymmetric r-z formulation:

    S = [S_rr, S_θθ, S_zz, γ_rz]

where γ_rz = 2S_rz = ∂u_r/∂z + ∂u_z/∂r.

Assumes u = [u_r, u_z], x = [r, z], and ∇u[i,j] = ∂u_i/∂x_j.
"""
function strain(::AxisymmetricRZ, u, ∇u, x)
    r = x[1]
    r > zero(r) || throw(ArgumentError("axisymmetric strain radius must be positive, got r=$r"))

    Sᵣᵣ = ∇u[1, 1]
    Sθθ = u[1] / r
    Szz = ∇u[2, 2]
    γᵣz = ∇u[1, 2] + ∇u[2, 1]

    return SVector(Sᵣᵣ, Sθθ, Szz, γᵣz)
end

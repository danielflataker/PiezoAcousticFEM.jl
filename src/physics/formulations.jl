# formulations.jl

"""
    AbstractFormulation

Supertype for feltformuleringer, for eksempel aksesymmetrisk r-z-formulering
eller senere full 3D-/akustisk formulering.
"""
abstract type AbstractFormulation end

"""
    AxisymmetricRZ

Aksesymmetrisk formulering i r-z-planet. En 2D mesh i koordinatene `x = [r, z]`
representerer et 3D legeme rotert rundt z-aksen.
"""
struct AxisymmetricRZ <: AbstractFormulation end

coordinate_measure(::AxisymmetricRZ, x) = x[1]
symmetry_measure(::AxisymmetricRZ) = 2π

"""
    integration_weight(::AxisymmetricRZ, x, dΩ)

Returnerer aksesymmetrisk integrasjonsvekt

    2π*r*dΩ

hvor `r = x[1]`. Her er `dΩ` arealelementet i r-z-planet fra Ferrite.
"""
function integration_weight(kin::AxisymmetricRZ, x, dΩ)
    r = coordinate_measure(kin, x)
    r > zero(r) || throw(ArgumentError("axisymmetric quadrature radius must be positive, got r=$r"))
    dΩ > zero(dΩ) || throw(ArgumentError("axisymmetric quadrature area measure must be positive, got dΩ=$dΩ"))

    return symmetry_measure(kin) * r * dΩ
end

"""
    strain(::AxisymmetricRZ, u, ∇u, x)

Returnerer Voigt-strain for aksesymmetrisk r-z-formulering:

    S = [S_rr, S_θθ, S_zz, γ_rz]

der γ_rz = 2S_rz = ∂u_r/∂z + ∂u_z/∂r.

Antar at u = [u_r, u_z], x = [r, z], og at ∇u[i,j] = ∂u_i/∂x_j.
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

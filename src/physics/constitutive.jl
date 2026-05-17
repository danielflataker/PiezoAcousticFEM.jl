# Constitutive helpers for the current axisymmetric piezoelectric model.

"""
    electric_field(::AxisymmetricRZ, ∇ϕ)

Return the electric field for the axisymmetric r-z formulation:

    E = [E_r, E_z] = -∇ϕ

Assumes ∇ϕ = [∂ϕ/∂r, ∂ϕ/∂z].
"""
function electric_field(::AxisymmetricRZ, ∇ϕ)
    Eᵣ = -∇ϕ[1]
    Ez = -∇ϕ[2]

    return SVector(Eᵣ, Ez)
end

"""
    stress(material::AxisymmetricRZPiezoMaterial, S, E)

Return the mechanical stress in Voigt form for a linear piezoelectric material:

    T = cᴱ*S - transpose(e)*E

Here `S = [S_rr, S_θθ, S_zz, γ_rz]` and `E = [E_r, E_z]`.
This uses `transpose(material.e)` instead of `material.e'`, so the expression
does not complex-conjugate the piezoelectric matrix if later loss models add
complex constants.
"""
function stress(material::AxisymmetricRZPiezoMaterial, S, E)
    return material.cᴱ * S - transpose(material.e) * E
end

"""
    electric_displacement(material::AxisymmetricRZPiezoMaterial, S, E)

Return the electric displacement for a linear piezoelectric material:

    D = e*S + εˢ*E

Here `S = [S_rr, S_θθ, S_zz, γ_rz]` and `E = [E_r, E_z]`.
"""
function electric_displacement(material::AxisymmetricRZPiezoMaterial, S, E)
    return material.e * S + material.εˢ * E
end

# Constitutive helpers for the current axisymmetric piezoelectric model.

"""
    electric_field(::AxisymmetricRZ, ∇ϕ)

Returnerer elektrisk felt for aksesymmetrisk r-z-formulering:

    E = [E_r, E_z] = -∇ϕ

Antar at ∇ϕ = [∂ϕ/∂r, ∂ϕ/∂z].
"""
function electric_field(::AxisymmetricRZ, ∇ϕ)
    Eᵣ = -∇ϕ[1]
    Ez = -∇ϕ[2]

    return SVector(Eᵣ, Ez)
end

"""
    stress(material::Piezo6mmAxi, S, E)

Returnerer mekanisk stress på Voigt-form for et lineært piezoelektrisk materiale:

    T = cᴱ*S - transpose(e)*E

Her er `S = [S_rr, S_θθ, S_zz, γ_rz]` og `E = [E_r, E_z]`.
Bruker `transpose(material.e)` i stedet for `material.e'`, slik at uttrykket ikke
komplekskonjugerer piezomatrisen dersom tapsmodeller med komplekse konstanter
legges til senere.
"""
function stress(material::Piezo6mmAxi, S, E)
    return material.cᴱ * S - transpose(material.e) * E
end

"""
    electric_displacement(material::Piezo6mmAxi, S, E)

Returnerer elektrisk flukstetthet for et lineært piezoelektrisk materiale:

    D = e*S + εˢ*E

Her er `S = [S_rr, S_θθ, S_zz, γ_rz]` og `E = [E_r, E_z]`.
"""
function electric_displacement(material::Piezo6mmAxi, S, E)
    return material.e * S + material.εˢ * E
end

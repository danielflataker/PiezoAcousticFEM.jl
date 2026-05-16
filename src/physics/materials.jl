# materials.jl

"""
    AbstractMaterial

Felles supertype for materialmodeller i PiezoAcousticFEM.
"""

abstract type AbstractMaterial end

"""
    AbstractPiezoMaterial <: AbstractMaterial

Supertype for lineære piezoelektriske materialmodeller.
"""
abstract type AbstractPiezoMaterial <: AbstractMaterial end

"""
    AbstractLossModel

Marker for material/loss policies selected at the problem layer.
"""
abstract type AbstractLossModel end

abstract type AbstractMaterialPolicy end
abstract type AbstractSystemDamping end

"""
    Lossless()

Material policy for the current lossless solvers. It reduces the as-given
material to the requested formulation and drops imaginary parts explicitly, so
future lossy policies can enter through the same problem-level API.
"""
struct Lossless <: AbstractLossModel end

"""
    RealMaterial()

Material policy that uses only the real part of the supplied constants.
"""
struct RealMaterial <: AbstractMaterialPolicy end

"""
    MaterialAsGiven()

Material policy that preserves the supplied constants after formulation
reduction. This is the entry point for complex material constants before any
extra damping model is added.
"""
struct MaterialAsGiven <: AbstractMaterialPolicy end

"""
    NoSystemDamping()

System-damping policy with no added damping matrix.
"""
struct NoSystemDamping <: AbstractSystemDamping end

"""
    PhysicalLoss(material_policy, system_damping)

Problem-level loss description split into material handling and added system
damping. This keeps future losses out of solver keywords.
"""
struct PhysicalLoss{MP<:AbstractMaterialPolicy,SD<:AbstractSystemDamping} <: AbstractLossModel
    material_policy::MP
    system_damping::SD
end

"""
    Piezo6mmConstants

Full 6mm piezoelektrisk materialmodell i Voigt-form.
Feltene kan ha uavhengige elementtyper, slik at elastisitet, piezokobling,
permittivitet og massetetthet ikke kunstig konverteres til samme scalar-type.
"""
struct Piezo6mmConstants{TC,TE,Tε,Tρ} <: AbstractPiezoMaterial
    cᴱ::SMatrix{6,6,TC,36}
    e::SMatrix{3,6,TE,18}
    εˢ::SMatrix{3,3,Tε,9}
    ρ::Tρ
end

"""
    Piezo6mmAxi

Piezoelektrisk 6mm-materiale eksplisitt redusert til aksesymmetrisk
r-z-formulering.

Feltene har dimensjoner

    cᴱ :: 4×4
    e  :: 2×4
    εˢ :: 2×2

og brukes sammen med

    S = [S_rr, S_θθ, S_zz, γ_rz]
    E = [E_r, E_z]

slik at

    T = cᴱ*S - transpose(e)*E
    D = e*S + εˢ*E
"""
struct Piezo6mmAxi{TC,TE,Tε,Tρ} <: AbstractPiezoMaterial
    cᴱ::SMatrix{4,4,TC,16}
    e::SMatrix{2,4,TE,8}
    εˢ::SMatrix{2,2,Tε,4}
    ρ::Tρ
end

"""
    PZT5A(; TC=Float64, TE=TC, Tε=TC, Tρ=TC)

Returnerer full 6mm PZT-5A-materialkonstanter. Typeparametrene er uavhengige:
`TC` for elastisitet, `TE` for piezokobling, `Tε` for permittivitet og `Tρ`
for massetetthet.

Konstantene er i SI-enheter og følger PZT-5A-tallene brukt i Kocbach/FEMP:
elastisitet i Pa, piezokobling i C/m^2, permittivitet i F/m og tetthet i
kg/m^3.
"""
function PZT5A(; T=Float64, TC=T, TE=TC, Tε=TC, Tρ=TC)
    c11 = TC(12.1e10)
    c12 = TC(7.54e10)
    c13 = TC(7.52e10)
    c33 = TC(11.1e10)
    c44 = TC(2.26e10)
    c66 = (c11 - c12) / TC(2)

    e31 = TE(-5.4)
    e33 = TE(15.8)
    e15 = TE(12.3)

    ε0 = Tε(8.8541878128e-12)
    ε11 = Tε(916) * ε0
    ε33 = Tε(830) * ε0

    ρ = Tρ(7750.0)

    cᴱ = @SMatrix [
        c11 c12 c13 TC(0) TC(0) TC(0)
        c12 c11 c13 TC(0) TC(0) TC(0)
        c13 c13 c33 TC(0) TC(0) TC(0)
        TC(0) TC(0) TC(0) c44 TC(0) TC(0)
        TC(0) TC(0) TC(0) TC(0) c44 TC(0)
        TC(0) TC(0) TC(0) TC(0) TC(0) c66
    ]

    e = @SMatrix [
        TE(0) TE(0) TE(0) TE(0) e15 TE(0)
        TE(0) TE(0) TE(0) e15 TE(0) TE(0)
        e31 e31 e33 TE(0) TE(0) TE(0)
    ]

    εˢ = @SMatrix [
        ε11 Tε(0) Tε(0)
        Tε(0) ε11 Tε(0)
        Tε(0) Tε(0) ε33
    ]

    return Piezo6mmConstants(cᴱ, e, εˢ, ρ)
end

"""
    reduce_material(material, formulation)

Eksplisitt reduksjon fra full materialmodell til formuleringens feltbasis.
"""
reduce_material(material::Piezo6mmAxi, ::AxisymmetricRZ) = material

function reduce_material(material::Piezo6mmConstants, ::AxisymmetricRZ)
    c = material.cᴱ
    e = material.e
    ε = material.εˢ

    c_axi = @SMatrix [
        c[1, 1] c[1, 2] c[1, 3] c[1, 5]
        c[2, 1] c[2, 2] c[2, 3] c[2, 5]
        c[3, 1] c[3, 2] c[3, 3] c[3, 5]
        c[5, 1] c[5, 2] c[5, 3] c[5, 5]
    ]

    e_axi = @SMatrix [
        e[1, 1] e[1, 2] e[1, 3] e[1, 5]
        e[3, 1] e[3, 2] e[3, 3] e[3, 5]
    ]

    ε_axi = @SMatrix [
        ε[1, 1] ε[1, 3]
        ε[3, 1] ε[3, 3]
    ]

    return Piezo6mmAxi(c_axi, e_axi, ε_axi, material.ρ)
end


"""
    effective_material(material, formulation, loss)

Return the formulation-specific material actually used during assembly.
"""
effective_material(material::AbstractPiezoMaterial, formulation::AbstractFormulation, loss::Lossless) =
    lossless_material(reduce_material(material, formulation))

effective_material(
    material::AbstractPiezoMaterial,
    formulation::AbstractFormulation,
    loss::PhysicalLoss,
) = apply_material_policy(reduce_material(material, formulation), loss.material_policy)


lossless_material(material::Piezo6mmConstants) =
    Piezo6mmConstants(real.(material.cᴱ), real.(material.e), real.(material.εˢ), real(material.ρ))

lossless_material(material::Piezo6mmAxi) =
    Piezo6mmAxi(real.(material.cᴱ), real.(material.e), real.(material.εˢ), real(material.ρ))

apply_material_policy(material::AbstractPiezoMaterial, ::RealMaterial) = lossless_material(material)
apply_material_policy(material::AbstractPiezoMaterial, ::MaterialAsGiven) = material

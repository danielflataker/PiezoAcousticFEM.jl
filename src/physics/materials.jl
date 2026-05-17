# materials.jl

"""
    AbstractMaterial

Common supertype for material models in PiezoAcousticFEM.
"""

abstract type AbstractMaterial end

"""
    AbstractPiezoMaterial <: AbstractMaterial

Supertype for linear piezoelectric material models in stress-charge form.
"""
abstract type AbstractPiezoMaterial <: AbstractMaterial end

"""
    AbstractReducedPiezoMaterial <: AbstractPiezoMaterial

Supertype for formulation-specific piezoelectric materials after reducing a
full Voigt source material to the field basis used by an FE formulation.
"""
abstract type AbstractReducedPiezoMaterial <: AbstractPiezoMaterial end

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
    PiezoComplexMaterialLoss(; Qm, tan_delta, Qe)

Simple material-loss policy for direct harmonic smoke tests. With the current
`exp(im*ω*t)` convention it scales stiffness by `1 + im / Qm`, piezoelectric
coupling by `1 + im / Qe`, and permittivity by `1 - im * tan_delta`.
"""
struct PiezoComplexMaterialLoss{QM,TD,QE} <: AbstractMaterialPolicy
    Qm::QM
    tan_delta::TD
    Qe::QE
end

function PiezoComplexMaterialLoss(; Qm, tan_delta, Qe)
    Qm > zero(Qm) || throw(ArgumentError("Qm must be positive"))
    tan_delta >= zero(tan_delta) || throw(ArgumentError("tan_delta must be nonnegative"))
    Qe > zero(Qe) || throw(ArgumentError("Qe must be positive"))

    return PiezoComplexMaterialLoss(Qm, tan_delta, Qe)
end

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
    PiezoMaterial

Full piezoelectric material in Voigt matrix form. This is the source-material
storage type used by built-in material constructors. It stores

    cᴱ :: 6×6
    e  :: 3×6
    εˢ :: 3×3

for the stress-charge constitutive equations

    T = cᴱ*S - transpose(e)*E
    D = e*S + εˢ*E

with Voigt ordering `[11, 22, 33, 23, 13, 12]` and engineering shear strains.
The type itself does not enforce symmetry or a 6mm material class; constructors
and validation policy are responsible for documenting such assumptions.
"""
struct PiezoMaterial{TC,TE,Tε,Tρ} <: AbstractPiezoMaterial
    cᴱ::SMatrix{6,6,TC,36}
    e::SMatrix{3,6,TE,18}
    εˢ::SMatrix{3,3,Tε,9}
    ρ::Tρ
end

"""
    AxisymmetricRZPiezoMaterial

Piezoelectric material reduced to the `AxisymmetricRZ` formulation.

The fields have dimensions

    cᴱ :: 4×4
    e  :: 2×4
    εˢ :: 2×2

and are used with

    S = [S_rr, S_θθ, S_zz, γ_rz]
    E = [E_r, E_z]

so that

    T = cᴱ*S - transpose(e)*E
    D = e*S + εˢ*E
"""
struct AxisymmetricRZPiezoMaterial{TC,TE,Tε,Tρ} <: AbstractReducedPiezoMaterial
    cᴱ::SMatrix{4,4,TC,16}
    e::SMatrix{2,4,TE,8}
    εˢ::SMatrix{2,2,Tε,4}
    ρ::Tρ
end

"""
    PZT5A(; TC=Float64, TE=TC, Tε=TC, Tρ=TC)

Return full-Voigt PZT-5A material constants. The type parameters are independent:
`TC` for stiffness, `TE` for piezoelectric coupling, `Tε` for permittivity, and
`Tρ` for mass density.

The constants use SI units and follow FEMP material `40`: stiffness in Pa,
piezoelectric coupling in C/m^2, permittivity in F/m, and density in kg/m^3.
"""
function PZT5A(; T=Float64, TC=T, TE=TC, Tε=TC, Tρ=TC)
    c11 = TC(12.1e10)
    c12 = TC(7.54e10)
    c13 = TC(7.52e10)
    c33 = TC(11.1e10)
    c44 = TC(2.11e10)
    c55 = TC(2.11e10)
    c66 = TC(2.26e10)

    e31 = TE(-5.4)
    e33 = TE(15.8)
    e15 = TE(12.3)

    ε11 = Tε(8.11026e-9)
    ε33 = Tε(7.34882e-9)

    ρ = Tρ(7750.0)

    cᴱ = @SMatrix [
        c11 c12 c13 TC(0) TC(0) TC(0)
        c12 c11 c13 TC(0) TC(0) TC(0)
        c13 c13 c33 TC(0) TC(0) TC(0)
        TC(0) TC(0) TC(0) c44 TC(0) TC(0)
        TC(0) TC(0) TC(0) TC(0) c55 TC(0)
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

    return PiezoMaterial(cᴱ, e, εˢ, ρ)
end

"""
    reduce_material(material, formulation)

Explicit reduction from a full material model to the formulation field basis.
For `AxisymmetricRZ`, the reduced strain basis is
`[S_rr, S_θθ, S_zz, γ_rz]`, corresponding to full Voigt indices
`[1, 2, 3, 5]` in `[11, 22, 33, 23, 13, 12]` ordering.
"""
reduce_material(material::AxisymmetricRZPiezoMaterial, ::AxisymmetricRZ) = material

function reduce_material(material::PiezoMaterial, ::AxisymmetricRZ)
    c = material.cᴱ
    e = material.e
    ε = material.εˢ
    basis = (1, 2, 3, 5)

    c_axi = @SMatrix [
        c[basis[1], basis[1]] c[basis[1], basis[2]] c[basis[1], basis[3]] c[basis[1], basis[4]]
        c[basis[2], basis[1]] c[basis[2], basis[2]] c[basis[2], basis[3]] c[basis[2], basis[4]]
        c[basis[3], basis[1]] c[basis[3], basis[2]] c[basis[3], basis[3]] c[basis[3], basis[4]]
        c[basis[4], basis[1]] c[basis[4], basis[2]] c[basis[4], basis[3]] c[basis[4], basis[4]]
    ]

    e_axi = @SMatrix [
        e[1, basis[1]] e[1, basis[2]] e[1, basis[3]] e[1, basis[4]]
        e[3, basis[1]] e[3, basis[2]] e[3, basis[3]] e[3, basis[4]]
    ]

    ε_axi = @SMatrix [
        ε[1, 1] ε[1, 3]
        ε[3, 1] ε[3, 3]
    ]

    return AxisymmetricRZPiezoMaterial(c_axi, e_axi, ε_axi, material.ρ)
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


lossless_material(material::PiezoMaterial) =
    PiezoMaterial(real.(material.cᴱ), real.(material.e), real.(material.εˢ), real(material.ρ))

lossless_material(material::AxisymmetricRZPiezoMaterial) =
    AxisymmetricRZPiezoMaterial(real.(material.cᴱ), real.(material.e), real.(material.εˢ), real(material.ρ))

apply_material_policy(material::AbstractPiezoMaterial, ::RealMaterial) = lossless_material(material)
apply_material_policy(material::AbstractPiezoMaterial, ::MaterialAsGiven) = material

function apply_material_policy(material::PiezoMaterial, loss::PiezoComplexMaterialLoss)
    return PiezoMaterial(
        material.cᴱ .* complex(one(loss.Qm), inv(loss.Qm)),
        material.e .* complex(one(loss.Qe), inv(loss.Qe)),
        material.εˢ .* complex(one(loss.tan_delta), -loss.tan_delta),
        material.ρ,
    )
end

function apply_material_policy(material::AxisymmetricRZPiezoMaterial, loss::PiezoComplexMaterialLoss)
    return AxisymmetricRZPiezoMaterial(
        material.cᴱ .* complex(one(loss.Qm), inv(loss.Qm)),
        material.e .* complex(one(loss.Qe), inv(loss.Qe)),
        material.εˢ .* complex(one(loss.tan_delta), -loss.tan_delta),
        material.ρ,
    )
end

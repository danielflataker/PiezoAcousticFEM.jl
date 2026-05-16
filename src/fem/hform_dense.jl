"""
    HFormSystem

H-form-system etter at indre potensial-DOF-er er kondensert ut og
den ikke-jordede elektroden er representert med ett skalarpotensial.

```text
-ω² [Muu 0] [u] + [Huu  Huϕ] [u] = [F]
    [0   0] [ϕ]   [Hϕu  Hϕϕ] [ϕ]   [-Q]
```
"""
struct HFormSystem{HUU,HUP,HPU,HPP,MUU}
    Huu::HUU
    Huϕ::HUP
    Hϕu::HPU
    Hϕϕ::HPP
    Muu::MUU
end


"""
    h_form_dense(system, partition)

Transformer globale tette K-form-blokker til Kocbachs H-form for små
verifikasjonsproblemer. Dette følger KLV 1999 Eq. (2.192)-(2.203) /
Kocbach 2000 Eq. (3.190)-(3.192):

- jordede potensial-DOF-er er Dirichlet-DOF-er og brukes ikke i reduksjonen,
- indre potensial-DOF-er kondenseres ut med et Schur-komplement,
- alle noder på den drivne elektroden bindes til ett skalarpotensial via
  `Ip = [1, 1, ..., 1]ᵀ`.

Implementasjonen bruker lineære solves (`Kii \\ ...`) i stedet for eksplisitt
invers, men uttrykkene er algebraisk de samme som hos Kocbach. Den materialiserer
tette matriser og skal derfor ikke brukes som produksjonsvei for store nett.
"""
function h_form_dense(system::KFormSystem, partition::HarmonicVoltageDofPartition)
    validate_partition(system, partition)

    i = partition.internal
    p = partition.driven
    Ip_uϕ = ones(eltype(system.Kuϕ), length(p))
    Ip_ϕu = ones(eltype(system.Kϕu), length(p))
    Ip_ϕϕ = ones(eltype(system.Kϕϕ), length(p))

    Kuu = system.Kuu
    Muu = system.Muu
    Kup = system.Kuϕ[:, p]
    Kpu = system.Kϕu[p, :]
    Kpp = system.Kϕϕ[p, p]

    if isempty(i)
        Huu = Matrix(Kuu)
        Huϕ = Kup * Ip_uϕ
        Hϕu = transpose(Ip_ϕu) * Kpu
        Hϕϕ = sum(Kpp * Ip_ϕϕ)
        return HFormSystem(Huu, Huϕ, Matrix(Hϕu), Hϕϕ, Matrix(Muu))
    end

    Kui = Matrix(system.Kuϕ[:, i])
    Kiu = Matrix(system.Kϕu[i, :])
    Kii = Matrix(system.Kϕϕ[i, i])
    Kip = Matrix(system.Kϕϕ[i, p])
    Kpi = Matrix(system.Kϕϕ[p, i])

    Kii_fact = factorize(Kii)
    Huu = Kuu - Kui * (Kii_fact \ Kiu)
    Huϕ = (Kup - Kui * (Kii_fact \ Kip)) * Ip_uϕ
    Hϕu = transpose(Ip_ϕu) * (Kpu - Kpi * (Kii_fact \ Kiu))
    Hϕϕ = sum((Kpp - Kpi * (Kii_fact \ Kip)) * Ip_ϕϕ)

    return HFormSystem(Matrix(Huu), Huϕ, Matrix(Hϕu), Hϕϕ, Matrix(Muu))
end


"""
    short_circuit_h_form_dense(system, partition)

Dense short-circuit mechanical operator for small modal reference problems.
All electrode potential DOFs are homogeneous Dirichlet DOFs, and internal
potential DOFs are condensed out.
"""
function short_circuit_h_form_dense(system::KFormSystem, partition::ShortCircuitDofPartition)
    validate_partition(system, partition)

    i = partition.internal
    Kuu = system.Kuu
    Muu = system.Muu

    if isempty(i)
        return (Kuu=Matrix(Kuu), Muu=Matrix(Muu))
    end

    Kui = Matrix(system.Kuϕ[:, i])
    Kiu = Matrix(system.Kϕu[i, :])
    Kii = Matrix(system.Kϕϕ[i, i])
    Kii_fact = factorize(Kii)
    Huu = Kuu - Kui * (Kii_fact \ Kiu)

    return (Kuu=Matrix(Huu), Muu=Matrix(Muu))
end

"""
    ElectrodeReducedKForm

K-form der jordede potensial-DOF-er er fjernet og alle DOF-er på den drivne
elektroden er kollapset til ett skalarpotensial `V`. Ukjente i den frie
delen er `u` og de interne potensialene `ϕᵢ`; elektrodekolonner og -rader er
lagret separat slik at en direkte spenningssolve kan sette `V = V0`.
"""
struct ElectrodeReducedKForm{
    KUU<:AbstractMatrix,
    KUI<:AbstractMatrix,
    KIU<:AbstractMatrix,
    KII<:AbstractMatrix,
    MUU<:AbstractMatrix,
    KUP<:AbstractVector,
    KIP<:AbstractVector,
    KPU<:AbstractVector,
    KPI<:AbstractVector,
    KPP,
}
    Kuu::KUU
    Kui::KUI
    Kiu::KIU
    Kii::KII
    Muu::MUU
    KuP::KUP
    KiP::KIP
    KPu::KPU
    KPi::KPI
    KPP::KPP
    partition::HarmonicVoltageDofPartition
end


"""
    electrode_reduced_k_form(system, partition)

Collapse the driven electrode in a `KFormSystem` without condensing internal
electrical DOFs. Matrix blocks may be dense reference matrices or sparse
production matrices.
"""
function electrode_reduced_k_form(system::KFormSystem, partition::HarmonicVoltageDofPartition)
    validate_partition(system, partition)

    i = partition.internal
    p = partition.driven
    Ip_uϕ = ones(eltype(system.Kuϕ), length(p))
    Ip_ϕu = ones(eltype(system.Kϕu), length(p))
    Ip_ϕϕ = ones(eltype(system.Kϕϕ), length(p))

    Kuu = system.Kuu
    Muu = system.Muu
    Kui = system.Kuϕ[:, i]
    Kiu = system.Kϕu[i, :]
    Kii = system.Kϕϕ[i, i]

    KuP = system.Kuϕ[:, p] * Ip_uϕ
    KiP = system.Kϕϕ[i, p] * Ip_ϕϕ
    KPu = vec(transpose(Ip_ϕu) * system.Kϕu[p, :])
    KPi = vec(transpose(Ip_ϕϕ) * system.Kϕϕ[p, i])
    KPP = sum(system.Kϕϕ[p, p] * Ip_ϕϕ)

    return ElectrodeReducedKForm(Kuu, Kui, Kiu, Kii, Muu, KuP, KiP, KPu, KPi, KPP, partition)
end

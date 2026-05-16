"""
    ElectrodeReducedKForm

K-form where grounded potential DOFs have been removed and all DOFs on the
driven electrode have been collapsed to one scalar potential `V`. The unknowns
in the free part are `u` and the internal potentials `ϕᵢ`; electrode columns
and rows are stored separately so a direct voltage solve can set `V = V0`.

`KuP` and `KiP` are the prescribed-voltage columns that move to the active
solve right-hand side. `KPu`, `KPi`, and `KPP` are the removed driven-electrode
row after electrode tying; after the active solve, that row is evaluated as the
reaction equation used to recover terminal charge.
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
production matrices. The returned reduction keeps active solve blocks and the
driven-electrode reaction row separate.
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

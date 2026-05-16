"""
    HarmonicVoltageDofPartition(nϕ, internal, driven, grounded)

Partition of the global electric potential DOFs for harmonic voltage drive:

- `internal`: nodes not on electrodes (`ϕᵢ` in Kocbach's notation).
- `driven`: nodes on the non-grounded electrode (`ϕₚ = Ip * ϕ`).
- `grounded`: nodes on the reference electrode, set to zero and therefore
  removed.

The indices are 1-based Julia indices into the global `ϕ` vector.
"""
struct HarmonicVoltageDofPartition
    nϕ::Int
    internal::Vector{Int}
    driven::Vector{Int}
    grounded::Vector{Int}

    function HarmonicVoltageDofPartition(nϕ::Integer, internal, driven, grounded)
        nϕ >= 0 || throw(ArgumentError("nϕ must be nonnegative"))
        i = collect(Int, internal)
        p = collect(Int, driven)
        g = collect(Int, grounded)
        isempty(p) && throw(ArgumentError("partition.driven must contain at least one electrode DOF"))
        validate_partition_coverage(Int(nϕ), i, p, g)

        return new(Int(nϕ), i, p, g)
    end
end

"""
    ShortCircuitDofPartition(nϕ, internal, grounded)

Partition of potential DOFs for short-circuit modal analysis. The electrodes
are homogeneous electric Dirichlet DOFs; only internal potential DOFs are
condensed.
"""
struct ShortCircuitDofPartition
    nϕ::Int
    internal::Vector{Int}
    grounded::Vector{Int}

    function ShortCircuitDofPartition(nϕ::Integer, internal, grounded)
        nϕ >= 0 || throw(ArgumentError("nϕ must be nonnegative"))
        i = collect(Int, internal)
        g = collect(Int, grounded)
        validate_partition_coverage(Int(nϕ), i, g)

        return new(Int(nϕ), i, g)
    end
end

"""
    OpenCircuitDofPartition(nϕ, internal, floating, reference)

Placeholder partition shape for later open-circuit reductions. It records the
distinct electrical roles without pretending they are voltage-driven DOFs.
"""
struct OpenCircuitDofPartition
    nϕ::Int
    internal::Vector{Int}
    floating::Vector{Int}
    reference::Vector{Int}

    function OpenCircuitDofPartition(nϕ::Integer, internal, floating, reference)
        nϕ >= 0 || throw(ArgumentError("nϕ must be nonnegative"))
        i = collect(Int, internal)
        f = collect(Int, floating)
        r = collect(Int, reference)
        isempty(f) && throw(ArgumentError("partition.floating must contain at least one electrode DOF"))
        isempty(r) && throw(ArgumentError("partition.reference must contain at least one electrode DOF"))
        validate_partition_coverage(Int(nϕ), i, f, r)

        return new(Int(nϕ), i, f, r)
    end
end

function validate_partition_coverage(nϕ::Int, parts...)
    all_phi = reduce(vcat, parts; init=Int[])
    length(all_phi) == nϕ ||
        throw(ArgumentError("potential DOF partition must cover every potential DOF exactly once"))
    all(1 .<= all_phi .<= nϕ) || throw(ArgumentError("potential DOF index outside 1:$nϕ"))
    length(unique(all_phi)) == length(all_phi) ||
        throw(ArgumentError("potential DOF partition contains duplicate indices"))
    sort(all_phi) == collect(1:nϕ) ||
        throw(ArgumentError("potential DOF partition must cover every potential DOF exactly once"))

    return nothing
end


function validate_system_potential_size(system::KFormSystem, expected_nϕ::Int)
    system_nϕ = size(system.Kϕϕ, 1)
    size(system.Kϕϕ, 2) == system_nϕ || throw(DimensionMismatch("Kϕϕ must be square"))
    size(system.Kuϕ, 2) == system_nϕ || throw(DimensionMismatch("Kuϕ columns must match Kϕϕ"))
    size(system.Kϕu, 1) == system_nϕ || throw(DimensionMismatch("Kϕu rows must match Kϕϕ"))
    size(system.Kuu, 1) == size(system.Kuu, 2) || throw(DimensionMismatch("Kuu must be square"))
    size(system.Muu) == size(system.Kuu) || throw(DimensionMismatch("Muu must match Kuu"))
    size(system.Kuϕ, 1) == size(system.Kuu, 1) || throw(DimensionMismatch("Kuϕ rows must match Kuu"))
    size(system.Kϕu, 2) == size(system.Kuu, 2) || throw(DimensionMismatch("Kϕu columns must match Kuu"))

    expected_nϕ == system_nϕ ||
        throw(ArgumentError("partition.nϕ=$expected_nϕ does not match system potential size $system_nϕ"))

    return nothing
end


function validate_partition(system::KFormSystem, partition::HarmonicVoltageDofPartition)
    nϕ = size(system.Kϕϕ, 1)
    validate_system_potential_size(system, partition.nϕ)
    isempty(partition.driven) &&
        throw(ArgumentError("partition.driven must contain at least one electrode DOF"))
    validate_partition_coverage(nϕ, partition.internal, partition.driven, partition.grounded)

    return nothing
end


function validate_partition(system::KFormSystem, partition::ShortCircuitDofPartition)
    nϕ = size(system.Kϕϕ, 1)
    validate_system_potential_size(system, partition.nϕ)
    validate_partition_coverage(nϕ, partition.internal, partition.grounded)

    return nothing
end

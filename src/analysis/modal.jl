"""
    ShortCircuitModalAnalysis(nev)

Lossless short-circuit modal reference analysis. This first implementation is
dense and intended for small verification meshes. Use `nothing` to keep all
modes.
"""
struct ShortCircuitModalAnalysis{N}
    nev::N
end


"""
    ShortCircuitModalReduction

Dense short-circuit H-form plus homogeneous mechanical constraints.
"""
struct ShortCircuitModalReduction{A,P,S,C,F}
    assembled::A
    partition::P
    system::S
    mechanical_dirichlet::C
    free_dofs::F
end


"""
    ShortCircuitModalResult

Eigenpairs and metadata for `solve(problem, ::ShortCircuitModalAnalysis)`.
Modes are stored as full compact displacement vectors. `normalization` records
the convention used for the columns of `modes`.
"""
struct ShortCircuitModalResult{P,A,AP,R,Λ,W,F,M,N}
    problem::P
    analysis::A
    assembled::AP
    reduction::R
    eigenvalues::Λ
    angular_frequencies::W
    frequencies::F
    modes::M
    normalization::N
end


function reduce(
    assembled::AssembledPiezoProblem,
    ::ShortCircuitModalAnalysis,
)
    problem = assembled.problem
    partition = short_circuit_partition(assembled.assembly, problem.electrodes)
    system = short_circuit_h_form_dense(assembled.assembly.system, partition)
    constraints = mechanical_dirichlet(assembled.assembly, problem.boundary_conditions)
    constrained = apply_modal_constraints(system.Kuu, system.Muu, constraints)

    return ShortCircuitModalReduction(
        assembled,
        partition,
        system,
        constraints,
        constrained.free,
    )
end


function solve(problem::PiezoProblem, analysis::ShortCircuitModalAnalysis)
    assembled = assemble(problem)
    reduction = reduce(assembled, analysis)
    constrained = apply_modal_constraints(
        reduction.system.Kuu,
        reduction.system.Muu,
        reduction.mechanical_dirichlet,
    )
    eig = eigen(Symmetric(constrained.K), Symmetric(constrained.M))
    order = sortperm(real.(eig.values))
    selected = select_modes(order, analysis.nev)
    λ = real.(eig.values[selected])
    free_modes = eig.vectors[:, selected]
    normalize_modes!(free_modes, constrained.M)

    modes = zeros(eltype(free_modes), size(reduction.system.Kuu, 1), size(free_modes, 2))
    modes[constrained.free, :] = free_modes
    ω = sqrt.(max.(λ, zero(eltype(λ))))
    f = ω ./ (2π)

    return ShortCircuitModalResult(problem, analysis, assembled, reduction, λ, ω, f, modes, :mass)
end


function apply_modal_constraints(K, M, mechanical_dirichlet)
    size(K, 1) == size(K, 2) || throw(DimensionMismatch("modal stiffness must be square"))
    size(M) == size(K) || throw(DimensionMismatch("modal mass must match stiffness"))
    n = size(K, 1)

    if mechanical_dirichlet === nothing
        free = collect(1:n)
        return (K=Matrix(K), M=Matrix(M), free=free)
    end

    all(iszero, mechanical_dirichlet.values) ||
        throw(ArgumentError("modal mechanical constraints must be homogeneous"))
    all(1 .<= mechanical_dirichlet.indices .<= n) ||
        throw(ArgumentError("modal mechanical constraint index outside 1:$n"))

    free = setdiff(collect(1:n), mechanical_dirichlet.indices)

    return (K=Matrix(K[free, free]), M=Matrix(M[free, free]), free=free)
end


function select_modes(order, nev)
    nev === nothing && return order
    0 <= nev <= length(order) || throw(ArgumentError("nev must be between 0 and $(length(order))"))

    return order[1:nev]
end


function normalize_modes!(modes, M)
    for j in axes(modes, 2)
        mode = view(modes, :, j)
        mnorm = sqrt(real(dot(mode, M * mode)))
        iszero(mnorm) && throw(ArgumentError("cannot mass-normalize a zero modal vector"))
        mode ./= mnorm
    end

    return modes
end

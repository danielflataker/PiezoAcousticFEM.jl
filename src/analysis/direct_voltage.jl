"""
    AbstractLinearSolverConfig

Configuration object for the linear solve used by voltage analyses.
"""
abstract type AbstractLinearSolverConfig end

"""
    BackslashSolver()

Use Julia's direct `A \\ b` solve for each linear system.
"""
struct BackslashSolver <: AbstractLinearSolverConfig end

"""
    FactorizedDirectSolver()

Factorize each reduced linear system with `factorize(A)` before solving. This
is currently a per-solve factorization hook; reusable caches are future work.
"""
struct FactorizedDirectSolver <: AbstractLinearSolverConfig end


"""
    DirectVoltageSolverInfo

Metadata for the linear solve part of a direct voltage analysis.
"""
struct DirectVoltageSolverInfo{N,R}
    method::Symbol
    matrix_size::Tuple{Int,Int}
    rhs_norm::N
    reduced_residual_norm::R
    reduced_relative_residual::R
end


"""
    DirectVoltageSolution

Result from a direct voltage solve on an `ElectrodeReducedKForm`.
`drive_terminal_charge` is the total charge on the drive terminal,
`drive_terminal_current` is `im * ω * drive_terminal_charge`, and
`drive_terminal_admittance` is `drive_terminal_current / voltage`.
"""
struct DirectVoltageSolution{
    U<:AbstractVector,
    PI<:AbstractVector,
    P<:AbstractVector,
    V,
    Q,
    I,
    Y,
    W,
    A,
    C,
    R<:AbstractVector,
    S,
}
    displacement::U
    internal_potential::PI
    potential::P
    voltage::V
    drive_terminal_charge::Q
    drive_terminal_current::I
    drive_terminal_admittance::Y
    ω::W
    analysis::A
    convention::C
    reduced_residual::R
    solver_info::S
end


"""
    solve_direct_voltage(reduced, ω, voltage; mechanical_dirichlet=nothing)

Solve the electrode-collapsed K-form with prescribed drive-terminal voltage.
The free unknowns are mechanical displacement DOFs and internal electrical
potential DOFs. `mechanical_dirichlet` may be a `DirichletDofs` object whose
indices refer to the compact displacement vector.
"""
function solve_direct_voltage(
    reduced::ElectrodeReducedKForm,
    ω,
    voltage;
    mechanical_dirichlet=nothing,
    analysis=nothing,
    convention=:exp_iomega_t,
    solver::AbstractLinearSolverConfig=BackslashSolver(),
)
    iszero(voltage) && throw(ArgumentError("voltage must be nonzero when computing admittance"))

    nᵤ = size(reduced.Kuu, 1)
    nᵢ = length(reduced.KiP)
    A, rhs = build_direct_voltage_system(reduced, ω, voltage)
    constrained = apply_direct_voltage_constraints(A, rhs, mechanical_dirichlet, nᵤ)
    x, reduced_residual, solver_info = linear_solve(constrained, solver)

    u = x[1:nᵤ]
    ϕᵢ = nᵢ == 0 ? Vector{eltype(x)}(undef, 0) : x[(nᵤ + 1):(nᵤ + nᵢ)]
    observables = direct_voltage_observables(reduced, u, ϕᵢ, voltage, ω)

    return DirectVoltageSolution(
        u,
        ϕᵢ,
        observables.potential,
        voltage,
        observables.drive_terminal_charge,
        observables.drive_terminal_current,
        observables.drive_terminal_admittance,
        ω,
        analysis,
        convention,
        reduced_residual,
        solver_info,
    )
end


function build_direct_voltage_system(reduced::ElectrodeReducedKForm, ω, voltage)
    Kdyn = reduced.Kuu - ω^2 * reduced.Muu

    if isempty(reduced.KiP)
        return Kdyn, -reduced.KuP * voltage
    end

    A = [
        Kdyn reduced.Kui
        reduced.Kiu reduced.Kii
    ]
    rhs = -vcat(reduced.KuP, reduced.KiP) * voltage

    return A, rhs
end


function apply_direct_voltage_constraints(A, rhs, mechanical_dirichlet, nᵤ)
    if mechanical_dirichlet === nothing
        return (A=A, rhs=rhs, reduction=nothing)
    end

    validate_mechanical_dirichlet(mechanical_dirichlet, nᵤ)

    return (A=A, rhs=rhs, reduction=apply_dirichlet(A, rhs, mechanical_dirichlet))
end


function linear_solve(system, solver::AbstractLinearSolverConfig)
    if system.reduction === nothing
        x = solve_linear_system(solver, system.A, system.rhs)
        residual = system.A * x - system.rhs
        solver_info = direct_voltage_solver_info(solver, system.A, system.rhs, residual)
        warn_large_direct_voltage_residual(solver_info)

        return x, residual, solver_info
    end

    reduction = system.reduction
    x_free = solve_linear_system(solver, reduction.A, reduction.b)
    residual = reduction.A * x_free - reduction.b
    solver_info = direct_voltage_solver_info(solver, reduction.A, reduction.b, residual)
    warn_large_direct_voltage_residual(solver_info)
    x = reconstruct_solution(reduction, x_free)

    return x, residual, solver_info
end


solve_linear_system(::BackslashSolver, A, rhs) = A \ rhs

function solve_linear_system(::FactorizedDirectSolver, A, rhs)
    F = factorize(A)

    return F \ rhs
end


function direct_voltage_observables(reduced::ElectrodeReducedKForm, u, ϕᵢ, voltage, ω)
    ϕ = reconstruct_potential(reduced, ϕᵢ, voltage)
    # Evaluate the removed drive-terminal row as a reaction equation. With
    # the Kocbach/KLV sign convention this reaction is minus the terminal charge.
    drive_terminal_charge = -(sum(reduced.KPu .* u) + sum(reduced.KPi .* ϕᵢ) + reduced.KPP * voltage)
    drive_terminal_current = im * ω * drive_terminal_charge
    drive_terminal_admittance = drive_terminal_current / voltage

    return (
        potential=ϕ,
        drive_terminal_charge=drive_terminal_charge,
        drive_terminal_current=drive_terminal_current,
        drive_terminal_admittance=drive_terminal_admittance,
    )
end


function validate_mechanical_dirichlet(dirichlet::DirichletDofs, nᵤ)
    all(1 .<= dirichlet.indices .<= nᵤ) ||
        throw(ArgumentError("mechanical Dirichlet indices must be displacement indices in 1:$nᵤ"))

    return nothing
end


"""
    reconstruct_potential(reduced, internal_potential, voltage)

Reconstruct the full compact electric potential vector from internal
potentials, prescribed drive-terminal voltage, and grounded electrode values.
"""
function reconstruct_potential(reduced::ElectrodeReducedKForm, internal_potential, voltage)
    partition = reduced.partition
    length(internal_potential) == length(partition.internal) ||
        throw(DimensionMismatch("internal_potential length must match partition.internal"))

    T = promote_type(eltype(internal_potential), typeof(voltage))
    ϕ = zeros(T, partition.nϕ)
    ϕ[partition.internal] = internal_potential
    ϕ[partition.driven] .= voltage
    ϕ[partition.grounded] .= zero(T)

    return ϕ
end


solver_method(::BackslashSolver) = :backslash
solver_method(::FactorizedDirectSolver) = :factorized_direct
function direct_voltage_solver_info(solver::AbstractLinearSolverConfig, A, rhs, residual)
    rhs_norm = norm(rhs)
    residual_norm = norm(residual)
    relative_residual = iszero(rhs_norm) ? residual_norm : residual_norm / rhs_norm

    return DirectVoltageSolverInfo(
        solver_method(solver),
        size(A),
        rhs_norm,
        residual_norm,
        relative_residual,
    )
end


function warn_large_direct_voltage_residual(info::DirectVoltageSolverInfo)
    T = typeof(float(real(info.reduced_relative_residual)))
    threshold = sqrt(eps(T))
    info.reduced_relative_residual <= threshold && return nothing

    method = info.method
    matrix_size = info.matrix_size
    relative_residual = info.reduced_relative_residual
    @warn "direct voltage solve returned a large reduced residual" method matrix_size relative_residual threshold

    return nothing
end

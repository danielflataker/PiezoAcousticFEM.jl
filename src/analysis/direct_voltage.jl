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
`charge` is the total charge on the driven electrode and `admittance` is
`im * ω * charge / voltage`.
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
    charge::Q
    current::I
    admittance::Y
    ω::W
    analysis::A
    convention::C
    reduced_residual::R
    solver_info::S
end


"""
    solve_direct_voltage(reduced, ω, voltage; mechanical_dirichlet=nothing)

Solve the electrode-collapsed K-form with prescribed driven-electrode voltage.
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
)
    iszero(voltage) && throw(ArgumentError("voltage must be nonzero when computing admittance"))

    nᵤ = size(reduced.Kuu, 1)
    nᵢ = length(reduced.KiP)
    A, rhs = build_direct_voltage_system(reduced, ω, voltage)
    constrained = apply_direct_voltage_constraints(A, rhs, mechanical_dirichlet, nᵤ)
    x, reduced_residual, solver_info = linear_solve(constrained)

    u = x[1:nᵤ]
    ϕᵢ = nᵢ == 0 ? Vector{eltype(x)}(undef, 0) : x[(nᵤ + 1):(nᵤ + nᵢ)]
    observables = direct_voltage_observables(reduced, u, ϕᵢ, voltage, ω)

    return DirectVoltageSolution(
        u,
        ϕᵢ,
        observables.potential,
        voltage,
        observables.charge,
        observables.current,
        observables.admittance,
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


function linear_solve(system)
    if system.reduction === nothing
        x = system.A \ system.rhs
        residual = system.A * x - system.rhs
        solver_info = direct_voltage_solver_info(:backslash, system.A, system.rhs, residual)

        return x, residual, solver_info
    end

    reduction = system.reduction
    x_free = reduction.A \ reduction.b
    residual = reduction.A * x_free - reduction.b
    solver_info = direct_voltage_solver_info(:backslash, reduction.A, reduction.b, residual)
    x = reconstruct_solution(reduction, x_free)

    return x, residual, solver_info
end


function direct_voltage_observables(reduced::ElectrodeReducedKForm, u, ϕᵢ, voltage, ω)
    ϕ = reconstruct_potential(reduced, ϕᵢ, voltage)
    charge = -(sum(reduced.KPu .* u) + sum(reduced.KPi .* ϕᵢ) + reduced.KPP * voltage)
    current = im * ω * charge
    admittance = current / voltage

    return (potential=ϕ, charge=charge, current=current, admittance=admittance)
end


function validate_mechanical_dirichlet(dirichlet::DirichletDofs, nᵤ)
    all(1 .<= dirichlet.indices .<= nᵤ) ||
        throw(ArgumentError("mechanical Dirichlet indices must be displacement indices in 1:$nᵤ"))

    return nothing
end


"""
    reconstruct_potential(reduced, internal_potential, voltage)

Reconstruct the full compact electric potential vector from internal
potentials, prescribed driven-electrode voltage, and grounded electrode values.
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


function direct_voltage_solver_info(method::Symbol, A, rhs, residual)
    rhs_norm = norm(rhs)
    residual_norm = norm(residual)
    relative_residual = iszero(rhs_norm) ? residual_norm : residual_norm / rhs_norm

    return DirectVoltageSolverInfo(method, size(A), rhs_norm, residual_norm, relative_residual)
end

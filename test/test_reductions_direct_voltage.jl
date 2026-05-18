@testset "global H-form Schur complement" begin
    Kuu = [10.0 1.0; 1.0 8.0]
    Kuϕ = [2.0 3.0 100.0; 4.0 5.0 200.0]
    Kϕu = Matrix(transpose(Kuϕ))
    Kϕϕ = [
        -6.0 -1.0 0.0
        -1.0 -7.0 0.0
        0.0 0.0 -9.0
    ]
    Muu = Matrix{Float64}(I, 2, 2)

    system = PiezoAcousticFEM.KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    partition = PiezoAcousticFEM.HarmonicVoltageDofPartition(3, [1], [2], [3])
    h = PiezoAcousticFEM.h_form_dense(system, partition)

    Kii = Kϕϕ[1:1, 1:1]
    @test h.Huu ≈ Kuu - Kuϕ[:, 1:1] * (Kii \ Kϕu[1:1, :])
    @test h.Huϕ ≈ vec(Kuϕ[:, 2:2] - Kuϕ[:, 1:1] * (Kii \ Kϕϕ[1:1, 2:2]))
    @test h.Hϕu ≈ Kϕu[2:2, :] - Kϕϕ[2:2, 1:1] * (Kii \ Kϕu[1:1, :])
    @test h.Hϕϕ ≈ only(Kϕϕ[2:2, 2:2] - Kϕϕ[2:2, 1:1] * (Kii \ Kϕϕ[1:1, 2:2]))
    @test h.Muu == Muu
end

@testset "K-form block types are independent" begin
    Kuu = ComplexF64[10.0 1.0; 1.0 8.0]
    Kuϕ = Float32[2.0 3.0 100.0; 4.0 5.0 200.0]
    Kϕu = Float32[
        2.0 4.0
        3.0 5.0
        100.0 200.0
    ]
    Kϕϕ = BigFloat[
        -6.0 -1.0 0.0
        -1.0 -7.0 0.0
        0.0 0.0 -9.0
    ]
    Muu = Float64[1.0 0.0; 0.0 1.0]

    system = PiezoAcousticFEM.KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    partition = PiezoAcousticFEM.HarmonicVoltageDofPartition(3, [1], [2], [3])
    reduced = PiezoAcousticFEM.electrode_reduced_k_form(system, partition)

    @test eltype(reduced.Kuu) == ComplexF64
    @test eltype(reduced.KuP) == Float32
    @test eltype(reduced.KPu) == Float32
    @test eltype(reduced.Kii) == BigFloat
    @test eltype(reduced.Muu) == Float64
    @test typeof(reduced.KPP) == BigFloat
end

@testset "dense electrode-reduced direct voltage solve" begin
    Kuu = [10.0 1.0; 1.0 8.0]
    Kuϕ = [2.0 3.0 100.0; 4.0 5.0 200.0]
    Kϕu = Matrix(transpose(Kuϕ))
    Kϕϕ = [
        -6.0 -1.0 0.0
        -1.0 -7.0 0.0
        0.0 0.0 -9.0
    ]
    Muu = Matrix{Float64}(I, 2, 2)
    system = PiezoAcousticFEM.KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    partition = PiezoAcousticFEM.HarmonicVoltageDofPartition(3, [1], [2], [3])

    reduced = PiezoAcousticFEM.electrode_reduced_k_form(system, partition)
    h = PiezoAcousticFEM.h_form_dense(system, partition)
    ω = 3.0
    V0 = 2.5

    direct = PiezoAcousticFEM.solve_direct_voltage(reduced, ω, V0)
    factorized_direct = PiezoAcousticFEM.solve_direct_voltage(
        reduced,
        ω,
        V0;
        solver=FactorizedDirectSolver(),
    )
    D = h.Huu - ω^2 * h.Muu
    u_h = D \ (-h.Huϕ * V0)
    removed_electrode_row_residual = only(h.Hϕu * u_h) + h.Hϕϕ * V0
    Q_h = -removed_electrode_row_residual
    Y_h = im * ω * Q_h / V0

    @test direct.displacement ≈ u_h
    @test direct.charge ≈ -removed_electrode_row_residual
    @test direct.charge ≈ Q_h
    @test direct.current ≈ im * ω * Q_h
    @test direct.admittance ≈ Y_h
    @test direct.ω == ω
    @test direct.analysis === nothing
    @test direct.convention == :exp_iomega_t
    @test direct.solver_info.method == :backslash
    @test direct.solver_info.matrix_size == (3, 3)
    @test direct.solver_info.reduced_relative_residual < 1.0e-12
    @test factorized_direct.displacement ≈ direct.displacement
    @test factorized_direct.admittance ≈ direct.admittance
    @test factorized_direct.solver_info.method == :factorized_direct
    @test direct.potential[partition.internal] ≈ direct.internal_potential
    @test direct.potential[partition.driven] == fill(V0, length(partition.driven))
    @test direct.potential[partition.grounded] == zeros(length(partition.grounded))
    @test PiezoAcousticFEM.solve_direct_voltage(reduced, ω, V0).admittance ≈ Y_h

    mechanical_dirichlet = PiezoAcousticFEM.DirichletDofs([1], [0.25])
    constrained_direct =
        PiezoAcousticFEM.solve_direct_voltage(reduced, ω, V0; mechanical_dirichlet)
    A = [
        Kuu - ω^2 * Muu Kuϕ[:, 1:1]
        Kϕu[1:1, :] Kϕϕ[1:1, 1:1]
    ]
    rhs = -vcat(Kuϕ[:, 2:2], Kϕϕ[1:1, 2:2])[:, 1] * V0
    manual_reduction = PiezoAcousticFEM.apply_dirichlet(A, rhs, mechanical_dirichlet)
    manual_x = PiezoAcousticFEM.reconstruct_solution(manual_reduction, manual_reduction.A \ manual_reduction.b)
    manual_u = manual_x[1:2]
    manual_ϕᵢ = manual_x[3:3]
    manual_reaction_residual = (
        sum(reduced.KPu .* manual_u) +
        sum(reduced.KPi .* manual_ϕᵢ) +
        reduced.KPP * V0
    )
    manual_charge = -manual_reaction_residual

    @test constrained_direct.displacement ≈ manual_u
    @test constrained_direct.internal_potential ≈ manual_ϕᵢ
    @test constrained_direct.displacement[1] == 0.25
    @test constrained_direct.charge ≈ -manual_reaction_residual
    @test constrained_direct.charge ≈ manual_charge
    @test constrained_direct.current ≈ im * ω * manual_charge
    @test constrained_direct.admittance ≈ im * ω * manual_charge / V0
    @test constrained_direct.solver_info.matrix_size == (2, 2)
    @test constrained_direct.solver_info.reduced_relative_residual < 1.0e-12
    @test PiezoAcousticFEM.solve_direct_voltage(reduced, ω, V0; mechanical_dirichlet).admittance ≈
        constrained_direct.admittance
    @test_throws ArgumentError PiezoAcousticFEM.solve_direct_voltage(
        reduced,
        ω,
        V0;
        mechanical_dirichlet=PiezoAcousticFEM.DirichletDofs([3], [0.0]),
    )
end

@testset "charge and current sign convention" begin
    C = 3.5e-9
    Kuu = ones(1, 1)
    Kuϕ = zeros(1, 1)
    Kϕu = zeros(1, 1)
    Kϕϕ = fill(-C, 1, 1)
    Muu = zeros(1, 1)
    system = PiezoAcousticFEM.KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    partition = PiezoAcousticFEM.HarmonicVoltageDofPartition(1, Int[], [1], Int[])
    reduced = PiezoAcousticFEM.electrode_reduced_k_form(system, partition)

    ω = 2π * 10_000.0
    V0 = 2.0
    solution = PiezoAcousticFEM.solve_direct_voltage(reduced, ω, V0)
    removed_electrode_row_residual = reduced.KPP * V0

    @test removed_electrode_row_residual ≈ -C * V0
    @test solution.charge ≈ C * V0
    @test solution.current ≈ im * ω * C * V0
    @test solution.admittance ≈ im * ω * C
end

@testset "FE direct admittance agrees with dense H-form" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (1, 2),
        Vec{2}((1.0e-3, 0.0)),
        Vec{2}((2.0e-3, 2.0e-3)),
    )
    dense_assembly = PiezoAcousticFEM.assemble_k_form_dense(grid, mat, kin, ip, qr)
    sparse_assembly = PiezoAcousticFEM.assemble_k_form_sparse(grid, mat, kin, ip, qr)
    partition = PiezoAcousticFEM.potential_partition(
        dense_assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )
    sparse_partition = PiezoAcousticFEM.potential_partition(
        sparse_assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )
    reduced = PiezoAcousticFEM.electrode_reduced_k_form(dense_assembly.system, partition)
    sparse_reduced = PiezoAcousticFEM.electrode_reduced_k_form(sparse_assembly.system, sparse_partition)
    h = PiezoAcousticFEM.h_form_dense(dense_assembly.system, partition)

    ω = 2π * 10_000.0
    V0 = 1.25
    direct = PiezoAcousticFEM.solve_direct_voltage(reduced, ω, V0)
    sparse_direct = PiezoAcousticFEM.solve_direct_voltage(sparse_reduced, ω, V0)
    sparse_A, _ = PiezoAcousticFEM.build_direct_voltage_system(sparse_reduced, ω, V0)
    D = h.Huu - ω^2 * h.Muu
    u_h = D \ (-h.Huϕ * V0)
    Q_h = -(only(h.Hϕu * u_h) + h.Hϕϕ * V0)
    Y_h = im * ω * Q_h / V0

    @test length(direct.internal_potential) == 2
    @test issparse(sparse_A)
    @test direct.displacement ≈ u_h
    @test direct.charge ≈ Q_h
    @test direct.admittance ≈ Y_h
    @test sparse_direct.displacement ≈ direct.displacement
    @test sparse_direct.internal_potential ≈ direct.internal_potential
    @test sparse_direct.potential ≈ direct.potential
    @test sparse_direct.charge ≈ direct.charge
    @test sparse_direct.current ≈ direct.current
    @test sparse_direct.admittance ≈ direct.admittance
end

@testset "FE direct solve with axis displacement constraint" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (1, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 2.0e-3)),
    )
    assembly = PiezoAcousticFEM.assemble_k_form_sparse(grid, mat, kin, ip, qr)
    partition = PiezoAcousticFEM.potential_partition(
        assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )
    reduced = PiezoAcousticFEM.electrode_reduced_k_form(assembly.system, partition)
    axis_constraint =
        PiezoAcousticFEM.axis_radial_displacement_constraint(assembly; axis=AxisBoundary(FacetBoundary("left")))

    result = PiezoAcousticFEM.solve_direct_voltage(
        reduced,
        2π * 10_000.0,
        1.0;
        mechanical_dirichlet=axis_constraint,
    )

    @test length(partition.internal) == 2
    @test !isempty(axis_constraint.indices)
    @test result.displacement[axis_constraint.indices] == axis_constraint.values
    @test all(isfinite, result.displacement)
    @test all(isfinite, result.internal_potential)
    @test result.potential[partition.internal] ≈ result.internal_potential
    @test result.potential[partition.driven] == ones(length(partition.driven))
    @test result.potential[partition.grounded] == zeros(length(partition.grounded))
    @test isfinite(result.charge)
    @test isfinite(result.current)
    @test isfinite(result.admittance)
end

@testset "harmonic solve accepts solver config" begin
    kin = AxisymmetricRZ()
    mat = PZT5A()
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (1, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 2.0e-3)),
    )
    problem = PiezoProblem(
        grid,
        mat,
        kin,
        ip,
        qr;
        electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
        boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
        loss=Lossless(),
    )
    analysis = HarmonicVoltageAnalysis(2π * 10_000.0, 1.0, :exp_iomega_t)

    factorized = solve(problem, analysis; solver=FactorizedDirectSolver())

    @test isfinite(factorized.solution.admittance)
    @test factorized.solution.solver_info.method == :factorized_direct
end

@testset "large direct solve residual warns" begin
    info = PiezoAcousticFEM.DirectVoltageSolverInfo(
        :backslash,
        (1, 1),
        1.0,
        1.0e-3,
        1.0e-3,
    )

    @test_logs (:warn, "direct voltage solve returned a large reduced residual") PiezoAcousticFEM.warn_large_direct_voltage_residual(info)
end

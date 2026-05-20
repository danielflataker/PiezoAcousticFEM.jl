mutable struct FailsSecondLinearSolve <: AbstractLinearSolverConfig
    calls::Int
end

FailsSecondLinearSolve() = FailsSecondLinearSolve(0)

PiezoAcousticFEM.solver_method(::FailsSecondLinearSolve) = :fails_second_linear_solve

function PiezoAcousticFEM.solve_linear_system(solver::FailsSecondLinearSolve, A, rhs)
    solver.calls += 1
    solver.calls == 2 && throw(ErrorException("intentional sweep failure"))

    return A \ rhs
end


@testset "axisymmetric direct voltage workflow smoke test" begin
    kin = AxisymmetricRZ()
    mat = PZT5A()
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (8, 12),
        Vec{2}((0.0, 0.0)),
        Vec{2}((4.0e-3, 6.0e-3)),
    )

    analysis = HarmonicVoltageAnalysis(2π * 10_000.0, 1.0, :exp_iomega_t)
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
    assembled = assemble(problem)
    reduction = prepare_analysis(assembled, analysis)
    result = solve(problem, analysis)
    prepared_result = solve(reduction, analysis)

    @test problem.material_source === mat
    @test problem.loss isa Lossless
    @test problem.electrodes.drive isa FacetElectrode
    @test problem.electrodes.reference isa FacetElectrode
    @test only(problem.boundary_conditions.mechanical) isa AxisBoundary
    @test assembled.problem === problem
    @test assembled.material isa AxisymmetricRZPiezoMaterial
    @test reduction.assembled === assembled
    @test reduction.reduced isa PiezoAcousticFEM.ElectrodeReducedKForm
    @test result isa HarmonicVoltageResult
    @test result.reduction.assembled.problem === problem
    @test issparse(result.assembled.assembly.system.Kuu)
    @test issparse(result.assembled.assembly.system.Kuϕ)
    @test issparse(result.assembled.assembly.system.Kϕu)
    @test issparse(result.assembled.assembly.system.Kϕϕ)
    @test issparse(result.assembled.assembly.system.Muu)
    @test length(result.reduction.partition.internal) > 0
    @test !isempty(result.reduction.partition.driven)
    @test !isempty(result.reduction.partition.grounded)
    @test isempty(intersect(result.reduction.partition.driven, result.reduction.partition.grounded))
    @test !isempty(result.reduction.mechanical_dirichlet.indices)
    @test result.solution.displacement[result.reduction.mechanical_dirichlet.indices] ==
        result.reduction.mechanical_dirichlet.values
    @test result.solution.potential[result.reduction.partition.driven] ==
        ones(length(result.reduction.partition.driven))
    @test result.solution.potential[result.reduction.partition.grounded] ==
        zeros(length(result.reduction.partition.grounded))
    @test all(isfinite, result.solution.displacement)
    @test all(isfinite, result.solution.potential)
    @test isfinite(result.solution.drive_terminal_charge)
    @test isfinite(result.solution.drive_terminal_current)
    @test isfinite(result.solution.drive_terminal_admittance)
    @test result.problem === problem
    @test result.analysis === analysis
    @test result.solution.analysis === analysis
    @test result.solution.convention == analysis.convention
    @test prepared_result.problem === problem
    @test prepared_result.assembled === assembled
    @test prepared_result.reduction === reduction
    @test prepared_result.solution.drive_terminal_admittance ≈ result.solution.drive_terminal_admittance
    sweep_analyses = [
        HarmonicVoltageAnalysis(2π * f, 1.0, :exp_iomega_t)
        for f in (10_000.0, 12_000.0)
    ]
    sweep = solve(reduction, sweep_analyses)
    @test sweep isa FrequencySweepResult
    @test sweep.problem === problem
    @test sweep.analyses === sweep_analyses
    @test sweep.assembled === assembled
    @test sweep.reduction === reduction
    @test length(sweep.results) == length(sweep_analyses)
    @test all(point -> point isa HarmonicSweepPointResult, sweep.results)
    @test all(point -> point.status == :success, sweep.results)
    @test all(point -> point.error === nothing, sweep.results)
    @test all(point -> point.result.reduction === reduction, sweep.results)
    @test sweep.reuse_level == :frequency_change
    @test_throws ArgumentError solve(reduction, HarmonicVoltageAnalysis[])
    voltage_sweep = solve(reduction, [
        HarmonicVoltageAnalysis(2π * 10_000.0, voltage, :exp_iomega_t)
        for voltage in (1.0, 2.0)
    ])
    @test voltage_sweep.reuse_level == :analysis_change

    failing_solver = FailsSecondLinearSolve()
    failed_sweep = solve(reduction, sweep_analyses; solver=failing_solver)
    @test failed_sweep.results[1].status == :success
    @test failed_sweep.results[1].result isa HarmonicVoltageResult
    @test failed_sweep.results[1].error === nothing
    @test failed_sweep.results[2].status == :failed
    @test failed_sweep.results[2].result === nothing
    @test failed_sweep.results[2].error isa ErrorException
    @test failed_sweep.results[2].analysis === sweep_analyses[2]

    explicit_solution = PiezoAcousticFEM.solve_direct_voltage(
        reduction.reduced,
        analysis.ω,
        analysis.voltage;
        mechanical_dirichlet=reduction.mechanical_dirichlet,
        analysis,
        convention=analysis.convention,
    )
    @test result.solution.drive_terminal_admittance ≈ explicit_solution.drive_terminal_admittance
end

@testset "direct voltage supports complex material loss" begin
    kin = AxisymmetricRZ()
    mat = PZT5A()
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (2, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    problem = PiezoProblem(
        grid,
        mat,
        kin,
        ip,
        qr;
        electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
        boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
        loss=PhysicalLoss(
            PiezoComplexMaterialLoss(; Qm=50.0, tan_delta=0.02, Qe=80.0),
            NoSystemDamping(),
        ),
    )
    result = solve(problem, HarmonicVoltageAnalysis(2π * 10_000.0, 1.0, :exp_iomega_t))

    @test eltype(result.assembled.material.cᴱ) <: Complex
    @test eltype(result.assembled.assembly.system.Kuu) <: Complex
    @test eltype(result.solution.displacement) <: Complex
    @test eltype(result.solution.potential) <: Complex
    @test isfinite(real(result.solution.drive_terminal_admittance))
    @test isfinite(imag(result.solution.drive_terminal_admittance))
end

@testset "short-circuit modal dense reference" begin
    kin = AxisymmetricRZ()
    mat = PZT5A()
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (2, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
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
    analysis = ShortCircuitModalAnalysis(4)
    result = solve(problem, analysis)

    K = result.reduction.system.Kuu
    M = result.reduction.system.Muu
    free = result.reduction.free_dofs

    @test result isa ShortCircuitModalResult
    @test result.assembled.problem === problem
    @test result.reduction.partition isa PiezoAcousticFEM.ShortCircuitDofPartition
    @test sort(vcat(result.reduction.partition.internal, result.reduction.partition.grounded)) ==
        collect(1:result.reduction.partition.nϕ)
    @test length(result.eigenvalues) == 4
    @test issorted(result.eigenvalues)
    @test minimum(result.eigenvalues) >= -1.0e-10 * maximum(abs, result.eigenvalues)
    @test result.normalization == :mass
    @test size(result.modes) == (length(result.assembled.assembly.dofmap.u_dofs), 4)
    @test result.modes[result.reduction.mechanical_dirichlet.indices, :] ≈
        zeros(length(result.reduction.mechanical_dirichlet.indices), 4)

    Kff = K[free, free]
    Mff = M[free, free]
    for j in axes(result.modes, 2)
        mode = result.modes[:, j]
        free_mode = mode[free]
        @test dot(mode, M * mode) ≈ 1.0 atol = 1.0e-8
        residual = Kff * free_mode - result.eigenvalues[j] * Mff * free_mode
        scale = opnorm(Kff, Inf) * norm(free_mode) +
            abs(result.eigenvalues[j]) * opnorm(Mff, Inf) * norm(free_mode)
        @test norm(residual) <= 1.0e-8 * max(scale, 1.0)
    end

    fields = reconstruct_fields(result, 1)
    @test fields isa NodalFieldOutput
    @test length(fields.displacement) == getnnodes(grid)
    @test length(fields.potential) == getnnodes(grid)
    @test fields.metadata.eigenvalue == result.eigenvalues[1]
    @test fields.metadata.angular_frequency == result.angular_frequencies[1]
    @test fields.metadata.frequency == result.frequencies[1]
    @test fields.metadata.normalization == :mass
    @test_throws ArgumentError reconstruct_fields(result, 0)

    complex_loss_problem = PiezoProblem(
        grid,
        mat,
        kin,
        ip,
        qr;
        electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
        boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
        loss=PhysicalLoss(
            PiezoComplexMaterialLoss(; Qm=50.0, tan_delta=0.02, Qe=80.0),
            NoSystemDamping(),
        ),
    )
    as_given_problem = PiezoProblem(
        grid,
        mat,
        kin,
        ip,
        qr;
        electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
        boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
        loss=PhysicalLoss(MaterialAsGiven(), NoSystemDamping()),
    )
    real_policy_problem = PiezoProblem(
        grid,
        mat,
        kin,
        ip,
        qr;
        electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
        boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
        loss=PhysicalLoss(RealMaterial(), NoSystemDamping()),
    )

    @test solve(real_policy_problem, analysis) isa ShortCircuitModalResult
    @test_throws ArgumentError solve(complex_loss_problem, analysis)
    @test_throws ArgumentError solve(as_given_problem, analysis)
end

@testset "short-circuit modal analysis validates constructor arguments" begin
    @test ShortCircuitModalAnalysis(nothing).nev === nothing
    @test ShortCircuitModalAnalysis(0).nev == 0
    @test_throws ArgumentError ShortCircuitModalAnalysis(-1)
    @test_throws ArgumentError ShortCircuitModalAnalysis(1.5)
end

@testset "VTK solution output" begin
    kin = AxisymmetricRZ()
    mat = PZT5A()
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (2, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
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
    run = solve(problem, HarmonicVoltageAnalysis(2π * 10_000.0, 1.0, :exp_iomega_t))
    basename = tempname()
    fields = reconstruct_fields(run)
    @test fields isa NodalFieldOutput
    @test fields.metadata.output_kind == :nodal_fields
    @test fields.metadata.analysis == :harmonic_voltage
    filename = write_vtk(basename, grid, fields)
    result_basename = tempname()
    result_filename = write_vtk(result_basename, run)

    @test isfile(filename)
    @test filesize(filename) > 0
    @test isfile(result_filename)
    @test filesize(result_filename) > 0
    @test occursin("piezoacousticfem_displacement_unit", read(filename, String))
    result_vtk_text = read(result_filename, String)
    @test occursin("piezoacousticfem_harmonic_convention", result_vtk_text)
    @test occursin("piezoacousticfem_quantity_interpretation", result_vtk_text)
    @test length(fields.displacement) == getnnodes(grid)
    @test length(fields.potential) == getnnodes(grid)
    @test length(fields.radius) == getnnodes(grid)
    rm(filename; force=true)
    rm(result_filename; force=true)
end

@testset "VTK modal output" begin
    kin = AxisymmetricRZ()
    mat = PZT5A()
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (2, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
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
    result = solve(problem, ShortCircuitModalAnalysis(2))
    basename = tempname()
    fields = reconstruct_fields(result, 1)
    @test fields isa NodalFieldOutput
    @test fields.metadata.output_kind == :nodal_fields
    @test fields.metadata.analysis == :short_circuit_modal
    @test fields.metadata.mode_index == 1
    filename = write_vtk(basename, grid, fields)
    result_basename = tempname()
    result_filename = write_vtk(result_basename, result; mode_index=1)

    @test isfile(filename)
    @test isfile(result_filename)
    @test occursin("piezoacousticfem_modal_normalization", read(result_filename, String))
    @test occursin("piezoacousticfem_mode_index", read(result_filename, String))
    @test length(fields.displacement) == getnnodes(grid)
    @test length(fields.potential) == getnnodes(grid)

    rm(filename; force=true)
    rm(result_filename; force=true)
end

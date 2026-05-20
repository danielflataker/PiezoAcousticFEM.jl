function axisymmetric_disk_problem(; cells=(96, 144))
    radius = 4.0e-3
    thickness = 6.0e-3

    grid = generate_grid(
        Quadrilateral,
        cells,
        Vec{2}((0.0, 0.0)),
        Vec{2}((radius, thickness)),
    )

    return PiezoProblem(
        grid,
        PZT5A(),
        AxisymmetricRZ(),
        Lagrange{RefQuadrilateral,1}(),
        QuadratureRule{RefQuadrilateral}(2);
        electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
        boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
        loss=Lossless(),
    )
end

const PIEZO_PIPELINE_CELLS = (96, 144)
const PIEZO_PIPELINE_ANALYSIS = HarmonicVoltageAnalysis(2π * 10_000.0, 1.0, :exp_iomega_t)

assemble_disk_benchmark() =
    @benchmarkable assemble(problem) setup=(problem = axisymmetric_disk_problem(; cells=PIEZO_PIPELINE_CELLS))

prepare_harmonic_reduction_benchmark() =
    @benchmarkable prepare_analysis(assembled, analysis) setup=(
        problem = axisymmetric_disk_problem(; cells=PIEZO_PIPELINE_CELLS);
        analysis = PIEZO_PIPELINE_ANALYSIS;
        assembled = assemble(problem)
    )

solve_prepared_harmonic_benchmark() =
    @benchmarkable solve(reduction, analysis) setup=(
        problem = axisymmetric_disk_problem(; cells=PIEZO_PIPELINE_CELLS);
        analysis = PIEZO_PIPELINE_ANALYSIS;
        reduction = prepare_analysis(assemble(problem), analysis)
    )

full_harmonic_pipeline_benchmark() =
    @benchmarkable solve(problem, analysis) setup=(
        problem = axisymmetric_disk_problem(; cells=PIEZO_PIPELINE_CELLS);
        analysis = PIEZO_PIPELINE_ANALYSIS
    )

reconstruct_fields_benchmark() =
    @benchmarkable reconstruct_fields(result) setup=(
        problem = axisymmetric_disk_problem(; cells=PIEZO_PIPELINE_CELLS);
        analysis = PIEZO_PIPELINE_ANALYSIS;
        result = solve(problem, analysis)
    )

register_benchmark_suite!(BenchmarkSuiteSpec(
    "piezo-pipeline",
    "piezo pipeline",
    [
        BenchmarkSpec(
            "assemble_disk_96x144",
            "Assemble sparse K-form for a 96x144 axisymmetric disk",
            assemble_disk_benchmark(),
        ),
        BenchmarkSpec(
            "prepare_harmonic_reduction_disk_96x144",
            "Prepare electrode and Dirichlet reductions for a 96x144 disk",
            prepare_harmonic_reduction_benchmark(),
        ),
        BenchmarkSpec(
            "solve_prepared_harmonic_disk_96x144",
            "Solve a prepared harmonic voltage problem for a 96x144 disk",
            solve_prepared_harmonic_benchmark(),
        ),
        BenchmarkSpec(
            "reconstruct_fields_disk_96x144",
            "Reconstruct nodal output fields for a solved 96x144 disk",
            reconstruct_fields_benchmark(),
        ),
        BenchmarkSpec(
            "full_harmonic_pipeline_disk_96x144",
            "Assemble, prepare, and solve a harmonic voltage problem for a 96x144 disk",
            full_harmonic_pipeline_benchmark(),
        ),
    ],
))

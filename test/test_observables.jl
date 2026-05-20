@testset "observables extract harmonic voltage quantities" begin
    grid = generate_grid(
        Quadrilateral,
        (2, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    problem = validation_test_problem(; grid, ip=Lagrange{RefQuadrilateral,1}(), qr=QuadratureRule{RefQuadrilateral}(2))
    result = solve(problem, HarmonicVoltageAnalysis(2π * 10_000.0, 1.0, :exp_iomega_t))

    @test evaluate(AdmittanceObservable(), result) == result.solution.drive_terminal_admittance
    @test evaluate(ChargeObservable(), result) == result.solution.drive_terminal_charge
    @test evaluate(CurrentObservable(), result) == result.solution.drive_terminal_current
end

@testset "observables extract modal quantities" begin
    grid = generate_grid(
        Quadrilateral,
        (2, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    problem = validation_test_problem(; grid, ip=Lagrange{RefQuadrilateral,1}(), qr=QuadratureRule{RefQuadrilateral}(2))
    result = solve(problem, ShortCircuitModalAnalysis(4))

    @test evaluate(ModalFrequenciesObservable(), result) == result.frequencies
    @test_throws MethodError evaluate(AdmittanceObservable(), result)
end

@testset "observables extract harmonic voltage quantities" begin
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
    result = solve(problem, HarmonicVoltageAnalysis(2π * 10_000.0, 1.0, :exp_iomega_t))

    @test evaluate(AdmittanceObservable(), result) == result.solution.admittance
    @test evaluate(ChargeObservable(), result) == result.solution.charge
    @test evaluate(CurrentObservable(), result) == result.solution.current
end

@testset "observables extract modal quantities" begin
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
    result = solve(problem, ShortCircuitModalAnalysis(4))

    @test evaluate(ModalFrequenciesObservable(), result) == result.frequencies
    @test_throws MethodError evaluate(AdmittanceObservable(), result)
end

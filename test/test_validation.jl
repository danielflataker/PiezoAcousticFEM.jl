function validation_test_problem(;
    grid,
    ip,
    qr,
    electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
    boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
)
    return PiezoProblem(
        grid,
        PZT5A(),
        AxisymmetricRZ(),
        ip,
        qr;
        electrodes,
        boundary_conditions,
        loss=Lossless(),
    )
end


@testset "problem validation" begin
    grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)

    problem = validation_test_problem(; grid, ip, qr)
    @test validate(problem) === nothing

    wrong_shape_problem = validation_test_problem(;
        grid,
        ip=Lagrange{RefTriangle,1}(),
        qr=QuadratureRule{RefTriangle}(2),
    )
    @test_throws ArgumentError validate(wrong_shape_problem)

    missing_electrode_problem = validation_test_problem(;
        grid,
        ip,
        qr,
        electrodes=TwoTerminalElectrodes(FacetElectrode("missing"), FacetElectrode("bottom")),
    )
    @test_throws ArgumentError validate(missing_electrode_problem)

    overlapping_electrodes_problem = validation_test_problem(;
        grid,
        ip,
        qr,
        electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("top")),
    )
    @test_throws ArgumentError validate(overlapping_electrodes_problem)

    off_axis_problem = validation_test_problem(;
        grid,
        ip,
        qr,
        boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("right")),)),
    )
    @test_throws ArgumentError validate(off_axis_problem)
end

function validation_test_problem(; grid, ip, qr)
    return PiezoProblem(
        grid,
        PZT5A(),
        AxisymmetricRZ(),
        ip,
        qr;
        electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
        boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
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

    quadratic_problem = validation_test_problem(;
        grid,
        ip=Lagrange{RefQuadrilateral,2}(),
        qr=QuadratureRule{RefQuadrilateral}(3),
    )
    @test_throws ArgumentError assemble(quadratic_problem)

    wrong_shape_problem = validation_test_problem(;
        grid,
        ip=Lagrange{RefTriangle,1}(),
        qr=QuadratureRule{RefTriangle}(2),
    )
    @test_throws ArgumentError validate(wrong_shape_problem)
end

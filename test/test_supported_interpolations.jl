@testset "supported interpolations" begin
    grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )

    quadratic_problem = validation_test_problem(;
        grid,
        ip=Lagrange{RefQuadrilateral,2}(),
        qr=QuadratureRule{RefQuadrilateral}(3),
    )
    @test validate(quadratic_problem) === nothing

    serendipity_problem = validation_test_problem(;
        grid,
        ip=Serendipity{RefQuadrilateral,2}(),
        qr=QuadratureRule{RefQuadrilateral}(3),
    )
    serendipity_result = solve(serendipity_problem, HarmonicVoltageAnalysis(2π * 10_000.0, 1.0, :exp_iomega_t))
    serendipity_fields = reconstruct_field_dofs(serendipity_result)
    @test length(serendipity_fields.displacement) == 16
    @test length(serendipity_fields.potential) == 8
    @test isfinite(serendipity_result.solution.drive_terminal_admittance)

    serendipity_modal = solve(serendipity_problem, ShortCircuitModalAnalysis(2))
    serendipity_mode_fields = reconstruct_field_dofs(serendipity_modal, 1)
    @test length(serendipity_mode_fields.displacement) == 16
    @test length(serendipity_mode_fields.potential) == 8
end

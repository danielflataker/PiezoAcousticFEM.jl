@testset "element-local derived fields" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((1.0, 0.0)),
        Vec{2}((2.0, 1.0)),
    )
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    assembly = PiezoAcousticFEM.assemble_k_form_sparse(grid, mat, kin, ip, qr)

    displacement = zeros(length(assembly.dofmap.u_dofs))
    potential = zeros(length(assembly.dofmap.ϕ_dofs))
    for nodeid in 1:getnnodes(grid)
        x = Ferrite.get_node_coordinate(grid, nodeid)
        displacement[PiezoAcousticFEM.compact_displacement_dof(assembly.dofmap, nodeid, 1)] = 2x[1] + 3x[2]
        displacement[PiezoAcousticFEM.compact_displacement_dof(assembly.dofmap, nodeid, 2)] = -x[1] + 4x[2]
        potential[PiezoAcousticFEM.compact_potential_node_dof(assembly.dofmap, nodeid)] = 5x[1] - 7x[2]
    end
    fields = CompactFieldDofs(displacement, potential; metadata=(analysis=:unit_test,))

    ξ = Vec{2}((0.25, -0.5))
    output = PiezoAcousticFEM.evaluate_derived_fields(
        assembly,
        mat,
        kin,
        ip,
        fields,
        1,
        [ξ],
    )

    @test output isa ElementDerivedFieldOutput
    @test output.cellid == 1
    @test output.reference_points == [ξ]
    @test output.metadata.output_kind == :element_derived_fields
    @test output.metadata.analysis == :unit_test
    @test output.metadata.electric_field_units == "V/m"

    x = only(output.coordinates)
    expected_u = Vec{2}((2x[1] + 3x[2], -x[1] + 4x[2]))
    expected_∇u = [2.0 3.0; -1.0 4.0]
    expected_∇ϕ = Vec{2}((5.0, -7.0))
    expected_E = Vec{2}((-5.0, 7.0))
    expected_S = PiezoAcousticFEM.strain(kin, expected_u, expected_∇u, x)

    @test only(output.displacement) ≈ expected_u
    @test only(output.potential) ≈ 5x[1] - 7x[2]
    @test only(output.displacement_gradient) ≈ expected_∇u
    @test only(output.potential_gradient) ≈ expected_∇ϕ
    @test only(output.strain) ≈ expected_S
    @test only(output.electric_field) ≈ expected_E
    @test only(output.stress) ≈ PiezoAcousticFEM.stress(mat, expected_S, expected_E)
    @test only(output.electric_displacement) ≈
          PiezoAcousticFEM.electric_displacement(mat, expected_S, expected_E)
    @test_throws ArgumentError PiezoAcousticFEM.evaluate_derived_fields(
        assembly,
        mat,
        kin,
        ip,
        fields,
        2,
        [ξ],
    )
end

@testset "per-element derived field sampling" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((1.0, 0.0)),
        Vec{2}((2.0, 1.0)),
    )
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    assembly = PiezoAcousticFEM.assemble_k_form_sparse(grid, mat, kin, ip, qr)

    displacement = zeros(length(assembly.dofmap.u_dofs))
    potential = zeros(length(assembly.dofmap.ϕ_dofs))
    for nodeid in 1:getnnodes(grid)
        x = Ferrite.get_node_coordinate(grid, nodeid)
        displacement[PiezoAcousticFEM.compact_displacement_dof(assembly.dofmap, nodeid, 1)] = 2x[1] + 3x[2]
        displacement[PiezoAcousticFEM.compact_displacement_dof(assembly.dofmap, nodeid, 2)] = -x[1] + 4x[2]
        potential[PiezoAcousticFEM.compact_potential_node_dof(assembly.dofmap, nodeid)] = 5x[1] - 7x[2]
    end
    fields = CompactFieldDofs(displacement, potential)

    output = sample_element_fields(
        assembly,
        mat,
        kin,
        ip,
        fields,
        1;
        samples_per_axis=(3, 2),
    )

    @test output isa ElementDerivedFieldOutput
    @test output.metadata.output_kind == :element_sampled_fields
    @test output.metadata.evaluation == :reference_lattice
    @test output.metadata.samples_per_axis == (3, 2)
    @test length(output.reference_points) == 6
    @test output.reference_points == [
        Vec{2}((-1.0, -1.0)),
        Vec{2}((0.0, -1.0)),
        Vec{2}((1.0, -1.0)),
        Vec{2}((-1.0, 1.0)),
        Vec{2}((0.0, 1.0)),
        Vec{2}((1.0, 1.0)),
    ]

    first_x = first(output.coordinates)
    @test first(output.displacement) ≈ Vec{2}((2first_x[1] + 3first_x[2], -first_x[1] + 4first_x[2]))
    @test first(output.potential) ≈ 5first_x[1] - 7first_x[2]

    @test_throws ArgumentError sample_element_fields(
        assembly,
        mat,
        kin,
        ip,
        fields,
        1;
        samples_per_axis=1,
    )
end

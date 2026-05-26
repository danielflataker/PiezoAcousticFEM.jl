function linear_derived_field_fixture()
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

    return (; kin, mat, assembly, ip, fields=CompactFieldDofs(displacement, potential; metadata=(analysis=:unit_test,)))
end

function linear_assembled_fixture()
    fixture = linear_derived_field_fixture()
    problem = PiezoProblem(
        Ferrite.get_grid(fixture.assembly.dofhandler),
        PZT5A(),
        fixture.kin,
        fixture.ip,
        QuadratureRule{RefQuadrilateral}(2);
        electrodes=TwoTerminalElectrodes(FacetElectrode("right"), FacetElectrode("left")),
        boundary_conditions=AxisymmetricBoundaryConditions(()),
        loss=Lossless(),
    )

    return (; fixture..., assembled=PiezoAcousticFEM.AssembledPiezoProblem(problem, fixture.mat, fixture.assembly))
end

@testset "element-local derived fields" begin
    fixture = linear_assembled_fixture()

    ξ = Vec{2}((0.25, -0.5))
    output = evaluate_derived_fields(
        fixture.assembled,
        fixture.fields,
        1,
        [ξ],
    )

    @test output isa PiezoAcousticFEM.ElementDerivedFieldOutput
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
    expected_S = PiezoAcousticFEM.strain(fixture.kin, expected_u, expected_∇u, x)

    @test only(output.displacement) ≈ expected_u
    @test only(output.potential) ≈ 5x[1] - 7x[2]
    @test only(output.displacement_gradient) ≈ expected_∇u
    @test only(output.potential_gradient) ≈ expected_∇ϕ
    @test only(output.strain) ≈ expected_S
    @test only(output.electric_field) ≈ expected_E
    @test only(output.stress) ≈ PiezoAcousticFEM.stress(fixture.mat, expected_S, expected_E)
    @test only(output.electric_displacement) ≈
          PiezoAcousticFEM.electric_displacement(fixture.mat, expected_S, expected_E)
    @test_throws ArgumentError evaluate_derived_fields(
        fixture.assembled,
        fixture.fields,
        2,
        [ξ],
    )
end

@testset "per-element derived field sampling" begin
    fixture = linear_assembled_fixture()

    output = sample_element_fields(
        fixture.assembled,
        fixture.fields,
        1;
        samples_per_axis=(3, 2),
    )

    @test output isa PiezoAcousticFEM.ElementDerivedFieldOutput
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
        fixture.assembled,
        fixture.fields,
        1;
        samples_per_axis=1,
    )
end

@testset "derived field VTK cell output" begin
    fixture = linear_assembled_fixture()
    grid = Ferrite.get_grid(fixture.assembly.dofhandler)
    output = sample_element_fields(
        fixture.assembled,
        fixture.fields,
        1;
        samples_per_axis=(2, 2),
    )

    filename = write_vtk(tempname(), grid, [output])
    vtk_text = read(filename, String)

    @test isfile(filename)
    @test occursin("<CellData", vtk_text)
    @test occursin("derived_potential", vtk_text)
    @test occursin("derived_displacement_r", vtk_text)
    @test occursin("derived_strain_rr", vtk_text)
    @test occursin("derived_electric_field_r", vtk_text)
    @test occursin("derived_stress_rr", vtk_text)
    @test occursin("piezoacousticfem_output_kind", vtk_text)
    @test occursin("piezoacousticfem_cell_data_policy", vtk_text)
    @test occursin("piezoacousticfem_evaluation", vtk_text)
    @test_throws ArgumentError write_vtk(tempname(), grid, typeof(output)[])

    rm(filename; force=true)
end

@testset "physical RZ-grid derived field sampling" begin
    fixture = linear_assembled_fixture()

    r_coordinates = [1.0, 1.5, 2.0]
    z_coordinates = [0.0, 0.25, 1.0]
    output = sample_physical_grid(
        fixture.assembled,
        fixture.fields,
        r_coordinates,
        z_coordinates,
    )

    @test output isa PiezoAcousticFEM.PhysicalGridDerivedFieldOutput
    @test output.r_coordinates == r_coordinates
    @test output.z_coordinates == z_coordinates
    @test output.metadata.output_kind == :physical_grid_derived_fields
    @test output.metadata.evaluation == :physical_grid
    @test output.metadata.sample_order == :z_major_r_fastest
    @test length(output.coordinates) == length(r_coordinates) * length(z_coordinates)
    @test output.coordinates == [Vec{2}((r, z)) for z in z_coordinates for r in r_coordinates]
    @test all(==(1), output.cellids)
    @test length(output.reference_points) == length(output.coordinates)

    for (i, x) in pairs(output.coordinates)
        @test output.displacement[i] ≈ Vec{2}((2x[1] + 3x[2], -x[1] + 4x[2]))
        @test output.potential[i] ≈ 5x[1] - 7x[2]
        @test output.potential_gradient[i] ≈ Vec{2}((5.0, -7.0))
    end

    @test_throws ArgumentError sample_physical_grid(
        fixture.assembled,
        fixture.fields,
        [0.5],
        [0.5],
    )
    @test_throws ArgumentError sample_physical_grid(
        fixture.assembled,
        fixture.fields,
        Float64[],
        [0.5],
    )
end

@testset "Dirichlet elimination" begin
    A = [
        4.0 1.0 2.0
        1.0 3.0 0.0
        2.0 0.0 5.0
    ]
    b = [1.0, 2.0, 3.0]
    constraints = PiezoAcousticFEM.DirichletDofs([2, 3], [10.0, -1.0])
    reduction = PiezoAcousticFEM.apply_dirichlet(A, b, constraints)

    @test reduction.free == [1]
    @test reduction.constrained == [2, 3]
    @test reduction.A == A[1:1, 1:1]
    @test reduction.b ≈ [1.0 - dot(A[1, [2, 3]], constraints.values)]

    x_free = reduction.A \ reduction.b
    x = PiezoAcousticFEM.reconstruct_solution(reduction, x_free)
    @test x[1] ≈ only(x_free)
    @test x[2] == 10.0
    @test x[3] == -1.0

    sparse_reduction = PiezoAcousticFEM.apply_dirichlet(sparse(A), b, constraints)
    @test issparse(sparse_reduction.A)
    @test Matrix(sparse_reduction.A) == reduction.A
    @test sparse_reduction.b ≈ reduction.b

    @test_throws ArgumentError PiezoAcousticFEM.DirichletDofs([1, 1], [0.0, 1.0])
    @test_throws DimensionMismatch PiezoAcousticFEM.DirichletDofs([1], [0.0, 1.0])
    @test_throws ArgumentError PiezoAcousticFEM.apply_dirichlet(A, b, PiezoAcousticFEM.DirichletDofs([4], [0.0]))
end

@testset "Ferrite adapter for electrode partitions" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)

    grid = generate_grid(
        Quadrilateral,
        (2, 1),
        Vec{2}((1.0e-3, 0.0)),
        Vec{2}((3.0e-3, 1.0e-3)),
    )

    assembly = PiezoAcousticFEM.assemble_k_form_sparse(grid, mat, kin, ip, qr)
    partition = PiezoAcousticFEM.potential_partition(
        assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )

    @test isempty(partition.internal)
    @test length(partition.driven) == 3
    @test length(partition.grounded) == 3
    @test sort(vcat(partition.driven, partition.grounded)) == collect(1:6)
    @test isempty(intersect(partition.driven, partition.grounded))

    h = PiezoAcousticFEM.h_form_dense(assembly.system, partition)
    @test size(h.Huu) == (12, 12)
    @test size(h.Huϕ) == (12,)
    @test size(h.Hϕu) == (1, 12)
    @test h.Huu ≈ assembly.system.Kuu

    through_thickness_grid = generate_grid(
        Quadrilateral,
        (1, 2),
        Vec{2}((1.0e-3, 0.0)),
        Vec{2}((2.0e-3, 2.0e-3)),
    )
    through_thickness_assembly = PiezoAcousticFEM.assemble_k_form_sparse(through_thickness_grid, mat, kin, ip, qr)
    through_thickness_partition = PiezoAcousticFEM.potential_partition(
        through_thickness_assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )

    @test length(through_thickness_partition.internal) == 2
    @test length(through_thickness_partition.driven) == 2
    @test length(through_thickness_partition.grounded) == 2
    @test sort(vcat(
        through_thickness_partition.internal,
        through_thickness_partition.driven,
        through_thickness_partition.grounded,
    )) == collect(1:6)
    @test through_thickness_partition.nϕ == 6

    explicit_partition = PiezoAcousticFEM.HarmonicVoltageDofPartition(6, [3, 4], [5, 6], [1, 2])
    @test explicit_partition.nϕ == 6
    @test_throws ArgumentError PiezoAcousticFEM.HarmonicVoltageDofPartition(6, [3], [5, 6], [1, 2])
    @test_throws ArgumentError PiezoAcousticFEM.HarmonicVoltageDofPartition(6, [3, 3], [5, 6], [1, 2])
    @test_throws ArgumentError PiezoAcousticFEM.HarmonicVoltageDofPartition(6, [3, 4], [6, 7], [1, 2])

    short_partition = PiezoAcousticFEM.short_circuit_partition(
        through_thickness_assembly,
        TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
    )
    @test length(short_partition.internal) == 2
    @test length(short_partition.grounded) == 4
    @test sort(vcat(short_partition.internal, short_partition.grounded)) == collect(1:6)
end

@testset "Ferrite adapter for displacement constraints" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    assembly = PiezoAcousticFEM.assemble_k_form_sparse(grid, mat, kin, ip, qr)

    radial = PiezoAcousticFEM.displacement_component_dofs_on_facets(assembly, getfacetset(grid, "left"), :r)
    axial = PiezoAcousticFEM.displacement_component_dofs_on_facets(assembly, getfacetset(grid, "left"), :z)
    axis_constraint =
        PiezoAcousticFEM.axis_radial_displacement_constraint(assembly; axis=AxisBoundary(FacetBoundary("left")))

    @test radial == [1, 7]
    @test axial == [2, 8]
    @test PiezoAcousticFEM.compact_displacement_dof(assembly.dofmap, 1, :r) == 1
    @test PiezoAcousticFEM.compact_displacement_dof(assembly.dofmap, 1, :z) == 2
    @test PiezoAcousticFEM.compact_potential_node_dof(assembly.dofmap, 1) == 1
    @test axis_constraint.indices == radial
    @test axis_constraint.values == zeros(length(radial))

    x_free = ones(length(assembly.dofmap.u_dofs) - length(radial))
    reduction = PiezoAcousticFEM.DirichletReduction(
        Matrix{Float64}(I, length(x_free), length(x_free)),
        ones(length(x_free)),
        setdiff(collect(1:length(assembly.dofmap.u_dofs)), radial),
        radial,
        axis_constraint.values,
        length(assembly.dofmap.u_dofs),
    )
    full_u = PiezoAcousticFEM.reconstruct_solution(reduction, x_free)
    @test full_u[radial] == axis_constraint.values
    @test all(==(1.0), full_u[axial])
end


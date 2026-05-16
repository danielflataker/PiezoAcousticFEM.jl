using Test
using LinearAlgebra
using SparseArrays
using Ferrite
using PiezoAcousticFEM

@testset "PZT5A constants" begin
    mat = effective_material(PZT5A(), AxisymmetricRZ(), Lossless())
    @test mat.cᴱ[4, 4] == 2.26e10
    @test mat.e[1, 4] == 12.3
    @test mat.e[2, 1] == -5.4
    @test mat.e[2, 3] == 15.8
    @test mat.ρ == 7750

    full = PZT5A(; TC=ComplexF64, TE=Float32, Tε=BigFloat, Tρ=Float64)
    reduced = reduce_material(full, AxisymmetricRZ())
    @test full isa Piezo6mmConstants
    @test reduced isa Piezo6mmAxi
    @test eltype(reduced.cᴱ) == ComplexF64
    @test eltype(reduced.e) == Float32
    @test eltype(reduced.εˢ) == BigFloat
    @test typeof(reduced.ρ) == Float64

    base = PZT5A()
    complex_full = Piezo6mmConstants(
        complex.(base.cᴱ, base.cᴱ .* 0.01),
        complex.(base.e, base.e .* 0.01),
        complex.(base.εˢ, base.εˢ .* 0.01),
        base.ρ,
    )
    effective = effective_material(complex_full, AxisymmetricRZ(), Lossless())
    @test effective isa Piezo6mmAxi
    @test !(eltype(effective.cᴱ) <: Complex)
    @test effective.cᴱ ≈ real.(reduce_material(complex_full, AxisymmetricRZ()).cᴱ)

    as_given = effective_material(
        complex_full,
        AxisymmetricRZ(),
        PhysicalLoss(MaterialAsGiven(), NoSystemDamping()),
    )
    real_policy = effective_material(
        complex_full,
        AxisymmetricRZ(),
        PhysicalLoss(RealMaterial(), NoSystemDamping()),
    )
    @test eltype(as_given.cᴱ) <: Complex
    @test as_given.cᴱ ≈ reduce_material(complex_full, AxisymmetricRZ()).cᴱ
    @test !(eltype(real_policy.cᴱ) <: Complex)
    @test real_policy.cᴱ ≈ effective.cᴱ
end

@testset "one element K-form blocks" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())

    grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((1.0e-3, 0.0)),
        Vec{2}((2.0e-3, 1.0e-3)),
    )

    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    cellvalues_u = CellValues(qr, ip)
    cellvalues_ϕ = CellValues(qr, ip)

    cell = first(CellIterator(grid))
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_ϕ, cell)

    elem = piezo_element_matrices(cellvalues_u, cellvalues_ϕ, getcoordinates(cell), mat, kin)

    @test size(elem.Kuu) == (8, 8)
    @test size(elem.Kuϕ) == (8, 4)
    @test size(elem.Kϕu) == (4, 8)
    @test size(elem.Kϕϕ) == (4, 4)
    @test size(elem.Muu) == (8, 8)
    @test elem.Kuu ≈ transpose(elem.Kuu)
    @test elem.Kϕϕ ≈ transpose(elem.Kϕϕ)
    @test elem.Kuϕ ≈ transpose(elem.Kϕu)
    @test maximum(eigvals(Symmetric(elem.Kϕϕ))) <= 1.0e-18
end

@testset "Dirichlet elimination" begin
    A = [
        4.0 1.0 2.0
        1.0 3.0 0.0
        2.0 0.0 5.0
    ]
    b = [1.0, 2.0, 3.0]
    constraints = DirichletDofs([2, 3], [10.0, -1.0])
    reduction = apply_dirichlet(A, b, constraints)

    @test reduction.free == [1]
    @test reduction.constrained == [2, 3]
    @test reduction.A == A[1:1, 1:1]
    @test reduction.b ≈ [1.0 - dot(A[1, [2, 3]], constraints.values)]

    x_free = reduction.A \ reduction.b
    x = reconstruct_solution(reduction, x_free)
    @test x[1] ≈ only(x_free)
    @test x[2] == 10.0
    @test x[3] == -1.0

    sparse_reduction = apply_dirichlet(sparse(A), b, constraints)
    @test issparse(sparse_reduction.A)
    @test Matrix(sparse_reduction.A) == reduction.A
    @test sparse_reduction.b ≈ reduction.b

    @test_throws ArgumentError DirichletDofs([1, 1], [0.0, 1.0])
    @test_throws DimensionMismatch DirichletDofs([1], [0.0, 1.0])
    @test_throws ArgumentError apply_dirichlet(A, b, DirichletDofs([4], [0.0]))
end

@testset "axisymmetric geometry guards" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)

    @test_throws ArgumentError integration_weight(kin, Vec{2}((0.0, 0.0)), 1.0)
    @test_throws ArgumentError integration_weight(kin, Vec{2}((1.0, 0.0)), 0.0)
    @test_throws ArgumentError strain(kin, Vec{2}((1.0, 0.0)), zeros(2, 2), Vec{2}((0.0, 0.0)))

    axis_touching_grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    axis_touching_assembly = assemble_k_form_sparse(axis_touching_grid, mat, kin, ip, qr)
    @test all(isfinite, axis_touching_assembly.system.Kuu)

    negative_radius_grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((-1.0e-3, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    @test_throws ArgumentError assemble_k_form_sparse(negative_radius_grid, mat, kin, ip, qr)
end

@testset "dense global K-form assembly" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)

    one_cell_grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((1.0e-3, 0.0)),
        Vec{2}((2.0e-3, 1.0e-3)),
    )

    assembly = assemble_k_form_dense(one_cell_grid, mat, kin, ip, qr)

    cellvalues_u = CellValues(qr, ip)
    cellvalues_ϕ = CellValues(qr, ip)
    cell = first(CellIterator(one_cell_grid))
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_ϕ, cell)
    elem = piezo_element_matrices(cellvalues_u, cellvalues_ϕ, getcoordinates(cell), mat, kin)

    @test assembly.system.Kuu ≈ elem.Kuu
    @test assembly.system.Kuϕ ≈ elem.Kuϕ
    @test assembly.system.Kϕu ≈ elem.Kϕu
    @test assembly.system.Kϕϕ ≈ elem.Kϕϕ
    @test assembly.system.Muu ≈ elem.Muu

    two_cell_grid = generate_grid(
        Quadrilateral,
        (2, 1),
        Vec{2}((1.0e-3, 0.0)),
        Vec{2}((3.0e-3, 1.0e-3)),
    )

    two_cell_assembly = assemble_k_form_dense(two_cell_grid, mat, kin, ip, qr)
    @test size(two_cell_assembly.system.Kuu) == (12, 12)
    @test size(two_cell_assembly.system.Kuϕ) == (12, 6)
    @test size(two_cell_assembly.system.Kϕu) == (6, 12)
    @test size(two_cell_assembly.system.Kϕϕ) == (6, 6)
    @test size(two_cell_assembly.system.Muu) == (12, 12)
    @test two_cell_assembly.system.Kuu ≈ transpose(two_cell_assembly.system.Kuu)
    @test two_cell_assembly.system.Kϕϕ ≈ transpose(two_cell_assembly.system.Kϕϕ)
    @test two_cell_assembly.system.Kuϕ ≈ transpose(two_cell_assembly.system.Kϕu)
end

@testset "sparse global K-form assembly matches dense reference" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (2, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((2.0e-3, 2.0e-3)),
    )

    dense = assemble_k_form_dense(grid, mat, kin, ip, qr)
    sparse_assembly = assemble_k_form_sparse(grid, mat, kin, ip, qr)

    @test issparse(sparse_assembly.system.Kuu)
    @test issparse(sparse_assembly.system.Kuϕ)
    @test issparse(sparse_assembly.system.Kϕu)
    @test issparse(sparse_assembly.system.Kϕϕ)
    @test issparse(sparse_assembly.system.Muu)
    @test sparse_assembly.dofmap.u_dofs == dense.dofmap.u_dofs
    @test sparse_assembly.dofmap.ϕ_dofs == dense.dofmap.ϕ_dofs
    @test Matrix(sparse_assembly.system.Kuu) ≈ dense.system.Kuu
    @test Matrix(sparse_assembly.system.Kuϕ) ≈ dense.system.Kuϕ
    @test Matrix(sparse_assembly.system.Kϕu) ≈ dense.system.Kϕu
    @test Matrix(sparse_assembly.system.Kϕϕ) ≈ dense.system.Kϕϕ
    @test Matrix(sparse_assembly.system.Muu) ≈ dense.system.Muu
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

    assembly = assemble_k_form_sparse(grid, mat, kin, ip, qr)
    partition = potential_partition(
        assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )

    @test isempty(partition.internal)
    @test length(partition.driven) == 3
    @test length(partition.grounded) == 3
    @test sort(vcat(partition.driven, partition.grounded)) == collect(1:6)
    @test isempty(intersect(partition.driven, partition.grounded))

    h = h_form_dense(assembly.system, partition)
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
    through_thickness_assembly = assemble_k_form_sparse(through_thickness_grid, mat, kin, ip, qr)
    through_thickness_partition = potential_partition(
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

    explicit_partition = HarmonicVoltageDofPartition(6, [3, 4], [5, 6], [1, 2])
    @test explicit_partition.nϕ == 6
    @test_throws ArgumentError HarmonicVoltageDofPartition(6, [3], [5, 6], [1, 2])
    @test_throws ArgumentError HarmonicVoltageDofPartition(6, [3, 3], [5, 6], [1, 2])
    @test_throws ArgumentError HarmonicVoltageDofPartition(6, [3, 4], [6, 7], [1, 2])

    short_partition = short_circuit_partition(
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
    assembly = assemble_k_form_sparse(grid, mat, kin, ip, qr)

    radial = displacement_component_dofs_on_facets(assembly, getfacetset(grid, "left"), :r)
    axial = displacement_component_dofs_on_facets(assembly, getfacetset(grid, "left"), :z)
    axis_constraint =
        axis_radial_displacement_constraint(assembly; axis=AxisBoundary(FacetBoundary("left")))

    @test radial == [1, 7]
    @test axial == [2, 8]
    @test compact_displacement_dof(assembly.dofmap, 1, :r) == 1
    @test compact_displacement_dof(assembly.dofmap, 1, :z) == 2
    @test compact_potential_node_dof(assembly.dofmap, 1) == 1
    @test axis_constraint.indices == radial
    @test axis_constraint.values == zeros(length(radial))

    x_free = ones(length(assembly.dofmap.u_dofs) - length(radial))
    reduction = DirichletReduction(
        Matrix{Float64}(I, length(x_free), length(x_free)),
        ones(length(x_free)),
        setdiff(collect(1:length(assembly.dofmap.u_dofs)), radial),
        radial,
        axis_constraint.values,
        length(assembly.dofmap.u_dofs),
    )
    full_u = reconstruct_solution(reduction, x_free)
    @test full_u[radial] == axis_constraint.values
    @test all(==(1.0), full_u[axial])
end

@testset "global H-form Schur complement" begin
    Kuu = [10.0 1.0; 1.0 8.0]
    Kuϕ = [2.0 3.0 100.0; 4.0 5.0 200.0]
    Kϕu = Matrix(transpose(Kuϕ))
    Kϕϕ = [
        -6.0 -1.0 0.0
        -1.0 -7.0 0.0
        0.0 0.0 -9.0
    ]
    Muu = Matrix{Float64}(I, 2, 2)

    system = KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    partition = HarmonicVoltageDofPartition(3, [1], [2], [3])
    h = h_form_dense(system, partition)

    Kii = Kϕϕ[1:1, 1:1]
    @test h.Huu ≈ Kuu - Kuϕ[:, 1:1] * (Kii \ Kϕu[1:1, :])
    @test h.Huϕ ≈ vec(Kuϕ[:, 2:2] - Kuϕ[:, 1:1] * (Kii \ Kϕϕ[1:1, 2:2]))
    @test h.Hϕu ≈ Kϕu[2:2, :] - Kϕϕ[2:2, 1:1] * (Kii \ Kϕu[1:1, :])
    @test h.Hϕϕ ≈ only(Kϕϕ[2:2, 2:2] - Kϕϕ[2:2, 1:1] * (Kii \ Kϕϕ[1:1, 2:2]))
    @test h.Muu == Muu
end

@testset "K-form block types are independent" begin
    Kuu = ComplexF64[10.0 1.0; 1.0 8.0]
    Kuϕ = Float32[2.0 3.0 100.0; 4.0 5.0 200.0]
    Kϕu = Float32[
        2.0 4.0
        3.0 5.0
        100.0 200.0
    ]
    Kϕϕ = BigFloat[
        -6.0 -1.0 0.0
        -1.0 -7.0 0.0
        0.0 0.0 -9.0
    ]
    Muu = Float64[1.0 0.0; 0.0 1.0]

    system = KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    partition = HarmonicVoltageDofPartition(3, [1], [2], [3])
    reduced = electrode_reduced_k_form(system, partition)

    @test eltype(reduced.Kuu) == ComplexF64
    @test eltype(reduced.KuP) == Float32
    @test eltype(reduced.KPu) == Float32
    @test eltype(reduced.Kii) == BigFloat
    @test eltype(reduced.Muu) == Float64
    @test typeof(reduced.KPP) == BigFloat
end

@testset "dense electrode-reduced direct voltage solve" begin
    Kuu = [10.0 1.0; 1.0 8.0]
    Kuϕ = [2.0 3.0 100.0; 4.0 5.0 200.0]
    Kϕu = Matrix(transpose(Kuϕ))
    Kϕϕ = [
        -6.0 -1.0 0.0
        -1.0 -7.0 0.0
        0.0 0.0 -9.0
    ]
    Muu = Matrix{Float64}(I, 2, 2)
    system = KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    partition = HarmonicVoltageDofPartition(3, [1], [2], [3])

    reduced = electrode_reduced_k_form(system, partition)
    h = h_form_dense(system, partition)
    ω = 3.0
    V0 = 2.5

    direct = solve_direct_voltage(reduced, ω, V0)
    D = h.Huu - ω^2 * h.Muu
    u_h = D \ (-h.Huϕ * V0)
    Q_h = -(only(h.Hϕu * u_h) + h.Hϕϕ * V0)
    Y_h = im * ω * Q_h / V0

    @test direct.displacement ≈ u_h
    @test direct.charge ≈ Q_h
    @test direct.current ≈ im * ω * Q_h
    @test direct.admittance ≈ Y_h
    @test direct.ω == ω
    @test direct.analysis === nothing
    @test direct.convention == :exp_iomega_t
    @test direct.solver_info.method == :backslash
    @test direct.solver_info.matrix_size == (3, 3)
    @test direct.solver_info.reduced_relative_residual < 1.0e-12
    @test direct.potential[partition.internal] ≈ direct.internal_potential
    @test direct.potential[partition.driven] == fill(V0, length(partition.driven))
    @test direct.potential[partition.grounded] == zeros(length(partition.grounded))
    @test solve_direct_voltage(reduced, ω, V0).admittance ≈ Y_h

    mechanical_dirichlet = DirichletDofs([1], [0.25])
    constrained_direct =
        solve_direct_voltage(reduced, ω, V0; mechanical_dirichlet)
    A = [
        Kuu - ω^2 * Muu Kuϕ[:, 1:1]
        Kϕu[1:1, :] Kϕϕ[1:1, 1:1]
    ]
    rhs = -vcat(Kuϕ[:, 2:2], Kϕϕ[1:1, 2:2])[:, 1] * V0
    manual_reduction = apply_dirichlet(A, rhs, mechanical_dirichlet)
    manual_x = reconstruct_solution(manual_reduction, manual_reduction.A \ manual_reduction.b)
    manual_u = manual_x[1:2]
    manual_ϕᵢ = manual_x[3:3]
    manual_charge = -(
        sum(reduced.KPu .* manual_u) +
        sum(reduced.KPi .* manual_ϕᵢ) +
        reduced.KPP * V0
    )

    @test constrained_direct.displacement ≈ manual_u
    @test constrained_direct.internal_potential ≈ manual_ϕᵢ
    @test constrained_direct.displacement[1] == 0.25
    @test constrained_direct.charge ≈ manual_charge
    @test constrained_direct.current ≈ im * ω * manual_charge
    @test constrained_direct.admittance ≈ im * ω * manual_charge / V0
    @test constrained_direct.solver_info.matrix_size == (2, 2)
    @test constrained_direct.solver_info.reduced_relative_residual < 1.0e-12
    @test solve_direct_voltage(reduced, ω, V0; mechanical_dirichlet).admittance ≈
        constrained_direct.admittance
    @test_throws ArgumentError solve_direct_voltage(
        reduced,
        ω,
        V0;
        mechanical_dirichlet=DirichletDofs([3], [0.0]),
    )
end

@testset "charge and current sign convention" begin
    C = 3.5e-9
    Kuu = ones(1, 1)
    Kuϕ = zeros(1, 1)
    Kϕu = zeros(1, 1)
    Kϕϕ = fill(-C, 1, 1)
    Muu = zeros(1, 1)
    system = KFormSystem(Kuu, Kuϕ, Kϕu, Kϕϕ, Muu)
    partition = HarmonicVoltageDofPartition(1, Int[], [1], Int[])
    reduced = electrode_reduced_k_form(system, partition)

    ω = 2π * 10_000.0
    V0 = 2.0
    solution = solve_direct_voltage(reduced, ω, V0)

    @test solution.charge ≈ C * V0
    @test solution.current ≈ im * ω * C * V0
    @test solution.admittance ≈ im * ω * C
end

@testset "FE direct admittance agrees with dense H-form" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (1, 2),
        Vec{2}((1.0e-3, 0.0)),
        Vec{2}((2.0e-3, 2.0e-3)),
    )
    dense_assembly = assemble_k_form_dense(grid, mat, kin, ip, qr)
    sparse_assembly = assemble_k_form_sparse(grid, mat, kin, ip, qr)
    partition = potential_partition(
        dense_assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )
    sparse_partition = potential_partition(
        sparse_assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )
    reduced = electrode_reduced_k_form(dense_assembly.system, partition)
    sparse_reduced = electrode_reduced_k_form(sparse_assembly.system, sparse_partition)
    h = h_form_dense(dense_assembly.system, partition)

    ω = 2π * 10_000.0
    V0 = 1.25
    direct = solve_direct_voltage(reduced, ω, V0)
    sparse_direct = solve_direct_voltage(sparse_reduced, ω, V0)
    sparse_A, _ = PiezoAcousticFEM.build_direct_voltage_system(sparse_reduced, ω, V0)
    D = h.Huu - ω^2 * h.Muu
    u_h = D \ (-h.Huϕ * V0)
    Q_h = -(only(h.Hϕu * u_h) + h.Hϕϕ * V0)
    Y_h = im * ω * Q_h / V0

    @test length(direct.internal_potential) == 2
    @test issparse(sparse_A)
    @test direct.displacement ≈ u_h
    @test direct.charge ≈ Q_h
    @test direct.admittance ≈ Y_h
    @test sparse_direct.displacement ≈ direct.displacement
    @test sparse_direct.internal_potential ≈ direct.internal_potential
    @test sparse_direct.potential ≈ direct.potential
    @test sparse_direct.charge ≈ direct.charge
    @test sparse_direct.current ≈ direct.current
    @test sparse_direct.admittance ≈ direct.admittance
end

@testset "FE direct solve with axis displacement constraint" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (1, 2),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 2.0e-3)),
    )
    assembly = assemble_k_form_sparse(grid, mat, kin, ip, qr)
    partition = potential_partition(
        assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )
    reduced = electrode_reduced_k_form(assembly.system, partition)
    axis_constraint =
        axis_radial_displacement_constraint(assembly; axis=AxisBoundary(FacetBoundary("left")))

    result = solve_direct_voltage(
        reduced,
        2π * 10_000.0,
        1.0;
        mechanical_dirichlet=axis_constraint,
    )

    @test length(partition.internal) == 2
    @test !isempty(axis_constraint.indices)
    @test result.displacement[axis_constraint.indices] == axis_constraint.values
    @test all(isfinite, result.displacement)
    @test all(isfinite, result.internal_potential)
    @test result.potential[partition.internal] ≈ result.internal_potential
    @test result.potential[partition.driven] == ones(length(partition.driven))
    @test result.potential[partition.grounded] == zeros(length(partition.grounded))
    @test isfinite(result.charge)
    @test isfinite(result.current)
    @test isfinite(result.admittance)
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
    reduction = reduce(assembled, analysis)
    result = solve(problem, analysis)

    @test problem.material_source === mat
    @test problem.loss isa Lossless
    @test problem.electrodes.signal isa FacetElectrode
    @test problem.electrodes.reference isa FacetElectrode
    @test only(problem.boundary_conditions.mechanical) isa AxisBoundary
    @test assembled.problem === problem
    @test assembled.material isa Piezo6mmAxi
    @test reduction.assembled === assembled
    @test reduction.reduced isa ElectrodeReducedKForm
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
    @test isfinite(result.solution.charge)
    @test isfinite(result.solution.current)
    @test isfinite(result.solution.admittance)
    @test result.problem === problem
    @test result.analysis === analysis
    @test result.solution.analysis === analysis
    @test result.solution.convention == analysis.convention
    explicit_solution = solve_direct_voltage(
        reduction.reduced,
        analysis.ω,
        analysis.voltage;
        mechanical_dirichlet=reduction.mechanical_dirichlet,
        analysis,
        convention=analysis.convention,
    )
    @test result.solution.admittance ≈ explicit_solution.admittance
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
    @test result.reduction.partition isa ShortCircuitDofPartition
    @test sort(vcat(result.reduction.partition.internal, result.reduction.partition.grounded)) ==
        collect(1:result.reduction.partition.nϕ)
    @test length(result.eigenvalues) == 4
    @test issorted(result.eigenvalues)
    @test minimum(result.eigenvalues) >= -1.0e-10 * maximum(abs, result.eigenvalues)
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
    filename = write_vtk(basename, grid, fields)

    @test isfile(filename)
    @test filesize(filename) > 0
    @test length(fields.displacement) == getnnodes(grid)
    @test length(fields.potential) == getnnodes(grid)
    @test length(fields.radius) == getnnodes(grid)
    rm(filename; force=true)
end

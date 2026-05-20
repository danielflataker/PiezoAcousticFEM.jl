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

    elem = PiezoAcousticFEM.piezo_element_matrices(cellvalues_u, cellvalues_ϕ, getcoordinates(cell), mat, kin)

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

@testset "axisymmetric geometry guards" begin
    kin = AxisymmetricRZ()
    mat = effective_material(PZT5A(), kin, Lossless())
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)

    @test_throws ArgumentError PiezoAcousticFEM.integration_weight(kin, Vec{2}((0.0, 0.0)), 1.0)
    @test_throws ArgumentError PiezoAcousticFEM.integration_weight(kin, Vec{2}((1.0, 0.0)), 0.0)
    @test_throws ArgumentError PiezoAcousticFEM.strain(kin, Vec{2}((1.0, 0.0)), zeros(2, 2), Vec{2}((0.0, 0.0)))

    axis_touching_grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((0.0, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    axis_touching_assembly = PiezoAcousticFEM.assemble_k_form_sparse(axis_touching_grid, mat, kin, ip, qr)
    @test all(isfinite, axis_touching_assembly.system.Kuu)

    negative_radius_grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((-1.0e-3, 0.0)),
        Vec{2}((1.0e-3, 1.0e-3)),
    )
    @test_throws ArgumentError PiezoAcousticFEM.assemble_k_form_sparse(negative_radius_grid, mat, kin, ip, qr)
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

    assembly = PiezoAcousticFEM.assemble_k_form_dense(one_cell_grid, mat, kin, ip, qr)

    cellvalues_u = CellValues(qr, ip)
    cellvalues_ϕ = CellValues(qr, ip)
    cell = first(CellIterator(one_cell_grid))
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_ϕ, cell)
    elem = PiezoAcousticFEM.piezo_element_matrices(cellvalues_u, cellvalues_ϕ, getcoordinates(cell), mat, kin)

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

    two_cell_assembly = PiezoAcousticFEM.assemble_k_form_dense(two_cell_grid, mat, kin, ip, qr)
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

    dense = PiezoAcousticFEM.assemble_k_form_dense(grid, mat, kin, ip, qr)
    sparse_assembly = PiezoAcousticFEM.assemble_k_form_sparse(grid, mat, kin, ip, qr)

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

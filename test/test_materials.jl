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

@testset "specialized scalar types survive assembly reductions" begin
    kin = AxisymmetricRZ()
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    grid = generate_grid(
        Quadrilateral,
        (1, 1),
        Vec{2}((1.0e-3, 0.0)),
        Vec{2}((2.0e-3, 1.0e-3)),
    )
    base = PZT5A()

    function scaled_material(scale)
        return Piezo6mmConstants(
            base.cᴱ .* scale,
            base.e .* scale,
            base.εˢ .* scale,
            base.ρ * scale,
        )
    end

    complex_material = scaled_material(1.0 + 0.02im)
    complex_assembly = assemble_k_form_dense(
        grid,
        effective_material(
            complex_material,
            kin,
            PhysicalLoss(MaterialAsGiven(), NoSystemDamping()),
        ),
        kin,
        ip,
        qr,
    )
    complex_partition = potential_partition(
        complex_assembly;
        driven=FacetElectrode("top"),
        grounded=FacetElectrode("bottom"),
    )
    complex_reduced = electrode_reduced_k_form(complex_assembly.system, complex_partition)

    @test eltype(complex_assembly.system.Kuu) <: Complex
    @test eltype(complex_reduced.KuP) <: Complex
    @test typeof(complex_reduced.KPP) <: Complex

    function reduced_trace(scale)
        material = effective_material(
            scaled_material(scale),
            kin,
            PhysicalLoss(MaterialAsGiven(), NoSystemDamping()),
        )
        assembly = assemble_k_form_dense(grid, material, kin, ip, qr)
        partition = potential_partition(
            assembly;
            driven=FacetElectrode("top"),
            grounded=FacetElectrode("bottom"),
        )
        reduced = electrode_reduced_k_form(assembly.system, partition)

        return sum(reduced.Kuu) + sum(reduced.KuP) + reduced.KPP + sum(reduced.Muu)
    end

    dual_value = reduced_trace(ForwardDiff.Dual(1.0, 1.0))
    @test dual_value isa ForwardDiff.Dual
    @test ForwardDiff.partials(dual_value, 1) != 0
    @test ForwardDiff.derivative(reduced_trace, 1.0) ≈ ForwardDiff.partials(dual_value, 1)
end


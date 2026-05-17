using StaticArrays: @SMatrix

@testset "PZT5A FEMP material 40 full Voigt contract" begin
    material = PZT5A()

    expected_c = @SMatrix [
        1.21000e11 7.54000e10 7.52000e10 0.0 0.0 0.0
        7.54000e10 1.21000e11 7.52000e10 0.0 0.0 0.0
        7.52000e10 7.52000e10 1.11000e11 0.0 0.0 0.0
        0.0 0.0 0.0 2.11000e10 0.0 0.0
        0.0 0.0 0.0 0.0 2.11000e10 0.0
        0.0 0.0 0.0 0.0 0.0 2.26000e10
    ]
    expected_e = @SMatrix [
        0.0 0.0 0.0 0.0 12.3 0.0
        0.0 0.0 0.0 12.3 0.0 0.0
        -5.4 -5.4 15.8 0.0 0.0 0.0
    ]
    expected_ε = @SMatrix [
        8.11026e-9 0.0 0.0
        0.0 8.11026e-9 0.0
        0.0 0.0 7.34882e-9
    ]

    @test material isa PiezoMaterial
    @test material isa AbstractPiezoMaterial
    @test material.cᴱ == expected_c
    @test material.e == expected_e
    @test material.εˢ == expected_ε
    @test material.ρ == 7750.0
    @test material.cᴱ[6, 6] == 2.26e10
    @test material.cᴱ[6, 6] != (material.cᴱ[1, 1] - material.cᴱ[1, 2]) / 2
end

@testset "AxisymmetricRZ material reduction basis" begin
    c = @SMatrix [
        11 12 13 14 15 16
        21 22 23 24 25 26
        31 32 33 34 35 36
        41 42 43 44 45 46
        51 52 53 54 55 56
        61 62 63 64 65 66
    ]
    e = @SMatrix [
        101 102 103 104 105 106
        201 202 203 204 205 206
        301 302 303 304 305 306
    ]
    ε = @SMatrix [
        1001 1002 1003
        2001 2002 2003
        3001 3002 3003
    ]
    material = PiezoMaterial(c, e, ε, 7000.0)

    reduced = reduce_material(material, AxisymmetricRZ())

    expected_c = @SMatrix [
        11 12 13 15
        21 22 23 25
        31 32 33 35
        51 52 53 55
    ]
    expected_e = @SMatrix [
        101 102 103 105
        301 302 303 305
    ]
    expected_ε = @SMatrix [
        1001 1003
        3001 3003
    ]

    @test reduced isa AxisymmetricRZPiezoMaterial
    @test reduced isa AbstractReducedPiezoMaterial
    @test reduced.cᴱ == expected_c
    @test reduced.e == expected_e
    @test reduced.εˢ == expected_ε
    @test reduced.ρ == material.ρ
    @test reduced.cᴱ[4, 4] == material.cᴱ[5, 5]
end

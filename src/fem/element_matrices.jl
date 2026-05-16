# Local element matrices for axisymmetric piezoelectric vacuum elements.

"""
    PiezoElementMatrices

Lokale elementmatriser for et aksesymmetrisk piezoelement.

DOF-rekkefølgen antas å være:
- mekanisk: `u = [u_r, u_z]`
- elektrisk: `ϕ`

Matrisene er blokkene

    Kuu, Kuϕ, Kϕu, Kϕϕ, Muu

før de eventuelt settes inn i en global matrise.
"""
struct PiezoElementMatrices{KUU,KUP,KPU,KPP,MUU}
    Kuu::KUU
    Kuϕ::KUP
    Kϕu::KPU
    Kϕϕ::KPP
    Muu::MUU
end


function zero_piezo_element_matrices(
    ::Type{TK},
    ::Type{TE},
    ::Type{Tε},
    ::Type{Tρ},
    nᵤ::Integer,
    nᵩ::Integer,
) where {TK,TE,Tε,Tρ}
    return PiezoElementMatrices(
        zeros(TK, 2nᵤ, 2nᵤ),
        zeros(TE, 2nᵤ, nᵩ),
        zeros(TE, nᵩ, 2nᵤ),
        zeros(Tε, nᵩ, nᵩ),
        zeros(Tρ, 2nᵤ, 2nᵤ),
    )
end


function displacement_shape_value(cellvalues, q, a)
    N = shape_value(cellvalues, q, a)

    return @SMatrix [
        N 0
        0 N
    ]
end


function displacement_shape_gradient(cellvalues, q, a)
    ∇N = shape_gradient(cellvalues, q, a)

    # ∇u[i,j] = ∂u_i/∂x_j
    #
    # For node a:
    # u_r = N_a * u_ra
    # u_z = N_a * u_za
    #
    # Bidrag fra u_ra:
    # ∇u = [∂N/∂r  ∂N/∂z
    #       0       0]
    #
    # Bidrag fra u_za:
    # ∇u = [0       0
    #       ∂N/∂r  ∂N/∂z]
    return ∇N
end


function Bu_matrix(kin::AxisymmetricRZ, cellvalues, q, x)
    nᵤ = getnbasefunctions(cellvalues)
    T = eltype(x)
    r = x[1]
    r > zero(r) || throw(ArgumentError("axisymmetric B_u radius must be positive, got r=$r"))

    Bu = zeros(T, 4, 2nᵤ)

    for a in 1:nᵤ
        N = shape_value(cellvalues, q, a)
        ∇N = shape_gradient(cellvalues, q, a)

        dNdr = ∇N[1]
        dNdz = ∇N[2]

        colᵣ = 2a - 1
        colz = 2a

        # S = [S_rr, S_θθ, S_zz, γ_rz]
        Bu[1, colᵣ] = dNdr
        Bu[2, colᵣ] = N / r
        Bu[4, colᵣ] = dNdz

        Bu[3, colz] = dNdz
        Bu[4, colz] = dNdr
    end

    return Bu
end


function Bϕ_matrix(kin::AxisymmetricRZ, cellvalues, q)
    nᵩ = getnbasefunctions(cellvalues)
    T = typeof(shape_value(cellvalues, q, 1))

    Bϕ = zeros(T, 2, nᵩ)

    for a in 1:nᵩ
        ∇N = shape_gradient(cellvalues, q, a)

        # E = -∇ϕ.
        #
        # Her lar vi Bϕ være gradientoperatoren ∇ϕ.
        # Fortegn håndteres i elementintegralene.
        Bϕ[1, a] = ∇N[1]
        Bϕ[2, a] = ∇N[2]
    end

    return Bϕ
end


function Nu_matrix(cellvalues, q)
    nᵤ = getnbasefunctions(cellvalues)
    T = typeof(shape_value(cellvalues, q, 1))

    Nu = zeros(T, 2, 2nᵤ)

    for a in 1:nᵤ
        N = shape_value(cellvalues, q, a)

        Nu[1, 2a-1] = N
        Nu[2, 2a] = N
    end

    return Nu
end


"""
    piezo_element_matrices(cellvalues_u, cellvalues_ϕ, xᵉ, material, formulation)

Lager lokale elementmatriser for ett aksesymmetrisk piezoelement.

Forutsetter at `reinit!` allerede er kalt på cellvalues-objektene.
`xᵉ` er elementets fysiske nodekoordinater, typisk fra `getcoordinates(cell)`.
"""
function piezo_element_matrices(cellvalues_u, cellvalues_ϕ, xᵉ,
    material::Piezo6mmAxi,
    formulation::AxisymmetricRZ)

    validate_element_geometry(formulation, xᵉ)

    nᵤ = getnbasefunctions(cellvalues_u)
    nᵩ = getnbasefunctions(cellvalues_ϕ)
    getnquadpoints(cellvalues_u) == getnquadpoints(cellvalues_ϕ) ||
        throw(DimensionMismatch("cellvalues_u and cellvalues_ϕ must use the same quadrature rule"))

    TK, TE, Tε, Tρ = k_form_block_types_from_coordinates(material, xᵉ)
    elem = zero_piezo_element_matrices(TK, TE, Tε, Tρ, nᵤ, nᵩ)

    for q in 1:getnquadpoints(cellvalues_u)
        x = spatial_coordinate(cellvalues_u, q, xᵉ)
        dV = integration_weight(formulation, x, getdetJdV(cellvalues_u, q))

        Bu = Bu_matrix(formulation, cellvalues_u, q, x)
        Bϕ = Bϕ_matrix(formulation, cellvalues_ϕ, q)
        Nu = Nu_matrix(cellvalues_u, q)

        elem.Kuu .+= transpose(Bu) * material.cᴱ * Bu * dV

        # E = -∇ϕ.
        # Fra T = cᴱS - eᵀE = cᴱS + eᵀ∇ϕ.
        # Dette gir fortegnkonvensjonen under for den mekaniske ligningen.
        elem.Kuϕ .+= transpose(Bu) * transpose(material.e) * Bϕ * dV

        # Fra D = eS + εˢE = eS - εˢ∇ϕ.
        # Dette matcher Kocbach/KLV sin K-form med negativ dielektrisk blokk.
        elem.Kϕu .+= transpose(Bϕ) * material.e * Bu * dV
        elem.Kϕϕ .-= transpose(Bϕ) * material.εˢ * Bϕ * dV

        elem.Muu .+= material.ρ * (transpose(Nu) * Nu) * dV
    end

    return elem
end


function k_form_block_types_from_coordinates(material::Piezo6mmAxi, xᵉ)
    Tx = coordinate_eltype(xᵉ)

    return (
        promote_type(eltype(material.cᴱ), Tx),
        promote_type(eltype(material.e), Tx),
        promote_type(eltype(material.εˢ), Tx),
        promote_type(typeof(material.ρ), Tx),
    )
end


coordinate_eltype(xᵉ) = eltype(first(xᵉ))


function validate_element_geometry(::AxisymmetricRZ, xᵉ)
    for x in xᵉ
        r = x[1]
        r >= zero(r) || throw(ArgumentError("axisymmetric element node radius must be nonnegative, got r=$r"))
    end

    return nothing
end

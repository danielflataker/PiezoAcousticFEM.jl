"""
    ElementDerivedFieldOutput(...)

Derived physical fields evaluated inside one finite element at explicitly
provided reference coordinates. These values are local postprocessing samples,
not nodal fields and not projected VTK fields.
"""
struct ElementDerivedFieldOutput{R,X,U,P,GU,GP,S,E,T,D,M}
    cellid::Int
    reference_points::R
    coordinates::X
    displacement::U
    potential::P
    displacement_gradient::GU
    potential_gradient::GP
    strain::S
    electric_field::E
    stress::T
    electric_displacement::D
    metadata::M
end

"""
    PhysicalGridDerivedFieldOutput(...)

Derived physical fields evaluated on a rectangular physical `(r, z)` sample
grid. Samples are stored in `z`-major order with `r` varying fastest, matching
the order produced by `sample_physical_grid`.
"""
struct PhysicalGridDerivedFieldOutput{R,Z,X,C,Ξ,U,P,GU,GP,S,E,T,D,M}
    r_coordinates::R
    z_coordinates::Z
    coordinates::X
    cellids::C
    reference_points::Ξ
    displacement::U
    potential::P
    displacement_gradient::GU
    potential_gradient::GP
    strain::S
    electric_field::E
    stress::T
    electric_displacement::D
    metadata::M
end

"""
    sample_physical_grid(assembled, fields, r_coordinates, z_coordinates)

Evaluate derived fields on a rectangular physical `(r, z)` grid across the
mesh. Every requested point must lie inside a cell; this function is a sampling
primitive and deliberately does not extrapolate, project, or choose
visualization phase transforms.
"""
function sample_physical_grid(
    assembled::AssembledPiezoProblem,
    fields::CompactFieldDofs,
    r_coordinates,
    z_coordinates,
)
    return _sample_physical_grid(
        assembled.assembly,
        assembled.material,
        assembled.problem.formulation,
        assembled.problem.interpolation,
        fields,
        r_coordinates,
        z_coordinates;
        metadata=merge_physical_grid_metadata(fields.metadata),
    )
end


function _sample_physical_grid(
    assembly::KFormAssembly,
    material::AxisymmetricRZPiezoMaterial,
    formulation::AxisymmetricRZ,
    interpolation,
    fields::CompactFieldDofs,
    r_coordinates,
    z_coordinates;
    metadata=merge_physical_grid_metadata(fields.metadata),
)
    r_values = validate_physical_grid_axis(r_coordinates, :r_coordinates)
    z_values = validate_physical_grid_axis(z_coordinates, :z_coordinates)
    points = [Vec{2}((r, z)) for z in z_values for r in r_values]
    grid = ferrite_grid(assembly)
    point_handler = PointEvalHandler(grid, points; warn=false)
    samples = map(points, point_handler.cells, point_handler.local_coords) do point, cellid, ξ
        cellid === nothing &&
            throw(ArgumentError("physical sample point $point is outside the mesh"))
        ξ === nothing &&
            throw(ArgumentError("physical sample point $point could not be mapped to a cell"))

        return first_element_derived_sample(
            evaluate_derived_fields(
                assembly,
                material,
                formulation,
                interpolation,
                fields,
                cellid,
                [ξ];
                metadata,
            )
        )
    end

    return PhysicalGridDerivedFieldOutput(
        r_values,
        z_values,
        points,
        Int.(point_handler.cells),
        point_handler.local_coords,
        [sample.displacement for sample in samples],
        [sample.potential for sample in samples],
        [sample.displacement_gradient for sample in samples],
        [sample.potential_gradient for sample in samples],
        [sample.strain for sample in samples],
        [sample.electric_field for sample in samples],
        [sample.stress for sample in samples],
        [sample.electric_displacement for sample in samples],
        metadata,
    )
end


function first_element_derived_sample(output::ElementDerivedFieldOutput)
    return (
        displacement=only(output.displacement),
        potential=only(output.potential),
        displacement_gradient=only(output.displacement_gradient),
        potential_gradient=only(output.potential_gradient),
        strain=only(output.strain),
        electric_field=only(output.electric_field),
        stress=only(output.stress),
        electric_displacement=only(output.electric_displacement),
    )
end

"""
    sample_element_fields(..., cellid; samples_per_axis=5)

Evaluate derived fields on a tensor-product lattice of reference points inside
one element. This is a per-element visualization/sample primitive; it does not
locate arbitrary physical points or project values to nodes.
"""
function sample_element_fields(
    assembled::AssembledPiezoProblem,
    fields::CompactFieldDofs,
    cellid::Integer;
    samples_per_axis=5,
)
    return sample_element_fields(
        assembled.assembly,
        assembled.material,
        assembled.problem.formulation,
        assembled.problem.interpolation,
        fields,
        cellid;
        samples_per_axis,
        metadata=merge_sampled_field_metadata(fields.metadata, samples_per_axis),
    )
end


function sample_element_fields(
    result::HarmonicVoltageResult,
    cellid::Integer;
    samples_per_axis=5,
)
    fields = reconstruct_field_dofs(result)

    return sample_element_fields(
        result.assembled.assembly,
        result.assembled.material,
        result.problem.formulation,
        result.problem.interpolation,
        fields,
        cellid;
        samples_per_axis,
        metadata=merge_sampled_field_metadata(
            merge(
                fields.metadata,
                (
                    analysis=:harmonic_voltage,
                    harmonic_convention=result.analysis.convention,
                    quantity_interpretation="complex amplitude",
                ),
            ),
            samples_per_axis,
        ),
    )
end


function sample_element_fields(
    assembly::KFormAssembly,
    material::AxisymmetricRZPiezoMaterial,
    formulation::AxisymmetricRZ,
    interpolation,
    fields::CompactFieldDofs,
    cellid::Integer;
    samples_per_axis=5,
    metadata=merge_sampled_field_metadata(fields.metadata, samples_per_axis),
)
    reference_points = reference_sample_points(formulation, samples_per_axis)

    return evaluate_derived_fields(
        assembly,
        material,
        formulation,
        interpolation,
        fields,
        cellid,
        reference_points;
        metadata,
    )
end

"""
    evaluate_derived_fields(assembled, fields, cellid, reference_points)
    evaluate_derived_fields(result, cellid, reference_points)

Evaluate axisymmetric derived fields inside one element at local reference
coordinates. The returned quantities include interpolated primary values,
gradients, strain, electric field, stress, and electric displacement.
"""
function evaluate_derived_fields(
    assembled::AssembledPiezoProblem,
    fields::CompactFieldDofs,
    cellid::Integer,
    reference_points,
)
    return evaluate_derived_fields(
        assembled.assembly,
        assembled.material,
        assembled.problem.formulation,
        assembled.problem.interpolation,
        fields,
        cellid,
        reference_points;
        metadata=merge_derived_field_metadata(fields.metadata),
    )
end


function evaluate_derived_fields(result::HarmonicVoltageResult, cellid::Integer, reference_points)
    fields = reconstruct_field_dofs(result)

    return evaluate_derived_fields(
        result.assembled.assembly,
        result.assembled.material,
        result.problem.formulation,
        result.problem.interpolation,
        fields,
        cellid,
        reference_points;
        metadata=merge_derived_field_metadata(
            merge(
                fields.metadata,
                (
                    analysis=:harmonic_voltage,
                    harmonic_convention=result.analysis.convention,
                    quantity_interpretation="complex amplitude",
                ),
            ),
        ),
    )
end


function evaluate_derived_fields(
    assembly::KFormAssembly,
    material::AxisymmetricRZPiezoMaterial,
    formulation::AxisymmetricRZ,
    interpolation,
    fields::CompactFieldDofs,
    cellid::Integer,
    reference_points;
    metadata=merge_derived_field_metadata(fields.metadata),
)
    grid = ferrite_grid(assembly)
    1 <= cellid <= getncells(grid) ||
        throw(ArgumentError("cellid must be in 1:$(getncells(grid)), got $cellid"))

    points = collect(reference_points)
    cell = CellCache(assembly.dofhandler)
    reinit!(cell, Int(cellid))
    coordinates = getcoordinates(cell)
    cell_dofs = celldofs(cell)
    u_range = piezo_field_dof_range(assembly.dofhandler, :u)
    ϕ_range = piezo_field_dof_range(assembly.dofhandler, :ϕ)

    u_local = fields.displacement[[
        compact_displacement_dof(assembly.dofmap, dof) for dof in cell_dofs[u_range]
    ]]
    ϕ_local = fields.potential[[
        compact_potential_dof(assembly.dofmap, dof) for dof in cell_dofs[ϕ_range]
    ]]

    pointvalues_u = PointValues(interpolation)
    pointvalues_ϕ = PointValues(interpolation)

    samples = [
        evaluate_derived_field_point(
            pointvalues_u,
            pointvalues_ϕ,
            coordinates,
            u_local,
            ϕ_local,
            material,
            formulation,
            ξ,
        )
        for ξ in points
    ]

    return ElementDerivedFieldOutput(
        Int(cellid),
        points,
        [sample.coordinate for sample in samples],
        [sample.displacement for sample in samples],
        [sample.potential for sample in samples],
        [sample.displacement_gradient for sample in samples],
        [sample.potential_gradient for sample in samples],
        [sample.strain for sample in samples],
        [sample.electric_field for sample in samples],
        [sample.stress for sample in samples],
        [sample.electric_displacement for sample in samples],
        metadata,
    )
end


function evaluate_derived_field_point(
    pointvalues_u,
    pointvalues_ϕ,
    coordinates,
    u_local,
    ϕ_local,
    material,
    formulation,
    ξ,
)
    reinit!(pointvalues_u, coordinates, ξ)
    reinit!(pointvalues_ϕ, coordinates, ξ)

    x = point_coordinate(pointvalues_u, coordinates)
    Nu = Nu_matrix(pointvalues_u, 1)
    Bϕ = Bϕ_matrix(formulation, pointvalues_ϕ, 1)

    u_values = Nu * u_local
    u = Vec{2}((u_values[1], u_values[2]))
    ∇u = displacement_gradient_from_pointvalues(pointvalues_u, u_local)
    ∇ϕ_values = Bϕ * ϕ_local
    ∇ϕ = Vec{2}((∇ϕ_values[1], ∇ϕ_values[2]))
    S = strain(formulation, u, ∇u, x)
    E = electric_field(formulation, ∇ϕ)

    return (
        coordinate=x,
        displacement=u,
        potential=dot_shape_values(pointvalues_ϕ, ϕ_local),
        displacement_gradient=∇u,
        potential_gradient=∇ϕ,
        strain=S,
        electric_field=E,
        stress=stress(material, S, E),
        electric_displacement=electric_displacement(material, S, E),
    )
end


function merge_derived_field_metadata(metadata)
    return merge(
        metadata,
        (
            output_kind=:element_derived_fields,
            evaluation=:reference_points,
            coordinate_units="m",
            potential_units="V",
            strain_units="1",
            electric_field_units="V/m",
            stress_units="Pa",
            electric_displacement_units="C/m^2",
        ),
    )
end


function merge_sampled_field_metadata(metadata, samples_per_axis)
    return merge(
        merge_derived_field_metadata(metadata),
        (
            output_kind=:element_sampled_fields,
            evaluation=:reference_lattice,
            samples_per_axis=normalize_samples_per_axis(samples_per_axis),
        ),
    )
end


function merge_physical_grid_metadata(metadata)
    return merge(
        merge_derived_field_metadata(metadata),
        (
            output_kind=:physical_grid_derived_fields,
            evaluation=:physical_grid,
            sample_order=:z_major_r_fastest,
        ),
    )
end


function reference_sample_points(::AxisymmetricRZ, samples_per_axis)
    nr, nz = normalize_samples_per_axis(samples_per_axis)
    ξr = range(-1.0, 1.0; length=nr)
    ξz = range(-1.0, 1.0; length=nz)

    return [Vec{2}((r, z)) for z in ξz for r in ξr]
end


function validate_physical_grid_axis(coordinates, name::Symbol)
    values = collect(coordinates)
    !isempty(values) || throw(ArgumentError("$name must contain at least one coordinate"))
    all(isfinite, values) || throw(ArgumentError("$name must contain only finite coordinates"))

    return values
end


function normalize_samples_per_axis(samples_per_axis::Integer)
    samples_per_axis >= 2 ||
        throw(ArgumentError("samples_per_axis must be at least 2, got $samples_per_axis"))

    return (Int(samples_per_axis), Int(samples_per_axis))
end


function normalize_samples_per_axis(samples_per_axis)
    length(samples_per_axis) == 2 ||
        throw(ArgumentError("samples_per_axis must be an integer or a length-2 tuple/vector"))
    nr, nz = samples_per_axis
    nr >= 2 && nz >= 2 ||
        throw(ArgumentError("samples_per_axis entries must be at least 2, got $samples_per_axis"))

    return (Int(nr), Int(nz))
end


function dot_shape_values(pointvalues, local_values)
    value = zero(eltype(local_values))
    for a in 1:getnbasefunctions(pointvalues)
        value += shape_value(pointvalues, 1, a) * local_values[a]
    end

    return value
end


function point_coordinate(pointvalues, coordinates)
    x = zero(first(coordinates))
    for a in 1:getnbasefunctions(pointvalues)
        x += shape_value(pointvalues, 1, a) * coordinates[a]
    end

    return x
end


function displacement_gradient_from_pointvalues(pointvalues, u_local)
    T = eltype(u_local)
    ∇u = zero_displacement_gradient(T)

    for a in 1:getnbasefunctions(pointvalues)
        ∇N = shape_gradient(pointvalues, 1, a)
        uᵣ = u_local[2a-1]
        uz = u_local[2a]
        ∇u += @SMatrix [
            uᵣ * ∇N[1] uᵣ * ∇N[2]
            uz * ∇N[1] uz * ∇N[2]
        ]
    end

    return ∇u
end


zero_displacement_gradient(::Type{T}) where {T} = @SMatrix zeros(T, 2, 2)

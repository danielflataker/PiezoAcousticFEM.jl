"""
    write_vtk(filename, grid, fields)

Write reconstructed nodal fields to a VTK `.vtu` file. Complex fields are
exported as explicit real/imag/abs/phase components. Field-data metadata
records basic units and output interpretation.
"""
function write_vtk(filename::AbstractString, grid, fields; metadata=nothing)
    VTKGridFile(String(filename), grid) do vtk
        write_displacement_data(vtk, fields.displacement)
        write_potential_data(vtk, fields.potential)
        write_node_data(vtk, fields.radius, "radius")
        write_vtk_metadata(vtk, vtk_metadata(fields; metadata))
    end

    return String(filename) * ".vtu"
end

"""
    write_vtk(filename, grid, element_outputs::AbstractVector{<:ElementDerivedFieldOutput})

Write element-local derived fields as VTK cell data. Multiple samples inside
each element are averaged to one value per cell; no nodal projection or
smoothing is performed.
"""
function write_vtk(
    filename::AbstractString,
    grid,
    element_outputs::AbstractVector{<:ElementDerivedFieldOutput};
    metadata=nothing,
)
    cell_outputs = derived_outputs_by_cell(grid, element_outputs)

    VTKGridFile(String(filename), grid) do vtk
        write_derived_cell_data(vtk, cell_outputs)
        write_vtk_metadata(vtk, vtk_derived_cell_metadata(element_outputs; metadata))
    end

    return String(filename) * ".vtu"
end

"""
    write_vtk(filename, result::HarmonicVoltageResult)

Reconstruct nodal fields from a harmonic voltage result and write them to a
VTK `.vtu` file.
"""
function write_vtk(filename::AbstractString, result::HarmonicVoltageResult)
    grid = result.problem.grid
    fields = reconstruct_fields(result)

    return write_vtk(filename, grid, fields)
end

"""
    write_vtk(filename, result::ShortCircuitModalResult; mode_index = 1)

Reconstruct nodal fields for one short-circuit mode and write them to a VTK
`.vtu` file.
"""
function write_vtk(filename::AbstractString, result::ShortCircuitModalResult; mode_index::Integer=1)
    grid = result.assembled.problem.grid
    fields = reconstruct_fields(result, mode_index)

    return write_vtk(filename, grid, fields)
end

function vtk_metadata(fields; metadata=nothing)
    entries = Dict{String,String}(
        "piezoacousticfem_output_kind" => "nodal_fields",
        "piezoacousticfem_displacement_unit" => "m",
        "piezoacousticfem_potential_unit" => "V",
        "piezoacousticfem_radius_unit" => "m",
        "piezoacousticfem_phase_unit" => "rad",
        "piezoacousticfem_component_convention" =>
            "complex fields are exported as real, imaginary, absolute-value, and phase nodal arrays",
    )

    if hasproperty(fields, :metadata)
        merge_vtk_metadata!(entries, fields.metadata)
    end

    metadata === nothing && return entries

    merge_vtk_metadata!(entries, metadata)

    return entries
end

function merge_vtk_metadata!(entries, metadata)
    for (key, value) in pairs(metadata)
        vtk_key = key === :normalization ? :modal_normalization : key
        entries["piezoacousticfem_$(vtk_key)"] = string(value)
    end

    return entries
end

function write_vtk_metadata(vtk, metadata)
    for name in sort!(collect(keys(metadata)))
        vtk.vtk[name, VTKFieldData()] = metadata[name]
    end

    return vtk
end

function derived_outputs_by_cell(grid, element_outputs)
    ncells = getncells(grid)
    length(element_outputs) == ncells ||
        throw(ArgumentError("derived VTK cell output requires one element output per grid cell"))

    by_cell = Vector{eltype(element_outputs)}(undef, ncells)
    seen = falses(ncells)
    for output in element_outputs
        1 <= output.cellid <= ncells ||
            throw(ArgumentError("derived output cellid must be in 1:$ncells, got $(output.cellid)"))
        !seen[output.cellid] ||
            throw(ArgumentError("derived output contains duplicate cellid $(output.cellid)"))

        by_cell[output.cellid] = output
        seen[output.cellid] = true
    end
    all(seen) || throw(ArgumentError("derived output must include every grid cell exactly once"))

    return by_cell
end

function vtk_derived_cell_metadata(element_outputs; metadata=nothing)
    entries = Dict{String,String}()

    !isempty(element_outputs) && hasproperty(first(element_outputs), :metadata) &&
        merge_vtk_metadata!(entries, first(element_outputs).metadata)
    metadata !== nothing && merge_vtk_metadata!(entries, metadata)

    entries["piezoacousticfem_output_kind"] = "derived_cell_fields"
    entries["piezoacousticfem_evaluation"] = "cell_average"
    entries["piezoacousticfem_cell_data_policy"] =
        "element-local samples are averaged per cell without nodal projection or smoothing"
    return entries
end

function write_derived_cell_data(vtk, cell_outputs)
    write_cell_components(vtk, averaged_cell_field(cell_outputs, :displacement), "derived_displacement", ("r", "z"))
    write_cell_scalar_data(vtk, averaged_cell_field(cell_outputs, :potential), "derived_potential")
    write_cell_components(vtk, averaged_cell_field(cell_outputs, :strain), "derived_strain", ("rr", "thetatheta", "zz", "rz"))
    write_cell_components(vtk, averaged_cell_field(cell_outputs, :electric_field), "derived_electric_field", ("r", "z"))
    write_cell_components(vtk, averaged_cell_field(cell_outputs, :stress), "derived_stress", ("rr", "thetatheta", "zz", "rz"))

    return vtk
end

averaged_cell_field(cell_outputs, field) = cell_average.(getproperty.(cell_outputs, field))

function cell_average(values)
    !isempty(values) || throw(ArgumentError("derived cell output must contain at least one sample"))

    total = sum(values)
    return total / length(values)
end

function write_cell_scalar_data(vtk, data, name)
    if eltype(data) <: Complex
        write_cell_data(vtk, real.(data), "$(name)_real")
        write_cell_data(vtk, imag.(data), "$(name)_imag")
        write_cell_data(vtk, abs.(data), "$(name)_abs")
        write_cell_data(vtk, angle.(data), "$(name)_phase")
    else
        write_cell_data(vtk, data, name)
    end

    return vtk
end

function write_cell_components(vtk, data, name, component_names)
    for (component, component_name) in enumerate(component_names)
        write_cell_scalar_data(vtk, [x[component] for x in data], "$(name)_$(component_name)")
    end

    return vtk
end

function write_displacement_data(vtk, displacement)
    if is_complex_vec_data(displacement)
        displacement_real = map_vec(real, displacement)
        displacement_imag = map_vec(imag, displacement)
        displacement_abs = map_vec(abs, displacement)
        u_r = [u[1] for u in displacement]
        u_z = [u[2] for u in displacement]

        write_node_data(vtk, displacement_real, "displacement_real")
        write_node_data(vtk, displacement_imag, "displacement_imag")
        write_node_data(vtk, displacement_abs, "displacement_abs")
        write_node_data(vtk, real.(u_r), "u_r_real")
        write_node_data(vtk, imag.(u_r), "u_r_imag")
        write_node_data(vtk, abs.(u_r), "u_r_abs")
        write_node_data(vtk, angle.(u_r), "u_r_phase")
        write_node_data(vtk, real.(u_z), "u_z_real")
        write_node_data(vtk, imag.(u_z), "u_z_imag")
        write_node_data(vtk, abs.(u_z), "u_z_abs")
        write_node_data(vtk, angle.(u_z), "u_z_phase")
        write_node_data(vtk, [norm(u) for u in displacement_abs], "u_magnitude_abs")
    else
        u_r = [u[1] for u in displacement]
        u_z = [u[2] for u in displacement]
        u_magnitude = [norm(u) for u in displacement]

        write_node_data(vtk, displacement, "displacement")
        write_node_data(vtk, u_r, "u_r")
        write_node_data(vtk, u_z, "u_z")
        write_node_data(vtk, u_magnitude, "u_magnitude")
    end

    return vtk
end


function write_potential_data(vtk, potential)
    if eltype(potential) <: Complex
        write_node_data(vtk, real.(potential), "phi_real")
        write_node_data(vtk, imag.(potential), "phi_imag")
        write_node_data(vtk, abs.(potential), "phi_abs")
        write_node_data(vtk, angle.(potential), "phi_phase")
    else
        write_node_data(vtk, potential, "phi")
    end

    return vtk
end


is_complex_vec_data(data) = !isempty(data) && eltype(first(data)) <: Complex

map_vec(f, data) = [Vec{2}((f(u[1]), f(u[2]))) for u in data]

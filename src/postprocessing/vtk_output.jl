"""
    write_vtk(filename, grid, fields)

Write reconstructed nodal fields to a VTK `.vtu` file. Complex fields are
exported as explicit real/imag/abs/phase components.
"""
function write_vtk(filename::AbstractString, grid, fields)
    VTKGridFile(String(filename), grid) do vtk
        write_displacement_data(vtk, fields.displacement)
        write_potential_data(vtk, fields.potential)
        write_node_data(vtk, fields.radius, "radius")
    end

    return String(filename) * ".vtu"
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

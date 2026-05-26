using Ferrite
using PiezoAcousticFEM

"""
Small axisymmetric PZT-5A disk harmonic-voltage example.

The script solves a voltage-driven vacuum disk model, evaluates the drive
terminal admittance, and writes nodal displacement/potential fields to VTK.
"""

radius = 4.0e-3
thickness = 6.0e-3
frequency = 10.0e3
voltage = 1.0
output_basename = joinpath(tempdir(), "disk_harmonic_voltage")

grid = generate_grid(
    Quadrilateral,
    (8, 12),
    Vec{2}((0.0, 0.0)),
    Vec{2}((radius, thickness)),
)

problem = PiezoProblem(
    grid,
    PZT5A(),
    AxisymmetricRZ(),
    Lagrange{RefQuadrilateral,1}(),
    QuadratureRule{RefQuadrilateral}(2);
    electrodes=TwoTerminalElectrodes(FacetElectrode("top"), FacetElectrode("bottom")),
    boundary_conditions=AxisymmetricBoundaryConditions((AxisBoundary(FacetBoundary("left")),)),
    loss=Lossless(),
)
analysis = HarmonicVoltageAnalysis(2π * frequency, voltage, :exp_iomega_t)

assembled = assemble(problem)
prepared = prepare_analysis(assembled, analysis)

result = solve(prepared, analysis)
admittance = result.solution.drive_terminal_admittance
vtk_file = write_vtk(output_basename, result)

if abspath(PROGRAM_FILE) == @__FILE__
    println("admittance = ", admittance)
    println("wrote ", vtk_file)
end

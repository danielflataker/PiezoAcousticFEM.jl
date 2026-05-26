# PiezoAcousticFEM.jl

module PiezoAcousticFEM

using Ferrite
using LinearAlgebra
using SparseArrays
using StaticArrays
using WriteVTK

include("physics/formulations.jl")
include("physics/materials.jl")
include("physics/constitutive.jl")

include("fem/element_matrices.jl")
include("fem/systems.jl")

include("boundary_conditions/electrodes.jl")
include("boundary_conditions/boundaries.jl")
include("boundary_conditions/electrical_partition.jl")
include("boundary_conditions/constraints.jl")

include("fem/ferrite_backend.jl")
include("fem/dof_map.jl")
include("fem/dofhandler.jl")
include("fem/assembly.jl")
include("fem/hform_dense.jl")

include("boundary_conditions/electrode_reduction.jl")

include("analysis/direct_voltage.jl")
include("analysis/harmonic_voltage.jl")
include("analysis/modal.jl")
include("analysis/validation.jl")

include("postprocessing/reconstruction.jl")
include("postprocessing/derived_fields.jl")
include("postprocessing/observables.jl")
include("postprocessing/vtk_output.jl")

export AbstractFormulation, AxisymmetricRZ
export AbstractMaterial, AbstractPiezoMaterial, AbstractReducedPiezoMaterial
export PiezoMaterial, AxisymmetricRZPiezoMaterial, PZT5A, reduce_material, effective_material
export Lossless, PhysicalLoss, RealMaterial, MaterialAsGiven, PiezoComplexMaterialLoss
export NoSystemDamping
export FacetElectrode, TwoTerminalElectrodes
export FacetBoundary, AxisBoundary, AxisymmetricBoundaryConditions
export PiezoProblem, HarmonicVoltageAnalysis, ShortCircuitModalAnalysis
export AbstractLinearSolverConfig, BackslashSolver, FactorizedDirectSolver
export solve, assemble, prepare_analysis, validate
export HarmonicVoltageResult, HarmonicSweepPointResult, FrequencySweepResult, ShortCircuitModalResult
export AbstractObservable, AdmittanceObservable, ChargeObservable, CurrentObservable
export ModalFrequenciesObservable, evaluate
export CompactFieldDofs, NodalFieldOutput, ElementDerivedFieldOutput, PhysicalGridDerivedFieldOutput
export reconstruct_field_dofs, reconstruct_fields, evaluate_derived_fields, sample_element_fields
export sample_physical_grid, write_vtk

end

# PiezoAcousticFEM.jl

module PiezoAcousticFEM

using Ferrite
using LinearAlgebra
using SparseArrays
using StaticArrays
using WriteVTK

import Base: reduce

include("physics/formulations.jl")
include("physics/materials.jl")
include("physics/constitutive.jl")

include("fem/element_matrices.jl")
include("fem/systems.jl")

include("boundary_conditions/electrodes.jl")
include("boundary_conditions/boundaries.jl")
include("boundary_conditions/electrical_partition.jl")
include("boundary_conditions/constraints.jl")

include("fem/dof_map.jl")
include("fem/dofhandler.jl")
include("fem/assembly.jl")
include("fem/hform_dense.jl")

include("boundary_conditions/electrode_reduction.jl")

include("analysis/direct_voltage.jl")
include("analysis/harmonic_voltage.jl")
include("analysis/modal.jl")

include("postprocessing/reconstruction.jl")
include("postprocessing/observables.jl")
include("postprocessing/vtk_output.jl")

export AbstractFormulation, AxisymmetricRZ
export Piezo6mmConstants, Piezo6mmAxi, PZT5A, reduce_material, effective_material
export Lossless, PhysicalLoss, RealMaterial, MaterialAsGiven, PiezoComplexMaterialLoss
export NoSystemDamping
export strain, electric_field, integration_weight
export stress, electric_displacement
export PiezoElementMatrices, KFormSystem, piezo_element_matrices
export Bu_matrix, Bϕ_matrix, Nu_matrix
export HarmonicVoltageDofPartition, ShortCircuitDofPartition, OpenCircuitDofPartition
export HFormSystem, h_form_dense, short_circuit_h_form_dense
export DirichletDofs, DirichletReduction, apply_dirichlet, reconstruct_solution
export ElectrodeReducedKForm, electrode_reduced_k_form
export DirectVoltageSolution, DirectVoltageSolverInfo, solve_direct_voltage, reconstruct_potential
export PiezoFieldDofMap, compact_displacement_dof, compact_potential_node_dof
export KFormAssembly, piezo_dofhandler, assemble_k_form_dense, assemble_k_form_sparse
export potential_dofs_on_facets, potential_partition, short_circuit_partition, mechanical_dirichlet
export displacement_component_dofs_on_facets, axis_radial_displacement_constraint
export FacetElectrode, TwoTerminalElectrodes
export FacetBoundary, AxisBoundary, AxisymmetricBoundaryConditions
export PiezoProblem, AssembledPiezoProblem, HarmonicVoltageAnalysis, ShortCircuitModalAnalysis
export solve, assemble, reduce
export HarmonicVoltageResult, ShortCircuitModalReduction, ShortCircuitModalResult
export AbstractObservable, AdmittanceObservable, ChargeObservable, CurrentObservable
export ModalFrequenciesObservable, evaluate
export reconstruct_fields, write_vtk

end

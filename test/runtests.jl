using Test
using LinearAlgebra
using SparseArrays
using Ferrite
using ForwardDiff
using PiezoAcousticFEM


include("test_materials.jl")
include("test_material_contract.jl")
include("test_validation.jl")
include("test_supported_interpolations.jl")
include("test_element_assembly.jl")
include("test_boundary_conditions.jl")
include("test_reductions_direct_voltage.jl")
include("test_workflows_modal_postprocessing.jl")
include("test_observables.jl")

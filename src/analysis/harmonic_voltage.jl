"""
    PiezoProblem(grid, material_source, formulation, interpolation, quadrature; electrodes, boundary_conditions, loss)

Problem definition for a piezoelectric structure. `material_source` is stored
as supplied; formulation- and loss-specific material is created during assembly
via [`effective_material`](@ref).
"""
struct PiezoProblem{G,M,F,I,Q,E,B,L}
    grid::G
    material_source::M
    formulation::F
    interpolation::I
    quadrature::Q
    electrodes::E
    boundary_conditions::B
    loss::L
end


function PiezoProblem(
    grid,
    material_source::AbstractPiezoMaterial,
    formulation::AbstractFormulation,
    interpolation,
    quadrature;
    electrodes,
    boundary_conditions,
    loss,
)
    return PiezoProblem(
        grid,
        material_source,
        formulation,
        interpolation,
        quadrature,
        electrodes,
        boundary_conditions,
        loss,
    )
end


"""
    AssembledPiezoProblem

Problem plus effective material and assembled K-form blocks.
"""
struct AssembledPiezoProblem{P,M,A}
    problem::P
    material::M
    assembly::A
end


"""
    HarmonicVoltageAnalysis(ω, voltage, convention)

Direct harmonic voltage analysis. `convention` describes the time convention
used for admittance and phase postprocessing.
"""
struct HarmonicVoltageAnalysis{W,V,C}
    ω::W
    voltage::V
    convention::C
end


"""
    HarmonicVoltageReduction

Analysis-specific reduction from assembled K-form to direct voltage solve data.
"""
struct HarmonicVoltageReduction{A,P,R,C}
    assembled::A
    partition::P
    reduced::R
    mechanical_dirichlet::C
end


"""
    HarmonicVoltageResult

Result container for `solve(problem, ::HarmonicVoltageAnalysis)`.
"""
struct HarmonicVoltageResult{P,A,AP,R,S}
    problem::P
    analysis::A
    assembled::AP
    reduction::R
    solution::S
end


"""
    assemble(problem)

Build the formulation/loss-specific material and assemble the sparse K-form.
"""
function assemble(problem::PiezoProblem)
    material = effective_material(problem.material_source, problem.formulation, problem.loss)
    assembly = assemble_k_form_sparse(
        problem.grid,
        material,
        problem.formulation,
        problem.interpolation,
        problem.quadrature,
    )

    return AssembledPiezoProblem(problem, material, assembly)
end


"""
    reduce(assembled, analysis)

Build the analysis-specific electrode and mechanical reductions.
"""
function reduce(
    assembled::AssembledPiezoProblem,
    analysis::HarmonicVoltageAnalysis,
)
    problem = assembled.problem
    partition = potential_partition(
        assembled.assembly;
        driven=problem.electrodes.signal,
        grounded=problem.electrodes.reference,
    )
    reduced = electrode_reduced_k_form(assembled.assembly.system, partition)
    constraints = mechanical_dirichlet(assembled.assembly, problem.boundary_conditions)

    return HarmonicVoltageReduction(assembled, partition, reduced, constraints)
end


"""
    solve(problem, analysis)

Assemble, reduce, and solve a `PiezoProblem` for the given analysis.
"""
function solve(problem::PiezoProblem, analysis::HarmonicVoltageAnalysis)
    assembled = assemble(problem)
    reduction = reduce(assembled, analysis)
    solution = solve_direct_voltage(
        reduction.reduced,
        analysis.ω,
        analysis.voltage;
        mechanical_dirichlet=reduction.mechanical_dirichlet,
        analysis,
        convention=analysis.convention,
    )

    return HarmonicVoltageResult(
        problem,
        analysis,
        assembled,
        reduction,
        solution,
    )
end

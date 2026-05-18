# PiezoAcousticFEM.jl

[![CI](https://github.com/danielflataker/PiezoAcousticFEM.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/danielflataker/PiezoAcousticFEM.jl/actions/workflows/CI.yml)
[![Periodic CI](https://github.com/danielflataker/PiezoAcousticFEM.jl/actions/workflows/PeriodicCI.yml/badge.svg)](https://github.com/danielflataker/PiezoAcousticFEM.jl/actions/workflows/PeriodicCI.yml)

PiezoAcousticFEM.jl is an early-stage Julia/Ferrite project for finite-element
models of piezoelectric transducers. The current implementation focuses on a
small, explicit core for axisymmetric piezoelectric solids in vacuum: assembly
of K-form matrices, voltage-driven harmonic response, and short-circuit modal
reference calculations.

The mathematical formulation follows standard piezoelectric finite-element
theory, with the main references being Jan Kocbach's 2000 dissertation and the
1999 FEMP technical report by Kocbach, Lunde, and Vestrheim. The implementation
in this repository is original Julia code.

This package aims to provide a clear, testable, package-oriented workflow:
explicit problem objects, assembly, analysis preparation, result objects,
observables, and reproducible validation tests.

All public numeric inputs and outputs use SI units unless a docstring explicitly
states otherwise.

## Current Scope

Implemented or in progress:

- axisymmetric `r-z` piezoelectric solid formulation in vacuum;
- displacement field `u = (u_r, u_z)` and scalar electric potential `phi`;
- sparse K-form assembly with dense reference paths for small systems;
- voltage-driven harmonic response;
- short-circuit modal reference analysis;
- nodal field reconstruction and VTK output for smoke testing.

Not yet implemented:

- acoustic/fluid domains, radiation, or fluid loading;
- passive elastic domains and full transducer stacks;
- broad mesh/import workflows comparable to mature FEM tools;
- generated documentation, benchmarks, or independent benchmark comparisons.

## Validation Status

The test suite currently checks material constants, element assembly, boundary
and electrode validation, direct voltage solves, modal reference behavior,
observables, VTK output, and the public example script. These tests are useful
regressions for the current axisymmetric piezoelectric-vacuum core. Independent
benchmark comparisons are still future work.

## References

- Jan Kocbach (2000). *Finite Element Modeling of Ultrasonic Piezoelectric
  Transducers. Influence of geometry and material parameters on vibration,
  response functions and radiated field*. Dr. scient. dissertation, University
  of Bergen. [PDF](https://kocbach.net/thesis/diss.pdf)
- Jan Kocbach, Per Lunde, and Magne Vestrheim (1999). *FEMP - Finite Element
  Modeling of Piezoelectric Structures. Theory and Verification for
  Piezoceramic Disks*. Scientific/Technical Report 1999-07, University of
  Bergen, Department of Physics. [ResearchGate](https://www.researchgate.net/publication/315788533_FEMP_-_Finite_element_modeling_of_piezoelectric_structures_Theory_and_verification_for_piezoceramic_disks)

## Status

This repository is not yet a stable public API. Expect names, exports, and
result types to change while the material contract, validation cases, and public
documentation are being tightened.

## Example

Run the current public harmonic-voltage workflow with:

```sh
julia --project=. examples/disk_harmonic_voltage.jl
```

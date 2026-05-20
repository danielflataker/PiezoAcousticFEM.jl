# piezo pipeline Benchmark

- Commit: `22f8466b92b1bff6ebf537a3e81094dea1063ded`
- Julia: `1.12.6`
- Runner: `Linux`

| Benchmark | Time | Memory | Allocations | Samples |
|---|---:|---:|---:|---:|
| `assemble_disk_96x144` | 539 ms | 453.96 MiB | 1963867 | 5 |
| `prepare_harmonic_reduction_disk_96x144` | 9 ms | 16.21 MiB | 3002 | 5 |
| `solve_prepared_harmonic_disk_96x144` | 373 ms | 261.17 MiB | 297 | 5 |
| `reconstruct_fields_disk_96x144` | 1 ms | 0.43 MiB | 9 | 5 |
| `full_harmonic_pipeline_disk_96x144` | 889 ms | 731.34 MiB | 1967167 | 5 |

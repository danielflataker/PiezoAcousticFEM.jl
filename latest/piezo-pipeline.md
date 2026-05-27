# piezo pipeline Benchmark

- Commit: `b8aa226e8a117852d78e6eeee156607f37d4d62e`
- Julia: `1.12.6`
- Runner: `Linux`

| Benchmark | Time | Memory | Allocations | Samples |
|---|---:|---:|---:|---:|
| `assemble_disk_96x144` | 567 ms | 453.96 MiB | 1963867 | 5 |
| `prepare_harmonic_reduction_disk_96x144` | 5 ms | 16.21 MiB | 3002 | 5 |
| `solve_prepared_harmonic_disk_96x144` | 383 ms | 261.17 MiB | 297 | 5 |
| `reconstruct_fields_disk_96x144` | 2 ms | 0.43 MiB | 9 | 5 |
| `full_harmonic_pipeline_disk_96x144` | 944 ms | 731.34 MiB | 1967167 | 5 |

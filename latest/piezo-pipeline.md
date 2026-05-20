# piezo pipeline Benchmark

- Commit: `c0a6bc11bb9d4cd0f1be47cf9bb7f3ea799cfa86`
- Julia: `1.12.6`
- Runner: `Linux`

| Benchmark | Time | Memory | Allocations | Samples |
|---|---:|---:|---:|---:|
| `assemble_disk_96x144` | 487 ms | 453.96 MiB | 1963867 | 5 |
| `prepare_harmonic_reduction_disk_96x144` | 5 ms | 16.21 MiB | 3002 | 5 |
| `solve_prepared_harmonic_disk_96x144` | 566 ms | 261.17 MiB | 298 | 5 |
| `reconstruct_fields_disk_96x144` | 1 ms | 0.43 MiB | 9 | 5 |
| `full_harmonic_pipeline_disk_96x144` | 823 ms | 731.34 MiB | 1967167 | 5 |

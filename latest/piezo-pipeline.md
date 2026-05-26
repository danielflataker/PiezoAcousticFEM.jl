# piezo pipeline Benchmark

- Commit: `75d12854105a49fb61d74a5b9b8eb333825edce3`
- Julia: `1.12.6`
- Runner: `Linux`

| Benchmark | Time | Memory | Allocations | Samples |
|---|---:|---:|---:|---:|
| `assemble_disk_96x144` | 502 ms | 453.96 MiB | 1963867 | 5 |
| `prepare_harmonic_reduction_disk_96x144` | 6 ms | 16.21 MiB | 3002 | 5 |
| `solve_prepared_harmonic_disk_96x144` | 645 ms | 261.17 MiB | 298 | 5 |
| `reconstruct_fields_disk_96x144` | 3 ms | 0.43 MiB | 9 | 5 |
| `full_harmonic_pipeline_disk_96x144` | 868 ms | 731.34 MiB | 1967167 | 5 |

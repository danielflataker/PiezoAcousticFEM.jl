# Benchmarks

Benchmarks are grouped into named suites under `benchmark/suites/`. The common
runner, `benchmark/runbenchmarks.jl`, loads every suite file, runs each
registered `BenchmarkSpec`, and writes three output surfaces:

- `latest/<suite>.json`: latest full result for humans and tools.
- `history/<suite>/<run>.jsonl`: one JSON record per benchmark result, suitable
  for later plotting.
- `badge/<suite>.json`: a Shields.io endpoint badge for the final benchmark in
  the suite.

The `benchmarks` Git branch is used only as a data branch for these generated
outputs. It intentionally stays separate from `main` and `gh-pages`, so source
history and future documentation pages are not churned by benchmark updates.
Each benchmark run writes a new history file; the workflow keeps old files when
publishing the branch.

## Current Suite

`piezo-pipeline` benchmarks a moderately sized axisymmetric PZT-5A disk with a
96x144 quadrilateral mesh. It records the main pieces of the current public
workflow:

- sparse K-form assembly;
- harmonic voltage reduction preparation;
- prepared harmonic voltage solve;
- full assemble-prepare-solve pipeline;
- nodal field reconstruction.

The badge reports the full harmonic pipeline time. GitHub-hosted runner timing
is noisy, so treat the badge as a trend indicator rather than an absolute
machine-independent performance claim.

## Extending

Add new benchmark families as separate files in `benchmark/suites/`. Keep suite
ids stable, because they become filenames in the `benchmarks` data branch.

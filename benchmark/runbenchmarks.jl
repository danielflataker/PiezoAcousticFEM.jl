using Pkg

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIR = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(REPO_ROOT, "benchmark-output")

Pkg.activate(; temp=true)
Pkg.develop(PackageSpec(path=REPO_ROOT))
Pkg.add([
    PackageSpec(name="BenchmarkTools"),
    PackageSpec(name="Ferrite"),
    PackageSpec(name="JSON"),
])

using BenchmarkTools
using Dates
using JSON
using Printf

using PiezoAcousticFEM
using Ferrite

struct BenchmarkSpec
    id::String
    label::String
    benchmarkable::Any
end

struct BenchmarkSuiteSpec
    id::String
    label::String
    benchmarks::Vector{BenchmarkSpec}
end

const REGISTERED_SUITES = BenchmarkSuiteSpec[]

function register_benchmark_suite!(suite::BenchmarkSuiteSpec)
    push!(REGISTERED_SUITES, suite)
    return suite
end

function write_json(path, value)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, value)
        println(io)
    end
    return path
end

function write_jsonl(path, values)
    mkpath(dirname(path))
    open(path, "w") do io
        for value in values
            JSON.print(io, value)
            println(io)
        end
    end
    return path
end

function git_output(args...)
    try
        return strip(read(`git -C $REPO_ROOT $args`, String))
    catch
        return ""
    end
end

function run_metadata()
    now = Dates.now(Dates.UTC)
    timestamp = Dates.format(now, dateformat"yyyy-mm-ddTHH:MM:SSZ")
    run_slug = get(ENV, "GITHUB_RUN_ID", "")
    isempty(run_slug) && (run_slug = Dates.format(now, dateformat"yyyymmddTHHMMSSZ"))

    return (;
        timestamp,
        run_slug,
        commit=get(ENV, "GITHUB_SHA", git_output("rev-parse", "HEAD")),
        branch=get(ENV, "GITHUB_REF_NAME", git_output("branch", "--show-current")),
        workflow=get(ENV, "GITHUB_WORKFLOW", ""),
        workflow_run_id=get(ENV, "GITHUB_RUN_ID", ""),
        workflow_run_attempt=get(ENV, "GITHUB_RUN_ATTEMPT", ""),
        runner=get(ENV, "RUNNER_OS", Sys.KERNEL),
        julia_version=string(VERSION),
    )
end

function color_for_time_ms(value_ms)
    value_ms < 1_000 && return "brightgreen"
    value_ms < 5_000 && return "green"
    value_ms < 15_000 && return "yellow"
    value_ms < 30_000 && return "orange"
    return "red"
end

function format_time_ms(value_ms)
    if value_ms < 1_000
        return @sprintf("%.0f ms", value_ms)
    elseif value_ms < 10_000
        return @sprintf("%.2f s", value_ms / 1_000)
    else
        return @sprintf("%.1f s", value_ms / 1_000)
    end
end

function run_benchmark(spec::BenchmarkSpec; samples::Int)
    trial = run(spec.benchmarkable; samples, evals=1)
    estimate = median(trial)
    return (;
        id=spec.id,
        label=spec.label,
        unit="ms",
        value=estimate.time / 1.0e6,
        time_ns=estimate.time,
        memory_bytes=estimate.memory,
        allocs=estimate.allocs,
        samples=length(trial.times),
    )
end

function write_summary(path, suite::BenchmarkSuiteSpec, metadata, results)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", suite.label, " Benchmark")
        println(io)
        println(io, "- Commit: `", metadata.commit, "`")
        println(io, "- Julia: `", metadata.julia_version, "`")
        println(io, "- Runner: `", metadata.runner, "`")
        println(io)
        println(io, "| Benchmark | Time | Memory | Allocations | Samples |")
        println(io, "|---|---:|---:|---:|---:|")
        for result in results
            memory_mib = result.memory_bytes / 1024^2
            println(
                io,
                "| `", result.id, "` | ",
                format_time_ms(result.value), " | ",
                @sprintf("%.2f MiB", memory_mib), " | ",
                result.allocs, " | ",
                result.samples, " |",
            )
        end
    end
    return path
end

function load_suites()
    suites_dir = joinpath(@__DIR__, "suites")
    for file in sort(readdir(suites_dir; join=true))
        endswith(file, ".jl") && include(file)
    end
    isempty(REGISTERED_SUITES) && error("no benchmark suites were registered from $suites_dir")
    return REGISTERED_SUITES
end

function main()
    samples = parse(Int, get(ENV, "PIEZO_BENCHMARK_SAMPLES", "5"))
    metadata = run_metadata()
    suites = load_suites()

    mkpath(OUTPUT_DIR)
    for suite in suites
        results = [run_benchmark(spec; samples) for spec in suite.benchmarks]
        badge_result = last(results)

        latest = (;
            schema_version=1,
            suite=suite.id,
            label=suite.label,
            metadata,
            benchmarks=results,
        )
        history_entries = [
            (;
                schema_version=1,
                suite=suite.id,
                suite_label=suite.label,
                benchmark=result.id,
                benchmark_label=result.label,
                unit=result.unit,
                value=result.value,
                time_ns=result.time_ns,
                memory_bytes=result.memory_bytes,
                allocs=result.allocs,
                samples=result.samples,
                metadata...,
            )
            for result in results
        ]
        badge = (;
            schemaVersion=1,
            label=suite.label,
            message=format_time_ms(badge_result.value),
            color=color_for_time_ms(badge_result.value),
        )

        write_json(joinpath(OUTPUT_DIR, "latest", suite.id * ".json"), latest)
        history_file = string(metadata.run_slug, "-", first(metadata.commit, 12), ".jsonl")
        write_jsonl(joinpath(OUTPUT_DIR, "history", suite.id, history_file), history_entries)
        write_json(joinpath(OUTPUT_DIR, "badge", suite.id * ".json"), badge)
        write_summary(joinpath(OUTPUT_DIR, "latest", suite.id * ".md"), suite, metadata, results)
    end
end

main()

#!/usr/bin/env julia
# Sequential H4 Nk sweep:
#   1. build `.dat-c` with NCTSSoS (Julia only),
#   2. choose CPU/GPU backend (default: auto probe, but --backend=cpu bypasses GPU),
#   3. solve Nk values sequentially with a mediocre time/iteration budget.

using Dates
using Printf

function parse_args(argv)
    opts = Dict{String,String}()
    flags = Set{String}()
    for arg in argv
        if startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "=", limit = 2)
            opts[String(key)] = String(value)
        elseif startswith(arg, "--")
            push!(flags, String(arg[3:end]))
        else
            error("unexpected positional argument: $arg")
        end
    end
    nks_raw = get(opts, "nks", "2,3,4,5,6,7,8")
    nks = parse.(Int, split(nks_raw, ","))
    backend_raw = Symbol(get(opts, "backend", "auto"))
    backend = backend_raw === :gpu ? :cuda : backend_raw
    backend in (:auto, :cpu, :cuda) || error("--backend must be auto, cpu, or cuda")
    return (; nks,
        backend,
        norb = parse(Int, get(opts, "norb", "4")),
        nelec_per_cell = parse(Int, get(opts, "nelec-per-cell", "4")),
        outroot = get(opts, "outroot", "output/h4_nk_sweep"),
        bpsdp_project = get(opts, "bpsdp-project", "/tmp/h4_bpsdp_jl_env"),
        bpsdp_path = get(opts, "bpsdp-path", "/home/ubuntu/BPSDP.jl"),
        threads = parse(Int, get(opts, "threads", "8")),
        export_timeout = parse(Int, get(opts, "export-timeout", "1800")),
        solve_seconds = parse(Int, get(opts, "solve-seconds", "900")),
        solve_max_iter = parse(Int, get(opts, "solve-max-iter", "2000")),
        probe_seconds = parse(Int, get(opts, "probe-seconds", "180")),
        probe_iter = parse(Int, get(opts, "probe-iter", "50")),
        speedup_threshold = parse(Float64, get(opts, "speedup-threshold", "1.25")),
        force = "force" in flags)
end

function shell_quote(s::AbstractString)
    return "'" * replace(s, "'" => "'\\''") * "'"
end

function run_logged(cmd::Cmd, log_path::AbstractString)
    mkpath(dirname(log_path))
    open(log_path, "w") do io
        println(io, "# started: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
        println(io, "# command: ", cmd)
        flush(io)
        try
            run(pipeline(cmd; stdout = io, stderr = io))
            println(io, "# finished: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
            println(io, "# exit_code: 0")
            return 0
        catch err
            code = err isa ProcessFailedException ? err.procs[1].exitcode : -1
            println(io, "# error: ", sprint(showerror, err))
            println(io, "# finished: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
            println(io, "# exit_code: ", code)
            return code
        end
    end
end

function ensure_bpsdp_env(args)
    project_file = joinpath(args.bpsdp_project, "Project.toml")
    if isfile(project_file)
        return 0
    end
    cmd = `julia --startup-file=no --project=$(args.bpsdp_project) -e $("using Pkg; Pkg.develop(path=\"$(args.bpsdp_path)\"); Pkg.add(name=\"CUDA\"); Pkg.instantiate()")`
    return run_logged(cmd, joinpath(args.outroot, "setup_bpsdp_env.log"))
end

function export_cmd(nk::Int, args)
    dat_dir = joinpath(args.outroot, "nk$(nk)")
    basename = "h4_nk$(nk)_mock_pqg"
    script = joinpath("probes", "h4_nk_periodic_moment_sos_export_libsdp.jl")
    julia_cmd = `julia --startup-file=no --project=. $script --nk=$(nk) --norb=$(args.norb) --nelec-per-cell=$(args.nelec_per_cell) --outdir=$dat_dir --basename=$basename --paper-spin --no-1d`
    return `timeout $(args.export_timeout) $julia_cmd`
end

function solve_cmd(dat_path::AbstractString, backend::Symbol, summary_path::AbstractString, args; probe::Bool = false)
    script = joinpath("probes", "solve_datc_bpsdp_primitive.jl")
    max_seconds = probe ? args.probe_seconds : args.solve_seconds
    max_iter = probe ? args.probe_iter : args.solve_max_iter
    monitor_period = probe ? 5 : 25
    return `env JULIA_NUM_THREADS=$(args.threads) julia --startup-file=no --project=$(args.bpsdp_project) $script --problem=$dat_path --backend=$(String(backend)) --max-iter=$max_iter --max-seconds=$max_seconds --monitor-period=$monitor_period --summary=$summary_path`
end

function parse_summary(path::AbstractString)
    out = Dict{String,String}()
    isfile(path) || return out
    for line in eachline(path)
        occursin(" = ", line) || continue
        key, value = split(line, " = ", limit = 2)
        out[String(strip(key))] = String(strip(value))
    end
    return out
end

function useful_probe(summary::Dict{String,String})
    get(summary, "status", "error") != "error" || return false
    parse(Int, get(summary, "outer_iterations", "0")) > 0 || return false
    parse(Float64, get(summary, "solve_seconds", "Inf")) > 0 || return false
    return true
end

function seconds_per_outer(summary::Dict{String,String})
    outer = parse(Int, get(summary, "outer_iterations", "0"))
    seconds = parse(Float64, get(summary, "solve_seconds", "Inf"))
    outer > 0 || return Inf
    return seconds / outer
end

function choose_backend(first_dat::AbstractString, args)
    probe_dir = joinpath(args.outroot, "backend_probe")
    mkpath(probe_dir)

    gpu_summary_path = joinpath(probe_dir, "gpu.summary")
    gpu_log = joinpath(probe_dir, "gpu.log")
    println("[", now(), "] GPU probe first")
    gpu_code = run_logged(solve_cmd(first_dat, :cuda, gpu_summary_path, args; probe = true), gpu_log)
    gpu_summary = parse_summary(gpu_summary_path)
    if gpu_code != 0 || !useful_probe(gpu_summary)
        println("[", now(), "] GPU probe failed/unsupported; using CPU")
        return :cpu
    end

    cpu_summary_path = joinpath(probe_dir, "cpu.summary")
    cpu_log = joinpath(probe_dir, "cpu.log")
    println("[", now(), "] CPU comparison probe")
    cpu_code = run_logged(solve_cmd(first_dat, :cpu, cpu_summary_path, args; probe = true), cpu_log)
    cpu_summary = parse_summary(cpu_summary_path)
    if cpu_code != 0 || !useful_probe(cpu_summary)
        println("[", now(), "] CPU probe failed but GPU ran; using GPU")
        return :cuda
    end

    gpu_rate = seconds_per_outer(gpu_summary)
    cpu_rate = seconds_per_outer(cpu_summary)
    speedup = cpu_rate / gpu_rate
    @printf("[%s] probe rates: CPU %.4f s/outer, GPU %.4f s/outer, speedup %.2fx\n",
        string(now()), cpu_rate, gpu_rate, speedup)
    return speedup >= args.speedup_threshold ? :cuda : :cpu
end

function append_summary_md(io, nk, export_code, solve_code, backend, dat_path, solve_summary_path)
    summary = parse_summary(solve_summary_path)
    println(io, "| ", nk,
        " | ", export_code,
        " | ", backend,
        " | ", solve_code,
        " | ", get(summary, "status", "missing"),
        " | ", get(summary, "stop_reason", "missing"),
        " | ", get(summary, "outer_iterations", ""),
        " | ", get(summary, "primal_objective", ""),
        " | ", get(summary, "primal_error", ""),
        " | ", get(summary, "dual_error", ""),
        " | ", get(summary, "solve_seconds", ""),
        " | ", dat_path,
        " |")
    flush(io)
end

function main(argv = ARGS)
    args = parse_args(argv)
    mkpath(args.outroot)
    println("== H4 Nk sequential NCTSSoS/BPSDP sweep ==")
    println("nks=", args.nks)
    println("norb=", args.norb, " nelec_per_cell=", args.nelec_per_cell)
    println("outroot=", args.outroot)
    println("backend=", args.backend)
    println("No Python is used. Integrals are deterministic Julia mock data unless --integrals is wired into the exporter.")

    setup_code = ensure_bpsdp_env(args)
    setup_code == 0 || error("failed to set up BPSDP env; see setup log")

    summary_path = joinpath(args.outroot, "summary.md")
    open(summary_path, "w") do io
        println(io, "# H4 Nk sequential NCTSSoS/BPSDP sweep")
        println(io)
        println(io, "- started: ", now())
        println(io, "- nks: ", join(args.nks, ", "))
        println(io, "- active spatial orbitals per k: ", args.norb)
        println(io, "- electrons per cell: ", args.nelec_per_cell)
        println(io, "- solve max seconds per Nk: ", args.solve_seconds)
        println(io, "- solve max iterations per Nk: ", args.solve_max_iter)
        println(io, "- backend policy: ", args.backend)
        println(io, "- GPU speedup threshold: ", args.speedup_threshold, "x")
        println(io)
        println(io, "| Nk | export exit | backend | solve exit | status | stop | outer | primal objective | primal error | dual error | solve seconds | dat-c |")
        println(io, "|---:|---:|:---:|---:|:---|:---|---:|---:|---:|---:|---:|:---|")
        flush(io)

        chosen_backend = args.backend === :auto ? nothing : args.backend
        first_dat = ""
        for nk in args.nks
            nk_dir = joinpath(args.outroot, "nk$(nk)")
            basename = "h4_nk$(nk)_mock_pqg"
            dat_path = joinpath(nk_dir, basename * ".dat-c")
            export_log = joinpath(nk_dir, "export.log")
            solve_log = joinpath(nk_dir, "solve.log")
            solve_summary = joinpath(nk_dir, "solve.summary")

            export_code = 0
            if args.force || !isfile(dat_path)
                println("[", now(), "] exporting Nk=", nk)
                export_code = run_logged(export_cmd(nk, args), export_log)
            else
                println("[", now(), "] export exists for Nk=", nk, "; skipping")
            end

            if export_code != 0 || !isfile(dat_path)
                append_summary_md(io, nk, export_code, -1, :none, dat_path, solve_summary)
                continue
            end

            if chosen_backend === nothing
                first_dat = dat_path
                chosen_backend = choose_backend(first_dat, args)
                println(io)
                println(io, "Chosen backend after probe: `", chosen_backend, "`.")
                println(io)
                flush(io)
            elseif args.backend !== :auto && nk == first(args.nks)
                println(io)
                println(io, "Forced backend: `", chosen_backend, "`; GPU probe skipped.")
                println(io)
                flush(io)
            end

            println("[", now(), "] solving Nk=", nk, " with ", chosen_backend)
            solve_code = run_logged(solve_cmd(dat_path, chosen_backend, solve_summary, args), solve_log)
            append_summary_md(io, nk, export_code, solve_code, chosen_backend, dat_path, solve_summary)
        end
        println(io)
        println(io, "- finished: ", now())
    end

    println("summary: ", summary_path)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

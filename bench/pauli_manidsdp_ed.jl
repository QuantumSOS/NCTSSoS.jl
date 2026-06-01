#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using NCTSSoS
using ManiDSDP

const CUDA_IMPORT_ERROR = Ref{Any}(nothing)
try
    import CUDA
catch err
    CUDA_IMPORT_ERROR[] = err
end

const MOI = NCTSSoS.MOI
const σX_ED = ComplexF64[0 1; 1 0]
const σY_ED = ComplexF64[0 -im; im 0]
const σZ_ED = ComplexF64[1 0; 0 -1]
const I2_ED = Matrix{ComplexF64}(I, 2, 2)

function parse_options(args)
    opts = Dict{Symbol,Any}(
        :sizes => collect(3:6),
        :order => 1,
        :rank => 4,
        :backend => :cpu,
        :gpu => 0,
        :max_iter => 200,
        :tol => 1e-7,
        :inner_max_iter => 1_000,
        :inner_tol => 1e-8,
        :seed => 1234,
        :ts => :none,
        :ed => :auto,
        :ed_max_n => 12,
        :fail_fast => false,
    )
    for arg in args
        startswith(arg, "--sizes=") && (opts[:sizes] = parse_sizes(split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--order=") && (opts[:order] = parse(Int, split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--rank=") && (opts[:rank] = parse(Int, split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--backend=") && (opts[:backend] = Symbol(lowercase(split(arg, "=", limit=2)[2])); continue)
        startswith(arg, "--gpu=") && (opts[:gpu] = parse(Int, split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--max_iter=") && (opts[:max_iter] = parse(Int, split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--tol=") && (opts[:tol] = parse(Float64, split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--inner_max_iter=") && (opts[:inner_max_iter] = parse(Int, split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--inner_tol=") && (opts[:inner_tol] = parse(Float64, split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--seed=") && (opts[:seed] = parse(Int, split(arg, "=", limit=2)[2]); continue)
        startswith(arg, "--ts=") && (opts[:ts] = Symbol(lowercase(split(arg, "=", limit=2)[2])); continue)
        startswith(arg, "--ed=") && (opts[:ed] = Symbol(lowercase(split(arg, "=", limit=2)[2])); continue)
        startswith(arg, "--ed_max_n=") && (opts[:ed_max_n] = parse(Int, split(arg, "=", limit=2)[2]); continue)
        arg == "--fail-fast" && (opts[:fail_fast] = true; continue)
        error("unknown option $arg")
    end
    opts[:backend] in (:cpu, :cuda) || error("--backend must be cpu or cuda")
    opts[:ed] in (:auto, :yes, :no) || error("--ed must be auto, yes, or no")
    return opts
end

function parse_sizes(raw::AbstractString)
    if occursin(':', raw)
        parts = parse.(Int, split(raw, ':'))
        length(parts) in (2, 3) || error("--sizes range must be lo:hi or lo:step:hi")
        return length(parts) == 2 ? collect(parts[1]:parts[2]) : collect(parts[1]:parts[2]:parts[3])
    end
    return parse.(Int, split(raw, ','))
end

function maybe_init_cuda!(opts)
    opts[:backend] === :cuda || return nothing
    isdefined(@__MODULE__, :CUDA) || error("CUDA backend requested, but CUDA.jl could not be imported: $(CUDA_IMPORT_ERROR[])")
    CUDA.functional() || error("CUDA backend requested but CUDA.functional() is false")
    CUDA.device!(Int(opts[:gpu]))
    println("CUDA device: ", CUDA.device())
    println("ManiDSDP CUDA extension loaded: ", Base.get_extension(ManiDSDP, :ManiDSDPCUDAExt) !== nothing)
    return nothing
end

function kron_all(ops::Vector{Matrix{ComplexF64}})
    out = ops[1]
    for k in 2:length(ops)
        out = kron(out, ops[k])
    end
    return out
end

function heisenberg_periodic_ed(n::Integer)
    n >= 2 || throw(ArgumentError("need at least two qubits"))
    dim = 1 << n
    H = zeros(ComplexF64, dim, dim)
    for i in 1:n
        j = i == n ? 1 : i + 1
        for σ in (σX_ED, σY_ED, σZ_ED)
            ops = [I2_ED for _ in 1:n]
            ops[i] = σ
            ops[j] = σ
            H .+= 0.25 .* kron_all(ops)
        end
    end
    return eigmin(Hermitian(H))
end

function maybe_ed(n::Integer, opts)
    run_ed = opts[:ed] === :yes || (opts[:ed] === :auto && n <= Int(opts[:ed_max_n]))
    run_ed || return (value=NaN, time=0.0)
    ed_time = @elapsed ed = heisenberg_periodic_ed(n)
    return (value=Float64(ed), time=ed_time)
end

function heisenberg_periodic_poly(n::Integer)
    registry, (σx, σy, σz) = create_pauli_variables(1:n)
    H = sum(ComplexF64(0.25) * op[i] * op[mod1(i + 1, n)] for op in (σx, σy, σz) for i in 1:n)
    return polyopt(H, registry)
end

function mani_optimizer(n::Integer, opts)
    rank = min(Int(opts[:rank]), 1 + 3n)
    return MOI.OptimizerWithAttributes(
        ManiDSDP.Optimizer,
        MOI.RawOptimizerAttribute("rank") => rank,
        MOI.RawOptimizerAttribute("backend") => Symbol(opts[:backend]),
        MOI.RawOptimizerAttribute("max_iter") => Int(opts[:max_iter]),
        MOI.RawOptimizerAttribute("tol") => Float64(opts[:tol]),
        MOI.RawOptimizerAttribute("inner_max_iter") => Int(opts[:inner_max_iter]),
        MOI.RawOptimizerAttribute("inner_tol") => Float64(opts[:inner_tol]),
        MOI.Silent() => true,
    )
end

function ts_algorithm(name::Symbol)
    name === :none && return NoElimination()
    name === :mmd && return MMD()
    name === :mf && return MF()
    error("unsupported --ts=$name; expected none, mmd, or mf")
end

function run_case(n::Integer, opts)
    Random.seed!(Int(opts[:seed]) + n)
    ed = maybe_ed(n, opts)
    pop_time = @elapsed pop = heisenberg_periodic_poly(n)
    cfg = SolverConfig(
        optimizer = mani_optimizer(n, opts),
        order = Int(opts[:order]),
        ts_algo = ts_algorithm(Symbol(opts[:ts])),
    )
    sparsity_time = @elapsed sparsity = compute_sparsity(pop, cfg)
    moment_time = @elapsed mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    solve_time = @elapsed solved = NCTSSoS.solve_sdp(
        mp,
        cfg.optimizer;
        dualize=true,
        formulation=:moment_variables,
        representation=:real,
        orphan_policy=:error,
    )
    bound = Float64(real(solved.objective))
    ed_value = Float64(ed.value)
    return (
        ok = true,
        error = "",
        n = n,
        order = Int(opts[:order]),
        rank = min(Int(opts[:rank]), 1 + 3n),
        backend = Symbol(opts[:backend]),
        ts = Symbol(opts[:ts]),
        status = solved.status,
        ed = ed_value,
        bound = bound,
        gap = isfinite(ed_value) ? bound - ed_value : NaN,
        ed_per_site = isfinite(ed_value) ? ed_value / n : NaN,
        bound_per_site = bound / n,
        pop_time = pop_time,
        sparsity_time = sparsity_time,
        moment_time = moment_time,
        solve_time = solve_time,
        total_time = pop_time + sparsity_time + moment_time + solve_time,
        ed_time = Float64(ed.time),
        block_sizes = NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities),
        n_unique = solved.n_unique_elements,
    )
end

function print_row(stats)
    @printf("| %d | %d | %s | %s | %d | %s | %s | %d | %.12g | %.12g | %.3e | %.12g | %.12g | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |\n",
            stats.n,
            stats.order,
            String(stats.backend),
            String(stats.ts),
            stats.rank,
            sprint(show, stats.status),
            sprint(show, stats.block_sizes),
            stats.n_unique,
            stats.bound,
            stats.ed,
            stats.gap,
            stats.bound_per_site,
            stats.ed_per_site,
            stats.pop_time,
            stats.sparsity_time,
            stats.moment_time,
            stats.solve_time,
            stats.total_time,
            stats.ed_time)
end

function print_error_row(n, opts, err, elapsed)
    msg = replace(sprint(showerror, err), '\n' => ' ')
    @printf("| %d | %d | %s | %s | %d | ERROR | %s | 0 | NaN | NaN | NaN | NaN | NaN | NaN | NaN | NaN | %.3f | %.3f | 0.000 |\n",
            n,
            Int(opts[:order]),
            String(Symbol(opts[:backend])),
            String(Symbol(opts[:ts])),
            min(Int(opts[:rank]), 1 + 3n),
            repr(first(msg, min(lastindex(msg), 120))),
            elapsed,
            elapsed)
end

function main(args=ARGS)
    opts = parse_options(args)
    maybe_init_cuda!(opts)
    println("Pauli Heisenberg periodic chain: ManiDSDP SOS-dual bound vs ED")
    println("Julia threads: ", Threads.nthreads())
    println("opts: ", opts)
    println()
    println("| n | order | backend | ts | rank | status | block sizes | nuniq | ManiDSDP objective | ED objective | gap | Mani/site | ED/site | poly s | sparsity s | moment s | solve s | total s | ED s |")
    println("|---:|---:|:---|:---|---:|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for n in opts[:sizes]
        t0 = time()
        try
            print_row(run_case(n, opts))
        catch err
            print_error_row(n, opts, err, time() - t0)
            opts[:fail_fast] && rethrow()
        end
    end
end

main()

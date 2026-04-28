# demos/h4_periodic_term_sparsity_cosmo_benchmark.jl
#
# Purpose
# -------
# Try to actually solve the H4 periodic (Nk = 2) order-2 V2RDM relaxation with
# term sparsity + COSMO, and benchmark each stage of the pipeline.
#
# Context
# -------
# The pinned expectation `correlative_sparsity_only_attempt` in
# `test/data/expectations/h4_periodic_v2rdm.toml` records that a CS-only
# (no term sparsity) 32-mode order-2 solve produces ~508k constraints and a
# 2081 x 2081 PSD block and OOMs on Mosek.  The existing size-probe demo
# `demos/h4_periodic_term_sparsity_size.jl` shows that TS with `MMD()`
# shatters that single dense block into many smaller blocks.
#
# This script goes one step further: it actually asks COSMO to solve the
# TS-reduced SDP and reports whether the solver terminates, how long each
# stage takes, the resulting energy vs. the paper's V2RDM[4,8] target
# (-2.188 Ha/cell), and what the block footprint looks like.
#
# Expect non-trivial wall time (minutes) and memory (~10 GiB transiently)
# for the symbolic term-sparsity graph.  The SDP solve itself may be the
# slower stage or the faster one — that is exactly what we want to measure.
#
# Usage
# -----
#   julia --project=./docs demos/h4_periodic_term_sparsity_cosmo_benchmark.jl
#   julia --project=./docs demos/h4_periodic_term_sparsity_cosmo_benchmark.jl \
#       demos/results/h4_periodic_term_sparsity_cosmo_benchmark.md
#
# Optional CLI flags (must come before the markdown output path):
#   --dualize=true|false   SOS dual (cs_nctssos default) vs. primal moment SDP
#                          (default: false here — see note below).
#   --max-iter=N           COSMO max iterations (default 3000).
#   --eps=TOL              COSMO eps_abs / eps_rel (default 1e-4).
#
# Why default `--dualize=false`:
#   A first run with `dualize=true` (cs_nctssos's default) built a JuMP SOS
#   dual with **2,006,058** variables and 172,195 constraints, and COSMO did
#   not terminate within 45 minutes.  The primal moment route keeps the same
#   313 PSD blocks (max 168 × 168), but the assembled Hermitian JuMP model is
#   171,882 scalar vars representing 26,817 unique symbolic moments.  That
#   scale is much friendlier for first-order ADMM.

using Printf
using LinearAlgebra
using NCTSSoS
using JuMP
using COSMO

const MOI = JuMP.MOI
BLAS.set_num_threads(1)

# ---------------------------------------------------------------------------
# stdout flushing helpers
# ---------------------------------------------------------------------------
# Julia's `println` is buffered when stdout is not a tty (e.g. piped to `tee`).
# If the process is killed mid-run, the log on disk can be 0 bytes even after
# 45 minutes of work.  `flushln` forces Julia's IO buffer *and* the underlying
# C stdio buffer (which COSMO's verbose printer goes through) to disk after
# every progress line.  The cost is a handful of extra syscalls per stage —
# negligible next to minutes of symbolic work and ADMM iterations.
function flushln(args...)
    println(args...)
    flush(stdout)
    Base.Libc.flush_cstdio()
    return nothing
end
# Flush once up-front so even the banner survives an immediate abort.
flush(stdout); Base.Libc.flush_cstdio()

# ---------------------------------------------------------------------------
# CLI parsing — flags first, at most one positional arg (md output path)
# ---------------------------------------------------------------------------

function parse_flags(args)
    dualize   = false
    max_iter  = 3_000
    eps       = 1e-4
    decompose = true           # COSMO default; TS already decomposed blocks so this is often redundant
    rest      = String[]
    for arg in args
        if startswith(arg, "--dualize=")
            v = lowercase(split(arg, "="; limit = 2)[2])
            dualize = v in ("1", "true", "yes")
        elseif startswith(arg, "--max-iter=")
            max_iter = parse(Int, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--eps=")
            eps = parse(Float64, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--decompose=")
            v = lowercase(split(arg, "="; limit = 2)[2])
            decompose = v in ("1", "true", "yes")
        elseif startswith(arg, "--")
            error("Unknown flag: $(arg)")
        else
            push!(rest, arg)
        end
    end
    return (; dualize, max_iter, eps, decompose, positional = rest)
end

const CLI = parse_flags(ARGS)

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: fermionic_order2_basis_size,
                         fermionic_order2_nuniq,
                         load_nk2_asset

# ---------------------------------------------------------------------------
# Hamiltonian builder (mirrors demos/h4_periodic_term_sparsity_size.jl)
# ---------------------------------------------------------------------------

function build_h4_nk2_hamiltonian(h1e, eri; nk::Int, norb::Int)
    registry, ((c_up_k0, c_up_k0_dag),
               (c_dn_k0, c_dn_k0_dag),
               (c_up_k1, c_up_k1_dag),
               (c_dn_k1, c_dn_k1_dag)) = create_fermionic_variables([
        ("c_up_k0", 1:norb),
        ("c_dn_k0", 1:norb),
        ("c_up_k1", 1:norb),
        ("c_dn_k1", 1:norb),
    ])

    ann = Dict(
        (0, :up) => c_up_k0, (0, :dn) => c_dn_k0,
        (1, :up) => c_up_k1, (1, :dn) => c_dn_k1,
    )
    dag = Dict(
        (0, :up) => c_up_k0_dag, (0, :dn) => c_dn_k0_dag,
        (1, :up) => c_up_k1_dag, (1, :dn) => c_dn_k1_dag,
    )
    spin_channels = ((:up, :up), (:up, :dn), (:dn, :up), (:dn, :dn))

    ham = (0.0 + 0.0im) * (c_up_k0_dag[1] * c_up_k0[1])

    # One-body: diagonal in k and spin, scaled by 1/Nk.
    for k in 0:nk-1
        h_block = h1e[k] / nk
        for p in 1:norb, s in 1:norb, sigma in (:up, :dn)
            coeff = h_block[p, s]
            iszero(coeff) && continue
            ham += coeff * dag[(k, sigma)][p] * ann[(k, sigma)][s]
        end
    end

    # Two-body: only the 8 momentum-conserving ERI blocks survive.
    for (k1, k2, k3, k4) in sort!(collect(keys(eri)))
        block = eri[(k1, k2, k3, k4)] / nk^2
        for p in 1:norb, r in 1:norb, q in 1:norb, s in 1:norb
            coeff = block[p, r, q, s]
            iszero(coeff) && continue
            for (sigma, tau) in spin_channels
                ham += 0.5 * coeff *
                    dag[(k1, sigma)][p] * dag[(k2, tau)][q] *
                    ann[(k4, tau)][s] * ann[(k3, sigma)][r]
            end
        end
    end

    # Scrub roundoff asymmetry before feeding the symbolic pipeline.
    return registry,
           ((c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
            (c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag)),
           0.5 * (ham + adjoint(ham))
end

# ---------------------------------------------------------------------------
# Helpers for block-footprint reporting
# ---------------------------------------------------------------------------

format_bytes(b) = b < 1024        ? @sprintf("%d B",   b) :
                  b < 1024^2      ? @sprintf("%.2f KiB", b / 1024) :
                  b < 1024^3      ? @sprintf("%.2f MiB", b / 1024^2) :
                                    @sprintf("%.2f GiB", b / 1024^3)

summarize_blocks(blocks) = length.(blocks)

function block_stats(sizes)
    isempty(sizes) && return (n_blocks=0, max_size=0, total_scalar_entries=0, sum_sq=0)
    return (
        n_blocks = length(sizes),
        max_size = maximum(sizes),
        total_scalar_entries = sum(s * (s + 1) ÷ 2 for s in sizes),
        sum_sq = sum(s^2 for s in sizes),
    )
end

# ---------------------------------------------------------------------------
# Build Hamiltonian + polyopt
# ---------------------------------------------------------------------------

flushln("="^78)
flushln("H4 periodic (Nk=2) — order-2 V2RDM — term sparsity (MMD) + COSMO solve")
flushln("="^78)

asset = load_nk2_asset()
flushln("Asset preflight: Nk = $(asset.nk), norb_per_k = $(asset.n_active_orb), ",
        "total spin-orbital modes = $(asset.total_spin_orbital_modes), ",
        "HF constant shift = $(round(asset.preflight["hf_constant_shift"]; digits=10)) Ha")

build_time = @timed begin
    registry, vars, ham = build_h4_nk2_hamiltonian(
        asset.h1e, asset.eri; nk = asset.nk, norb = asset.n_active_orb)
    (c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
    (c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag) = vars

    n_up = 1.0 * sum(
        c_up_k0_dag[i] * c_up_k0[i] + c_up_k1_dag[i] * c_up_k1[i]
        for i in 1:asset.n_active_orb)
    n_dn = 1.0 * sum(
        c_dn_k0_dag[i] * c_dn_k0[i] + c_dn_k1_dag[i] * c_dn_k1[i]
        for i in 1:asset.n_active_orb)

    polyopt(ham, registry;
            moment_eq_constraints = [n_up - 4.0 * one(ham),
                                     n_dn - 4.0 * one(ham)])
end
pop = build_time.value
flushln(@sprintf("[stage] build polyopt        : %.2f s  |  alloc %s",
                 build_time.time, format_bytes(build_time.bytes)))

flushln(@sprintf("\nSolver config: dualize = %s, COSMO max_iter = %d, eps_abs/rel = %.1e",
                 CLI.dualize, CLI.max_iter, CLI.eps))

flushln(@sprintf("Extra COSMO flag: decompose = %s", CLI.decompose))

config = SolverConfig(
    optimizer = optimizer_with_attributes(
        COSMO.Optimizer,
        MOI.Silent()        => false,
        "max_iter"          => CLI.max_iter,
        "eps_abs"           => CLI.eps,
        "eps_rel"           => CLI.eps,
        "eps_prim_inf"      => 1e-5,
        "eps_dual_inf"      => 1e-5,
        "check_termination" => 200,
        "verbose"           => true,    # show ADMM progress so partial runs are diagnosable
        "verbose_timing"    => true,
        "decompose"         => CLI.decompose,
    ),
    order    = 2,
    ts_algo  = MMD(),
)

# ---------------------------------------------------------------------------
# Stage 1: compute term-sparsity structure
# ---------------------------------------------------------------------------

flushln("\n[stage] compute_sparsity     : running order-2 TS with MMD() ...")
sparsity_time = @timed compute_sparsity(pop, config)
sparsity = sparsity_time.value
flushln(@sprintf("[stage] compute_sparsity     : %.2f s  |  alloc %s",
                 sparsity_time.time, format_bytes(sparsity_time.bytes)))

blocks = sparsity.cliques_term_sparsities[1][1].block_bases
sizes = summarize_blocks(blocks)
stats = block_stats(sizes)

dense_block_dim = Int(fermionic_order2_basis_size(asset.total_spin_orbital_modes))
flushln("\n-- block footprint at order 2 --")
flushln(@sprintf("  dense baseline PSD block          : %d × %d", dense_block_dim, dense_block_dim))
flushln(@sprintf("  correlative cliques               : %d", length(sparsity.corr_sparsity.cliques)))
flushln(@sprintf("  TS PSD blocks                     : %d", stats.n_blocks))
flushln(@sprintf("  TS largest block                  : %d × %d", stats.max_size, stats.max_size))
flushln(@sprintf("  sum(block_size^2) proxy           : %d", stats.sum_sq))
flushln(@sprintf("  sum of upper-tri PSD scalar slots : %d", stats.total_scalar_entries))

# ---------------------------------------------------------------------------
# Stage 2: moment relaxation (symbolic)
# ---------------------------------------------------------------------------

flushln("\n[stage] moment_relax         : assembling symbolic moment problem ...")
mp_time = @timed NCTSSoS.moment_relax(pop, sparsity.corr_sparsity,
                                      sparsity.cliques_term_sparsities)
moment_problem = mp_time.value
flushln(@sprintf("[stage] moment_relax         : %.2f s  |  alloc %s",
                 mp_time.time, format_bytes(mp_time.bytes)))
flushln(@sprintf("  unique moment variables           : %d",
                 moment_problem.n_unique_moment_matrix_elements))

# ---------------------------------------------------------------------------
# Stage 3: build the JuMP SDP — either SOS dual or primal moment
# ---------------------------------------------------------------------------

# Two routes:
#   - dualize = true  : sos_dualize first (build JuMP SOS dual), then optimize!
#   - dualize = false : solve_moment_problem builds the primal moment JuMP
#                       model and optimizes it in one shot. We label the
#                       combined wall time under `optimize!` and leave the
#                       symbolic-build column at 0 since the public API does
#                       not expose them as separate steps.

local jump_model
local jump_kind
if CLI.dualize
    jump_kind = "SOS dual"
    flushln("\n[stage] sos_dualize          : building JuMP SOS dual ...")
    dualize_time = @timed NCTSSoS.sos_dualize(moment_problem)
    sos_problem = dualize_time.value
    set_optimizer(sos_problem.model, config.optimizer)
    jump_model = sos_problem.model
    flushln(@sprintf("[stage] sos_dualize          : %.2f s  |  alloc %s",
                     dualize_time.time, format_bytes(dualize_time.bytes)))
    flushln(@sprintf("  JuMP num_variables                : %d", num_variables(jump_model)))
    flushln(@sprintf("  JuMP num_constraints (all types)  : %d",
                     sum(num_constraints(jump_model, F, S)
                         for (F, S) in list_of_constraint_types(jump_model); init = 0)))

    flushln("\n[stage] optimize! (COSMO)    : running SDP solver ...")
    flushln("---------------- begin COSMO output ----------------")
    solve_wall = @timed optimize!(jump_model)
    flushln("---------------- end   COSMO output ----------------")
    flushln(@sprintf("[stage] optimize! (COSMO)    : %.2f s (wall)", solve_wall.time))
else
    jump_kind = "primal moment"
    # Split the primal path explicitly into JuMP build vs COSMO optimize so we
    # can tell whether the wall time is in model assembly or in ADMM.  This
    # inlines the logic from `NCTSSoS._solve_complex_moment_problem` (H4 is a
    # fermionic/complex algebra), minus the combined `optimize!` call.
    flushln("\n[stage] primal moment build  : assembling JuMP primal moment SDP ...")
    build_primal = @timed begin
        Cr = real(eltype(coefficients(moment_problem.objective)))
        model = GenericModel{Cr}()
        Mty = eltype(moment_problem.total_basis)
        basis = [symmetric_canon(NCTSSoS.expval(m)) for m in moment_problem.total_basis]
        NCTSSoS.sorted_unique!(basis)
        n_basis = length(basis)
        @variable(model, y_re[1:n_basis], set_string_name = false)
        @variable(model, y_im[1:n_basis], set_string_name = false)
        one_sym = symmetric_canon(NCTSSoS.expval(one(Mty)))
        idx_one = findfirst(==(one_sym), basis)
        idx_one === nothing && error("Expected identity moment to be present in basis")
        @constraint(model, y_re[idx_one] == 1)
        @constraint(model, y_im[idx_one] == 0)
        basis_to_idx = Dict(m => i for (i, m) in enumerate(basis))
        flushln(@sprintf("  basis length / unique moments     : %d", n_basis))
        n_constraint_blocks = length(moment_problem.constraints)
        flushln(@sprintf("  constraint blocks to assemble     : %d", n_constraint_blocks))
        # Walk the symbolic constraints, convert them to JuMP expressions,
        # stage-print every 25%% block so a hang is visible in the log.
        report_every = max(1, cld(n_constraint_blocks, 4))
        for (k, (cone, mat)) in enumerate(moment_problem.constraints)
            dim = size(mat, 1)
            mat_re = Matrix{Any}(undef, dim, dim)
            mat_im = Matrix{Any}(undef, dim, dim)
            for i in 1:dim, j in 1:dim
                re_expr, im_expr = NCTSSoS._substitute_complex_poly(
                    mat[i,j], basis_to_idx, y_re, y_im)
                mat_re[i,j] = re_expr
                mat_im[i,j] = im_expr
            end
            if cone == :Zero
                @constraint(model, [mat_re[i,j] for i in 1:dim, j in 1:dim] .== 0)
                @constraint(model, [mat_im[i,j] for i in 1:dim, j in 1:dim] .== 0)
            elseif cone == :HPSD
                embedded = [
                    [mat_re[i,j] for i in 1:dim, j in 1:dim]   [-mat_im[i,j] for i in 1:dim, j in 1:dim]
                    [mat_im[i,j] for i in 1:dim, j in 1:dim]   [mat_re[i,j] for i in 1:dim, j in 1:dim]
                ]
                @constraint(model, embedded in PSDCone())
            else
                error("Unexpected cone type $cone for complex problem")
            end
            if k % report_every == 0 || k == n_constraint_blocks
                flushln(@sprintf("    ... constraint block %d / %d", k, n_constraint_blocks))
            end
        end
        obj_re, _ = NCTSSoS._substitute_complex_poly(moment_problem.objective,
                                                      basis_to_idx, y_re, y_im)
        @objective(model, Min, obj_re)
        set_optimizer(model, config.optimizer)
        # silent = false: leave COSMO verbose output on
        model
    end
    jump_model   = build_primal.value
    dualize_time = (time = build_primal.time, bytes = build_primal.bytes)
    flushln(@sprintf("[stage] primal moment build  : %.2f s  |  alloc %s",
                     build_primal.time, format_bytes(build_primal.bytes)))
    flushln(@sprintf("  JuMP num_variables                : %d", num_variables(jump_model)))
    flushln(@sprintf("  JuMP num_constraints (all types)  : %d",
                     sum(num_constraints(jump_model, F, S)
                         for (F, S) in list_of_constraint_types(jump_model); init = 0)))

    flushln("\n[stage] optimize! (COSMO)    : running SDP solver ...")
    flushln("---------------- begin COSMO output ----------------")
    solve_wall = @timed optimize!(jump_model)
    flushln("---------------- end   COSMO output ----------------")
    flushln(@sprintf("[stage] optimize! (COSMO)    : %.2f s (wall)  |  alloc %s",
                     solve_wall.time, format_bytes(solve_wall.bytes)))
end

term_status  = termination_status(jump_model)
primal_stat  = primal_status(jump_model)
dual_stat    = dual_status(jump_model)
moi_solve_t  = try solve_time(jump_model) catch; NaN end
sol_summary  = solution_summary(jump_model)

obj_active = try objective_value(jump_model) catch; NaN end
hf_shift   = asset.preflight["hf_constant_shift"]
e_total    = obj_active + hf_shift

figure_v2rdm = asset.preflight["figure_v2rdm_nk2"]
figure_hf    = asset.preflight["figure_hf_nk2"]
figure_unc   = asset.preflight["figure_uncertainty_ha"]

flushln("\n-- COSMO termination --")
flushln("  termination_status : ", term_status)
flushln("  primal_status      : ", primal_stat)
flushln("  dual_status        : ", dual_stat)
flushln(@sprintf("  MOI solve_time     : %.2f s", moi_solve_t))
flushln(sol_summary)

flushln("\n-- Energy (Ha / unit cell) --")
flushln(@sprintf("  active-space electronic (SDP obj.)   : %.8f", obj_active))
flushln(@sprintf("  + HF constant shift                  : %.8f", hf_shift))
flushln(@sprintf("  => total energy (NCTSSoS + shift)    : %.8f", e_total))
flushln(@sprintf("  paper figure V2RDM[4,8]              : %.3f  (± %.3f)", figure_v2rdm, figure_unc))
flushln(@sprintf("  paper figure HF                      : %.3f", figure_hf))
flushln(@sprintf("  Δ(SDP total − paper V2RDM)           : %+.4f  Ha", e_total - figure_v2rdm))
flushln(@sprintf("  Δ(SDP total − paper HF)              : %+.4f  Ha", e_total - figure_hf))

total_wall = build_time.time + sparsity_time.time + mp_time.time + dualize_time.time + solve_wall.time
total_alloc = build_time.bytes + sparsity_time.bytes + mp_time.bytes + dualize_time.bytes + solve_wall.bytes
flushln("\n-- Aggregate --")
flushln(@sprintf("  total wall time (all 5 stages) : %.2f s", total_wall))
flushln(@sprintf("  total allocated bytes          : %s", format_bytes(total_alloc)))

# ---------------------------------------------------------------------------
# Optional markdown dump
# ---------------------------------------------------------------------------

function markdown_report(;
    asset, stats, dense_block_dim, moment_problem, jump_model, jump_kind,
    build_time, sparsity_time, mp_time, dualize_time, solve_wall,
    term_status, primal_stat, dual_stat, moi_solve_t,
    obj_active, hf_shift, e_total,
    figure_v2rdm, figure_hf, figure_unc,
    total_wall, total_alloc, cli,
)
    io = IOBuffer()
    println(io, "# H4 periodic (Nk=2) — order-2 V2RDM with term sparsity + COSMO")
    println(io)
    println(io, "Full 32-mode Bloch Hamiltonian from the vendored asset, ",
                "`ts_algo = MMD()`, `cs_algo = NoElimination()` (the default CS here ",
                "collapses to one clique of 32 on the complete spatial graph anyway, ",
                "see `spin_orbital_order2_blocker` in the expectations file).")
    println(io)
    println(io, @sprintf("Solver: COSMO with `dualize = %s`, `eps_abs = eps_rel = %.1e`, `max_iter = %d`.",
                         cli.dualize, cli.eps, cli.max_iter))
    println(io)

    println(io, "## Block footprint at order 2")
    println(io)
    println(io, "| Quantity | Value |")
    println(io, "| :--- | ---: |")
    println(io, "| dense baseline PSD block | $(dense_block_dim) × $(dense_block_dim) |")
    println(io, "| TS PSD blocks | $(stats.n_blocks) |")
    println(io, "| TS largest block | $(stats.max_size) × $(stats.max_size) |")
    println(io, "| sum(block_size²) | $(stats.sum_sq) |")
    println(io, "| sum of upper-tri PSD scalar slots | $(stats.total_scalar_entries) |")
    println(io, "| unique moment variables | $(moment_problem.n_unique_moment_matrix_elements) |")
    println(io, "| JuMP num_variables ($jump_kind) | $(num_variables(jump_model)) |")
    println(io)

    println(io, "## Stage timings")
    println(io)
    println(io, "| Stage | wall (s) | allocated |")
    println(io, "| :--- | ---: | ---: |")
    println(io, @sprintf("| build polyopt | %.2f | %s |", build_time.time, format_bytes(build_time.bytes)))
    println(io, @sprintf("| compute_sparsity (MMD) | %.2f | %s |", sparsity_time.time, format_bytes(sparsity_time.bytes)))
    println(io, @sprintf("| moment_relax | %.2f | %s |", mp_time.time, format_bytes(mp_time.bytes)))
    build_label = cli.dualize ? "sos_dualize + JuMP build" : "primal moment JuMP build"
    println(io, @sprintf("| %s | %.2f | %s |", build_label, dualize_time.time, format_bytes(dualize_time.bytes)))
    solve_label = cli.dualize ? "optimize! (COSMO)" : "optimize! (COSMO, primal moment)"
    println(io, @sprintf("| %s | %.2f | %s |", solve_label, solve_wall.time, format_bytes(solve_wall.bytes)))
    println(io, @sprintf("| **total** | **%.2f** | **%s** |", total_wall, format_bytes(total_alloc)))
    println(io)

    println(io, "## COSMO termination")
    println(io)
    println(io, "- `termination_status`: `", term_status, "`")
    println(io, "- `primal_status`: `", primal_stat, "`")
    println(io, "- `dual_status`: `", dual_stat, "`")
    println(io, @sprintf("- MOI `solve_time`: %.2f s", moi_solve_t))
    println(io)

    println(io, "## Energy (Ha / unit cell)")
    println(io)
    println(io, "| Quantity | Value |")
    println(io, "| :--- | ---: |")
    println(io, @sprintf("| active-space electronic (SDP objective) | %.8f |", obj_active))
    println(io, @sprintf("| HF constant shift | %.8f |", hf_shift))
    println(io, @sprintf("| **total energy (SDP + shift)** | **%.8f** |", e_total))
    println(io, @sprintf("| paper figure V2RDM[4,8] | %.3f ± %.3f |", figure_v2rdm, figure_unc))
    println(io, @sprintf("| paper figure HF | %.3f |", figure_hf))
    println(io, @sprintf("| Δ(SDP total − paper V2RDM) | %+.4f |", e_total - figure_v2rdm))
    println(io, @sprintf("| Δ(SDP total − paper HF) | %+.4f |", e_total - figure_hf))
    println(io)

    return String(take!(io))
end

if !isempty(CLI.positional)
    out_path = first(CLI.positional)
    report = markdown_report(;
        asset, stats, dense_block_dim, moment_problem, jump_model, jump_kind,
        build_time, sparsity_time, mp_time, dualize_time, solve_wall,
        term_status, primal_stat, dual_stat, moi_solve_t,
        obj_active, hf_shift, e_total,
        figure_v2rdm, figure_hf, figure_unc,
        total_wall, total_alloc, cli = CLI,
    )
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        write(io, report)
    end
    flushln("\nReport written to $(out_path)")
end

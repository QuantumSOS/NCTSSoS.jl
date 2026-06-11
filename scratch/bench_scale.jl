# N-site periodic Heisenberg ring, relaxation order 2, Mosek.
# Symmetry-reduced solve (sympleq_symmetry_spec) first, dense baseline second,
# so partial results survive if the dense solve becomes intractable.
#
# Warm-up at N=6 (both code paths) is discarded: JIT compilation is
# per-method-signature, not per-problem-size, and a warm-up dense solve at the
# target N would double the cost of large runs.
#
# Run on HAI: julia --project=docs scratch/bench_scale.jl <N> [--no-dense]

using NCTSSoS, MosekTools
import JuMP

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)

N = parse(Int, ARGS[1])
run_dense = !("--no-dense" in ARGS)

function build(n)
    registry, (σx, σy, σz) = create_pauli_variables(1:n)
    H = sum(
        ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, n)]
        for op in (σx, σy, σz) for i in 1:n
    )
    return polyopt(H, registry), H
end

dense_cfg() = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
)

sym_cfg(spec) = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = spec,
)

println("=== N=$N Heisenberg ring, order 2 (dense=$(run_dense)) ===")
println()
println("=== Warm-up at N=6 (discarded) ===")
flush(stdout)

pop6, H6 = build(6)
t_w1 = @elapsed cs_nctssos(pop6, dense_cfg())
spec6 = sympleq_symmetry_spec(H6)
t_w2 = @elapsed cs_nctssos(pop6, sym_cfg(spec6))
println("warm-up done (dense $(round(t_w1; digits=1))s, sym $(round(t_w2; digits=1))s)")
flush(stdout)
pop6 = H6 = spec6 = nothing
GC.gc()

pop, H = build(N)

println()
println("=== SympleQ detection + symmetry-reduced solve (timed) ===")
flush(stdout)

GC.gc()
t_detect = @elapsed spec = sympleq_symmetry_spec(H)
println("detect time           = $(round(t_detect; digits=3))s")
println("generators            = $(length(spec.clifford_generators))")
flush(stdout)

GC.gc()
t_sym = @elapsed result = cs_nctssos(pop, sym_cfg(spec))
t_sym_mosek = JuMP.solve_time(result.model)

report = result.symmetry
sym_obj    = result.objective
sym_vars   = result.n_unique_moment_matrix_elements
sym_blocks = report.psd_block_sizes
sym_status = JuMP.termination_status(result.model)

println("sym objective         = $sym_obj")
println("sym status            = $sym_status")
println("group order           = $(report.group_order)")
println("psd block sizes       = $sym_blocks")
println("n blocks / max block  = $(length(sym_blocks)) / $(maximum(sym_blocks))")
println("sym moment vars       = $sym_vars")
println("sym total time        = $(round(t_sym; digits=2))s")
println("  of which Mosek      = $(round(t_sym_mosek; digits=3))s")
println("  construct (approx)  = $(round(t_sym - t_sym_mosek; digits=2))s")
flush(stdout)

group_order = report.group_order
n_gens = length(spec.clifford_generators)
result = nothing
report = nothing
GC.gc()

if run_dense
    println()
    println("=== Dense baseline (timed) ===")
    flush(stdout)

    GC.gc()
    t_dense = @elapsed dense_result = cs_nctssos(pop, dense_cfg())
    t_dense_mosek = JuMP.solve_time(dense_result.model)

    dense_obj    = dense_result.objective
    dense_vars   = dense_result.n_unique_moment_matrix_elements
    dense_sizes  = dense_result.moment_matrix_sizes
    dense_status = JuMP.termination_status(dense_result.model)

    println("dense objective       = $dense_obj")
    println("dense status          = $dense_status")
    println("dense moment sizes    = $dense_sizes")
    println("dense moment vars     = $dense_vars")
    println("dense total time      = $(round(t_dense; digits=2))s")
    println("  of which Mosek      = $(round(t_dense_mosek; digits=2))s")
    flush(stdout)

    println()
    println("=== Comparison (N=$N, steady-state; warm-up at N=6 discarded) ===")
    println("run       total        of which Mosek   moment vars   PSD blocks")
    println("dense     $(round(t_dense; digits=1))s      $(round(t_dense_mosek; digits=1))s        $dense_vars      $dense_sizes")
    println("symmetry  $(round(t_detect + t_sym; digits=1))s       $(round(t_sym_mosek; digits=2))s        $sym_vars       $(length(sym_blocks)) blocks, max $(maximum(sym_blocks))")
    println("  symmetry time = detect $(round(t_detect; digits=2))s + solve $(round(t_sym; digits=1))s")
    println("End-to-end speedup: $(round(t_dense / (t_detect + t_sym); digits=1))x")
    println("Solver-only speedup: $(round(t_dense_mosek / t_sym_mosek; digits=0))x")
    println("Moment variable reduction: $dense_vars -> $sym_vars")
    println("Objective difference: $(round(abs(dense_obj - sym_obj); sigdigits=3))")
    println("Group order: $group_order, generators: $n_gens")
end

println()
println("BENCH-OK N=$N")
flush(stdout)

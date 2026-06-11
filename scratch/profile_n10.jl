# Stage-by-stage profiling of the N=10 order-2 Heisenberg benchmark.
# Splits each pipeline (dense / auto symmetry / auto ∪ D10 symmetry) into:
#   1. compute_sparsity          (basis generation)
#   2. relaxation construction   (moment_relax / moment_relax_symmetric —
#                                 the latter includes symmetry-adapted basis construction)
#   3. sos_dualize               (SOS dual JuMP model build)
#   4. optimize!                 (wall) + JuMP.solve_time (raw Mosek inner time)
#
# Warm-up pass per path discarded.
#
# Run on HAI: julia --project=docs scratch/profile_n10.jl

using NCTSSoS, MosekTools, JuMP

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)

N = 10
registry, (σx, σy, σz) = create_pauli_variables(1:N)

H = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
    for op in (σx, σy, σz) for i in 1:N
)

pop = polyopt(H, registry)

M = typeof(σx[1])
translation = CliffordSymmetry(Dict{M,Tuple{Int,M}}(
    op[i] => (1, op[mod1(i + 1, N)]) for op in (σx, σy, σz), i in 1:N
))
reflection = CliffordSymmetry(Dict{M,Tuple{Int,M}}(
    op[i] => (1, op[N + 1 - i]) for op in (σx, σy, σz), i in 1:N
))

base_config = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
)

function profile_path(name, pop, config, spec)
    println("=== $name ===")
    flush(stdout)

    GC.gc()
    t_sparsity = @elapsed sparsity = compute_sparsity(pop, config)

    GC.gc()
    local mp, t_relax, report
    if spec === nothing
        t_relax = @elapsed mp = NCTSSoS.moment_relax(
            pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        report = nothing
    else
        t_relax = @elapsed (mp, report) = NCTSSoS.moment_relax_symmetric(
            pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities, spec)
    end

    GC.gc()
    t_dualize = @elapsed sos = NCTSSoS.sos_dualize(mp)

    set_optimizer(sos.model, SILENT_MOSEK)
    GC.gc()
    t_opt = @elapsed optimize!(sos.model)
    t_mosek = JuMP.solve_time(sos.model)

    obj = objective_value(sos.model)
    total = t_sparsity + t_relax + t_dualize + t_opt

    println("  objective            = $obj")
    if report !== nothing
        println("  group order          = $(report.group_order)")
        println("  psd block sizes      = $(report.psd_block_sizes)")
    end
    println("  moment vars          = $(mp.n_unique_moment_matrix_elements)")
    println("  compute_sparsity     = $(round(t_sparsity; digits=2))s")
    println("  relax construction   = $(round(t_relax; digits=2))s")
    println("  sos_dualize (JuMP)   = $(round(t_dualize; digits=2))s")
    println("  optimize! (wall)     = $(round(t_opt; digits=2))s")
    println("    of which Mosek     = $(round(t_mosek; digits=2))s")
    println("  TOTAL                = $(round(total; digits=2))s")
    println()
    flush(stdout)
    return total
end

println("=== Warm-up (discarded) ===")
flush(stdout)

auto_spec = sympleq_symmetry_spec(H)
union_spec = SymmetrySpec(auto_spec.clifford_generators..., translation, reflection)

t = @elapsed profile_path("warmup dense", pop, base_config, nothing)
t = @elapsed profile_path("warmup auto", pop, base_config, auto_spec)
t = @elapsed profile_path("warmup union", pop, base_config, union_spec)

println("=== Timed runs ===")
println()
flush(stdout)

GC.gc()
t_detect = @elapsed auto_spec2 = sympleq_symmetry_spec(H)
println("sympleq detection      = $(round(t_detect; digits=2))s")
println()
flush(stdout)

union_spec2 = SymmetrySpec(auto_spec2.clifford_generators..., translation, reflection)

total_dense = profile_path("dense", pop, base_config, nothing)
total_auto  = profile_path("auto symmetry (order 16)", pop, base_config, auto_spec2)
total_union = profile_path("auto ∪ D10 (order 160)", pop, base_config, union_spec2)

println("=== Summary ===")
println("dense total            = $(round(total_dense; digits=1))s")
println("auto total + detect    = $(round(total_auto + t_detect; digits=1))s  (speedup $(round(total_dense/(total_auto+t_detect); digits=1))x)")
println("union total + detect   = $(round(total_union + t_detect; digits=1))s  (speedup $(round(total_dense/(total_union+t_detect); digits=1))x)")

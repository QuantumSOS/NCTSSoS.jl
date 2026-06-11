# N=10 Heisenberg ring, order 2, Mosek: auto SympleQ spec UNION manual D10
# (translation + reflection). Compares against the auto-only result from
# scratch/bench_n10.jl (dense baseline: 129.5s, 436x436, 20686 moment vars).
#
# Run on HAI: julia --project=docs scratch/bench_n10_union.jl

using NCTSSoS, MosekTools

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

println("=== Warm-up (discarded) ===")
flush(stdout)

warm_auto = sympleq_symmetry_spec(H)
warm_union_spec = SymmetrySpec(warm_auto.clifford_generators..., translation, reflection)
warm_config = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = warm_union_spec,
)
t_warm = @elapsed cs_nctssos(pop, warm_config)
println("warm-up union done in $(round(t_warm; digits=1))s")
flush(stdout)
GC.gc()

println()
println("=== Detection (timed) ===")
GC.gc()
t_detect = @elapsed auto_spec = sympleq_symmetry_spec(H)
println("detect time           = $(round(t_detect; digits=3))s")
flush(stdout)

union_spec = SymmetrySpec(auto_spec.clifford_generators..., translation, reflection)

config = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = union_spec,
)

println()
println("=== auto ∪ D10 symmetry-reduced solve (timed) ===")
flush(stdout)

GC.gc()
t_union = @elapsed result = cs_nctssos(pop, config)

report = result.symmetry
println("union objective       = $(result.objective)")
println("group order           = $(report.group_order)")
println("psd block sizes       = $(report.psd_block_sizes)")
println("union moment vars     = $(result.n_unique_moment_matrix_elements)")
println("union solve time      = $(round(t_union; digits=2))s")
println()
println("total (detect+solve)  = $(round(t_detect + t_union; digits=1))s")
println("vs dense 129.5s       => speedup $(round(129.5 / (t_detect + t_union); digits=1))x")
println("vs auto-only 54.6s    => $(round(54.6 / (t_detect + t_union); digits=1))x better")
flush(stdout)

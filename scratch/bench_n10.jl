# N=10 site periodic Heisenberg ring, relaxation order 2, Mosek.
# Dense baseline vs sympleq_symmetry_spec symmetry-reduced solve.
# Warm-up runs discarded; detection and solve timed separately.
#
# Run on HAI: julia --project=docs scratch/bench_n10.jl

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

dense_config = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
)

println("=== Warm-up (discarded) ===")
flush(stdout)

warm_spec = sympleq_symmetry_spec(H)
warm_config = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = warm_spec,
)

t_warm_dense = @elapsed cs_nctssos(pop, dense_config)
println("warm-up dense done in $(round(t_warm_dense; digits=1))s")
flush(stdout)
GC.gc()

t_warm_sym = @elapsed cs_nctssos(pop, warm_config)
println("warm-up symmetry done in $(round(t_warm_sym; digits=1))s")
flush(stdout)
GC.gc()

println()
println("=== Dense baseline (timed) ===")
flush(stdout)

GC.gc()
t_dense = @elapsed dense_result = cs_nctssos(pop, dense_config)

println("dense objective       = $(dense_result.objective)")
println("dense moment sizes    = $(dense_result.moment_matrix_sizes)")
println("dense moment vars     = $(dense_result.n_unique_moment_matrix_elements)")
println("dense time            = $(round(t_dense; digits=2))s")
flush(stdout)

println()
println("=== SympleQ detection + symmetry-reduced solve (timed) ===")
flush(stdout)

GC.gc()
t_detect = @elapsed auto_spec = sympleq_symmetry_spec(H)
println("detect time           = $(round(t_detect; digits=3))s")
println("generators            = $(length(auto_spec.clifford_generators))")
flush(stdout)

config = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = auto_spec,
)

GC.gc()
t_sym = @elapsed result = cs_nctssos(pop, config)

report = result.symmetry
println("sym objective         = $(result.objective)")
println("group order           = $(report.group_order)")
println("psd block sizes       = $(report.psd_block_sizes)")
println("sym moment vars       = $(result.n_unique_moment_matrix_elements)")
println("sym solve time        = $(round(t_sym; digits=2))s")
flush(stdout)

println()
println("=== Comparison (steady-state; warm-up discarded) ===")
println("run       time         moment vars   PSD blocks")
println("dense     $(round(t_dense; digits=1))s      $(dense_result.n_unique_moment_matrix_elements)        $(dense_result.moment_matrix_sizes)")
println("symmetry  $(round(t_detect + t_sym; digits=1))s      $(result.n_unique_moment_matrix_elements)         $(report.psd_block_sizes)")
println("  symmetry time = detect $(round(t_detect; digits=2))s + solve $(round(t_sym; digits=1))s")
println("Speedup: $(round(t_dense / (t_detect + t_sym); digits=1))x")
println("Moment variable reduction: $(dense_result.n_unique_moment_matrix_elements) -> $(result.n_unique_moment_matrix_elements)")
println("Objective difference: $(round(abs(dense_result.objective - result.objective); sigdigits=3))")
println("Group order: $(report.group_order), generators: $(length(auto_spec.clifford_generators))")

using NCTSSoS, COSMO, LinearAlgebra
const MOI = NCTSSoS.MOI
const SILENT_SOLVER = MOI.OptimizerWithAttributes(COSMO.Optimizer, MOI.Silent() => true)

ED_GS = -3.6510934089371814  # exact diag ground state for N=8 periodic Heisenberg

N = 8
registry, (σx, σy, σz) = create_pauli_variables(1:N)
H = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)] for op in (σx, σy, σz) for i in 1:N)
pop = polyopt(H, registry)

println("=" ^ 60)
println("8-site Heisenberg ring — orders 1, 2, 3")
println("ED ground state: $ED_GS")
println("Solver: COSMO")
println("=" ^ 60)

# JIT warm-up
begin
    cfg0 = SolverConfig(optimizer=SILENT_SOLVER, order=1, cs_algo=NoElimination(), ts_algo=NoElimination())
    cs_nctssos(pop, cfg0)
    sp0 = sympleq_symmetry_spec(H)
    cfg0s = SolverConfig(optimizer=SILENT_SOLVER, order=1, cs_algo=NoElimination(), ts_algo=NoElimination(), symmetry=sp0)
    cs_nctssos(pop, cfg0s)
    println("JIT warm-up done\n")
end
GC.gc()

for d in [1, 2, 3]
    println("--- Order $d ---")

    # Dense (skip for order >= 3 — too large)
    if d <= 2
        dense_cfg = SolverConfig(optimizer=SILENT_SOLVER, order=d, cs_algo=NoElimination(), ts_algo=NoElimination())
        t_dense = @elapsed dense_result = cs_nctssos(pop, dense_cfg)
        dense_obj = real(dense_result.objective)
        println("  Dense:  $(round(t_dense;digits=2))s | blocks=$(dense_result.moment_matrix_sizes) | obj=$(round(dense_obj;digits=8)) | gap=$(round(ED_GS - dense_obj;digits=8))")
        GC.gc()
    else
        println("  Dense:  SKIPPED (1789×1789 too large for dense solve)")
    end

    # SympleQ + symmetry
    t_detect = @elapsed auto_spec = sympleq_symmetry_spec(H)
    sym_cfg = SolverConfig(optimizer=SILENT_SOLVER, order=d, cs_algo=NoElimination(), ts_algo=NoElimination(), symmetry=auto_spec)
    t_sym = @elapsed sym_result = cs_nctssos(pop, sym_cfg)
    report = sym_result.symmetry
    sym_obj = real(sym_result.objective)
    t_total = t_detect + t_sym
    println("  Sym:    $(round(t_total;digits=2))s (detect=$(round(t_detect;digits=2))s + solve=$(round(t_sym;digits=2))s)")
    println("          blocks=$(report.psd_block_sizes)")
    println("          max_block=$(maximum(report.psd_block_sizes)) | obj=$(round(sym_obj;digits=8)) | gap=$(round(ED_GS - sym_obj;digits=8))")
    if d <= 2
        println("  Speedup: $(round(t_dense/t_total;digits=2))x")
    end
    println("  Group order: $(report.group_order)")
    println()
    GC.gc()
end

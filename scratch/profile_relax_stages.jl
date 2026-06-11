# Stage-by-stage instrumentation of moment_relax_symmetric's Pauli/Clifford
# path for the N=10, order-2 Heisenberg benchmark (auto SympleQ spec, group
# order 16). Replays the internal pipeline with @elapsed around each phase:
#
#   1. compute_sparsity + total basis assembly
#   2. CliffordSymmetryGroup enumeration
#   3. _check_symmetry_invariance (objective invariance under all 16 elements)
#   4. _check_basis_closure (436-monomial basis closed under group)
#   5. _build_orbit_reducer (orbits of ~20k total-basis monomials)
#   6. _build_constraint_matrix (436x436 polynomial moment matrix)
#   7. _sw_decompose_half_basis (SymbolicWedderburn.symmetry_adapted_basis)
#   8a. diagonal congruence transforms  U_i' M U_i  (the blocks we keep)
#   8b. off-diagonal zero-verification  U_i' M U_j, i != j  (assert-only work)
#
# Two passes; pass 1 is JIT warm-up, pass 2 is reported.
#
# Run on HAI: julia --project=docs scratch/profile_relax_stages.jl

using NCTSSoS

const NC = NCTSSoS

N = 10
registry, (σx, σy, σz) = create_pauli_variables(1:N)

H = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
    for op in (σx, σy, σz) for i in 1:N
)

pop = polyopt(H, registry)

config = SolverConfig(
    optimizer = nothing,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
)

M = typeof(σx[1])
T = M.parameters[2]   # integer index type
A = PauliAlgebra
MP_C = NC._moment_problem_coeff_type(A, ComplexF64)
MP_P = NC.Polynomial{A,T,MP_C}

function run_pass(pass)
    println("===== pass $pass $(pass == 1 ? "(JIT warm-up)" : "(steady-state)") =====")
    flush(stdout)

    GC.gc()
    t_sparsity = @elapsed begin
        sparsity = compute_sparsity(pop, config)
        corr = sparsity.corr_sparsity
        cts = sparsity.cliques_term_sparsities
        total_basis, _, _ = NC._polynomial_total_basis(pop, corr, cts)
    end
    corr = sparsity.corr_sparsity
    cts = sparsity.cliques_term_sparsities
    println("  1. sparsity + total basis      = $(round(t_sparsity; digits=2))s  (total basis: $(length(total_basis)) monomials)")
    flush(stdout)

    spec = sympleq_symmetry_spec(H)

    GC.gc()
    t_group = @elapsed begin
        support_domain = NC._symmetry_domain(pop, corr, cts)
        nqubits = max(
            NC._clifford_nqubits_from_domain(support_domain),
            NC._clifford_max_generator_nqubits(spec.clifford_generators),
        )
        sw_group = CliffordSymmetryGroup(spec.clifford_generators; nqubits, integer_type=T, domain=support_domain)
        group = NC.CliffordSymmetry{T}[NC._clifford_group_value(el) for el in sw_group]
    end
    println("  2. group enumeration           = $(round(t_group; digits=2))s  (order $(length(group)))")
    flush(stdout)

    GC.gc()
    t_invar = @elapsed NC._check_symmetry_invariance(pop, group)
    println("  3. invariance check            = $(round(t_invar; digits=2))s")
    flush(stdout)

    GC.gc()
    t_closure = @elapsed begin
        for basis in corr.clq_mom_mtx_bases
            NC._check_basis_closure("moment basis", basis, group)
        end
        for term_sparsities in cts, ts in term_sparsities, basis in ts.block_bases
            NC._check_basis_closure("block basis", basis, group)
        end
    end
    println("  4. basis closure checks        = $(round(t_closure; digits=2))s")
    flush(stdout)

    GC.gc()
    t_reducer = @elapsed reducer = NC._build_orbit_reducer(total_basis, group)
    println("  5. orbit reducer               = $(round(t_reducer; digits=2))s")
    flush(stdout)

    basis = only(corr.clq_mom_mtx_bases)
    println("     (moment basis size: $(length(basis)))")

    GC.gc()
    t_mat = @elapsed begin
        _, raw_mat = NC._build_constraint_matrix(one(pop.objective), basis, :HPSD)
        mat = NC._convert_polynomial_matrix(MP_P, raw_mat)
    end
    println("  6. moment matrix build         = $(round(t_mat; digits=2))s  ($(size(mat,1))x$(size(mat,2)) polynomial matrix)")
    flush(stdout)

    GC.gc()
    t_sw = @elapsed blocks = NC._sw_decompose_half_basis(basis, sw_group)
    row_bases = [Matrix(b) for b in blocks]
    println("  7. SymbolicWedderburn decompose = $(round(t_sw; digits=2))s  (block rows: $(map(b -> size(b,1), row_bases)))")
    flush(stdout)

    GC.gc()
    t_diag = @elapsed diag_blocks = NC._diagonal_transformed_blocks(mat, row_bases, reducer)
    println("  8a. diagonal transforms        = $(round(t_diag; digits=2))s  (kept blocks: $(map(b -> size(b,1), diag_blocks)))")
    flush(stdout)

    GC.gc()
    t_full = @elapsed NC._reduce_transformed_blocks(mat, row_bases, reducer)
    println("  8b. off-diag zero verification = $(round(t_full - t_diag; digits=2))s  (= full reduce $(round(t_full; digits=2))s - diagonal $(round(t_diag; digits=2))s)")
    flush(stdout)

    total = t_sparsity + t_group + t_invar + t_closure + t_reducer + t_mat + t_sw + t_full
    println("  ---")
    println("  accounted total                = $(round(total; digits=2))s")
    println()
    flush(stdout)
end

run_pass(1)
run_pass(2)

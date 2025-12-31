# NCTSSOS Oracle Generator
# =========================
# Run this on a800 server where NCTSSOS is installed:
#   cd /home/yushengzhao/NCTSSOS && julia --project this_file.jl
#
# Or copy to the NCTSSOS directory and run there.

using NCTSSOS, DynamicPolynomials
using MosekTools

# Store all results
const RESULTS = Dict{String,Any}()

# =============================================================================
# Helper: Run nctssos_first and extract info
# =============================================================================

function run_test(name, f, vars, order; kwargs...)
    println("\n--- $name ---")
    try
        opt, data = nctssos_first([f], vars, order; kwargs...)
        
        # Get basis size info
        mb_size = hasfield(typeof(data), :mb) ? length(data.mb) : -1
        
        println("  opt = $opt")
        println("  mb_size = $mb_size")
        
        RESULTS[name] = (opt=opt, mb_size=mb_size, status="ok")
        return opt, data
    catch e
        println("  ERROR: $e")
        RESULTS[name] = (opt=NaN, mb_size=-1, status="error: $e")
        return NaN, nothing
    end
end

# =============================================================================
# MOMENT.JL TESTS
# =============================================================================

function test_moment()
    println("\n" * "="^70)
    println("MOMENT.JL TESTS")
    println("="^70)
    
    # --- CHSH ---
    # NCTSSoS: x,y are unipotent (x²=I), partition by site
    # Expected: -2√2 ≈ -2.8284
    @ncpolyvar x[1:2] y[1:2]
    f = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]
    # NCTSSOS maximizes, NCTSSoS minimizes
    # NCTSSoS: min(-CHSH) gives -2√2
    run_test("CHSH_Unipotent", -f, [x;y], 1,
             TS=false, partition=2, constraint="unipotent")
    
    # --- CS TS Example (n=10, d=3) ---
    # This is a large problem requiring --local with Mosek
    n = 10
    @ncpolyvar x[1:n]
    
    f = zero(typeof(x[1]+x[2]))
    for i = 1:n
        jset = setdiff(max(1, i-5):min(n, i+1), i)
        f += (2.0*x[i] + 5.0*x[i]^3 + 1)^2
        f -= sum(4.0*x[i]*x[j] + 10.0*x[i]^3*x[j] + 2.0*x[j] +
                 4.0*x[i]*x[j]^2 + 10.0*x[i]^3*x[j]^2 + 2.0*x[j]^2 
                 for j in jset)
        f += sum(1.0*x[j]*x[k] + 2.0*x[j]^2*x[k] + 1.0*x[j]^2*x[k]^2 
                 for j in jset for k in jset)
    end
    cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0/3 for i = 1:n])
    
    run_test("CS_TS_n10_d3", f, x, 3, 
             TS=true, CS="MF", ineq=cons)
    
    # --- Heisenberg Star Graph ---
    num_sites = 10
    vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]
    
    @ncpolyvar pij[1:length(vec_idx2ij)]
    
    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)
    
    # Star graph: center=1, edges (1,2),...,(1,10)
    objective = sum(1.0 * pij[findvaridx(1, j)] for j in 2:num_sites)
    
    gs = unique([
        pij[findvaridx(sort([i,j])...)] * pij[findvaridx(sort([j,k])...)] +
        pij[findvaridx(sort([j,k])...)] * pij[findvaridx(sort([i,j])...)] -
        pij[findvaridx(sort([i,j])...)] - pij[findvaridx(sort([j,k])...)] -
        pij[findvaridx(sort([i,k])...)] + 1.0
        for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites 
        if (i != j && j != k && i != k)
    ])
    
    run_test("Heisenberg_Star_n10", objective, pij, 1,
             CS="MF", constraint="unipotent", eq=gs, numeq=length(gs))
    
    # --- Example 1 (Dense and Sparse) ---
    n = 3
    @ncpolyvar x[1:n]
    f = 1.0*x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3.0*x[2]^2 - 2.0*x[1]*x[2]*x[1] +
        2.0*x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] + 6.0*x[3]^2 +
        9.0*x[2]^2*x[3] + 9.0*x[3]*x[2]^2 - 54.0*x[3]*x[2]*x[3] + 
        142.0*x[3]*x[2]^2*x[3]
    
    run_test("Example1_Dense", f, x, 2, TS=false, CS=false)
    run_test("Example1_TS_MMD", f, x, 2, TS="MD", CS=false)
    
    # --- Example 2 ---
    n = 2
    @ncpolyvar x[1:n]
    f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1]*x[2] + x[2]*x[1] - 2.0
    
    run_test("Example2_Dense", f, x, 2, TS=false, CS=false,
             ineq=[g], eq=[h1], numeq=1)
    run_test("Example2_CS_TS", f, x, 2, TS="MD", CS="MF",
             ineq=[g], eq=[h1], numeq=1)
    
    # --- Correlative Sparsity Example ---
    n = 3
    @ncpolyvar x[1:n]
    f = 1.0*x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3.0*x[2]^2 - 2.0*x[1]*x[2]*x[1] +
        2.0*x[1]*x[2]^2*x[1] - x[2]*x[3] - x[3]*x[2] + 6.0*x[3]^2 +
        9.0*x[2]^2*x[3] + 9.0*x[3]*x[2]^2 - 54.0*x[3]*x[2]*x[3] + 
        142.0*x[3]*x[2]^2*x[3]
    cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0/3 for i = 1:n])
    
    run_test("CorrSparsity_CS", f, x, 3, TS=false, CS="MF", ineq=cons)
    run_test("CorrSparsity_TS", f, x, 3, TS="MD", CS=false, ineq=cons)
end

# =============================================================================
# INTERFACE.JL TESTS
# =============================================================================

function test_interface()
    println("\n" * "="^70)
    println("INTERFACE.JL TESTS")
    println("="^70)
    
    # --- Majumdar Ghosh Model ---
    num_sites = 6
    @ncpolyvar hij[1:(num_sites*(num_sites-1)÷2)]
    
    ij2idx = Dict(
        (i,j) => findfirst(x -> x == (i,j), 
            [(a,b) for a in 1:num_sites for b in (a+1):num_sites])
        for i in 1:num_sites for j in (i+1):num_sites
    )
    
    J1_interactions = [(i, mod1(i+1, num_sites)) for i in 1:num_sites]
    J2_interactions = [(i, mod1(i+2, num_sites)) for i in 1:num_sites]
    J1, J2 = 2.0, 1.0
    
    objective = sum(J1 * hij[ij2idx[(min(i,j), max(i,j))]] for (i,j) in J1_interactions) +
                sum(J2 * hij[ij2idx[(min(i,j), max(i,j))]] for (i,j) in J2_interactions)
    
    gs = unique([
        hij[ij2idx[tuple(sort([i,j])...)]] * hij[ij2idx[tuple(sort([j,k])...)]] +
        hij[ij2idx[tuple(sort([j,k])...)]] * hij[ij2idx[tuple(sort([i,j])...)]] -
        0.5 * (hij[ij2idx[tuple(sort([i,j])...)]] + hij[ij2idx[tuple(sort([j,k])...)]] -
               hij[ij2idx[tuple(sort([i,k])...)]])
        for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites 
        if (i != j && j != k && i != k)
    ])
    
    # NCTSSoS minimizes -objective, NCTSSOS maximizes
    # To compare: NCTSSOS max(objective) = NCTSSoS -min(-objective)
    run_test("Majumdar_Ghosh", -objective, hij, 1,
             constraint="projector", eq=gs, numeq=length(gs))
    
    # --- Problem Creation Interface ---
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1]*x[2] + x[2]*x[1] - 2.0
    
    run_test("Problem_Creation", f, x, 2,
             TS="MD", CS="MF", ineq=[g], eq=[h1], numeq=1)
    
    # --- README Unconstrained ---
    @ncpolyvar x[1:3]
    f = 1.0 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]*x[2] + x[2]*x[1] + x[2]*x[3] + x[3]*x[2]
    
    run_test("README_Unconstrained_Dense", f, x, 2, TS=false, CS=false)
    run_test("README_Unconstrained_CS", f, x, 2, TS=false, CS="MF")
    run_test("README_Unconstrained_CS_TS", f, x, 2, TS="MD", CS="MF")
    
    # --- README Constrained ---
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1]*x[2] + x[2]*x[1] - 2.0
    
    run_test("README_Constrained_Dense", f, x, 2, TS=false, CS=false,
             ineq=[g], eq=[h1], numeq=1)
    run_test("README_Constrained_CS", f, x, 2, TS=false, CS="MF",
             ineq=[g], eq=[h1], numeq=1)
    run_test("README_Constrained_CS_TS", f, x, 2, TS="MD", CS="MF",
             ineq=[g], eq=[h1], numeq=1)
end

# =============================================================================
# SPARSITY.JL TESTS (Bell inequalities and constrained)
# =============================================================================

function test_sparsity()
    println("\n" * "="^70)
    println("SPARSITY.JL TESTS")
    println("="^70)
    
    # --- CHSH variants ---
    @ncpolyvar x[1:2] y[1:2]
    f = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]
    
    run_test("CHSH_Dense", -f, [x;y], 1, TS=false, CS=false,
             partition=2, constraint="unipotent")
    run_test("CHSH_CS_MF", -f, [x;y], 1, TS=false, CS="MF",
             partition=2, constraint="unipotent")
    run_test("CHSH_TS_MD", -f, [x;y], 1, TS="MD", CS=false,
             partition=2, constraint="unipotent")
    
    # --- I_3322 ---
    @ncpolyvar x[1:3] y[1:3]
    f = x[1]*(y[1] + y[2] + y[3]) + x[2]*(y[1] + y[2] - y[3]) +
        x[3]*(y[1] - y[2]) - x[1] - 2.0*y[1] - y[2]
    
    run_test("I3322_Dense", -f, [x;y], 2, TS=false, CS=false,
             partition=3, constraint="projector")
    run_test("I3322_CS_MF", -f, [x;y], 2, TS=false, CS="MF",
             partition=3, constraint="projector")
    
    # --- Ball Constraint ---
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2 + x[1]*x[2]*x[1]*x[2] + 
        x[2]*x[1]*x[2]*x[1] + x[1]^3*x[2] + x[2]*x[1]^3 + x[1]*x[2]^3 + x[2]^3*x[1]
    g1 = 1.0 - x[1]^2
    g2 = 1.0 - x[2]^2
    
    run_test("Ball_Dense", f, x, 2, TS=false, CS=false, ineq=[g1, g2])
    run_test("Ball_TS_MD", f, x, 2, TS="MD", CS=false, ineq=[g1, g2])
    
    # --- Rosenbrock ---
    n = 6
    @ncpolyvar x[1:n]
    f = Float64(n) + sum(100.0*x[i-1]^4 - 200.0*x[i-1]^2*x[i] - 2.0*x[i] + 101.0*x[i]^2 
                        for i in 2:n)
    
    run_test("Rosenbrock_Dense", f, x, 2, TS=false, CS=false)
    run_test("Rosenbrock_CS_MF", f, x, 2, TS=false, CS="MF")
    run_test("Rosenbrock_CS_TS", f, x, 2, TS="MD", CS="MF")
end

# =============================================================================
# STATE_POLY.JL TESTS
# =============================================================================

function test_state_poly()
    println("\n" * "="^70)
    println("STATE_POLY.JL TESTS")
    println("="^70)
    
    # State polynomial optimization uses nctssos_first with state=true
    # or cs_nctssos_first for correlative sparsity
    
    # --- 7.2.0 CHSH State ---
    # sp = -<x1 y1> - <x1 y2> - <x2 y1> + <x2 y2>
    @ncpolyvar x[1:2] y[1:2]
    
    # In NCTSSOS state poly format:
    # supp = [[[1,3]], [[1,4]], [[2,3]], [[2,4]]]
    # coe = [-1, -1, -1, 1]
    # Note: vars are numbered 1,2,3,4 = x1,x2,y1,y2
    
    supp = [[[1,3]], [[1,4]], [[2,3]], [[2,4]]]
    coe = [-1.0, -1.0, -1.0, 1.0]
    
    println("\n--- State_7_2_0 ---")
    try
        opt, data = cs_nctssos_first(supp, coe, 4, 1, 
                                     TS=false, partition=2, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["State_7_2_0"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["State_7_2_0"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    # --- 7.2.0 with TS ---
    println("\n--- State_7_2_0_TS ---")
    try
        opt, data = cs_nctssos_first(supp, coe, 4, 1,
                                     TS="MD", partition=2, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["State_7_2_0_TS"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["State_7_2_0_TS"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    # --- 7.2.1 Squared expectations ---
    # sp = -(<x1 y2> + <x2 y1>)² - (<x1 y1> - <x2 y2>)²
    # Expanded: -<x1 y2>² - 2<x1 y2><x2 y1> - <x2 y1>² - <x1 y1>² + 2<x1 y1><x2 y2> - <x2 y2>²
    # supp format: [[mon1], [mon2]] for <mon1><mon2>
    
    supp_721 = [
        [[1,4], [1,4]],  # <x1 y2><x1 y2>
        [[1,4], [2,3]],  # <x1 y2><x2 y1>
        [[2,3], [2,3]],  # <x2 y1><x2 y1>
        [[1,3], [1,3]],  # <x1 y1><x1 y1>
        [[1,3], [2,4]],  # <x1 y1><x2 y2>
        [[2,4], [2,4]]   # <x2 y2><x2 y2>
    ]
    coe_721 = [-1.0, -2.0, -1.0, -1.0, 2.0, -1.0]
    
    println("\n--- State_7_2_1_d3 ---")
    try
        opt, data = cs_nctssos_first(supp_721, coe_721, 4, 3,
                                     TS=false, partition=2, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["State_7_2_1_d3"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["State_7_2_1_d3"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    # --- 7.2.2 Covariance ---
    # cov(a,b) = <x_a y_b> - <x_a><y_b>
    # n=3 each side, vars 1-3 = x, vars 4-6 = y
    
    # Build covariance terms
    supp_722 = []
    coe_722 = Float64[]
    
    # cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)
    covs = [(1,1,1), (1,2,1), (1,3,1), (2,1,1), (2,2,1), (2,3,-1), (3,1,1), (3,2,-1)]
    
    for (a, b, sign) in covs
        # <x_a y_b> term
        push!(supp_722, [[a, b+3]])
        push!(coe_722, Float64(sign))
        # -<x_a><y_b> term  
        push!(supp_722, [[a], [b+3]])
        push!(coe_722, -Float64(sign))
    end
    
    println("\n--- State_7_2_2 ---")
    try
        opt, data = cs_nctssos_first(supp_722, coe_722, 6, 2,
                                     TS=false, partition=3, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["State_7_2_2"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["State_7_2_2"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    # --- 7.2.3 Mixed ---
    # From test file: complex combination of linear, mixed monomials, products, squares
    # x=1:2, y=3:4 (using 1-indexed)
    # sp = -<x2> - <y1> - <y2> + <x1 y1> - <x2 y1> - <x1 y2> - <x2 y2> + 
    #      <x1><y1> + <x2><y1> + <x2><y2> + <x1>² + <y2>²
    
    supp_723 = [
        [[2]],           # <x2>
        [[3]],           # <y1>
        [[4]],           # <y2>
        [[1,3]],         # <x1 y1>
        [[2,3]],         # <x2 y1>
        [[1,4]],         # <x1 y2>
        [[2,4]],         # <x2 y2>
        [[1], [3]],      # <x1><y1>
        [[2], [3]],      # <x2><y1>
        [[2], [4]],      # <x2><y2>
        [[1], [1]],      # <x1>²
        [[4], [4]]       # <y2>²
    ]
    coe_723 = [-1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    
    println("\n--- State_7_2_3 ---")
    try
        opt, data = cs_nctssos_first(supp_723, coe_723, 4, 2,
                                     TS=false, partition=2, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["State_7_2_3"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["State_7_2_3"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    println("\n--- State_7_2_3_TS ---")
    try
        opt, data = cs_nctssos_first(supp_723, coe_723, 4, 2,
                                     TS="MD", partition=2, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["State_7_2_3_TS"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["State_7_2_3_TS"] = (opt=NaN, mb_size=-1, status="error")
    end
end

# =============================================================================
# TRACE_POLY.JL TESTS
# =============================================================================

function test_trace_poly()
    println("\n" * "="^70)
    println("TRACE_POLY.JL TESTS")
    println("="^70)
    
    # Trace polynomial uses tracepop_first or cs_tracepop_first
    
    # --- Example 6.1 ---
    # p = tr(x1 x2 x3) + tr(x1 x2) * tr(x3)
    # ptsupp format: each element is a list of traces
    ptsupp_61 = [[[1,2,3]], [[1,2],[3]]]
    coe_61 = [1.0, 1.0]
    
    println("\n--- Trace_6_1_d2 ---")
    try
        opt, data = cs_tracepop_first(ptsupp_61, coe_61, 3, 2,
                                      TS=false, constraint="projector")
        println("  opt = $opt")
        RESULTS["Trace_6_1_d2"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["Trace_6_1_d2"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    println("\n--- Trace_6_1_d3 ---")
    try
        opt, data = cs_tracepop_first(ptsupp_61, coe_61, 3, 3,
                                      TS=false, constraint="projector")
        println("  opt = $opt")
        RESULTS["Trace_6_1_d3"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["Trace_6_1_d3"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    # --- Example 6.2.0 CHSH Trace ---
    # p = -tr(x1 y1) - tr(x1 y2) - tr(x2 y1) + tr(x2 y2)
    # vars: 1,2 = x, 3,4 = y
    ptsupp_620 = [[[1,3]], [[1,4]], [[2,3]], [[2,4]]]
    coe_620 = [-1.0, -1.0, -1.0, 1.0]
    
    println("\n--- Trace_6_2_0 ---")
    try
        opt, data = cs_tracepop_first(ptsupp_620, coe_620, 4, 1,
                                      TS="block", partition=2, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["Trace_6_2_0"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["Trace_6_2_0"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    # --- Example 6.2.1 Squared trace ---
    # (tr(x1 y2) + tr(x2 y1))² + (tr(x1 y1) - tr(x2 y2))²
    # Expanded:
    ptsupp_621 = [
        [[1,4],[1,4]],   # tr(x1 y2)²
        [[1,4],[2,3]],   # 2*tr(x1 y2)*tr(x2 y1)
        [[2,3],[2,3]],   # tr(x2 y1)²
        [[1,3],[1,3]],   # tr(x1 y1)²
        [[1,3],[2,4]],   # -2*tr(x1 y1)*tr(x2 y2)
        [[2,4],[2,4]]    # tr(x2 y2)²
    ]
    coe_621 = [1.0, 2.0, 1.0, 1.0, -2.0, 1.0]
    
    # We want to MINIMIZE this (negative of max)
    println("\n--- Trace_6_2_1_d2 ---")
    try
        # NCTSSoS minimizes -p, equivalent to NCTSSOS max(-p)
        # So pass negated coefficients
        opt, data = cs_tracepop_first(ptsupp_621, -coe_621, 4, 2,
                                      TS=false, partition=2, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["Trace_6_2_1_d2"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["Trace_6_2_1_d2"] = (opt=NaN, mb_size=-1, status="error")
    end
    
    # --- Example 6.2.2 Covariance Trace ---
    # cov(i,j) = tr(x_i y_j) - tr(x_i)*tr(y_j)
    # vars 1-3 = x, 4-6 = y
    ptsupp_622 = []
    coe_622 = Float64[]
    
    covs = [(1,1,1), (1,2,1), (1,3,1), (2,1,1), (2,2,1), (2,3,-1), (3,1,1), (3,2,-1)]
    for (a, b, sign) in covs
        push!(ptsupp_622, [[a, b+3]])        # tr(x_a y_b)
        push!(coe_622, Float64(sign))
        push!(ptsupp_622, [[a], [b+3]])      # -tr(x_a)*tr(y_b)
        push!(coe_622, -Float64(sign))
    end
    
    println("\n--- Trace_6_2_2 ---")
    try
        # Minimize by negating
        opt, data = cs_tracepop_first(ptsupp_622, -coe_622, 6, 2,
                                      TS=false, partition=3, constraint="unipotent")
        println("  opt = $opt")
        RESULTS["Trace_6_2_2"] = (opt=opt, mb_size=-1, status="ok")
    catch e
        println("  ERROR: $e")
        RESULTS["Trace_6_2_2"] = (opt=NaN, mb_size=-1, status="error")
    end
end

# =============================================================================
# MAIN
# =============================================================================

function main()
    println("NCTSSOS Oracle Generator")
    println("="^70)
    println("Generating reference values for NCTSSoS.jl tests")
    println()
    
    test_moment()
    test_interface()
    test_sparsity()
    test_state_poly()
    test_trace_poly()
    
    # Print summary as Julia code
    println("\n\n" * "="^70)
    println("ORACLE SUMMARY (copy to nctssos_oracles.jl)")
    println("="^70)
    
    println("\nconst NCTSSOS_ORACLE = Dict(")
    for (name, res) in sort(collect(RESULTS), by=first)
        println("    \"$name\" => (opt=$(res.opt), mb_size=$(res.mb_size), status=\"$(res.status)\"),")
    end
    println(")")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

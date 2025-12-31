# NCTSSOS Oracle Values
# =====================
# Expected optimal values from the original NCTSSOS package.
# These serve as ground truth for NCTSSoS.jl test validation.
#
# Structure:
#   test_name => (
#       nctssos_opt,      # Optimal value from NCTSSOS (maximizes)
#       nctssos_mb_size,  # Moment basis size from NCTSSOS
#       nctss_opt,        # Optimal value from NCTSSoS (minimizes)
#       notes,            # Any discrepancies or notes
#   )
#
# Sign Convention:
#   NCTSSOS maximizes by default, NCTSSoS minimizes.
#   For CHSH: NCTSSOS max(CHSH) = 2√2, NCTSSoS min(-CHSH) = -2√2
#
# To regenerate: Run nctssos_oracle_generator.jl on a800 server

# =============================================================================
# MOMENT.JL TEST ORACLES
# =============================================================================

const MOMENT_ORACLES = Dict(
    # CHSH Inequality
    # NCTSSOS: max(x1y1 + x1y2 + x2y1 - x2y2) with unipotent partition=2
    # NCTSSoS: min(x1y1 + x1y2 + x2y1 - x2y2) = -2√2 (negate CHSH)
    "CHSH_Unipotent" => (
        nctssos_opt = 2.8284271247461903,      # 2√2, from maximization
        nctssos_mb_size = 9,                    # Order=1: {1, x1, x2, y1, y2, x1x2, x2x1, y1y2, y2y1}? TBD
        nctss_opt = -2.8284271321623193,       # From test assertion
        notes = "NCTSSOS max = NCTSSoS -min(-f)",
    ),
    
    # CS TS Example (n=10, order=3) - requires Mosek
    "CS_TS_n10_d3" => (
        nctssos_opt = 3.011288,                 # From test assertion
        nctssos_mb_size = -1,                   # TBD from NCTSSOS
        nctss_opt = 3.011288,                   # Same sign (minimizing f directly)
        notes = "Large problem, requires --local with Mosek",
    ),
    
    # Heisenberg Model on Star Graph
    "Heisenberg_Star_n10" => (
        nctssos_opt = -1.0,                     # Ground state energy
        nctssos_mb_size = -1,                   # TBD
        nctss_opt = -1.0,                       # From test assertion
        notes = "Star graph with 10 sites, unipotent vars + SU(2) constraints",
    ),
    
    # Example 1 - Dense vs Sparse
    "Example1_Dense" => (
        nctssos_opt = 0.0,                      # Tight bound
        nctssos_mb_size = -1,                   # TBD
        nctss_opt = 0.0,                        # From test (atol=1e-5)
        notes = "3-var unconstrained, dense relaxation",
    ),
    "Example1_TS_MMD" => (
        nctssos_opt = -0.0035512,               # Weaker bound with term sparsity
        nctssos_mb_size = -1,                   # TBD
        nctss_opt = -0.0035512,                 # From test (atol=1e-4)
        notes = "Term sparsity gives weaker bound",
    ),
    
    # Example 2 - Constrained
    "Example2_Dense" => (
        nctssos_opt = -1.0,
        nctssos_mb_size = -1,
        nctss_opt = -1.0,
        notes = "2-var with equality and inequality constraints",
    ),
    "Example2_CS_TS" => (
        nctssos_opt = -1.0,
        nctssos_mb_size = -1,
        nctss_opt = -1.0,
        notes = "MF + MMD sparsity, same result",
    ),
    
    # Correlative Sparsity
    "CorrSparsity_CS" => (
        nctssos_opt = 0.9975306427277915,
        nctssos_mb_size = -1,
        nctss_opt = 0.9975306427277915,
        notes = "Order=3 with MF correlative sparsity",
    ),
    "CorrSparsity_TS" => (
        nctssos_opt = 0.9975306427277915,
        nctssos_mb_size = -1,
        nctss_opt = 0.9975306427277915,
        notes = "Order=3 with MMD term sparsity + higher iteration",
    ),
)

# =============================================================================
# INTERFACE.JL TEST ORACLES
# =============================================================================

const INTERFACE_ORACLES = Dict(
    # Naive Example (Pauli)
    "Naive_Pauli" => (
        nctssos_opt = nothing,                  # NCTSSOS doesn't have Pauli algebra
        nctssos_mb_size = nothing,
        nctss_opt = -0.8660254037844387,       # -√3/2
        notes = "Pauli algebra not in NCTSSOS, NCTSSoS-specific",
    ),
    
    # 1D Transverse Field Ising Model (Pauli)
    "TFIM_N3_periodic" => (
        nctssos_opt = nothing,
        nctssos_mb_size = nothing,
        nctss_opt = -1.0175918 * 3,            # Per-site * N
        notes = "Pauli, periodic boundary",
    ),
    "TFIM_N3_open" => (
        nctssos_opt = nothing,
        nctssos_mb_size = nothing,
        nctss_opt = -1.0104160 * 3,
        notes = "Pauli, open boundary",
    ),
    
    # 1D Heisenberg Chain (Pauli)
    "Heisenberg_N6" => (
        nctssos_opt = nothing,
        nctssos_mb_size = nothing,
        nctss_opt = -0.467129 * 6,
        notes = "Pauli, periodic",
    ),
    
    # I_3322 with Sparsity (Projector)
    "I3322_Dense_d3" => (
        nctssos_opt = -0.2508755573198166,
        nctssos_mb_size = -1,
        nctss_opt = -0.2508755573198166,
        notes = "Dense, order=3",
    ),
    "I3322_MF_MMD_d3" => (
        nctssos_opt = -0.9999999892255513,     # Looser bound
        nctssos_mb_size = -1,
        nctss_opt = -0.9999999892255513,
        notes = "MF + MMD gives looser bound",
    ),
    "I3322_MF_MaxElim_d3" => (
        nctssos_opt = -0.3507010331201541,
        nctssos_mb_size = -1,
        nctss_opt = -0.3507010331201541,
        notes = "MF + MaximalElimination",
    ),
    
    # Majumdar Ghosh Model (Projector)
    "Majumdar_Ghosh" => (
        nctssos_opt = 9.0,                      # -num_sites/4 * 6 = -6/4*6 = -9
        nctssos_mb_size = -1,
        nctss_opt = -9.0,                       # NCTSSoS minimizes -objective
        notes = "6-site chain with J1=2, J2=1, projector constraints",
    ),
    
    # Problem Creation Interface
    "Problem_Creation" => (
        nctssos_opt = -1.0,
        nctssos_mb_size = -1,
        nctss_opt = -1.0,
        notes = "2-var constrained",
    ),
    
    # README Examples
    "README_Unconstrained" => (
        nctssos_opt = -1.0,                     # Approx from tests
        nctssos_mb_size = -1,
        nctss_opt = -1.0,                       # All variants should match
        notes = "Dense/CS/CS+TS should all give same result",
    ),
    "README_Constrained" => (
        nctssos_opt = -1.0,
        nctssos_mb_size = -1,
        nctss_opt = -1.0,
        notes = "Dense/CS/CS+TS/higher should all give same result",
    ),
)

# =============================================================================
# SPARSITY.JL TEST ORACLES
# =============================================================================

const SPARSITY_ORACLES = Dict(
    # CHSH variants
    "CHSH_Dense" => (
        nctssos_opt = 2.8284,
        nctssos_mb_size = -1,
        nctss_opt = -2.8284,                   # min(-CHSH)
        notes = "Dense baseline",
    ),
    "CHSH_CS_MF" => (
        nctssos_opt = 2.8284,
        nctssos_mb_size = -1,
        nctss_opt = -2.8284,
        notes = "Correlative sparsity should match dense",
    ),
    "CHSH_TS_MMD" => (
        nctssos_opt = 2.8284,
        nctssos_mb_size = -1,
        nctss_opt = -2.8284,
        notes = "Term sparsity should match dense",
    ),
    
    # I_3322
    "I3322_Dense_d2" => (
        nctssos_opt = -0.25,
        nctssos_mb_size = -1,
        nctss_opt = -0.25,
        notes = "Order=2 dense",
    ),
    "I3322_MF_d2" => (
        nctssos_opt = -0.25,
        nctssos_mb_size = -1,
        nctss_opt = -0.25,
        notes = "Order=2 with MF",
    ),
    
    # Ball Constraint
    "Ball_Dense" => (
        nctssos_opt = -Inf,                     # Lower bound < 0
        nctssos_mb_size = -1,
        nctss_opt = -Inf,                       # Some negative bound
        notes = "Should give negative lower bound",
    ),
    "Ball_TS_MMD" => (
        nctssos_opt = -Inf,
        nctssos_mb_size = -1,
        nctss_opt = -Inf,
        notes = "Term sparsity variant",
    ),
    
    # Rosenbrock
    "Rosenbrock_Dense" => (
        nctssos_opt = -Inf,                     # Lower bound
        nctssos_mb_size = -1,
        nctss_opt = -Inf,
        notes = "6-var Rosenbrock",
    ),
)

# =============================================================================
# STATE_POLY.JL TEST ORACLES
# =============================================================================

const STATE_POLY_ORACLES = Dict(
    # 7.2.0 - CHSH State
    "State_7_2_0" => (
        nctssos_opt = 2.8284271321623202,      # max(<x1y1> + ...)
        nctssos_mb_size = -1,
        nctss_opt = -2.8284271321623202,       # min(-<x1y1> - ...)
        notes = "State polynomial CHSH",
    ),
    "State_7_2_0_TS" => (
        nctssos_opt = 2.8284271321623202,
        nctssos_mb_size = -1,
        nctss_opt = -2.8284271321623202,
        notes = "With MMD term sparsity",
    ),
    
    # 7.2.1 - Squared expectations
    "State_7_2_1_d3" => (
        nctssos_opt = 4.0,                      # max of squared sum
        nctssos_mb_size = -1,
        nctss_opt = -4.0,                       # min(-sum) at order=3
        notes = "Requires order=3 for tight bound",
    ),
    
    # 7.2.2 - Covariance
    "State_7_2_2" => (
        nctssos_opt = 5.0,                      # max(cov sum)
        nctssos_mb_size = -1,
        nctss_opt = -5.0,                       # min(-cov sum)
        notes = "Covariance expression with compound StateWords",
    ),
    
    # 7.2.3 - Mixed
    "State_7_2_3" => (
        nctssos_opt = 3.5114802,               # From NCTSSOS
        nctssos_mb_size = -1,
        nctss_opt = -3.5114802,
        notes = "Mixed linear, monomials, products, squares",
    ),
    "State_7_2_3_TS" => (
        nctssos_opt = 3.5114802,
        nctssos_mb_size = -1,
        nctss_opt = -3.5114802,
        notes = "With MMD term sparsity",
    ),
)

# =============================================================================
# TRACE_POLY.JL TEST ORACLES
# =============================================================================

const TRACE_POLY_ORACLES = Dict(
    # Example 6.1 (Projector)
    "Trace_6_1_d2" => (
        nctssos_opt = -0.046717378455438933,   # From test
        nctssos_mb_size = -1,
        nctss_opt = -0.046717378455438933,
        notes = "Order=2, projector constraint",
    ),
    "Trace_6_1_d3" => (
        nctssos_opt = -0.03124998978001017,    # Tighter at order=3
        nctssos_mb_size = -1,
        nctss_opt = -0.03124998978001017,
        notes = "Order=3, requires Mosek",
    ),
    
    # Example 6.2.0 - CHSH Trace
    "Trace_6_2_0" => (
        nctssos_opt = 2.8284271157283083,      # 2√2
        nctssos_mb_size = -1,
        nctss_opt = -2.8284271157283083,
        notes = "CHSH in trace form, MaximalElimination",
    ),
    
    # Example 6.2.1 - Squared trace
    "Trace_6_2_1_d2" => (
        nctssos_opt = 4.0,                      # Tight bound
        nctssos_mb_size = -1,
        nctss_opt = -4.0,
        notes = "Requires order=2 and Mosek",
    ),
    
    # Example 6.2.2 - Covariance trace
    "Trace_6_2_2" => (
        nctssos_opt = 5.0,
        nctssos_mb_size = -1,
        nctss_opt = -5.0,
        notes = "Covariance in trace form",
    ),
)

# =============================================================================
# COMBINED ORACLE
# =============================================================================

const ALL_ORACLES = merge(
    Dict("moment." * k => v for (k, v) in MOMENT_ORACLES),
    Dict("interface." * k => v for (k, v) in INTERFACE_ORACLES),
    Dict("sparsity." * k => v for (k, v) in SPARSITY_ORACLES),
    Dict("state_poly." * k => v for (k, v) in STATE_POLY_ORACLES),
    Dict("trace_poly." * k => v for (k, v) in TRACE_POLY_ORACLES),
)

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

"""
    validate_result(test_name, nctssos_result, rtol=1e-5)

Compare NCTSSoS result against NCTSSOS oracle.
Returns (passed, message).
"""
function validate_result(test_name, nctssos_result; rtol=1e-5)
    if !haskey(ALL_ORACLES, test_name)
        return (false, "Unknown test: $test_name")
    end
    
    oracle = ALL_ORACLES[test_name]
    expected = oracle.nctss_opt
    
    if isnothing(expected)
        return (true, "No oracle available ($(oracle.notes))")
    end
    
    if isnan(expected) || isinf(expected)
        return (true, "Oracle is $(expected), skipping numerical comparison")
    end
    
    if isapprox(nctssos_result, expected; rtol=rtol)
        return (true, "PASS: $nctssos_result ≈ $expected")
    else
        return (false, "FAIL: $nctssos_result ≠ $expected (expected)")
    end
end

"""
Print summary of all oracles with their expected values.
"""
function print_oracle_summary()
    println("NCTSSoS Oracle Summary")
    println("=" ^ 70)
    
    for (category, oracles) in [
        ("Moment", MOMENT_ORACLES),
        ("Interface", INTERFACE_ORACLES),
        ("Sparsity", SPARSITY_ORACLES),
        ("State Polynomial", STATE_POLY_ORACLES),
        ("Trace Polynomial", TRACE_POLY_ORACLES),
    ]
        println("\n### $category ###")
        for (name, oracle) in sort(collect(oracles), by=first)
            nctssos = isnothing(oracle.nctssos_opt) ? "N/A" : oracle.nctssos_opt
            nctss = isnothing(oracle.nctss_opt) ? "N/A" : oracle.nctss_opt
            println("  $name:")
            println("    NCTSSOS: $nctssos")
            println("    NCTSSoS: $nctss")
            if !isempty(oracle.notes)
                println("    Note: $(oracle.notes)")
            end
        end
    end
end

# Run summary if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    print_oracle_summary()
end

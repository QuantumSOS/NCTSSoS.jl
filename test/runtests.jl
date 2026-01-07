# =============================================================================
# NCTSSoS.jl Test Suite
# =============================================================================
#
# Quick Reference:
# ----------------
# | Command              | What                        | Time   |
# |----------------------|-----------------------------|--------|
# | make test-minimal    | Minimal suite (5 paths)     | ~25s   |
# | make test-polynomials| Algebra only (no solver)    | ~10s   |
# | make test-ci         | CI default (minimal + poly) | ~35s   |
# | make test            | Full suite (Mosek)          | ~5min  |
#
# Test Structure:
# ---------------
# test/
# ├── runtests.jl        - This file (entry point, parses flags)
# ├── test_minimal.jl    - Minimal test suite (includes from problems/)
# ├── TestUtils.jl       - SOLVER config + helpers
# ├── polynomials/       - --polynomials (no solver needed)
# ├── relaxations/       - --relaxations
# │   └── dualization.jl - SOS ≈ Moment equivalence
# ├── problems/          - --problems
# │   ├── bell_inequalities/
# │   │   ├── chsh_simple.jl      - Dense, CS, TS (order=1)
# │   │   ├── chsh_high_order.jl  - CS+TS (order=2)
# │   │   ├── chsh_state.jl       - State polynomial
# │   │   └── chsh_trace.jl       - Trace polynomial
# │   ├── nc_polynomial/
# │   │   ├── nc_example1.jl      - Unconstrained (3 vars)
# │   │   ├── nc_example2.jl      - Constrained (2 vars)
# │   │   └── nc_correlative.jl   - Correlative sparsity
# │   └── ...
# └── quality/           - --quality (Aqua, ExplicitImports, Doctest)
#
# Flags:
# ------
# --minimal      Fast correctness check (5 algorithm paths)
# --polynomials  Core algebra (no solver)
# --relaxations  SOS, sparsity components
# --problems     Problem-based tests
# --quality      Code quality checks
# --local        Use Mosek instead of COSMO
#
# CI Default: --minimal + --polynomials (no flags = full suite)
# =============================================================================

using NCTSSoS, Test

# =============================================================================
# Parse test arguments
# =============================================================================
const RUN_POLYNOMIALS = "--polynomials" in ARGS
const RUN_QUALITY = "--quality" in ARGS
const RUN_RELAXATIONS = "--relaxations" in ARGS
const RUN_PROBLEMS = "--problems" in ARGS
const RUN_MINIMAL = "--minimal" in ARGS
const USE_LOCAL = "--local" in ARGS

# If no specific test group flags, run all (some problem tests only with --local)
const RUN_ALL = !(RUN_POLYNOMIALS || RUN_QUALITY || RUN_RELAXATIONS || RUN_PROBLEMS || RUN_MINIMAL)

# Determine which test groups to run
const SHOULD_RUN_POLYNOMIALS = RUN_ALL || RUN_POLYNOMIALS
const SHOULD_RUN_QUALITY = RUN_ALL || RUN_QUALITY
const SHOULD_RUN_RELAXATIONS = RUN_ALL || RUN_RELAXATIONS
const SHOULD_RUN_PROBLEMS = RUN_ALL || RUN_PROBLEMS

# Print test plan
@info "Test configuration" USE_LOCAL ARGS
@info "Running tests" polynomials=SHOULD_RUN_POLYNOMIALS quality=SHOULD_RUN_QUALITY relaxations=SHOULD_RUN_RELAXATIONS problems=SHOULD_RUN_PROBLEMS minimal=RUN_MINIMAL

# =============================================================================
# Run selected test groups
# =============================================================================
@testset "NCTSSoS.jl" begin

    # =========================================================================
    # 1. Polynomial Algebra Tests (no solver needed, no JuMP imports)
    # =========================================================================
    if SHOULD_RUN_POLYNOMIALS
        @testset "Polynomials" begin
            include("polynomials/runtests.jl")
        end
    end

    # =========================================================================
    # 2. Code Quality Checks
    # =========================================================================
    if SHOULD_RUN_QUALITY
        @testset "Quality" begin
            include("quality/runtests.jl")
        end
    end

    # =========================================================================
    # Load solver configuration BEFORE relaxation/problem tests
    # (after polynomial tests to avoid JuMP's simplify shadowing NCTSSoS.simplify)
    # =========================================================================
    if SHOULD_RUN_RELAXATIONS || SHOULD_RUN_PROBLEMS
        include("TestUtils.jl")
    end

    # =========================================================================
    # 3. Relaxation Component Tests
    # =========================================================================
    if SHOULD_RUN_RELAXATIONS
        @testset "Relaxations" begin
            include("relaxations/runtests.jl")
        end
    end

    # =========================================================================
    # 4. Problem-Based Tests (--local recommended for full suite)
    # =========================================================================
    if SHOULD_RUN_PROBLEMS
        @testset "Problems" begin
            include("problems/runtests.jl")
        end
    end

    # =========================================================================
    # 5. Minimal Suite (fast smoke test for critical algorithm paths)
    # =========================================================================
    # Covers: Dense, CS, TS, constrained, dualization - DRY via includes
    # =========================================================================
    if RUN_MINIMAL
        include("TestUtils.jl")
        include("test_minimal.jl")
    end
end

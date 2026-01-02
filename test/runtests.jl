# =============================================================================
# NCTSSoS.jl Test Suite
# =============================================================================
#
# Test Structure:
# ---------------
# test/
# ├── oracles/         - Oracle scripts for generating test expectations (NCTSSOS)
# │   ├── scripts/     - NCTSSOS scripts
# │   └── results/     - Generated oracle values
# ├── polynomials/     - Core polynomial algebra (types, arithmetic, simplification)
# │   └── runtests.jl
# ├── relaxations/     - Relaxation algorithm components (SOS, sparsity, GNS)
# │   └── runtests.jl
# ├── problems/        - Problem-based tests organized by domain
# │   ├── bell_inequalities/   - CHSH, I_3322, Bell inequalities
# │   ├── condensed_matter/    - Heisenberg, Ising, XY, Bose-Hubbard, PXP
# │   ├── nc_polynomial/       - NC polynomial examples
# │   ├── state_polynomial/    - State polynomial examples (7.2.x)
# │   ├── trace_polynomial/    - Trace polynomial examples (6.x)
# │   ├── quantum_networks/    - Bilocal networks
# │   ├── fermionic/           - Fermionic systems
# │   ├── benchmarks/          - Classical optimization benchmarks
# │   └── runtests.jl
# ├── quality/         - Code quality checks (Aqua, ExplicitImports, Doctest)
# │   └── runtests.jl
# ├── setup.jl         - Shared solver configuration
# └── runtests.jl      - This file (main entry point)
#
# Solver Configuration:
# ---------------------
# - --local  → Mosek (commercial, fast, required for physics/condensed matter tests)
# - Default  → COSMO (open-source, sufficient for basic tests)
#
# Run commands:
# -------------
# CI suite (default):   Pkg.test("NCTSSoS")
# Full suite:           Pkg.test("NCTSSoS"; test_args=["--local"])
# Subset tests:         Pkg.test("NCTSSoS"; test_args=["--polynomials"])
#                       Pkg.test("NCTSSoS"; test_args=["--relaxations"])
#                       Pkg.test("NCTSSoS"; test_args=["--problems"])
#                       Pkg.test("NCTSSoS"; test_args=["--problems", "--local"])
#                       Pkg.test("NCTSSoS"; test_args=["--quality"])
# Multiple subsets:     Pkg.test("NCTSSoS"; test_args=["--polynomials", "--relaxations"])
#
# Makefile shortcuts:
# -------------------
# make test                  - Full suite with Mosek (--local)
# make test-ci               - CI suite (no condensed matter/fermionic, uses COSMO)
# make test-polynomials      - Polynomial tests only
# make test-relaxations      - Relaxation component tests only
# make test-problems         - Problem-based tests only
# make test-quality          - Code quality checks only
# =============================================================================

using NCTSSoS, Test

# =============================================================================
# Parse test arguments
# =============================================================================
const RUN_POLYNOMIALS = "--polynomials" in ARGS
const RUN_QUALITY = "--quality" in ARGS
const RUN_RELAXATIONS = "--relaxations" in ARGS
const RUN_PROBLEMS = "--problems" in ARGS
const USE_LOCAL = "--local" in ARGS

# If no specific test group flags, run all (some problem tests only with --local)
const RUN_ALL = !(RUN_POLYNOMIALS || RUN_QUALITY || RUN_RELAXATIONS || RUN_PROBLEMS)

# Determine which test groups to run
const SHOULD_RUN_POLYNOMIALS = RUN_ALL || RUN_POLYNOMIALS
const SHOULD_RUN_QUALITY = RUN_ALL || RUN_QUALITY
const SHOULD_RUN_RELAXATIONS = RUN_ALL || RUN_RELAXATIONS
const SHOULD_RUN_PROBLEMS = RUN_ALL || RUN_PROBLEMS

# Print test plan
@info "Test configuration" USE_LOCAL ARGS
@info "Running tests" polynomials=SHOULD_RUN_POLYNOMIALS quality=SHOULD_RUN_QUALITY relaxations=SHOULD_RUN_RELAXATIONS problems=SHOULD_RUN_PROBLEMS

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
        include("setup.jl")
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
end

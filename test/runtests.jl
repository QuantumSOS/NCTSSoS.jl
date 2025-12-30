# =============================================================================
# NCTSSoS.jl Test Suite
# =============================================================================
#
# Test Structure:
# ---------------
# test/
# ├── polynomials/     - Core polynomial algebra (types, arithmetic, simplification)
# │   └── runtests.jl
# ├── quality/         - Code quality checks (Aqua, ExplicitImports, Doctest)
# │   └── runtests.jl
# ├── solvers/         - SDP solver integration (moment, SOS, sparsity, GNS)
# │   └── runtests.jl
# ├── physics/         - Physics models (Heisenberg, XY, etc.) [--local only]
# │   └── runtests.jl
# ├── setup.jl         - Shared solver configuration
# └── runtests.jl      - This file (main entry point)
#
# Solver Configuration:
# ---------------------
# - --local  → Mosek (commercial, fast, required for physics tests)
# - Default  → COSMO (open-source, sufficient for basic tests)
#
# Run commands:
# -------------
# CI suite (default):   Pkg.test("NCTSSoS")
# Full suite:           Pkg.test("NCTSSoS"; test_args=["--local"])
# Subset tests:         Pkg.test("NCTSSoS"; test_args=["--polynomials"])
#                       Pkg.test("NCTSSoS"; test_args=["--solvers"])
#                       Pkg.test("NCTSSoS"; test_args=["--physics", "--local"])
#                       Pkg.test("NCTSSoS"; test_args=["--quality"])
# Multiple subsets:     Pkg.test("NCTSSoS"; test_args=["--polynomials", "--solvers"])
#
# Makefile shortcuts:
# -------------------
# make test                  - Full suite with Mosek (--local)
# make test-ci               - CI suite (no physics, uses COSMO)
# make test-polynomials      - Polynomial tests only
# make test-solvers          - Solver tests only
# make test-physics          - Physics tests only (requires Mosek)
# make test-quality          - Code quality checks only
# =============================================================================

using NCTSSoS, Test

# =============================================================================
# Parse test arguments
# =============================================================================
const RUN_POLYNOMIALS = "--polynomials" in ARGS
const RUN_QUALITY = "--quality" in ARGS
const RUN_SOLVERS = "--solvers" in ARGS
const RUN_PHYSICS = "--physics" in ARGS
const USE_LOCAL = "--local" in ARGS

# If no specific test group flags, run all (physics only with --local)
const RUN_ALL = !(RUN_POLYNOMIALS || RUN_QUALITY || RUN_SOLVERS || RUN_PHYSICS)

# Determine which test groups to run
const SHOULD_RUN_POLYNOMIALS = RUN_ALL || RUN_POLYNOMIALS
const SHOULD_RUN_QUALITY = RUN_ALL || RUN_QUALITY
const SHOULD_RUN_SOLVERS = RUN_ALL || RUN_SOLVERS

# Physics requires --local (for Mosek) unless explicitly requested
# If explicitly requested without --local, we'll warn but still try
const SHOULD_RUN_PHYSICS = RUN_PHYSICS || (RUN_ALL && USE_LOCAL)

if RUN_PHYSICS && !USE_LOCAL
    @warn "Running physics tests without --local. These may fail without Mosek."
end

# Print test plan
@info "Test configuration" USE_LOCAL ARGS
@info "Running tests" polynomials=SHOULD_RUN_POLYNOMIALS quality=SHOULD_RUN_QUALITY solvers=SHOULD_RUN_SOLVERS physics=SHOULD_RUN_PHYSICS

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
    # Load solver configuration BEFORE solver/physics tests
    # (after polynomial tests to avoid JuMP's simplify shadowing NCTSSoS.simplify)
    # =========================================================================
    if SHOULD_RUN_SOLVERS || SHOULD_RUN_PHYSICS
        include("setup.jl")
    end

    # =========================================================================
    # 3. Solver Integration Tests
    # =========================================================================
    if SHOULD_RUN_SOLVERS
        @testset "Solvers" begin
            include("solvers/runtests.jl")
        end
    end

    # =========================================================================
    # 4. Physics Model Tests (--local recommended - require Mosek)
    # =========================================================================
    if SHOULD_RUN_PHYSICS
        @testset "Physics" begin
            include("physics/runtests.jl")
        end
    end
end

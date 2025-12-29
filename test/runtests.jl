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
# ├── physics/         - Physics models (Heisenberg, XY, etc.) [LOCAL_TESTING only]
# │   └── runtests.jl
# ├── setup.jl         - Shared solver configuration
# └── runtests.jl      - This file (main entry point)
#
# Solver Configuration:
# ---------------------
# - LOCAL_TESTING=true  → Mosek (commercial, fast, required for large problems)
# - Default             → COSMO (open-source, sufficient for basic tests)
#
# Run commands:
# -------------
# Full suite:    julia --project -e 'using Pkg; Pkg.test()'
# Local suite:   LOCAL_TESTING=true julia --project -e 'using Pkg; Pkg.test()'
# Single file:   julia --project -e 'include("test/solvers/moment.jl")'
# Single folder: julia --project -e 'include("test/solvers/runtests.jl")'
# =============================================================================

using NCTSSoS, Test

@testset "NCTSSoS.jl" begin
    # =========================================================================
    # 1. Polynomial Algebra Tests (no solver needed, no JuMP imports)
    # =========================================================================
    @testset "Polynomials" begin
        include("polynomials/runtests.jl")
    end

    # =========================================================================
    # 2. Code Quality Checks
    # =========================================================================
    @testset "Quality" begin
        include("quality/runtests.jl")
    end

    # =========================================================================
    # Load solver configuration AFTER polynomial tests to avoid JuMP's 
    # simplify function shadowing NCTSSoS.simplify
    # =========================================================================
    include("setup.jl")

    # =========================================================================
    # 3. Solver Integration Tests
    # =========================================================================
    @testset "Solvers" begin
        include("solvers/runtests.jl")
    end

    # =========================================================================
    # 4. Physics Model Tests (LOCAL_TESTING only - require Mosek)
    # =========================================================================
    if LOCAL_TESTING
        @testset "Physics" begin
            include("physics/runtests.jl")
        end
    end
end

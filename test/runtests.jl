# =============================================================================
# NCTSSoS.jl Test Suite
# =============================================================================
#
# Test Categories:
# ----------------
# 1. polynomials/     - Core polynomial algebra (types, arithmetic, simplification)
# 2. quality/         - Code quality checks (Aqua, ExplicitImports)
# 3. solvers/         - SDP solver integration (moment, SOS, sparsity)
# 4. physics/         - Physics models (Heisenberg, XY, etc.) [LOCAL_TESTING only]
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
    # 2. Code Quality Checks (each defines its own testset)
    # =========================================================================
    include("quality/Aqua.jl")
    include("quality/ExplicitImports.jl")

    # =========================================================================
    # Load solver configuration AFTER polynomial tests to avoid JuMP's 
    # simplify function shadowing NCTSSoS.simplify
    # =========================================================================
    include("setup.jl")

    # =========================================================================
    # 3. Sparsity Algorithm Tests (unit tests, no solver needed)
    # =========================================================================
    @testset "Sparsity Algorithms" begin
        include("pop.jl")
        include("sparse.jl")
    end

    # =========================================================================
    # 4. SDP Solver Integration Tests
    # =========================================================================
    @testset "Solver Integration" begin
        include("solvers/moment.jl")
        include("solvers/sos.jl")
        include("solvers/interface.jl")
        include("solvers/sparsity.jl")
        include("solvers/state_poly.jl")
        include("solvers/trace_poly.jl")
    end

    # =========================================================================
    # 5. Physics Model Tests (LOCAL_TESTING only - require Mosek)
    # =========================================================================
    if LOCAL_TESTING
        @testset "Physics Models" begin
            include("physics/heisenberg.jl")
            include("physics/xy_model.jl")
            include("physics/bose_hubbard.jl")
            include("physics/bell_inequalities.jl")
            include("physics/fermionic.jl")
        end
    end
end

# NCTSSoS.jl test entrypoint. Canonical commands/flags: see TESTING.md.

using NCTSSoS, Test

# Parse test arguments
const IS_CI = lowercase(get(ENV, "CI", "")) in ("1", "true", "yes", "on")

const RUN_POLYNOMIALS = "--polynomials" in ARGS
const RUN_QUALITY = "--quality" in ARGS
const RUN_RELAXATIONS = "--relaxations" in ARGS
const RUN_PROBLEMS = "--problems" in ARGS
const RUN_MINIMAL = "--minimal" in ARGS
const USE_LOCAL = "--local" in ARGS

const HAS_ANY_SELECTOR = RUN_POLYNOMIALS || RUN_QUALITY || RUN_RELAXATIONS || RUN_PROBLEMS || RUN_MINIMAL

# If no specific test group flags:
# - Local: run all (some problem tests only with --local)
# - CI: run polynomials + relaxations + minimal for fast feedback
const RUN_ALL = !HAS_ANY_SELECTOR && !IS_CI
const USE_CI_DEFAULT = !HAS_ANY_SELECTOR && IS_CI

# Determine which test groups to run
const SHOULD_RUN_POLYNOMIALS = RUN_ALL || RUN_POLYNOMIALS || USE_CI_DEFAULT
const SHOULD_RUN_QUALITY = RUN_ALL || RUN_QUALITY
const SHOULD_RUN_RELAXATIONS = RUN_ALL || RUN_RELAXATIONS || USE_CI_DEFAULT
const SHOULD_RUN_PROBLEMS = RUN_ALL || RUN_PROBLEMS
const SHOULD_RUN_MINIMAL = RUN_MINIMAL || USE_CI_DEFAULT

# Print test plan
@info "Test configuration" IS_CI USE_LOCAL ARGS
@info "Running tests" polynomials=SHOULD_RUN_POLYNOMIALS quality=SHOULD_RUN_QUALITY relaxations=SHOULD_RUN_RELAXATIONS problems=SHOULD_RUN_PROBLEMS minimal=SHOULD_RUN_MINIMAL

# Run selected test groups
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
    # Coverage:
    # | # | Test                  | Algorithm Path              | Source File          |
    # |---|-----------------------|-----------------------------|----------------------|
    # | 1 | Dense baseline        | No sparsity (baseline)      | chsh_simple.jl       |
    # | 2 | Correlative sparsity  | MF clique decomposition     | chsh_simple.jl       |
    # | 3 | Term sparsity         | MMD block reduction         | nc_example1.jl       |
    # | 4 | Combined CS+TS        | Both + constraints          | nc_example2.jl       |
    # | 5 | Dualization           | SOS ≈ Moment equivalence    | dualization.jl       |
    # =========================================================================
    if SHOULD_RUN_MINIMAL
        include("TestUtils.jl")
        @testset "Minimal Suite" begin
            include("problems/bell_inequalities/chsh_simple.jl")
            include("problems/nc_polynomial/nc_example1.jl")
            include("problems/nc_polynomial/nc_example2.jl")
            include("relaxations/dualization.jl")
        end
    end
end

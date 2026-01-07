# =============================================================================
# NCTSSOS Oracle Utilities
# =============================================================================
# Shared helpers for oracle generation scripts.
#
# ENVIRONMENT SETUP:
#   These scripts must be run from the NCTSSOS repository (not NCTSSoS.jl).
#   Set NCTSSOS_PATH environment variable or use default locations:
#     - Local (macOS): /Users/yushengzhao/projects/NCTSSOS
#     - a800 server:   /home/yushengzhao/NCTSSOS
#
# RUNNING ORACLE SCRIPTS:
#   Option 1: Using make (from NCTSSoS.jl repo)
#     make oracle-chsh                         # default: mosek
#     make oracle-chsh SOLVER=cosmo            # use COSMO solver
#     make oracle-i3322 SOLVER=mosek           # use Mosek solver
#     NCTSSOS_PATH=/custom/path make oracle-chsh
#
#   Option 2: Manual (from NCTSSOS repo)
#     cd /path/to/NCTSSOS
#     julia --project /path/to/NCTSSoS.jl/test/oracles/scripts/nctssos_chsh.jl
#     ORACLE_SOLVER=cosmo julia --project /path/to/NCTSSoS.jl/test/oracles/scripts/nctssos_chsh.jl
#
# ORACLE FORMAT:
#   (opt, sides, nuniq) where:
#     opt   = optimal objective value (minimization)
#     sides = vector of moment matrix block sizes
#     nuniq = unique moment indices (affine constraints count)
#
# USAGE IN SCRIPTS:
#   include("oracle_utils.jl")
#   ... define problem-specific obj, vars, constraints ...
#   for v in VARIANTS
#       result = run_oracle(...)
#   end
# =============================================================================

# Detect NCTSSOS path for informational purposes
const NCTSSOS_PATH = get(ENV, "NCTSSOS_PATH",
    isdir("/Users/yushengzhao/projects/NCTSSOS") ? "/Users/yushengzhao/projects/NCTSSOS" :
    isdir("/home/yushengzhao/NCTSSOS") ? "/home/yushengzhao/NCTSSOS" :
    "(not detected - set NCTSSOS_PATH)"
)

# Solver selection via ORACLE_SOLVER env var (default: mosek)
const ORACLE_SOLVER_NAME = Symbol(get(ENV, "ORACLE_SOLVER", "mosek"))

using NCTSSOS, DynamicPolynomials, JuMP

if ORACLE_SOLVER_NAME == :cosmo
    using COSMO
    const ORACLE_SOLVER = optimizer_with_attributes(
        COSMO.Optimizer,
        "verbose" => false,
        "eps_abs" => 1e-8,
        "eps_rel" => 1e-8
    )
else
    using MosekTools
    const ORACLE_SOLVER = Mosek.Optimizer
end

@info "Oracle utilities loaded" NCTSSOS_PATH ORACLE_SOLVER_NAME

# Extract oracle info from NCTSSOS result
function extract_oracle(name, opt, data; use_cs=false)
    sides = if use_cs
        [size(M, 1) for clique in data.moment for M in clique]
    else
        [size(M, 1) for M in data.moment]
    end
    nuniq = length(data.ksupp)
    (name=name, opt=opt, sides=sides, nuniq=nuniq)
end

# Print oracle result in Dict format
function print_oracle(result)
    println("    \"$(result.name)\" => (opt=$(result.opt), sides=$(result.sides), nuniq=$(result.nuniq)),")
end

# Print summary Dict
function print_summary(prefix, results)
    println()
    println("=" ^ 70)
    println("$(uppercase(prefix)) ORACLE SUMMARY (solver: $ORACLE_SOLVER_NAME)")
    println("=" ^ 70)
    println()
    println("const $(uppercase(prefix))_ORACLES = Dict(")
    for res in sort(results, by=r -> r.name)
        print_oracle(res)
    end
    println(")")
end

# Print header
function print_header(title)
    println("=" ^ 70)
    println("NCTSSOS Oracle: $title (solver: $ORACLE_SOLVER_NAME)")
    println("=" ^ 70)
    println()
end

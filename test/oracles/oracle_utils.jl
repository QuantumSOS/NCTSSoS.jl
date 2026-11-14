# NCTSSOS Oracle Utilities
# Shared helpers for oracle generation scripts.
#
# ENVIRONMENT SETUP:
#   These scripts must be run from the NCTSSOS repository (not NCTSSoS.jl).
#   Set NCTSSOS_PATH to the legacy repository; there is no default location.
#
# RUNNING ORACLE SCRIPTS:
#   Option 1: Using make (from NCTSSoS.jl repo)
#     make oracle-chsh                         # default: mosek
#     NCTSSOS_PATH=/custom/path make oracle-chsh
#
#   Option 2: Manual (from NCTSSOS repo)
#     cd /path/to/NCTSSOS
#     julia --project /path/to/NCTSSoS.jl/test/oracles/nctssos_chsh.jl
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

# Read NCTSSOS path for informational purposes
const NCTSSOS_PATH = get(ENV, "NCTSSOS_PATH", "(not set - set NCTSSOS_PATH)")

# Generators use the legacy NCTSSOS solver default (Mosek).

using NCTSSOS, DynamicPolynomials

@info "Oracle utilities loaded (legacy solver default: Mosek)" NCTSSOS_PATH

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
    println("$(uppercase(prefix)) ORACLE SUMMARY (legacy solver default: Mosek)")
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
    println("NCTSSOS Oracle: $title (legacy solver default: Mosek)")
    println("=" ^ 70)
    println()
end

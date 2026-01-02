# NCTSSOS Oracle Utilities
# ========================
# Shared helpers for oracle generation scripts.
# Run scripts on server with NCTSSOS + MosekTools.
#
# Usage in scripts:
#   include("oracle_utils.jl")
#   include(joinpath(@__DIR__, "..", "problems", "my_problem.jl"))  # get VARIANTS
#   ... define problem-specific obj, vars, constraints ...
#   for v in MY_PROBLEM_VARIANTS
#       result = run_oracle(...)
#   end

using NCTSSOS, DynamicPolynomials
using MosekTools

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
    println("$(uppercase(prefix)) ORACLE SUMMARY")
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
    println("NCTSSOS Oracle: $title")
    println("=" ^ 70)
    println()
end

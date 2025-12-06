# Setup file for FastPolynomials tests
# This loads FastPolynomials directly without going through NCTSSoS
# (NCTSSoS source migration is Phase 3 - these tests are for Phase 2)

using Test

# Load FastPolynomials directly
const FASTPOLY_ROOT = joinpath(@__DIR__, "../../src/FastPolynomials/src")
include(joinpath(FASTPOLY_ROOT, "FastPolynomials.jl"))
using .FastPolynomials

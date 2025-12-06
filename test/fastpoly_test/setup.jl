# Setup file for FastPolynomials tests
# Uses NCTSSoS.FastPolynomials to ensure type compatibility with NCTSSoS tests
# that run afterward in the test suite.

using Test

# Use FastPolynomials from NCTSSoS (loaded by runtests.jl via `using NCTSSoS`)
using NCTSSoS.FastPolynomials

# Import @ncpolyvar explicitly to avoid ambiguity
import NCTSSoS.FastPolynomials: @ncpolyvar

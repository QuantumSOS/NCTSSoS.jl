# =============================================================================
# Solver Tests Runner
# =============================================================================
# Runs all solver-related tests including:
#   - Moment matrix construction
#   - SOS dualization
#   - High-level interface (polyopt, cs_nctssos)
#   - Sparsity algorithms (correlative + term)
#   - State polynomial optimization
#   - Trace polynomial optimization
#   - GNS reconstruction
#   - NC polynomial optimization benchmarks
#
# Prerequisites:
#   - SOLVER must be defined (from test/setup.jl or parent runtests.jl)
#   - LOCAL_TESTING determines whether Mosek or COSMO is used
# =============================================================================

using Test

@testset "Solvers" begin
    include("moment.jl")
    include("sos.jl")
    include("interface.jl")
    include("sparsity.jl")
    include("state_poly.jl")
    include("trace_poly.jl")
    include("gns.jl")
    include("ncpop_benchmarks.jl")
end

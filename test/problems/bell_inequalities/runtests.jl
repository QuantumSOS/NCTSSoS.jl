# =============================================================================
# test/problems/bell_inequalities/runtests.jl
# =============================================================================
# Tests: Bell inequality optimization problems
# Dependencies: SOLVER
# Requires --local: no
#
# Coverage:
# - chsh_simple.jl: CHSH Dense, CS, TS (order=1)
# - chsh_high_order.jl: CHSH CS+TS (order=2)
# - chsh_state.jl: CHSH State Polynomial formulation
# - chsh_trace.jl: CHSH Trace Polynomial formulation
# - i3322.jl: I_3322 inequality
# - bell_inequalities.jl: General Bell inequalities
# =============================================================================

using Test

@testset "Bell Inequalities" begin
    include("chsh_simple.jl")
    include("chsh_high_order.jl")
    include("chsh_state.jl")
    include("chsh_trace.jl")
    include("i3322.jl")
    include("bell_inequalities.jl")
end

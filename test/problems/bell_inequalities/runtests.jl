# =============================================================================
# Bell Inequality Tests Runner
# =============================================================================
# Tests for Bell inequality optimization problems:
#   - CHSH inequality
#   - I_3322 inequality
#   - General Bell inequalities
# =============================================================================

using Test

@testset "Bell Inequalities" begin
    include("chsh.jl")
    include("i3322.jl")
    include("bell_inequalities.jl")
end

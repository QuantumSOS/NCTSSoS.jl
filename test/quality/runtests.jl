# =============================================================================
# Code Quality Tests Runner
# =============================================================================
# Runs code quality checks:
#   - Aqua.jl: ambiguities, unbound args, piracy, etc.
#   - ExplicitImports.jl: import hygiene
#   - Doctest: docstring example verification
# =============================================================================

using Test

@testset "Code Quality" begin
    include("Aqua.jl")
    include("ExplicitImports.jl")
    include("Doctest.jl")
end

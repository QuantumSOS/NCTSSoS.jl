# =============================================================================
# test/quality.jl - Code Quality Checks
# =============================================================================
# Run: julia --project -e 'include("test/quality.jl")'
#
# Runs code quality checks:
#   - Aqua.jl: ambiguities, unbound args, piracy, etc.
#   - ExplicitImports.jl: import hygiene
#   - Doctest: docstring example verification
# =============================================================================

using Test

@testset "Code Quality" begin
    include("quality/Aqua.jl")
    include("quality/ExplicitImports.jl")
    include("quality/Doctest.jl")
end

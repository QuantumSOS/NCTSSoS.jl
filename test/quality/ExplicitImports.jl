# ExplicitImports.jl checks for clean imports
# ExplicitImports is a test dependency - only loaded through Pkg.test()
using Test

@testset "ExplicitImports" begin
    using NCTSSoS, ExplicitImports
    check_no_stale_explicit_imports(NCTSSoS)
    check_all_qualified_accesses_via_owners(NCTSSoS)
    check_no_self_qualified_accesses(NCTSSoS)
    check_all_qualified_accesses_are_public(NCTSSoS, ignore=(
        # Julia Base / LinearAlgebra internals used by dependencies
        :Zeros,
        :PositiveSemidefiniteConeSquare,
        :power_by_squaring,
        :show_default,
        :HasLength,
        # MathOptInterface termination codes (not marked public as of MOI v1.48)
        :OPTIMAL,
        :ALMOST_OPTIMAL,
        :LOCALLY_SOLVED,
        :ALMOST_LOCALLY_SOLVED,
    ))
end

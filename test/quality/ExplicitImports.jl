# ExplicitImports.jl checks for clean imports
using Test

@testset "ExplicitImports" begin
    using NCTSSoS, ExplicitImports
    check_no_stale_explicit_imports(NCTSSoS)
    check_all_qualified_accesses_via_owners(NCTSSoS)
    check_no_self_qualified_accesses(NCTSSoS)
    check_all_qualified_accesses_are_public(NCTSSoS, ignore=(
        # MathOptInterface enum values - correct API but not in names(MOI)
        # These are TerminationStatusCode and ResultStatusCode enum instances
        :OPTIMAL,
        :ALMOST_OPTIMAL,
        :LOCALLY_SOLVED,
        :ALMOST_LOCALLY_SOLVED,
        :ITERATION_LIMIT,
        :FEASIBLE_POINT,
        :NEARLY_FEASIBLE_POINT,
    ))
end

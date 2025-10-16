using Test, NCTSSoS, ExplicitImports

@testset "Stale Imports" begin
    check_no_stale_explicit_imports(NCTSSoS)
    check_all_qualified_accesses_via_owners(NCTSSoS)
    check_no_self_qualified_accesses(NCTSSoS, ignore=(:FastPolynomials,))
    check_all_qualified_accesses_are_public(NCTSSoS, ignore=(:Zeros, :PositiveSemidefiniteConeSquare, :power_by_squaring, :show_default))
end
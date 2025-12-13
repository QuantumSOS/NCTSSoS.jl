using Test, Documenter, NCTSSoS

@testset "DocTest" begin
    # Doctests need access to both NCTSSoS and FastPolynomials exports
    DocMeta.setdocmeta!(NCTSSoS, :DocTestSetup,
        :(using NCTSSoS, NCTSSoS.FastPolynomials, Clarabel, COSMO);
        recursive=true)
    doctest(NCTSSoS; manual=false, fix=false)
end
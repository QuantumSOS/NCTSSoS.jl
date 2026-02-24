using Test, Documenter, NCTSSoS
using Documenter: DocMeta

@testset "DocTest" begin
    DocMeta.setdocmeta!(NCTSSoS, :DocTestSetup, :(using NCTSSoS); recursive=true)
    doctest(NCTSSoS; manual=false)
end

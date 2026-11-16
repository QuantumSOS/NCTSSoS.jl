using Test, NCTSSoS

@testset "NCTSSoS public API documentation" begin
    for name in names(NCTSSoS; all=true)
        Base.ispublic(NCTSSoS, name) || continue
        @testset "$name" begin
            @test Docs.hasdoc(NCTSSoS, name)
        end
    end
end

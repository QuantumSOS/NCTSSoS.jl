using Test, NCTSSoS

@testset "Next public API documentation" begin
    for name in names(NCTSSoS.Next; all=true)
        Base.ispublic(NCTSSoS.Next, name) || continue
        @testset "$name" begin
            @test Docs.hasdoc(NCTSSoS.Next, name)
        end
    end
end

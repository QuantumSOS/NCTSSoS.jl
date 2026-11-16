function trace_literature(sos_solver)
    @testset "T3 NCSOStools dense trace polynomial" begin
        _,(a,)=site_variables(Free(),1,2;functional=NCTSSoS.RealTrace())
        x,y=a
        h=trace_symbol(1+2x+x^2+x*y^2+2y^2+y^2*x+y*x^2*y+y^4)
        built=NCTSSoS.build(Problem(h),Relaxation(order=2))
        oracle=fixture("t3_ncsostools_demo_polynomial","dense_order2")
        @test NCTSSoS.block_sides(built)==oracle["sides"]
        @test built.sizes.moments==oracle["nuniq"]
        @test solve(built;solver=sos_solver).bound ≈ oracle["objective"] atol=1e-6 skip=isnothing(sos_solver)
    end
    @testset "T6 projector trace dimensions" begin
        _,(p,)=site_variables(Projector(),1,3;functional=NCTSSoS.RealTrace())
        h=trace_symbol(p[1]*p[2]*p[3])+trace_symbol(p[1]*p[2])*trace_symbol(p[3])
        # The bounds at both orders are already covered by
        # FunctionalModelChecks.check(...;references=true).
        for order in 2:3
            built=NCTSSoS.build(Problem(h),Relaxation(order=order))
            oracle=fixture("t6_projection_optimization","dense_order$order")
            @test NCTSSoS.block_sides(built)==oracle["sides"]
            @test built.sizes.moments==oracle["nuniq"]
        end
    end
end

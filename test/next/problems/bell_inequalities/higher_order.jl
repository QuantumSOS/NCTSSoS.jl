function higher_bell(sos_solver)
    @testset "CHSH combined sparsity retains its loose bound" begin
        _,(a,b)=site_variables(Unipotent(),2,2)
        h=-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
        built=NCTSSoS.build(Problem(h),Relaxation(order=2,
            grouping=Correlative(algorithm=MF()),split=TermSplit(algorithm=MMD())))
        # The old test pins the count, not individual sides.
        @test length(NCTSSoS.block_sides(built))==12
        @test built.sizes.moments==5
        @test solve(built;solver=sos_solver).bound ≈ -3.999999999803662 atol=1e-6 skip=isnothing(sos_solver)
    end
    @testset "I3322 dense and term-sparse bounds" begin
        _,(a,b)=site_variables(Projector(),2,3)
        h=-(a[1]*(b[1]+b[2]+b[3])+a[2]*(b[1]+b[2]-b[3])+
            a[3]*(b[1]-b[2])-a[1]-2b[1]-b[2])
        # Reviewed Mosek targets and block multiplicities from the old i3322.jl.
        for (relaxation,target,sides,moments) in (
            (Relaxation(order=2),-0.25093979763955865,[28],154),
            (Relaxation(order=3),-0.25087575215587177,[88],868),
            (Relaxation(order=3,split=TermSplit(algorithm=MMD())),-0.25147090316786014,
                [fill(2,12);fill(4,36);fill(6,12);fill(9,12);fill(12,3);13],121))
            built=NCTSSoS.build(Problem(h),relaxation)
            @test sort(NCTSSoS.block_sides(built))==sides
            @test built.sizes.moments==moments
            @test solve(built;solver=sos_solver).bound ≈ target atol=1e-6 skip=isnothing(sos_solver)
        end
    end
end

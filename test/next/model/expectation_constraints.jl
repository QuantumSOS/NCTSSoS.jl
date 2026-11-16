module ExpectationConstraintChecks
using Test, LinearAlgebra, NCTSSoS
import NCTSSoS as N

function check(solver)
    @testset "Scalar expectation constraints" begin
        _,(letters,)=N.site_variables(N.Unipotent(),1,2)
        x,y=letters; objective=y*x*y; relaxation=N.Relaxation(order=3)
        unconstrained=N.build(N.Problem(objective),relaxation)
        problem=N.Problem(objective;constraints=(N.ExpectationPositive(x),))
        built=N.build(problem,relaxation)
        @test built.sizes.moments==unconstrained.sizes.moments
        @test built.sizes.rows==unconstrained.sizes.rows
        @test isempty(built.sizes.localizing_rows)
        @test built.sizes.blocks==unconstrained.sizes.blocks+1
        scalar=N.solve(built;solver)
        operator=N.solve(N.Problem(objective;constraints=(N.Positive(x),)),relaxation;solver)
        # Independent matrices: x=Z, y=X, state |0> has <x>=1 and <yxy>=-1.
        Z=diagm([1.,-1.]); X=[0. 1.;1. 0.]; state=[1.,0.]
        @test dot(state,Z*state)==1 && dot(state,X*Z*X*state)==-1
        @test scalar.bound ≈ -1 atol=1e-6
        @test operator.bound ≈ 1 atol=1e-6
        @test N.expectation(scalar,x)>=-1e-6
        @test occursin("scalar expectation constraints",sprint(show,problem))
        @test occursin("≥ 0 (scalar)",sprint(show,scalar.valid_for))
        @test occursin("⪰ 0 (operator)",sprint(show,operator.valid_for))

        for constructor in (N.ExpectationPositive,N.ExpectationZero)
            @test_throws r"not Hermitian.*self-adjoint" constructor(x*y)
            for F in (N.RealState(),N.RealTrace())
                _,(a,)=N.site_variables(N.Free(),1,1;functional=F)
                @test_throws r"Positive\(state_symbol\(p\)\).*outer averaged" constructor(a[1])
            end
        end
        _,(other,)=N.site_variables(N.Unipotent(),1,2)
        @test_throws r"different registry" N.Problem(x;constraints=(N.ExpectationZero(other[1]),))
        @test_throws r"uncovered word s1_2.*add rows" N.build(
            N.Problem(x;constraints=(N.ExpectationZero(y),)),N.Relaxation(rows=[one(x),x]))
        @test_throws r"uncovered word.*increase the order" N.build(
            N.Problem(x;constraints=(N.ExpectationZero(x*y*x),)),N.Relaxation(order=1))

        # Custom rows and a nonzero RHS: a single scalar equality, no localizer.
        equal=N.solve(N.Problem(x;constraints=(N.ExpectationZero(x-.25),)),
            N.Relaxation(rows=[one(x),x]);solver)
        @test equal.bound ≈ .25 atol=1e-6
        @test isempty(equal.sizes.localizing_rows)
        @test occursin("= 0 (scalar)",sprint(show,N.ExpectationZero(x)))

        _,(a,ad)=N.bosonic_variables(1)
        number=ad[1]*a[1]
        boson=N.solve(N.Problem(-number;constraints=(N.ExpectationZero(number-.25),)),N.Relaxation(order=1);solver)
        @test boson.bound ≈ -.25 atol=1e-6

        _,(px,_,pz)=N.pauli_variables(2)
        swap=N.SitePermutation([2,1])
        @test_throws r"does not preserve constraint 1" N.build(N.Problem(-px[1]-px[2];
            constraints=(N.ExpectationZero(pz[1]),),symmetry=(swap,)),N.Relaxation(order=1))
        invariant=N.Problem(-px[1]-px[2];constraints=(N.ExpectationZero(pz[1]+pz[2]),))
        for r in (N.Relaxation(order=1,grouping=N.Correlative()),
                  N.Relaxation(order=1,split=N.TermSplit()),
                  N.Relaxation(family=[[one(px[1]),px[1],px[2],pz[1],pz[2]]]))
            built=N.build(invariant,r)
            @test isempty(built.sizes.localizing_rows)
            @test last(built.program.blocks).kind==N.ZERO
        end
        @test_throws r"uncovered word Z1.*add rows.*increase the order" N.build(invariant,
            N.Relaxation(family=[[one(px[1]),px[1],px[2]]]))

        _,(px,_,pz)=N.pauli_variables(1)
        for constructor in (N.ExpectationPositive,N.ExpectationZero)
            result=N.solve(N.Problem(-px[1];constraints=(constructor(pz[1]),)),N.Relaxation(order=1);solver)
            @test result.bound ≈ -1 atol=1e-6
            g=N.gns(result;atol=1e-5,check_atol=1e-4)
            @test g.verdict in (:moments_reproduced,:optimizer_certified)
            @test N.expectation(g,pz[1]) ≈ N.expectation(result,pz[1]) atol=1e-4
            @test constructor===N.ExpectationPositive ? N.expectation(g,pz[1])>=-1e-4 : abs(N.expectation(g,pz[1]))<=1e-4
        end
        # Appending a scalar condition must not change the original cached build.
        original=N.build(N.Problem(-pz[1]),N.Relaxation(order=1))
        initial=N.solve(original;solver)
        updated=N.solve(initial,N.ExpectationZero(pz[1]-.25))
        @test updated.bound ≈ -.25 atol=1e-6
        @test N.solve(original;solver).bound ≈ -1 atol=1e-6

        # This scalar equality needs imaginary moments that the original real
        # relaxation did not have. Appending must rebuild that representation.
        _,(xy,)=N.site_variables(N.Free(),1,2)
        x,y=xy
        base=N.Problem(zero(x);constraints=(N.Positive(1-x^2),N.Positive(1-y^2)))
        initial=N.solve(base,N.Relaxation(order=1);solver)
        commutator=im*(x*y-y*x)
        complex=N.solve(initial,N.ExpectationZero(commutator-.2))
        @test complex.bound ≈ 0 atol=1e-6
        @test N.expectation(complex,commutator) ≈ .2 atol=1e-6
        symmetrized=N.solve(complex,N.ExpectationZero(x*y+y*x))
        @test symmetrized.bound ≈ 0 atol=1e-6
        @test N.expectation(symmetrized,x*y+y*x) ≈ 0 atol=1e-6
    end
    check_blocks(solver)
end

function check_blocks(solver)
    @testset "Scalar constraints on closed block families" begin
        reg,(x,y,z)=N.pauli_variables(2)
        total=sum(z); h=-total
        rows=[one(h);x;y;z]
        full=N.Relaxation(family=[rows])
        swap=N.SitePermutation([2,1])
        rotation=N.SpinRotation(reg,:z)
        for (constraint,kind) in ((N.ExpectationZero(total-.5),N.ZERO),
                                  (N.ExpectationPositive(.5-total),N.NONNEGATIVE))
            plain=N.solve(N.Problem(h;constraints=(constraint,)),N.Relaxation(order=1);solver)
            @test plain.bound≈-.5 atol=1e-6
            for symmetry in ((),(swap,),(rotation,),(swap,rotation))
                built=N.build(N.Problem(h;constraints=(constraint,),symmetry),full)
                @test built.sizes.rows==length(rows)
                @test isempty(built.sizes.localizing_rows)
                @test count(b->b.kind==kind,built.program.blocks)==1
                @test last(built.program.blocks).side==1
                @test N.solve(built;solver).bound≈plain.bound atol=1e-6
            end
            # Different program: separate site groups, with one global scalar sum.
            sparse=N.build(N.Problem(h;constraints=(constraint,)),
                N.Relaxation(order=1,grouping=N.Correlative()))
            @test length(sparse.sizes.closure.family)==2
            @test sparse.sizes.largest_block==4
            @test N.solve(sparse;solver).bound<=plain.bound+1e-6
        end
        # Each scalar restriction must be invariant, not just the constraint set.
        @test_throws r"does not preserve constraint 1" N.build(N.Problem(h;
            constraints=(N.ExpectationZero(z[1]),N.ExpectationZero(z[2])),symmetry=(swap,)),full)
        @test_throws r"constraint 1" N.build(N.Problem(h;
            constraints=(N.ExpectationZero(x[1]),),symmetry=(rotation,)),full)
        # Existing rows must cover every scalar term, even without a localizer.
        @test_throws r"uncovered word Z1 Z2.*add rows.*increase the order" N.build(N.Problem(-z[1];
            constraints=(N.ExpectationZero(z[1]*z[2]-.25),)),
            N.Relaxation(family=[[one(h),z[1]]]))
        for level in (1,2)
            p=N.Problem(-z[1];constraints=(N.ExpectationZero(x[1]*x[2]-.25),))
            built=N.build(p,N.Relaxation(order=1,split=N.TermSplit(;level)))
            @test N.solve(built;solver).bound<=N.solve(p,N.Relaxation(order=1);solver).bound+1e-6
        end
        # Sector elimination precedes insertion; the scalar row uses its final map.
        sector=N.Sector(z[1]*z[2],1)
        sector_bounds=Float64[]
        for symmetry in ((),(swap,rotation))
            p=N.Problem(h;constraints=(N.ExpectationZero(total-.5),),symmetry)
            built=N.build(p,N.Relaxation(order=1,strengthen=sector))
            @test last(built.program.blocks).kind==N.ZERO
            @test last(built.program.blocks).side==1
            push!(sector_bounds,N.solve(built;solver).bound)
        end
        @test all(b->isapprox(b,-.5;atol=1e-6),sector_bounds)
        # Appending to a solved block problem rebuilds without mutating the input.
        initial=N.solve(N.Problem(h;symmetry=(swap,rotation)),full;solver)
        updated=N.solve(initial,N.ExpectationZero(total-.25))
        @test updated.bound≈-.25 atol=1e-6
        @test initial.bound≈-2 atol=1e-6

        reg,(a,ad)=N.fermionic_variables(2)
        n=ad.*a; number=sum(n); pair=n[1]*n[2]
        constraint=N.ExpectationZero(pair-.25)
        @test_throws r"uncovered word.*add rows.*increase the order" N.build(
            N.Problem(number;constraints=(constraint,)),N.BlockRelaxation(order=1))
        # P/Q/G supplies the quartic moment absent from the order-one family.
        for symmetry in ((),(N.NumberRotation(reg),))
            built=N.build(N.Problem(number;constraints=(constraint,),symmetry),
                N.Relaxation(order=1,strengthen=N.FermionicConditions(1:2)))
            @test isempty(built.sizes.localizing_rows)
            @test last(built.program.blocks).kind==N.ZERO
            @test last(built.program.blocks).side==1
            result=N.solve(built;solver)
            # A mixture with probability 1/4 in |11> and 3/4 in |00> attains 1/2.
            @test result.bound≈.5 atol=1e-6
            @test N.expectation(result,pair)≈.25 atol=1e-6
        end
    end
end
end

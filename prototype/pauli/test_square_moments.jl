# Standalone with NCTSSoS on LOAD_PATH; no optimizer needed.
using Test
using Random: MersenneTwister, randperm
import NCTSSoS as NC
include("PauliWorkflowRedesign.jl")
const Workflow = PauliWorkflowRedesign

function square_oracle_minimum(W,L,letters)
    function image(swap,sx,sy,dx,dy)
        w = one(W)
        for (site,axis) in letters
            y,x = divrem(site-1,L)
            u,v = swap ? (y,x) : (x,y)
            mapped = 1+mod(sx*u+dx,L)+L*mod(sy*v+dy,L)
            w = last(NC.mul(w,NC.pauli_letter(W,mapped,axis)))
        end
        w
    end
    minimum(image(swap,sx,sy,dx,dy) for swap in (false,true) for sx in (-1,1) for sy in (-1,1)
        for dx in 0:L-1 for dy in 0:L-1)
end

function square_form_delta(a,at,b,bt)
    labels(t) = [(l.word,l.part) for l in t.labels]
    ak,bk = labels(at),labels(bt)
    @test Set(ak)==Set(bk)
    permutation = indexin(ak[2:end],bk[2:end])
    @test all(!isnothing,permutation)
    x,y = NC.moment_standard_form(a),NC.moment_standard_form(b)
    @test x.cones==y.cones
    ix = Int[i for i in permutation]
    @test x.b ≈ y.b atol=1e-12 rtol=0
    @test x.offset ≈ y.offset atol=1e-12 rtol=0
    @test x.c ≈ y.c[ix] atol=1e-12 rtol=0
    maximum(abs,x.A-y.A[:,ix];init=0.0)
end

@testset "Compact square orbit moments" begin
    @testset "Exact canonical order, including limb boundaries" begin
        rng = MersenneTwister(2014)
        for L in (2,6,8,9,16)
            W = NC.word_type(NC.Pauli(),L^2,0)
            t = Workflow.SquareOrbitMoments(W,L,true)
            for j in 1:24
                sites = randperm(rng,L^2)[1:min(L^2,1+mod(j,12))]
                axes = j%3==0 ? fill(2,length(sites)) : rand(rng,0:2,length(sites))
                letters = collect(zip(sites,axes))
                word = foldl(letters;init=one(W)) do w,(i,a)
                    last(NC.mul(w,NC.pauli_letter(W,i,a)))
                end
                @test Workflow._square_canonical(t,word)==square_oracle_minimum(W,L,letters)
            end
        end
    end

    @testset "Nonmutating alias, zero, missing and outside-domain lookups" begin
        W = NC.PauliWord{1}
        t = Workflow.SquareOrbitMoments(W,3,true)
        x(i) = NC.pauli_letter(W,i,0)
        y(i) = NC.pauli_letter(W,i,1)
        z(i) = NC.pauli_letter(W,i,2)
        alias = NC.orbit_coordinates!(t,x(8)).re
        @test NC.moment_terms(t,x(3))==((alias,1),)
        @test NC.moment_terms(t,z(8))===nothing
        @test NC.moment_terms(t,y(64))===nothing
        @test_throws ArgumentError NC.orbit_coordinates!(t,y(64))
        @test Workflow._pricing_native_table(t)
        @test Workflow._pricing_existing_terms(t,x(3),true)==((alias,1),)
        @test Workflow._pricing_existing_terms(t,x(3),false)===nothing
        @test Workflow._pricing_existing_terms(t,y(4),false)==()
        @test NC.moment_terms(t,y(4))==()
        @test isempty(t.zero_words)
        @test length(t.labels)==2
        @test !NC.moment_covered(t,alias)
        NC.orbit_coordinates!(t,y(4))
        NC.orbit_coordinates!(t,y(9))
        @test length(t.zero_words)==1
        _,(vx,_,_) = NC.pauli_variables(9)
        @test Set(NC._ward_seed_words(t,NC.spin_rotations(vx[1].registry)))==Set([one(W),x(1),y(1)])
        @test_throws ArgumentError NC._ward_seed_words(t,NC.spin_rotations(vx[1].registry;sites=1:3))

        unitary = Workflow.SquareOrbitMoments(W,3,false)
        id = NC.orbit_coordinates!(unitary,y(4)).re
        @test NC.moment_terms(unitary,y(9))==((id,1),)
    end

    @testset "Equivalent emitted forms before and after Ward elimination" begin
        L = 3
        _,xyz = NC.pauli_variables(L^2)
        site(x,y) = mod(x,L)+1+L*mod(y,L)
        edges = unique(minmax(site(x,y),site(x+dx,y+dy)) for x in 0:L-1 for y in 0:L-1 for (dx,dy) in ((1,0),(0,1)))
        H = sum(a[i]*a[j] for a in xyz for (i,j) in edges)
        word(p) = last(only(NC.terms(p)))
        rows = word.([one(H);reduce(vcat,xyz);[a[i]*b[j] for (i,j) in edges for a in xyz for b in xyz]])
        odd = word(xyz[2][1]*prod(xyz[3][i] for i in 2:5))
        for spin in (false,true)
            a,at,_ = Workflow.workflow_square_symmetry_program(H,rows,L;spin,coordinate_words=(odd,),compact_moments=false)
            b,bt,_ = Workflow.workflow_square_symmetry_program(H,rows,L;spin,coordinate_words=(odd,))
            @test square_form_delta(a,at,b,bt)<=1e-12
            @test isempty(NC.moment_terms(bt,odd))
            @test (spin ? bt.base : bt) isa Workflow.SquareOrbitMoments
        end
    end
end

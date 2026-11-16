isdefined(@__MODULE__, :SOLVER) || include("../../TestUtils.jl")

module PolynomialPQGChecks
using Test, LinearAlgebra
using NCTSSoS
import ..SOLVER

function model()
    r,(a,ad)=fermionic_variables(4)
    pairs=[(i,j) for i in 1:4 for j in i+1:4]
    rows=[one(a[1]);[a[i]*a[j] for (i,j) in pairs];
        [ad[i]*ad[j] for (i,j) in pairs];[ad[i]*a[j] for i in 1:4 for j in 1:4]]
    # As in the old monomials extraction, discard the CAR reordering phase.
    rows=[p/first(only(NCTSSoS.terms(p))) for p in rows]
    h=sum((-1+0.4k)/2*ad[2k+s]*a[2k+s] for k in 0:1 for s in 1:2)
    eri=Dict((k,l,k,l)=>ComplexF64(0.3+0.1k+0.2l) for k in 0:1 for l in 0:1)
    eri[(0,0,1,1)]=0.8im; eri[(1,1,0,0)]=-0.8im
    for ((k,l,p,t),v) in eri,s in 1:2,u in 1:2
        h+=v/8*ad[2k+s]*ad[2l+u]*a[2t+u]*a[2p+s]
    end
    r,h,rows,sum(ad[i]*a[i] for i in 1:4)
end

function check(;solver=JuMPBackend(SOLVER),sos_solver=nothing,references=false)
    @testset "Polynomial PQG custom rows retain fixed number" begin
        r,h,rows,number=model()
        expected=eigmin(Hermitian([-0.925 0.2im;-0.2im -0.45]))
        # Custom rows plus correlative grouping remain deferred after main.
        @test_throws ArgumentError Relaxation(;rows,grouping=Correlative(algorithm=MMD()),
            split=TermSplit(algorithm=MMD()),strengthen=Sector(number,2))
        for split in (nothing,TermSplit(algorithm=MMD()))
            relaxation=Relaxation(;rows,split,strengthen=Sector(number,2))
            for backend in (isnothing(sos_solver) ? (solver,) : (solver,sos_solver))
                plain=solve(Problem(h),relaxation;solver=backend)
                reduced=solve(Problem(h;symmetry=(NumberRotation(r),)),relaxation;solver=backend)
                @test reduced.bound ≈ plain.bound atol=2e-6
                @test reduced.bound ≈ expected atol=2e-6
            end
        end
    end
    @testset "Separate fermion cliques retain charged-input guard" begin
        r,(a,ad)=fermionic_variables(4)
        h=-(ad[1]*a[2]+ad[2]*a[1])-2im*ad[3]*a[4]+2im*ad[4]*a[3]+0.3ad[3]*a[3]
        relaxation=Relaxation(order=1,grouping=Correlative(algorithm=MMD()),split=TermSplit(algorithm=MMD()))
        @test_throws ArgumentError NCTSSoS.build(Problem(h+a[1]+ad[1];symmetry=(NumberRotation(r),)),relaxation)
        built=NCTSSoS.build(Problem(h;symmetry=(NumberRotation(r),)),relaxation)
        # The energy and charged-row Sector control are shared with the
        # relaxation batch, Number rotations with sparse complex hopping.
        @test maximum(NCTSSoS.block_sides(built))==2
        supports=Set(Tuple(sort!(unique!([s for w in rows for s in NCTSSoS.sites(w)])))
            for rows in NCTSSoS.block_rows(built))
        # An identity-only block is not an additional orbital clique.
        @test filter(!isempty,supports)==Set(((1,2),(3,4)))
    end
end
end

PolynomialPQGChecks.check()

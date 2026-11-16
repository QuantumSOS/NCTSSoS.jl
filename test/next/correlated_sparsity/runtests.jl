isdefined(@__MODULE__, :SOLVER) || include("../../TestUtils.jl")

module PublicCorrelativeChecks
using Test, TOML
using NCTSSoS.Next
import JuMP
import ..SOLVER

function broyden()
    _,(x,)=site_variables(Free(),1,5)
    h=zero(x[1])
    for i in 1:5
        neighbors=setdiff(max(1,i-5):min(5,i+1),i)
        residual=2x[i]+5x[i]^3+1-sum(x[j]+x[j]^2 for j in neighbors)
        h+=residual'*residual
    end
    constraints=Tuple(Positive(g) for g in [1 .- x.^2; x .- 1/3])
    Problem(h;constraints)
end

function check(;solver=JuMPBackend(SOLVER),sos_solver=nothing,references=false)
    @testset "Joint custom Pauli rows" begin
        r,ops=pauli_variables(3)
        h=sum(a[i]*a[i+1] for a in ops for i in 1:2)/4
        # The retired physics selector produced the complete degree-two family
        # on this three-site path. Keep its joint PSD constraint, not pair cones.
        result=solve(Problem(h),Relaxation(rows=collect(Next.basis(r,2)));solver)
        @test result.bound ≈ -1 atol=2e-5
    end
    @testset "Broyden banded polydisc eigenvalue" begin
        data=TOML.parsefile(joinpath(@__DIR__,"../../data/expectations/csts_broyden_banded_polydisc.toml"))
        expected=only(filter(c->c["id"]=="eigenvalue_n5",data["cases"]))["expected"]
        problem=broyden()
        dense=Next.build(problem,Relaxation(order=3))
        sparse=Next.build(problem,Relaxation(order=3,grouping=Correlative(algorithm=MF()),
            split=TermSplit(algorithm=MMD())))
        @test maximum(Next.block_sides(dense))==expected["dense_mb"]
        @test maximum(Next.block_sides(sparse))==expected["sparse_mb"]
        record=only(TOML.parsefile(joinpath(@__DIR__,
            "../../data/expectations/next_correlative_legacy.toml"))["cases"])["expected"]
        @test dense.program.nvars==record["dense_variables"]
        @test sparse.program.nvars==record["sparse_variables"]
        @test maximum(Next.block_sides(dense))==record["dense_max_side"]
        @test maximum(Next.block_sides(sparse))==record["sparse_max_side"]
        if references
            # Preserve the old case's COSMO iteration limit.
            backend=JuMPBackend(JuMP.optimizer_with_attributes(SOLVER.optimizer_constructor,
                SOLVER.params...,"max_iter"=>100_000))
            result=solve(sparse;solver=backend)
            @test result.bound ≈ expected["sparse_opt"] atol=1e-3
        end
    end
end
end

PublicCorrelativeChecks.check()

# Standalone native-writer checks. Construction never requests a license.
using Test, NCTSSoS, Mosek, JuMP, SparseArrays, LinearAlgebra
const N=NCTSSoS.Next

function mosek_forms()
    hermitian=N.MomentStandardForm(sparse([2,3],[1,2],[-sqrt(2),-sqrt(2)],4,2),
        [1.,0,0,1],[1.,2],0.3,[N.DataCone(:hpsd,2,1,4,0)])
    mixed=N.MomentStandardForm(sparse([1,4,5,8],[2,2,1,3],[-1.,-1,-1,1],8,3),
        [0.,sqrt(2),1,-2,0,3,4,0.25],[1.,1,2],0.,
        [N.DataCone(:psd,2,1,3,0),N.DataCone(:nonnegative,1,4,1,0),
         N.DataCone(:soc,3,5,3,0),N.DataCone(:zero,1,8,1,0)])
    linear=N.MomentStandardForm(sparse([1,2],[1,2],[-1.,1.],2,2),[-2.,0.25],
        [1.,2.],0.3,[N.DataCone(:nonnegative,1,1,1,0),N.DataCone(:zero,1,2,1,0)])
    constant=N.MomentStandardForm(spzeros(1,0),[1.],Float64[],2.5,[N.DataCone(:nonnegative,1,1,1,0)])
    return ((hermitian,0.3-sqrt(5),[-1,-2]/sqrt(5)),(mixed,7.5,[5,2,0.25]),
        (linear,2.8,[2.,0.25]),(constant,2.5,Float64[]))
end

function mosek_models()
    _,(x,y,z)=N.pauli_variables(4)
    h=sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3)
    _,(a,ad)=N.fermionic_variables(4)
    hub=-(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    _,(a,b)=N.site_variables(N.Unipotent(),2,2)
    chsh=-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
    return ((N.Problem(h),N.Relaxation(order=1),-9.),
        (N.Problem(hub),N.Relaxation(order=2),1-sqrt(5)),
        (N.Problem(chsh),N.Relaxation(order=1),-2sqrt(2)))
end

@testset "Mosek construction without a license" begin
    @test Base.get_extension(NCTSSoS,:NCTSSoSMosekExt)!==nothing
    @test_throws ArgumentError N.MosekBackend(form=:invalid)
    for (f,_,_) in mosek_forms(),form in (:auto,:moment,:sos)
        handle=N.prepare_standard(f,N.MosekBackend(;form))
        @test handle.form===f
        @test Mosek.getnumvar(handle.task)>=0
        Mosek.deletetask(handle.task)
    end
    @test_throws r"MSK_IPAR_" N.prepare_standard(first(mosek_forms())[1],N.MosekBackend(bad_option=1))
end

licensed=isfile(expanduser("~/mosek/mosek.lic")) || haskey(ENV,"MOSEKLM_LICENSE_FILE")
if !licensed
    @info "Skipping licensed Mosek solves: no ~/mosek/mosek.lic or MOSEKLM_LICENSE_FILE found. Construction tests ran."
else
    import MosekTools
    factory=optimizer_with_attributes(MosekTools.Optimizer,"MSK_IPAR_NUM_THREADS"=>1,
        "QUIET"=>true,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-10,"MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>1e-10,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>1e-10)
    @testset "Mosek native moment and SOS solutions" begin
        for form in (:moment,:sos,:auto)
            backend=N.MosekBackend(;form,MSK_IPAR_NUM_THREADS=1,MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-10,
                MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-10,MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-10)
            for (f,bound,values) in mosek_forms()
                ref=N.solve_standard(f,N.JuMPBackend(factory))
                result=N.solve_standard(f,backend)
                @test result.status==JuMP.MOI.OPTIMAL
                @test result.bound≈bound atol=1e-7
                @test result.bound≈ref.bound atol=1e-7
                @test result.values≈values atol=1e-7
            end
            for (p,r,expected) in mosek_models()
                result=N.solve(p,r;solver=backend)
                ref=N.solve(p,r;solver=N.JuMPBackend(factory))
                @test result.status==ref.status==JuMP.MOI.OPTIMAL
                @test result.bound≈ref.bound atol=1e-7
                @test result.bound≈expected atol=1e-7
                @test N.expectation(result,p.objective)≈expected atol=1e-7
            end
            limited=N.solve_standard(first(mosek_forms())[1],N.MosekBackend(;form,MSK_IPAR_INTPNT_MAX_ITERATIONS=0))
            @test limited.status==JuMP.MOI.ITERATION_LIMIT
            @test limited.bound===nothing
            timed=N.solve_standard(first(mosek_forms())[1],N.MosekBackend(;form,MSK_DPAR_OPTIMIZER_MAX_TIME=0.))
            @test timed.status==JuMP.MOI.TIME_LIMIT
            @test timed.bound===nothing
            _,ops=N.site_variables(N.Free(),1,1)
            x=only(ops[1])
            p=N.Problem(x;constraints=(N.Positive(1-x^2),))
            result=N.solve(p,N.Relaxation(order=2);solver=backend)
            @test result.bound≈-1 atol=1e-7
            @test N.expectation(result,x)≈-1 atol=1e-7
            cut=N.solve(result,N.Zero(x-0.6))
            @test cut.bound≈0.6 atol=1e-7
            @test N.expectation(cut,x)≈0.6 atol=1e-7
            @test !cut.solution.diagnostics.warm_started
            failed=N.solve(cut,N.Positive(x-0.8))
            @test failed.status==JuMP.MOI.INFEASIBLE
            @test failed.bound===nothing
            @test_throws ArgumentError N.expectation(failed,x)
        end
    end
    isdefined(@__MODULE__,:check_gns_backend) || include("../model/gns.jl")
    check_gns_backend(N.MosekBackend())
    isdefined(@__MODULE__,:FunctionalModelChecks) || include("../model/functional.jl")
    FunctionalModelChecks.check(N.MosekBackend())
    isdefined(@__MODULE__,:StagedSolveChecks) || include("../model/staged.jl")
    StagedSolveChecks.check(N.MosekBackend(MSK_IPAR_NUM_THREADS=1))
    @testset "One build, two native forms" begin
        _,(X,_,Z)=N.pauli_variables(1)
        problem=N.Problem(-X[1]-.3Z[1]); relaxation=N.Relaxation(order=1)
        built=N.build(problem,relaxation)
        for form in (:moment,:sos)
            solver=N.MosekBackend(;form,MSK_IPAR_NUM_THREADS=1)
            prepared=N.prepare(built;solver)
            result=N.solve(built;solver)
            direct=N.solve(problem,relaxation;solver)
            @test result.bound==direct.bound
            @test result.bound ≈ -sqrt(1.09) atol=1e-7
            @test result.solution.diagnostics.representation==form
            @test N.prepare(built;solver)===prepared
            @test StagedSolveChecks.preparation_bytes(built,solver)==0
        end
    end
end

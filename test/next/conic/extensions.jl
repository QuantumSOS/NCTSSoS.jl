# Standalone entry point: see the solver-extensions CI job for its pinned environment.
using Test, NCTSSoS, JuMP, LinearAlgebra, SparseArrays, SCS, Mosek, BPSDP
const N = NCTSSoS.Next
const LICENSED = get(ENV,"NCTSSOS_TEST_MOSEK","false")=="true"

@testset "Native extension loading and complex coordinates" begin
    for name in (:NCTSSoSSCSExt,:NCTSSoSMosekExt,:NCTSSoSBPSDPExt)
        @test Base.get_extension(NCTSSoS,name)!==nothing
    end
    M = ComplexF64[3 1+2im;1-2im 4]
    packed = [3,sqrt(2),-2sqrt(2),4]
    f = N.MomentStandardForm(sparse(reshape([1.,2,3,4],4,1)),packed,[0.7],0.3,[N.DataCone(:hpsd,2,1,4,0)])
    realified = N.real_psd_form(f)
    R = [real.(M) -imag.(M);imag.(M) real.(M)]
    svec(A) = [A[i,j]*(i==j ? 1 : sqrt(2)) for j in axes(A,2) for i in j:size(A,1)]
    @test realified.b ≈ svec(R)
    D = ComplexF64[1 (2-3im)/sqrt(2);(2+3im)/sqrt(2) 4]
    @test vec(Matrix(realified.A)) ≈ svec([real.(D) -imag.(D);imag.(D) real.(D)])
    mh = N.prepare_standard(f,N.MosekBackend())
    @test Mosek.getnumvar(mh.task)==1
    @test Mosek.getnumacc(mh.task)==1
    @test Mosek.getnumafe(mh.task)==10
    bh = N.prepare_standard(f,N.BPSDPBackend(print_level=0))
    X = vec(M)
    out = zeros(1)
    BPSDP.mul_A!(out,bh.solver.problem.operator,X)
    @test out ≈ transpose(f.A)*packed
    @test real(dot(bh.solver.problem.c,X)) ≈ dot(f.b,packed)
end

function extension_models()
    r,(x,y,z) = N.pauli_variables(4)
    h = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3)
    r,(a,ad) = N.fermionic_variables(4)
    hub = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    r,(x,y) = N.site_variables(N.Unipotent(),2,2)
    chsh = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
    return (("Heisenberg4",N.Problem(h),N.Relaxation(order=1),-9.),
        ("Hubbard2",N.Problem(hub),N.Relaxation(order=2),1-sqrt(5)),
        ("CHSH",N.Problem(chsh),N.Relaxation(order=1),-2sqrt(2)))
end

@testset "Native writers equal JuMP and known bounds" begin
    backends = Any[N.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7),
        N.BPSDPBackend(print_level=0,sdp_objective_convergence_abs=1e-7,sdp_objective_convergence_rel=1e-7,
            sdp_error_convergence_abs=1e-7,sdp_error_convergence_rel=1e-7)]
    LICENSED && push!(backends,N.MosekBackend())
    for (name,p,r,expected) in extension_models()
        ref = N.solve(p,r;solver=N.JuMPBackend(optimizer_with_attributes(SCS.Optimizer,"verbose"=>0,"eps_abs"=>1e-7,"eps_rel"=>1e-7)))
        @test ref.bound ≈ expected atol=2e-5
        for backend in backends
            sol = N.solve(p,r;solver=backend)
            @test sol.status==JuMP.MOI.OPTIMAL
            @test sol.bound ≈ ref.bound atol=2e-5
            f = N.moment_standard_form(N.build(p,r).program)
            values = sol.solution.values
            @test dot(f.c,values)+f.offset ≈ sol.bound atol=2e-5
            coordinates = f.b-f.A*values
            for cone in f.cones
                matrix = zeros(ComplexF64,cone.side,cone.side)
                k = cone.first
                for j in 1:cone.side,i in j:cone.side
                    z = coordinates[k]; k+=1
                    if i!=j
                        z = (z+im*coordinates[k])/sqrt(2); k+=1
                    end
                    matrix[i,j] = z; matrix[j,i] = conj(z)
                end
                @test eigmin(Hermitian(matrix)) >= -2e-5
            end
            @info "Native writer" name backend=typeof(backend) bound=sol.bound reference=ref.bound
        end
    end
    for backend in backends
        r,ops = N.site_variables(N.Free(),1,1)
        x = ops[1][1]
        sol = N.solve(N.Problem(x;constraints=(N.Positive(1-x^2),)),N.Relaxation(order=2);solver=backend)
        @test sol.bound ≈ -1 atol=2e-5
        limited = N.solve(sol,N.Positive(x-0.25))
        @test limited.bound ≈ 0.25 atol=2e-5
        @test limited.solution.handle===sol.solution.handle
        fixed = N.solve(limited,N.Zero(x-0.6))
        @test fixed.bound ≈ 0.6 atol=2e-5
    end
    _,p,r,expected = last(extension_models())
    f = N.moment_standard_form(N.build(p,r).program)
    sos = N.sos_standard_form(f)
    sol = N.solve_standard(sos,backends[2])
    @test sol.bound ≈ expected atol=2e-5
    @test norm(sos.A*sol.values-sos.b,Inf) <= 2e-5
    for backend in (N.SCSBackend(verbose=0,max_iters=1),N.BPSDPBackend(print_level=0,max_iter=0))
        failed = N.solve(p,r;solver=backend)
        @test failed.bound===nothing
        @test failed.bound_kind==:unavailable
    end
end

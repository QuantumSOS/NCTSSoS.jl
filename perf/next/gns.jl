# Run in an environment containing NCTSSoS, SCS and Mosek, with a valid Mosek
# license. Solver accuracy sweep plus warmed reconstruction on identical moments.
using NCTSSoS, SCS, Mosek, LinearAlgebra, Test
const N=NCTSSoS.Next
include(joinpath(@__DIR__,"../../test/next/model/gns.jl"))

fill_bytes(H,rows,sink)=@allocated N._gns_fill!(H,rows,sink)

function timed_samples(call;samples=5)
    call() # compile before measuring
    rows=NamedTuple[]
    for _ in 1:samples
        GC.gc()
        t=@timed call()
        push!(rows,(seconds=t.time,bytes=t.bytes))
    end
    return rows
end

function legacy_gns_comparison(result)
    g=N.gns(result)
    registry,_=create_pauli_variables(1:4)
    TI=eltype(indices(registry))
    W=eltype(g.basis)
    packed_values=Dict(w=>ComplexF64(N._gns_value(result,N._gns_coordinates(result.table,(w,),(1,)))) for w in N._gns_domain(result.table))
    monomap=Dict(NCTSSoS._gns_moment_key(NormalMonomial{PauliAlgebra,TI}(TI.(N.legacy_letters(w))))=>v for (w,v) in packed_values)
    old=gns_reconstruct(monomap,registry,g.full_degree;hankel_deg=g.principal_degree,atol=1e-6)
    @test old.rank==g.rank && old.full_rank==g.full_rank
    @test old.singular_values ≈ g.singular_values atol=1e-6
    old_error=maximum(abs(dot(old.xi,gns_matrix_word(w,old.matrices)*old.xi)-v) for (w,v) in packed_values)
    @test old_error<=1e-6
    println("LEGACY ranks=",(old.rank,old.full_rank)," all_moment_error=",old_error,
        " singular_value_error=",maximum(abs.(old.singular_values-g.singular_values)))
    H=zeros(ComplexF64,length(g.full_basis),length(g.full_basis))
    sink=N.GNSProductValue(packed_values,0.0im)
    N._gns_fill!(H,g.full_basis,sink)
    bytes=@allocated N._gns_fill!(H,g.full_basis,sink)
    @test bytes==0
    @test real(LinearAlgebra.tr(H)) ≈ length(g.full_basis) atol=1e-6
    println("FILL side=",size(H,1)," bytes=",bytes)
    for (name,call) in (("next_total",()->N.gns(result)),
        ("legacy_total",()->gns_reconstruct(monomap,registry,g.full_degree;hankel_deg=g.principal_degree,atol=1e-6)),
        ("next_fill",()->N._gns_fill!(H,g.full_basis,sink)),
        ("legacy_fill",()->NCTSSoS.hankel_matrix(monomap,registry,g.full_degree)))
        println("SAMPLES ",name," ",timed_samples(call))
    end
end

function main()
    BLAS.set_num_threads(1)
    println("ENV Julia=",VERSION," CPU=",Sys.CPU_NAME," Julia_threads=",Threads.nthreads()," BLAS_threads=",BLAS.get_num_threads())
    println("ENV SCS=",pkgversion(SCS)," Mosek=",pkgversion(Mosek))
    Mosek.maketask() do task
        Mosek.appendvars(task,1)
        Mosek.putvarbound(task,1,Mosek.MSK_BK_LO,0.,Inf)
        Mosek.putcj(task,1,1.)
        Mosek.optimize(task)
        @assert Mosek.getsolsta(task,Mosek.MSK_SOL_ITR)==Mosek.MSK_SOL_STA_OPTIMAL
    end
    comparison=nothing
    for (name,p,d) in gns_cases(), (tag,backend) in ((:scs8,N.SCSBackend(verbose=0,eps_abs=1e-8,eps_rel=1e-8,max_iters=200000)),
        (:scs6,N.SCSBackend(verbose=0,eps_abs=1e-6,eps_rel=1e-6,max_iters=200000)),(:mosek,N.MosekBackend()))
        r=N.solve(p,N.Relaxation(order=d);solver=backend)
        g=name===:chsh ? check_gns_chsh(r) : name===:free ? check_gns_free(r) : check_gns_xxz(r)
        println("ROW ",name," ",tag," status=",r.status," bound=",r.bound," verdict=",g.verdict,
            " pair=",(g.full_degree,g.principal_degree)," ranks=",(g.rank,g.full_rank))
        println("SPECTRA principal=",g.singular_values," full=",g.full_singular_values)
        println("FLATNESS ",g.flatness)
        println("CHECKS ",g.checks)
        name===:xxz && tag===:scs8 && (comparison=r)
        flush(stdout)
    end
    legacy_gns_comparison(comparison)
    _,(a,ad)=N.fermionic_variables(2)
    result=N.solve(N.Problem(-(ad[1]*a[2]+ad[2]*a[1])),N.Relaxation(order=2);
        solver=N.SCSBackend(verbose=0,eps_abs=1e-8,eps_rel=1e-8))
    g=N.gns(result)
    @test g.verdict===:moments_reproduced
    W=eltype(g.full_basis)
    values=Dict{W,ComplexF64}()
    @test isnothing(N._gns_required!(values,result,g.full_basis,N.PolynomialBuffer{W,ComplexF64}()))
    H=zeros(ComplexF64,length(g.full_basis),length(g.full_basis))
    sink=N.GNSProductValue(values,0.0im)
    N._gns_fill!(H,g.full_basis,sink)
    fill_bytes(H,g.full_basis,sink)
    bytes=fill_bytes(H,g.full_basis,sink)
    @test bytes==0
    println("FERMI_FILL side=",size(H,1)," bytes=",bytes," verdict=",g.verdict)
end

abspath(PROGRAM_FILE)==(@__FILE__) && main()

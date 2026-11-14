# Fresh-process diagnostic. CodeInstance inference times are stored as Float16
# seconds by Julia 1.12. Their sum is approximate; compile minus inference also
# contains compiler bookkeeping, not only LLVM code generation.
# Run each case in a separate process, three times, with one Julia/BLAS thread:
# julia --project=<SCS environment> --startup-file=no test/next/benchmark_first_call.jl <case>
# Cases: chain, legacy, xxz, split, hubbard, ball, density, soc, chsh.
const inferred = Any[]
function measure(label,f)
    empty!(inferred)
    ccall(:jl_set_newly_inferred,Cvoid,(Any,),inferred)
    t = @timed f()
    ccall(:jl_set_newly_inferred,Cvoid,(Any,),nothing)
    cis = unique(filter(x->x isa Core.CodeInstance,inferred))
    inference = sum(Float64(reinterpret(Float16,ci.time_infer_self)) for ci in cis;init=0.)
    println("PHASE ",label," wall=",t.time," compile=",t.compile_time,
        " inference_approx=",inference," codegen_and_other_approx=",t.compile_time-inference,
        " bytes=",t.bytes," inferred_instances=",length(cis))
    if label=="call"
        sort!(cis;by=ci->reinterpret(Float16,ci.time_infer_self),rev=true)
        for ci in Iterators.take(cis,12)
            println("INFERENCE ",Float64(reinterpret(Float16,ci.time_infer_self))," ",ci.def)
        end
    end
    flush(stdout)
    return t.value
end
measure("load",()->Core.eval(Main,:(using NCTSSoS)))
const N = NCTSSoS.Next
mode = only(ARGS)
mode in ("chain","legacy","xxz","split","hubbard","ball","density","soc","chsh") || error("Unknown case")
if mode ∉ ("chain","legacy")
    measure("solver_load",()->Core.eval(Main,:(using SCS)))
end
if mode in ("chain","legacy")
    include(joinpath(@__DIR__,"pipeline","benchmark_family.jl"))
end
function input(mode)
    if mode=="chain"
        return model_chain_case(16)
    elseif mode=="legacy"
        registry,ops=NCTSSoS.create_pauli_variables(1:16)
        return NCTSSoS.polyopt(NCTSSoS.heisenberg_chain_hamiltonian(ops;periodic=true,coupling=1),registry),ops
    elseif mode=="xxz"
        n=6
        _,(x,y,z)=N.pauli_variables(n)
        h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+0.7z[i]*z[mod1(i+1,n)] for i in 1:n)
        return N.Problem(h),N.BlockRelaxation(order=2,grouping=N.Correlative(),symmetry=N.FamilySymmetry([N.SitePermutation([2:n;1])]))
    elseif mode=="split"
        n=12
        _,(x,y,z)=N.pauli_variables(n)
        h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+0.7z[i]*z[mod1(i+1,n)] for i in 1:n)
        return N.Problem(h),N.BlockRelaxation(order=2,grouping=N.Correlative(),split=N.TermSplit(),
            symmetry=N.FamilySymmetry([N.SitePermutation([2:n;1]),N.SitePermutation([mod1(2-i,n) for i in 1:n])]))
    elseif mode=="hubbard"
        _,(a,ad)=N.fermionic_variables(4)
        h=-(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
        return N.Problem(h),N.Relaxation(order=2)
    elseif mode=="ball"
        _,ops=N.site_variables(N.Free(),1,2)
        a,b=ops[1]
        return N.Problem(a*b+b*a;constraints=(N.Positive(1-a^2-b^2),)),N.Relaxation(order=2)
    elseif mode in ("density","soc")
        registry,(x,y,z)=N.pauli_variables(3)
        h=-(x[1]*y[2]+x[1]*y[3]+x[2]*x[3]+y[2]*y[3]+z[2]*z[3])
        W=N.word_type(registry,4)
        symmetry=N.FamilySymmetry([N.SitePermutation([1,3,2])];moment_only=(N.AntiunitaryPauli(N.pauli_letter(W,1,1)),))
        kind=mode=="density" ? N.DensityMatrix(1:3) : N.InvariantCone(1:3,"YII")
        return N.Problem(h),N.BlockRelaxation(order=1,strengthen=(kind,),symmetry=symmetry)
    else
        _,(x,y)=N.site_variables(N.Unipotent(),2,2)
        return N.Problem(-(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])),N.Relaxation(order=1)
    end
end
args=measure("inputs",()->input(mode))
result=measure("call",()->mode=="legacy" ? NCTSSoS.pauli_translation_invariant_moment_relaxation(args...,4;axis_permutation_symmetry=true) :
    mode=="chain" ? N.build(args...) : N.solve(args...;solver=N.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7)))
if mode ∉ ("chain","legacy")
    @assert result.bound!==nothing
    expected=mode=="chsh" ? -2sqrt(2) : mode=="hubbard" ? 1-sqrt(5) :
        mode=="ball" ? -1.0 : mode=="density" ? -3.0 : mode=="soc" ? -sqrt(15) : nothing
    isnothing(expected) || @assert isapprox(result.bound,expected;atol=2e-6)
    println("RESULT bound=",result.bound," sizes=",result.sizes," residuals=",result.solution.diagnostics)
elseif mode=="chain"
    println("RESULT sizes=",result.sizes)
else
    println("RESULT type=",typeof(result))
end
println("ENV Julia=",VERSION," CPU=",Sys.cpu_info()[1].model," threads=",Threads.nthreads())

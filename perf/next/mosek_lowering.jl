# Bell n=4/order=3 stage profile. Run in an environment with MosekTools and a
# valid license. Julia allocations exclude native solver memory; use /usr/bin/time
# for whole-process peak RSS. This is a diagnostic, not a solve-time benchmark.
include(joinpath(@__DIR__,"probe.jl"))
import Mosek, MosekTools

function attached_model(f)
    model,_=N.to_jump(f)
    factory=JuMP.optimizer_with_attributes(MosekTools.Optimizer,
        "MSK_IPAR_NUM_THREADS"=>1,"MSK_DPAR_OPTIMIZER_MAX_TIME"=>120.,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-8)
    JuMP.set_optimizer(model,factory)
    JuMP.MOI.Utilities.attach_optimizer(JuMP.backend(model))
    return model
end

function profile_lowering()
    BLAS.set_num_threads(1)
    Mosek.maketask() do task
        Mosek.appendvars(task,1)
        Mosek.putvarbound(task,1,Mosek.MSK_BK_LO,0.,Inf)
        Mosek.putcj(task,1,1.)
        Mosek.optimize(task)
        @assert Mosek.getsolsta(task,Mosek.MSK_SOL_ITR)==Mosek.MSK_SOL_STA_OPTIMAL
    end
    f=N.moment_standard_form(input_case("bell","next",4,3,"dense","")().program)
    println("FORM ",size(f.A)," NNZ ",nnz(f.A))
    for (name,call) in (("LOWER",()->N.to_jump(f)),("LOWER_AND_ATTACH",()->attached_model(f)))
        call()
        for _ in 1:3
            GC.gc()
            t=@timed call()
            println(name," seconds=",t.time," bytes=",t.bytes)
        end
        GC.gc(); Profile.clear(); @profile call()
        Profile.print(format=:flat,mincount=5,sortedby=:count)
    end
    model=attached_model(f)
    GC.gc(); Profile.clear()
    println("ENTER_OPTIMIZE peak_rss=",Sys.maxrss()); flush(stdout)
    t=@timed @profile JuMP.optimize!(model)
    println("EXIT_OPTIMIZE status=",JuMP.termination_status(model),
        " seconds=",t.time," julia_bytes=",t.bytes," peak_rss=",Sys.maxrss())
    Profile.print(format=:flat,C=true,mincount=10,sortedby=:count)
end

abspath(PROGRAM_FILE)==(@__FILE__) && profile_lowering()

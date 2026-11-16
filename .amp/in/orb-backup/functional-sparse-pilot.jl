using Test, NCTSSoS, Mosek, JuMP, MosekTools
const N=NCTSSoS.Next
const solver=N.MosekBackend(MSK_IPAR_NUM_THREADS=1)
for trace in (false,true)
    F=trace ? N.RealTrace() : N.RealState()
    symbol=trace ? N.trace_symbol : N.state_symbol
    r,ops=N.site_variables(N.Unipotent(),trace ? 1 : 2,trace ? 4 : 2;functional=F)
    a,b=trace ? (ops[1][1:2],ops[1][3:4]) : ops
    h=-symbol(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
    p=N.Problem(h)
    for (name,rel) in (("plain",N.Relaxation(order=1)),("block",N.BlockRelaxation(order=1)),
        ("CS",N.Relaxation(order=1,grouping=N.Correlative())),
        ("TS",N.Relaxation(order=1,split=N.TermSplit(algorithm=MMD()))),
        ("CSTS",N.Relaxation(order=1,grouping=N.Correlative(),split=N.TermSplit(algorithm=MMD()))))
        built=N.build(p,rel)
        result=N.solve(built;solver)
        @test result.bound <= -2sqrt(2)+1e-6
        name in ("plain","block","TS") && @test result.bound ≈ -2sqrt(2) atol=1e-6
        println((;trace,name,bound=result.bound,sides=[b.side for b in built.program.blocks],kinds=[b.kind for b in built.program.blocks],vars=built.sizes.real_variables))
        flush(stdout)
    end
end
include(joinpath(pkgdir(NCTSSoS),"test/next/model/functional.jl"))
FunctionalModelChecks.check(solver)

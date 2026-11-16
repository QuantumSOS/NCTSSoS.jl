# Fresh process: julia --project=<solver env> backend_first_call.jl <state|pauli> <scs|both>
using NCTSSoS, SCS, LinearAlgebra
length(ARGS)==2 && ARGS[2] in ("scs","both") || error("Use <state|pauli> <scs|both>")
ARGS[2]=="both" && (@eval using Mosek)
BLAS.set_num_threads(1)
function input(domain)
    if domain=="state"
        _,(a,b)=site_variables(Unipotent(),2,3;functional=NCTSSoS.RealState())
        signs=[1 1 1;1 1 -1;1 -1 0]
        h=-sum(signs[i,j]*(state_symbol(a[i]*b[j])-state_symbol(a[i])*state_symbol(b[j]))
            for i in 1:3,j in 1:3 if signs[i,j]!=0)
    elseif domain=="pauli"
        _,(x,y,z)=pauli_variables(4)
        h=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+z[i]*z[mod1(i+1,4)] for i in 1:4)
    else
        error("Use state or pauli")
    end
    Problem(h)
end
p=input(ARGS[1]); relaxation=Relaxation(order=2)
backend=SCSBackend(verbose=0,eps_abs=1e-6,eps_rel=1e-6)
for call in 1:2
    t=@timed Base.invokelatest(solve,p,relaxation;solver=backend)
    @assert isapprox(t.value.bound,ARGS[1]=="state" ? -5 : -8;atol=1e-5)
    println("domain=",ARGS[1]," loaded=",ARGS[2]," call=",call," seconds=",t.time,
        " compile_seconds=",t.compile_time," bytes=",t.bytes," bound=",t.value.bound)
end

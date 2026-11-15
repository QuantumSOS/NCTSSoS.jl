# Measurement-only paths; these calls are deliberately absent from workloads.
t_load=@elapsed begin
    using NCTSSoS, SCS
end
const N=NCTSSoS.Next
const detection_seconds=Ref(0.0)

function heisenberg_ring(n)
    registry,(x,y,z)=N.pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    return registry,h
end
function inputs(::Val{:heisenberg8continuous})
    registry,h=heisenberg_ring(8)
    symmetry=(N.SitePermutation([2:8;1]),N.SitePermutation(collect(8:-1:1)),N.spin_rotations(registry)...)
    return N.Problem(h;symmetry),N.Relaxation(order=2)
end
function inputs(::Val{:fermionic4number})
    n=4
    registry,(a,ad)=N.fermionic_variables(n)
    h=-sum(ad[i]*a[mod1(i+1,n)]+ad[mod1(i+1,n)]*a[i] for i in 1:n)
    number=sum(ad[i]*a[i] for i in 1:n)
    return N.Problem(h;symmetry=(N.NumberRotation(registry),)),
        N.Relaxation(order=1,strengthen=(N.FermionicConditions(1:n;conditions=(:P,:Q,:G)),N.Sector(number,2)))
end
function inputs(::Val{:heisenberg8detect})
    _,h=heisenberg_ring(8)
    detection_seconds[]=@elapsed symmetry=N.detect_symmetry(h)
    return N.Problem(h;symmetry),N.Relaxation(order=1)
end
function inputs(::Val{:xxz6density})
    n=6; _,(x,y,z)=N.pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+0.7z[i]*z[mod1(i+1,n)] for i in 1:n)
    return N.Problem(h;symmetry=(N.SitePermutation([2:n;1]),)),
        N.Relaxation(order=1,strengthen=(N.DensityMatrix(1:3),))
end

case=only(ARGS)
backend=N.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7)
t_inputs=@elapsed args=inputs(Val(Symbol(case)))
t_call=@elapsed result=N.solve(args...;solver=backend)
t_second=@elapsed N.solve(args...;solver=backend)
println("PROBE case=",case," load=",round(t_load;digits=3),
    " inputs=",round(t_inputs-detection_seconds[];digits=3)," detection=",round(detection_seconds[];digits=3),
    " first_call=",round(t_call;digits=3)," second_call=",round(t_second;digits=4),
    " bound=",result.bound," status=",result.status)

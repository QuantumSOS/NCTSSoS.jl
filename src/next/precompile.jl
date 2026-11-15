using PrecompileTools: @compile_workload

# Common plain and translated-family routes. Less common options compile on use.
# Capture runtime-dispatched pipeline callees as well as inferred callers.
@compile_workload let
    n=4
    _,(x,y,z)=pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+0.7z[i]*z[mod1(i+1,n)] for i in 1:n)
    moment_standard_form(build(Problem(h),Relaxation(order=2)).program)
    translation=SitePermutation([2:n;1])
    built=build(Problem(h;symmetry=(translation,)),Relaxation(order=2,grouping=Correlative()))
    moment_standard_form(built.program)
    _,(a,b)=site_variables(Unipotent(),2,2)
    bell=Problem(-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2]))
    moment_standard_form(build(bell,Relaxation(order=1)).program)
end

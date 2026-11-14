using PrecompileTools: @compile_workload

# Small representatives of the one-call routes. No solver is loaded here.
# Capture runtime-dispatched pipeline callees as well as inferred callers.
@compile_workload let
    n=4
    registry,(x,y,z)=pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    problem=Problem(h)
    W=word_type(registry,4)
    translation=SitePermutation([2:n;1])
    point=AbstractAction[SitePermutation(collect(n:-1:1)),
        CliffordMap([pauli_letter(W,i,a) for a in (1,0) for i in 1:n]),
        CliffordMap([pauli_letter(W,i,a) for a in (1,2) for i in 1:n])]
    build(problem,BlockRelaxation(order=2,row_rule=Contiguous(periodic=true),
        symmetry=FamilySymmetry([translation;point];moment_only=(AntiunitaryPauli(one(W)),))))
    built=build(problem,BlockRelaxation(order=2,grouping=Correlative(),symmetry=FamilySymmetry([translation])))
    moment_standard_form(built.program)
    _,(a,b)=site_variables(Unipotent(),2,2)
    bell=Problem(-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2]))
    moment_standard_form(build(bell,Relaxation(order=1)).program)
end

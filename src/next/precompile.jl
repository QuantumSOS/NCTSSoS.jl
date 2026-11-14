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
    built=build(problem,BlockRelaxation(order=2,grouping=Correlative(),split=TermSplit(),
        symmetry=FamilySymmetry([translation,SitePermutation([mod1(2-i,n) for i in 1:n])])))
    moment_standard_form(built.program)
    _,(a,b)=site_variables(Unipotent(),2,2)
    bell=Problem(-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2]))
    moment_standard_form(build(bell,Relaxation(order=1)).program)
end

# Cover construction arithmetic and plain localizing forms as well as builds.
@compile_workload let
    _,(a,ad)=fermionic_variables(4)
    h=-(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    moment_standard_form(build(Problem(h),Relaxation(order=2)).program)
    _,ops=site_variables(Free(),1,2)
    x,y=ops[1]
    p=Problem(x*y+y*x;constraints=(Positive(1-x^2-y^2),))
    moment_standard_form(build(p,Relaxation(order=2)).program)
end

@compile_workload let
    registry,(x,y,z)=pauli_variables(3)
    h=-(x[1]*y[2]+x[1]*y[3]+x[2]*x[3]+y[2]*y[3]+z[2]*z[3])
    W=word_type(registry,4)
    symmetry=FamilySymmetry([SitePermutation([1,3,2])];moment_only=(AntiunitaryPauli(pauli_letter(W,1,1)),))
    for kind in (DensityMatrix(1:3),InvariantCone(1:3,"YII"))
        built=build(Problem(h),BlockRelaxation(order=1,strengthen=(kind,),symmetry=symmetry))
        moment_standard_form(built.program)
    end
end

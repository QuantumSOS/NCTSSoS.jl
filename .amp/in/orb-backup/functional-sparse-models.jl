module SparseFunctionalModels
using NCTSSoS
import NCTSSoS: Next as N

# Read the unchanged reviewed model definition, without executing its testsets.
function reference_function(path,name)
    function find_definition(ex)
        ex isa Expr || return nothing
        if ex.head==:function && ex.args[1] isa Expr && ex.args[1].head==:call && ex.args[1].args[1]==name
            return ex
        end
        for arg in ex.args
            found=find_definition(arg)
            isnothing(found) || return found
        end
        nothing
    end
    definition=find_definition(Meta.parseall(read(joinpath(pkgdir(NCTSSoS),path),String)))
    isnothing(definition) && error("Missing reviewed model $name")
    Core.eval(@__MODULE__,definition)
end
reference_function("test/problems/quantum_networks/bilocal_networks.jl",:create_bilocal_problem)
reference_function("test/problems/trace_polynomial/t2_broyden_tridiagonal_trace.jl",:broyden_tridiagonal_trace_problem)

function convert_problem(old,rule,nsites,nletters,trace)
    r,ops=N.site_variables(rule,nsites,nletters;functional=trace ? N.RealTrace() : N.RealState())
    variables=collect(Iterators.flatten(ops))
    index=Dict(zip(sort!(collect(indices(old.registry))),variables))
    identity=one(first(variables))
    operator(word)=foldl(*,(index[i] for i in word);init=identity)
    symbol=trace ? N.trace_symbol : N.state_symbol
    function convert_word(w)
        result=operator(w.nc_word.word)
        for s in w.sw.state_syms
            result=result*symbol(operator(s.mono))
        end
        result
    end
    polynomial(p)=N.polynomial(r,(c*convert_word(w) for (c,w) in zip(coefficients(p),monomials(p))))
    new=N.Problem(polynomial(old.objective);constraints=(map(p->N.Zero(polynomial(p)),old.eq_constraints)...,
        map(p->N.Positive(polynomial(p)),old.ineq_constraints)...))
    (;old,new,registry=r,convert_word)
end

function chsh(trace)
    r,ops=create_unipotent_variables(trace ? [("v",1:4)] : [("a",1:2),("b",1:2)])
    a,b=trace ? (ops[1][1:2],ops[1][3:4]) : ops
    symbol=trace ? NCTSSoS.tr : ς
    h=-symbol(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
    old=polyopt(h*one(typeof(a[1])),r)
    convert_problem(old,N.Unipotent(),trace ? 1 : 2,trace ? 4 : 2,trace)
end
bilocal()=convert_problem(first(create_bilocal_problem()),N.Unipotent(),3,2,false)
t2(n)=convert_problem(last(broyden_tridiagonal_trace_problem(n)),N.Free(),1,n,true)

function old_build(p,order,cs,ts)
    sparsity=compute_sparsity(p,SolverConfig(optimizer=nothing,order=order,cs_algo=cs,ts_algo=ts))
    moment=NCTSSoS.moment_relax(p,sparsity.corr_sparsity,sparsity.cliques_term_sparsities)
    (;sparsity,moment)
end
end

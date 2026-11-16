using NCTSSoS, JuMP, MosekTools
const solver=optimizer_with_attributes(MosekTools.Optimizer,"QUIET"=>true)
for trace in (false,true)
    r,ops=create_unipotent_variables(trace ? [("v",1:4)] : [("a",1:2),("b",1:2)])
    a,b=trace ? (ops[1][1:2],ops[1][3:4]) : ops
    symbol=trace ? NCTSSoS.tr : ς
    h=-symbol(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
    p=polyopt(h*one(typeof(a[1])),r)
    for (name,cs,ts) in (("dense",NoElimination(),NoElimination()),("CS",MF(),NoElimination()),
        ("TS",NoElimination(),MMD()),("CSTS",MF(),MMD()))
        config=SolverConfig(optimizer=solver,order=1,cs_algo=cs,ts_algo=ts)
        result=cs_nctssos(p,config;dualize=false)
        println((;trace,name,bound=result.objective,sides=reduce(vcat,result.moment_matrix_sizes),moments=result.n_unique_moment_matrix_elements))
        if name=="TS"
            display(result.sparsity)
            println("TS properties ",propertynames(first(result.sparsity.cliques_term_sparsities)[1]))
        end
        flush(stdout)
    end
end

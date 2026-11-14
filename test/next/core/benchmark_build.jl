# Standalone build-only comparison. No solver is invoked.
using NCTSSoS, COSMO, Statistics, Test
const Nx = NCTSSoS.Next

function measure_build(name,f::F,args::Vararg{Any,N}) where {F,N}
    f(args...)
    times = Float64[]
    for _ in 1:7
        GC.gc()
        push!(times,@elapsed f(args...))
    end
    GC.gc()
    result = @timed f(args...)
    println(name," median_ms=",round(1000median(times);digits=3),
        " q25_ms=",round(1000quantile(times,.25);digits=3)," q75_ms=",round(1000quantile(times,.75);digits=3),
        " bytes=",result.bytes," allocations=",Base.gc_alloc_count(result.gcstats))
    flush(stdout)
end
function new_build(algebra,W,n,order,objective)
    rows = Nx.basis_rows(W,algebra,1:n,order)
    return Nx.plain_program(rows,objective)
end
function new_standard(algebra,W,n,order,objective)
    p,t = new_build(algebra,W,n,order,objective)
    return Nx.moment_standard_form(p),t
end
old_build(pop,sp) = NCTSSoS.moment_relax(pop,sp.corr_sparsity,sp.cliques_term_sparsities)
function old_pipeline(pop,cfg)
    return old_build(pop,compute_sparsity(pop,cfg))
end
old_jump_pipeline(pop,cfg) = build_jump_model(old_pipeline(pop,cfg))
function compare_build(name,pop,algebra,W,n,order,convert_word)
    cfg = SolverConfig(optimizer=COSMO.Optimizer,order=order)
    sp = compute_sparsity(pop,cfg)
    old = old_build(pop,sp)
    objective = [(c,convert_word(m)) for (c,m) in terms(pop.objective)]
    p,t = new_build(algebra,W,n,order,objective)
    oldkeys = Set(Nx.moment_key(convert_word(m))[1] for m in values(old.linear.key_to_monomial))
    newkeys = Set(label.word for label in t.labels)
    @test oldkeys==newkeys
    println(name," rows_old=",length(sp.corr_sparsity.clq_mom_mtx_bases[1])," rows_new=",p.blocks[1].side,
        " keys=",length(newkeys)," real_variables=",p.nvars," entries=",length(p.forms)," terms=",length(p.forms.ids))
    measure_build("$name legacy moment_relax",old_build,pop,sp)
    measure_build("$name legacy rows+program",old_pipeline,pop,cfg)
    measure_build("$name legacy rows+program+JuMP",old_jump_pipeline,pop,cfg)
    measure_build("$name packed rows+program",new_build,algebra,W,n,order,objective)
    measure_build("$name packed rows+program+standard",new_standard,algebra,W,n,order,objective)
end
function main()
    println("Julia=",VERSION," CPU=",Sys.CPU_NAME," threads=",Threads.nthreads())
    reg,(x,y,z) = create_pauli_variables(1:10)
    ham = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:9)
    W = Nx.PauliWord{1}
    compare_build("Heisenberg10/d2",polyopt(ham,reg),Nx.Pauli(),W,10,2,m->Nx.word_from_legacy(W,m.word))
    reg,(a,ad) = create_fermionic_variables(1:4)
    ham = sum(-(ad[i]*a[i+1]+ad[i+1]*a[i]) for i in 1:3)+1.3ad[1]*a[1]*ad[2]*a[2]
    W = Nx.FermiWord{1}
    compare_build("Fermion4/d2",polyopt(ham,reg),Nx.Fermionic(),W,4,2,m->Nx.word_from_legacy(W,m.word))
    reg,(x,y) = create_unipotent_variables([("x",1:2),("y",1:2)])
    ham = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
    algebra = Nx.SiteAlgebra(Nx.Unipotent(),2)
    W = Nx.word_type(algebra,2,4)
    function convert_site(m)
        letters = [(mod1(NCTSSoS.decode_operator_id(l),2)<<8)|NCTSSoS.decode_site(l) for l in m.word]
        return Nx.word_from_legacy(W,letters)
    end
    compare_build("CHSH/d2",polyopt(ham,reg),algebra,W,2,2,convert_site)
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()

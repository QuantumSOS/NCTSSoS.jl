# Standalone part-1 validation/measurement, not part of runtests.jl.
# Needs MosekTools on LOAD_PATH and a working license. Run with one BLAS thread:
# OPENBLAS_NUM_THREADS=1 julia --project test/next/pipeline/benchmark_torus.jl
# The full degree-1 comparison took about 7 minutes in the orb.
# Timed builds exclude row/action construction, solver conversion and solves.
using NCTSSoS, LinearAlgebra, JuMP, MosekTools, Statistics
const Nx=NCTSSoS.Next
const SOLVER=optimizer_with_attributes(MosekTools.Optimizer,"QUIET"=>true,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-9,"MSK_IPAR_NUM_THREADS"=>1)
function bound(p)
    m,y=Nx.to_jump(p); set_optimizer(m,SOLVER); optimize!(m)
    @assert is_solved_and_feasible(m;allow_almost=false)
    objective_value(m)
end
function torus_case(nx,ny)
    n=nx*ny; W=Nx.word_type(Nx.Pauli(),n,2)
    site(x,y)=mod(x,nx)+nx*mod(y,ny)+1
    perm(f)=[site(f(x,y)...) for y in 0:ny-1 for x in 0:nx-1]
    tx=Nx.SitePermutation(perm((x,y)->(x+1,y))); ty=Nx.SitePermutation(perm((x,y)->(x,y+1)))
    point=Nx.AbstractAction[Nx.SitePermutation(perm((x,y)->(x,-y))),Nx.SitePermutation(perm((x,y)->(-y,x))),
        Nx.CliffordMap([Nx.pauli_letter(W,i,a) for a in (1,0) for i in 1:n]),
        Nx.CliffordMap([Nx.pauli_letter(W,i,a) for a in (1,2) for i in 1:n])]
    rows=Nx.basis_rows(W,Nx.Pauli(),1:n,1)
    objective=[(1.,Nx.mul(Nx.pauli_letter(W,site(x,y),a),Nx.pauli_letter(W,site(x+dx,y+dy),a))[2])
        for y in 0:ny-1 for x in 0:nx-1 for (dx,dy) in ((1,0),(0,1)) for a in 0:2]
    tr=Nx.TranslationGroup([tx,ty],[nx,ny])
    rows,objective,tr,point
end
println("Julia=",VERSION," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads()," CPU=",Sys.cpu_info()[1].model)
license=Model(SOLVER); @variable(license,x>=1); @objective(license,Min,x); optimize!(license)
@assert termination_status(license)==MOI.OPTIMAL && objective_value(license)≈1
println("license=OPTIMAL objective=",objective_value(license))
rows,objective,tr,point=torus_case(6,6)
build()=Nx.abelian_lattice_program(rows,objective,tr,point;moment_generators=(Nx.AntiunitaryPauli(one(eltype(rows))),))
cold=@timed build(); p=cold.value[1]
println("cold seconds=",cold.time," bytes=",cold.bytes," rows=",length(rows))
flush(stdout)
reduced=bound(p); plain,_=Nx.plain_program(rows,objective); dense=bound(plain)
@assert abs(reduced-dense)<1e-6
@assert abs(reduced+216)<1e-6
println("bound reduced=",reduced," plain=",dense," error=",abs(reduced-dense),
    " largest side=",maximum(b.side for b in p.blocks)," blocks=",length(p.blocks)," vars=",p.nvars)
@assert all(b.kind==Nx.PSD for b in p.blocks)
atoms=[[Nx.pauli_letter(eltype(rows),i,a) for a in 0:2] for i in 1:36]
discovered,_,report=Nx.generators_program(rows,objective,Nx.AbstractAction[tr.generators;point;Nx.AntiunitaryPauli(one(eltype(rows)))],atoms)
@assert abs(bound(discovered)-reduced)<1e-6
println("discovery order=",report.translation_order," source=",report.discovery.source)
flush(stdout)
samples=map(1:5) do _
    GC.gc(); m=@timed build(); (m.time,m.bytes)
end
println("warm samples=",samples)
println("warm median seconds=",median(first.(samples))," IQR=",quantile(first.(samples),[.25,.75])," bytes=",median(last.(samples)))

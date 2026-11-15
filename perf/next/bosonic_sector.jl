using NCTSSoS, Mosek, MosekTools, JuMP, LinearAlgebra, Test
const N=NCTSSoS.Next
const backend=N.MosekBackend(MSK_IPAR_NUM_THREADS=1,MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-10,
    MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-10,MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-10)
const factory=optimizer_with_attributes(MosekTools.Optimizer,"QUIET"=>true,"MSK_IPAR_NUM_THREADS"=>1,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-8,"MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>1e-8,
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>1e-8)
function exact_sector(L,edges)
    occupations=filter(v->sum(v)==L,vec(collect(Iterators.product(ntuple(_->0:L,L)...))))
    index=Dict(v=>i for (i,v) in enumerate(occupations))
    h=zeros(length(occupations),length(occupations))
    for (col,v) in enumerate(occupations)
        h[col,col]=sum(k*(k-1) for k in v)
        for (i,j) in edges, (to,from) in ((i,j),(j,i))
            v[from]==0 && continue
            u=ntuple(k->v[k]+(k==to)-(k==from),L)
            h[index[u],col]-=sqrt(v[from]*(v[to]+1))
        end
    end
    @assert h≈h'
    eigmin(Symmetric(h))
end
ham(b,bd,edges)=-sum(bd[i]*b[j]+bd[j]*b[i] for (i,j) in edges)+sum(bd[i]^2*b[i]^2 for i in eachindex(b))
function main()
    _,(_,_,z)=N.pauli_variables(1)
    @test N.solve(N.Problem(z[1]),N.Relaxation(order=1);solver=backend).bound≈-1 atol=1e-7
    open(isempty(ARGS) ? joinpath(@__DIR__,"bosonic_sector.tsv") : only(ARGS),"w") do io
        println(io,"sites\tperiodic\torder\tsymmetry\texact\tlegacy\tbound\tdifference\tvariables\telapsed")
        for (L,periodic,order) in ((2,false,2),(2,false,3),(3,false,2),(3,true,2),(4,false,2),(4,true,2))
            edges=[(i,i+1) for i in 1:L-1]
            periodic && push!(edges,(L,1))
            exact=exact_sector(L,edges)
            reg,(b,bd)=N.bosonic_variables(L)
            oldreg,(oldb,oldbd)=create_bosonic_variables(1:L)
            old=cs_nctssos(polyopt(ham(oldb,oldbd,edges),oldreg;moment_eq_constraints=[sum(oldbd.*oldb)-L]),
                SolverConfig(optimizer=factory,order=order,cs_algo=NoElimination(),ts_algo=NoElimination());dualize=false)
            @test isfinite(old.objective)
            h=ham(b,bd,edges); relax=N.Relaxation(;order,strengthen=N.Sector(sum(bd.*b),L))
            rotation=N.NumberRotation(reg)
            t=N.SitePermutation([2:L;1]); reflection=N.SitePermutation([1;L:-1:2])
            modes=periodic ? (:none,:number,:translation,:dihedral,:both,:moment_only) : (:none,:number)
            base=NaN; basevars=0
            for mode in modes
                symmetry=mode==:none ? N.Symmetry() : mode==:number ? N.Symmetry(;continuous=(rotation,)) :
                    mode==:translation ? N.Symmetry(t) : mode==:dihedral ? N.Symmetry(t,reflection) :
                    mode==:both ? N.Symmetry(t,reflection;continuous=(rotation,)) : N.Symmetry(;moment_only=(t,))
                sample=@timed N.solve(N.Problem(h;symmetry),relax;solver=backend)
                result=sample.value; vars=result.sizes.program.real_variables
                if mode==:none
                    base=result.bound; basevars=vars
                else
                    @test result.bound≈base atol=1e-6
                    @test vars<basevars
                end
                @test result.bound<=exact+1e-6
                @test result.bound>=old.objective-1e-6
                if L==2
                    @test result.bound≈1-sqrt(5) atol=1e-6
                    @test result.bound≈old.objective atol=1e-6
                end
                row=join((L,periodic,order,mode,exact,old.objective,result.bound,result.bound-old.objective,vars,sample.time),'\t')
                println(io,row); flush(io); println(row); flush(stdout)
            end
        end
    end
end
main()

using Test, NCTSSoS, LinearAlgebra, JuMP, MosekTools, Graphs, TOML
include("experiment.jl")
const E = SymmetrySparsityExperiment
const N = NCTSSoS
BLAS.set_num_threads(4)
const OPT = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>1e-8,
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>1e-8,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-8,
    "MSK_IPAR_NUM_THREADS"=>4)
const results = Dict{String,Any}[]
const X=ComplexF64[0 1;1 0]
const Y=ComplexF64[0 -im;im 0]
const Z=ComplexF64[1 0;0 -1]
const ID=Matrix{ComplexF64}(I,2,2)
function matrixword(w,n)
    factors=[ID for _ in 1:n]
    for k in w.word
        factors[N._pauli_site(k)]=(X,Y,Z)[N._pauli_type(k)+1]
    end
    reduce(kron,factors)
end
function enumerate_actions(actions,n)
    queue=[Tuple(1:n)]; seen=Set(queue); i=1
    while i<=length(queue)
        for a in actions
            b=Tuple(E.compose(a,queue[i]))
            if b ∉ seen
                push!(seen,b);push!(queue,b)
            end
        end
        i+=1
    end
    seen
end
function solvecheck(H,rows,groups,gs;reduced,dual_form=true)
    println("ASSEMBLE reduced=",reduced," groups=",length(groups)," max=",maximum(length,groups));flush(stdout)
    built=E.assemble(H,rows,groups,gs,OPT;reduced,dual_form)
    println("SOLVE reduced=",reduced," sizes=",built.sizes);flush(stdout)
    optimize!(built.model)
    @test termination_status(built.model)==MOI.OPTIMAL
    moment_values=built.moment_values()
    valueword(w)=begin
        i,s=built.coordinates.coordinates[w]
        i==0 ? 0.0 : s*moment_values[i]
    end
    mineig=Inf
    for C in groups
        M=[begin
            phase,w=E.product(rows[i],rows[j]); phase*valueword(w)
        end for i in C,j in C]
        @test norm(M-M') < 1e-10
        mineig=min(mineig,eigmin(Hermitian(M)))
    end
    @test mineig > -2e-6
    @test valueword(one(first(rows))) ≈ 1 atol=1e-8
    primal=sum(real(c)*valueword(w) for (c,w) in terms(H))
    dual=dual_form ? objective_value(built.model) : dual_objective_value(built.model)
    @test abs(primal-dual)<2e-6
    return (;primal,dual,mineig,sizes=built.sizes,variables=built.coordinates.nvars)
end

@testset "Small-system symmetry-sparsity validation" begin
    # A real solve, not an import, checks the license before relying on Mosek.
    license=Model(OPT); set_silent(license)
    @variable(license,t>=1); @objective(license,Min,t); optimize!(license)
    @test termination_status(license)==MOI.OPTIMAL
    @test value(t)≈1
    println("MOSEK_LICENSE_SOLVE_PASS");flush(stdout)

    @testset "Independent two-qubit group and stabilizers" begin
        reg,(x,y,z)=create_pauli_variables(1:2)
        rows=get_ncbasis(reg,2); mats=matrixword.(rows,2)
        h=(X+Z)/sqrt(2); s=ComplexF64[1 0;0 im]
        cnot=ComplexF64[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]
        signature(U)=[begin
            overlaps=[real(LinearAlgebra.tr(P'*U*Q*U'))/4 for P in mats]
            i=argmax(abs.(overlaps));round(Int,overlaps[i])*i
        end for Q in mats]
        allgroup=enumerate_actions(signature.([kron(h,ID),kron(ID,h),kron(s,ID),kron(ID,s),cnot]),length(rows))
        @test length(allgroup)==11520
        H=x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
        lookup=Dict(w=>i for (i,w) in enumerate(rows))
        cs=zeros(length(rows))
        for (c,w) in terms(H);cs[lookup[w]]=c;end
        exact=Set(a for a in allgroup if all(cs[abs(a[i])]==sign(a[i])*cs[i] for i in eachindex(rows)))
        @test length(exact)==192
        gs=sympleq_symmetry_spec(H).clifford_generators
        actions=E.actions_on(rows,gs)
        @test enumerate_actions(actions,length(rows))==exact
        seedgroups=E.window_groups(rows,rows,2,1)
        groups=E.close_groups(seedgroups,actions)
        @test Set(E.close_groups(seedgroups,reverse(actions)))==Set(groups)
        for orbit in E.block_orbits(groups,actions,length(rows))
            C=orbit.representative; loc=Dict(v=>i for (i,v) in enumerate(C))
            expected=Set(Tuple(sign(a[v])*loc[abs(a[v])] for v in C)
                for a in exact if E.image(collect(a),C)==C)
            actual=enumerate_actions(orbit.stabilizers,length(C))
            @test actual==expected
        end
        @test E.product(x[1],y[1])==(1im,z[1])
        @test E.product(y[1],x[1])==(-1im,z[1])
        for i in eachindex(rows), j in eachindex(rows)
            phase,w=E.product(rows[i],rows[j])
            @test phase*matrixword(w,2)≈mats[i]'*mats[j] atol=1e-12
        end
        for a in actions
            @test E.compose(a,E.inverse(a))==collect(eachindex(rows))
            for (i,j) in ((2,5),(3,8),(7,14))
                p,w=E.product(rows[i],rows[j])
                q,v=E.product(rows[abs(a[i])],rows[abs(a[j])])
                ai=a[lookup[w]]
                @test (p*sign(ai),rows[abs(ai)])==(sign(a[i])*sign(a[j])*q,v)
            end
        end
        coords=E.moment_coordinates(H,rows,[Tuple(eachindex(rows))],gs)
        @test coords.coordinates[x[1]][1]==0
        @test coords.coordinates[y[1]][1]==0
        @test_throws E.LimitExceeded E.close_groups(seedgroups,actions;budget=E.Budget(max_items=1))
        @test_throws E.LimitExceeded E.close_rows([one(x[1]),x[1]],gs;budget=E.Budget(max_items=2))
        @test_throws E.LimitExceeded E.close_groups(seedgroups,actions;budget=E.Budget(seconds=-1))
        @test E.decompose_signed(Vector{Int}[],3)==[Matrix{ComplexF64}(I,3,3)]
        @test_throws ErrorException E.moment_coordinates(H,rows,[(lookup[one(x[1])],)],gs)
    end

    @testset "Asymmetric complex objective, trivial group" begin
        reg,(x,y,z)=create_pauli_variables(1:2)
        rows=get_ncbasis(reg,2)
        H=0.37x[1]+0.83y[1]+1.19z[2]+0.41y[1]*x[2]+0.23z[1]*z[2]
        exact=eigmin(Hermitian(sum(c*matrixword(w,2) for (c,w) in terms(H))))
        for reduced in (false,true), dual_form in (false,true)
            result=solvecheck(H,rows,[Tuple(eachindex(rows))],CliffordSymmetry[];reduced,dual_form)
            @test result.primal≈exact atol=2e-6
        end
    end

    @testset "Completion must preserve a four-cycle rotation" begin
        C=[1,2,3,4]; G=cycle_graph(4); rotation=[2,3,4,1]
        firstfill=copy(G)
        for clique in N.clique_decomp(G,N.NaturalElimination())
            for j in eachindex(clique), i in 1:j-1
                add_edge!(firstfill,clique[i],clique[j])
            end
        end
        @test ne(firstfill)==5
        @test any(e->!has_edge(firstfill,rotation[src(e)],rotation[dst(e)]),edges(firstfill))
        graphs=Dict(C=>copy(G))
        result=E.complete_graphs!(graphs,[C],[rotation])
        @test result.blocks==[[1,2,3,4]]
        @test ne(graphs[C])==6
        @test first(N.CheckChordal(graphs[C]))
        @test_throws E.LimitExceeded E.complete_graphs!(Dict(C=>copy(G)),[C],[rotation];
            budget=E.Budget(max_entries=4))
    end

    @testset "Route E and level-1 D/F solve checks" begin
        for (n,degree,width) in ((2,2,2),(4,2,2),(4,4,4))
            println("CASE n=",n," degree=",degree," window=",width);flush(stdout)
            reg,ops=create_pauli_variables(1:n)
            H=n==2 ? sum(v[1]*v[2] for v in ops)/4 :
                heisenberg_chain_hamiltonian(ops;coupling=1/4,periodic=true)
            seed=pauli_contiguous_chain_basis(reg,degree)
            gs=sympleq_symmetry_spec(H).clifford_generators
            @test all(g->N._act_polynomial(g,H)==H,gs)
            rows=E.close_rows(seed,gs); actions=E.actions_on(rows,gs)
            println("ROWS ",length(rows));flush(stdout)
            parents=E.close_groups(E.window_groups(seed,rows,n,width),actions)
            println("PARENTS ",length(parents));flush(stdout)
            @test Set(v for C in parents for v in C)==Set(eachindex(rows))
            @test Set(E.close_groups(parents,actions))==Set(parents)
            exact=eigmin(Hermitian(sum(c*matrixword(w,n) for (c,w) in terms(H))))
            S=Set([one(first(rows));monomials(H)])
            println("TERM D");flush(stdout)
            D=E.term_blocks([Tuple(eachindex(rows))],rows,actions,S)
            println("TERM F");flush(stdout)
            F=E.term_blocks(parents,rows,actions,S)
            @test all(C->any(P->issubset(Set(C),Set(P)),parents),F.blocks)
            bounds=Dict{String,Float64}()
            for (route,groups) in (("E",parents),("D",D.blocks),("F",F.blocks))
                reduced=solvecheck(H,rows,groups,gs;reduced=true)
                reference=solvecheck(H,rows,groups,gs;reduced=false)
                @test reduced.primal≈reference.primal atol=2e-6
                @test reduced.primal <= exact+2e-6
                if n==2
                    @test reduced.primal≈-0.75 atol=2e-6
                    direct=solvecheck(H,rows,groups,gs;reduced=false,dual_form=false)
                    @test reduced.primal≈direct.primal atol=2e-6
                end
                bounds[route]=reduced.primal
                item=Dict{String,Any}("n"=>n,"degree"=>degree,"window"=>width,"route"=>route,
                    "rows"=>length(rows),"blocks"=>length(groups),"reduced_sizes"=>reduced.sizes,
                    "reduced_bound"=>reduced.primal,"unreduced_bound"=>reference.primal,
                    "exact_energy"=>exact,"minimum_full_block_eigenvalue"=>reduced.mineig,
                    "dual_bound"=>reduced.dual,"variables"=>reduced.variables)
                push!(results,item); println("PASS ",item);flush(stdout)
            end
            @test bounds["F"]<=bounds["E"]+2e-6
        end
    end
end
open(joinpath(@__DIR__,"validation.toml"),"w") do io
    TOML.print(io,Dict("results"=>results,"julia"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads()))
end
println("ALL_VALIDATION_PASSED")

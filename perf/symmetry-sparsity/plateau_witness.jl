using NCTSSoS, LinearAlgebra, TOML
include("experiment.jl")
include("plateau_check.jl")
const E=SymmetrySparsityExperiment
const N=NCTSSoS

function expectation(w,psi)
    value=0.0im
    for state in 0:length(psi)-1
        target=state;phase=1.0+0im
        for k in w.word
            site=N._pauli_site(k)-1;axis=N._pauli_type(k)
            @assert site<6
            bit=(state>>site)&1
            if axis==0
                target ⊻= 1<<site
            elseif axis==1
                target ⊻= 1<<site
                phase*=bit==0 ? im : -im
            else
                phase*=bit==0 ? 1 : -1
            end
        end
        value+=conj(psi[target+1])*phase*psi[state+1]
    end
    @assert abs(imag(value))<1e-12
    real(value)
end

function run()
    energy6,psi=ring_ground(6)
    reg,ops=create_pauli_variables(1:20)
    H=heisenberg_chain_hamiltonian(ops;coupling=1/4,periodic=true)
    seed=pauli_contiguous_chain_basis(reg,4)
    gs=sympleq_symmetry_spec(H).clifford_generators
    @assert all(g->N._act_polynomial(g,H)==H,gs)
    rows=E.close_rows(seed,gs);actions=E.actions_on(rows,gs)
    groups=E.close_groups(E.window_groups(seed,rows,20,6),actions)
    println("BUILD_COORDINATES rows=",length(rows)," blocks=",length(groups));flush(stdout)
    coords=E.moment_coordinates(H,rows,groups,gs;budget=E.Budget(max_items=2000000))
    lookup=Dict(w=>i for (i,w) in enumerate(rows))
    C=[lookup[w] for w in seed if all(k->N._pauli_site(k)<=6,w.word)]
    @assert length(C)==415
    localwords=Set(last(E.product(rows[i],rows[j])) for i in C,j in C)
    values=fill(NaN,coords.nvars); consistency=0.0
    for w in localwords
        v=expectation(w,psi);k,s=coords.coordinates[w]
        if k==0
            consistency=max(consistency,abs(v))
        elseif isnan(values[k])
            values[k]=s*v
        else
            consistency=max(consistency,abs(values[k]-s*v))
        end
    end
    println("LOCAL_ASSIGNMENT assigned=",count(!isnan,values)," variables=",length(values),
        " max_conflict=",consistency);flush(stdout)
    @assert all(!isnan,values)
    @assert consistency<1e-10
    valueword(w)=begin
        k,s=coords.coordinates[w];k==0 ? 0.0 : s*values[k]
    end
    symmetry_residual=0.0
    for w in keys(coords.coordinates),g in gs
        s,v=N._act_monomial(g,w)
        symmetry_residual=max(symmetry_residual,abs(valueword(w)-s*valueword(v)))
    end
    @assert symmetry_residual<1e-10
    orbits=E.block_orbits(groups,actions,length(rows))
    mineig=Inf
    for o in orbits
        M=[begin
            phase,w=E.product(rows[i],rows[j]);phase*valueword(w)
        end for i in o.representative,j in o.representative]
        @assert norm(M-M')<1e-10
        mineig=min(mineig,eigmin(Hermitian(M)))
    end
    @assert mineig>=-1e-10
    objective=sum(real(c)*valueword(w) for (c,w) in terms(H))
    @assert abs(valueword(one(first(rows)))-1)<1e-12
    @assert abs(objective-20energy6/6)<1e-10
    result=Dict("objective"=>objective,"energy6"=>energy6,
        "rows"=>length(rows),"generators"=>length(gs),"blocks"=>length(groups),
        "block_orbits"=>length(orbits),"moment_variables"=>coords.nvars,
        "local_assignment_conflict"=>consistency,"symmetry_residual"=>symmetry_residual,
        "minimum_original_block_eigenvalue"=>mineig)
    open(joinpath(@__DIR__,"plateau-witness.toml"),"w") do io
        TOML.print(io,result)
    end
    println("WITNESS_PASS ",result)
end
run()

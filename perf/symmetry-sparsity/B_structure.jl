using NCTSSoS, LinearAlgebra, JuMP, MosekTools, TOML
include("experiment.jl")
include("plateau_check.jl")
const E=SymmetrySparsityExperiment
const N=NCTSSoS
const OPT=optimizer_with_attributes(Mosek.Optimizer,"MSK_IPAR_NUM_THREADS"=>4)

function run()
    license=Model(OPT);set_silent(license)
    @variable(license,t>=1);@objective(license,Min,t);optimize!(license)
    @assert termination_status(license)==MOI.OPTIMAL && abs(value(t)-1)<1e-8
    energy,psi=ring_ground(6)
    correlations=[sum(psi[s+1]*psi[(s ⊻ 1 ⊻ (1<<d))+1] for s in 0:63) for d in 1:5]
    println("SIX_SITE_XX ",correlations);flush(stdout)
    completions=[]
    for known in (4,5)
        model=Model(OPT);set_silent(model)
        @variable(model,c[1:10]);@variable(model,margin)
        for d in 1:known
            fix(c[d],correlations[d];force=true)
        end
        lambda=[@expression(model,1+2sum(c[d]*cos(2pi*k*d/20) for d in 1:9)+c[10]*(-1)^k)
            for k in 0:19]
        constraints=[@constraint(model,a>=margin) for a in lambda]
        @objective(model,Max,margin);optimize!(model)
        @assert termination_status(model)==MOI.OPTIMAL
        cv=value.(c); lambdas=value.(lambda)
        C=[i==j ? 1.0 : cv[min(mod(i-j,20),mod(j-i,20))] for i in 1:20,j in 1:20]
        @assert abs(eigmin(Symmetric(C))-minimum(lambdas))<1e-10
        # A convex combination of Fourier eigenvalues gives an independent
        # upper bound on their minimum. Unknown-distance coefficients cancel.
        weights=dual.(constraints)
        println("LP_DIAGNOSTIC known=",known," margin=",objective_value(model),
            " weight_min=",minimum(weights)," weight_max=",maximum(weights)," weight_sum=",sum(weights));flush(stdout)
        # Make the weights nonnegative explicitly. Bound the small remaining
        # coefficients using |c[d]| <= 1, necessary for a PSD unit diagonal.
        weights=max.(weights,0);weights/=sum(weights)
        @assert minimum(weights)>=0 && abs(sum(weights)-1)<1e-10
        coeff=[sum(weights[k+1]*(d==10 ? (-1)^k : 2cos(2pi*k*d/20)) for k in 0:19) for d in 1:10]
        residual=maximum(abs,coeff[known+1:end])
        upper=1+sum(coeff[d]*correlations[d] for d in 1:known)
        @assert abs(upper-objective_value(model))<1e-8
        certified_upper=upper+sum(abs,coeff[known+1:end])
        result=Dict("known_distances"=>known,"max_min_eigenvalue"=>objective_value(model),
            "weighted_eigenvalue_upper_if_PSD"=>certified_upper,"free_coefficient_residual"=>residual,
            "weights"=>weights,"known_coefficients"=>coeff[1:known],
            "completion"=>cv,"eigenvalues"=>lambdas)
        push!(completions,result)
        println("GLOBAL_XX_COMPLETION ",result);flush(stdout)
    end

    reg,ops=create_pauli_variables(1:20);x,y,z=ops
    H=heisenberg_chain_hamiltonian(ops;coupling=1/4,periodic=true)
    seed=pauli_contiguous_chain_basis(reg,4)
    gs=sympleq_symmetry_spec(H).clifford_generators
    rows=E.close_rows(seed,gs);actions=E.actions_on(rows,gs)
    groups=E.close_groups(E.window_groups(seed,rows,20,5),actions)
    lookup=Dict(w=>i for (i,w) in enumerate(rows))
    pair=[lookup[x[1]],lookup[x[11]]]
    pair_present=any(C->issubset(pair,C),groups)
    original=Set(lookup[w] for w in seed)
    densepairs=Set{Tuple{Int,Int}}()
    for C in groups
        retained=filter(i->i in original,C)
        for j in eachindex(retained),i in 1:j
            push!(densepairs,minmax(retained[i],retained[j]))
        end
    end
    global_x=sort!([lookup[one(first(rows))];[lookup[w] for w in x]])
    extra=E.close_groups([global_x],actions)
    out=Dict("xx_correlations6"=>correlations,"completions"=>completions,
        "B_rows"=>length(seed),"E5_rows"=>length(rows),"E5_blocks"=>length(groups),
        "X1_X11_pair_in_E5"=>pair_present,"B_upper_triangle_pairs"=>length(seed)*(length(seed)+1)÷2,
        "B_pairs_in_E5"=>length(densepairs),"global_X_block_size"=>length(global_x),
        "global_X_block_orbit_size"=>length(extra))
    open(joinpath(@__DIR__,"B-structure.toml"),"w") do io
        TOML.print(io,out)
    end
    println("STRUCTURE_PASS ",filter(p->first(p)!="completions",out))
end
run()

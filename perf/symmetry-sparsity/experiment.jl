module SymmetrySparsityExperiment
using NCTSSoS, LinearAlgebra, JuMP, Graphs, SparseArrays
const N = NCTSSoS

struct LimitExceeded <: Exception
    message::String
end
Base.showerror(io::IO, e::LimitExceeded) = print(io,e.message)
mutable struct Budget
    deadline::Float64
    max_items::Int
    max_entries::Int
end
Budget(; seconds=1200, max_items=100_000, max_entries=5_000_000) =
    Budget(time()+seconds,max_items,max_entries)
function check(b::Budget, items=0, entries=0)
    time() <= b.deadline || throw(LimitExceeded("time limit; no partial result"))
    items <= b.max_items || throw(LimitExceeded("item limit: $items > $(b.max_items); no partial result"))
    entries <= b.max_entries || throw(LimitExceeded("entry limit: $entries > $(b.max_entries); no partial result"))
end

# Signed permutations: a[i] = sign * destination; compose(a,b) applies b first.
compose(a,b) = [sign(j)*a[abs(j)] for j in b]
function inverse(a)
    result = similar(a)
    for (i,j) in enumerate(a)
        result[abs(j)] = sign(j)*i
    end
    result
end
# Canonical block IDs are sorted vectors, never mutated after insertion in a
# set/dictionary. Long tuples encode basis size in the type and overwhelm JIT.
image(a,C) = sort!(abs.(a[collect(C)]))
function actions_on(rows, gs)
    lookup = Dict(w=>i for (i,w) in enumerate(rows))
    actions = map(gs) do g
        map(rows) do w
            s,v = N._act_monomial(g,w)
            haskey(lookup,v) || error("row action is not closed")
            s*lookup[v]
        end
    end
    @assert all(a->sort(abs.(a))==collect(eachindex(rows)),actions)
    actions
end
function close_rows(seed,gs; budget=Budget())
    rows=copy(seed); seen=Set(rows); i=1
    while i<=length(rows)
        check(budget,length(rows))
        for g in gs
            _,v=N._act_monomial(g,rows[i])
            if v ∉ seen
                check(budget,length(rows)+1)
                push!(seen,v); push!(rows,v)
            end
        end
        i+=1
    end
    sort!(rows)
end
function window_groups(seed,rows,n,width)
    lookup=Dict(w=>i for (i,w) in enumerate(rows))
    groups=Vector{Int}[]
    for i in 1:n
        sites=Set(mod1.(i:i+width-1,n))
        push!(groups,sort!([lookup[w] for w in seed
            if all(k->N._pauli_site(k) in sites,w.word)]))
    end
    sort!(unique(groups))
end
function close_groups(initial,actions; budget=Budget())
    groups=sort!(unique(collect.(initial))); seen=Set(groups); i=1
    entries=sum(length,groups)
    check(budget,length(groups),entries)
    while i<=length(groups)
        check(budget,length(groups),entries)
        for a in actions
            C=image(a,groups[i])
            if C ∉ seen
                check(budget,length(groups)+1,entries+length(C))
                push!(seen,C); push!(groups,C); entries+=length(C)
            end
        end
        i+=1
    end
    sort!(groups)
end

function product(u,v)
    word,phase=N.simplify(PauliAlgebra,vcat(u.word,v.word))
    return (1im)^Int(phase),typeof(u)(word)
end

# Return signed stabilizer actions restricted to each representative block.
# Only the orbit of a BLOCK is enumerated, never the Clifford group.
function block_orbits(groups,actions,nrows; budget=Budget())
    groups=collect.(groups)
    available=Set(groups); covered=Set{Vector{Int}}(); results=[]
    for C in groups
        C in covered && continue
        orbit=[C]; transport=Dict(C=>collect(1:nrows)); localgens=Vector{Int}[]
        seen=Set{Vector{Int}}(); i=1; lookup=Dict(v=>k for (k,v) in enumerate(C))
        while i<=length(orbit)
            check(budget,length(orbit),length(orbit)*nrows)
            D=orbit[i]
            for a in actions
                target=image(a,D)
                target in available || error("block collection is not closed")
                path=compose(a,transport[D])
                if !haskey(transport,target)
                    check(budget,length(orbit)+1,(length(orbit)+1)*nrows)
                    transport[target]=path; push!(orbit,target)
                end
                returnmap=compose(inverse(transport[target]),path)
                @assert image(returnmap,C)==C
                restricted=[sign(returnmap[v])*lookup[abs(returnmap[v])] for v in C]
                key=restricted
                if restricted != collect(1:length(C)) && key ∉ seen
                    push!(seen,key); push!(localgens,restricted)
                end
            end
            i+=1
        end
        union!(covered,orbit)
        push!(results,(representative=C,orbit=orbit,stabilizers=localgens))
    end
    @assert covered==available
    results
end
function decompose_signed(actions,n)
    isempty(actions) && return [Matrix{ComplexF64}(I,n,n)]
    signed=[(abs.(a),sign.(a)) for a in actions]
    labels,signs,count=N._signed_pair_orbits(n,signed)
    for trial in 1:3
        blocks=N._generator_blocks(n,signed,labels,signs,count,trial)
        isnothing(blocks) || return blocks
    end
    error("numerical stabilizer decomposition failed")
end

# The same S is used throughout: fill never advances the sparsity level.
function term_blocks(parents,rows,actions,S; budget=Budget())
    parents=collect.(parents)
    graphs=Dict(C=>SimpleGraph(length(C)) for C in parents)
    for C in parents
        check(budget)
        G=graphs[C]
        for j in eachindex(C), i in 1:j-1
            _,w=product(rows[C[i]],rows[C[j]])
            w in S && add_edge!(G,i,j)
        end
    end
    complete_graphs!(graphs,parents,actions;budget)
end

function complete_graphs!(graphs,parents,actions;budget=Budget())
    indices=Dict(C=>Dict(v=>i for (i,v) in enumerate(C)) for C in parents)
    rounds=0
    while true
        rounds+=1; check(budget)
        before=sum(ne,values(graphs))
        for G in values(graphs)
            for clique in N.clique_decomp(G,N.NaturalElimination())
                for j in eachindex(clique), i in 1:j-1
                    add_edge!(G,clique[i],clique[j])
                end
            end
        end
        # Iterate propagation to closure. Snapshot edges because a target may
        # also be the current graph; never mutate an active edge iterator.
        changed=true
        while changed
            changed=false
            for C in parents, a in actions
                check(budget,length(parents),sum(ne,values(graphs)))
                D=image(a,C); haskey(graphs,D) || error("parent blocks are not closed")
                for e in collect(edges(graphs[C]))
                    u=indices[D][abs(a[C[src(e)]])]
                    v=indices[D][abs(a[C[dst(e)]])]
                    changed |= add_edge!(graphs[D],u,v)
                end
            end
        end
        after=sum(ne,values(graphs)); check(budget,length(parents),after)
        if after==before
            @assert all(G->first(N.CheckChordal(G)),values(graphs))
            blocks=sort!(unique([sort!(C[c]) for C in parents
                for c in maximal_cliques(graphs[C])]))
            @assert Set(close_groups(blocks,actions;budget))==Set(blocks)
            return (;blocks,rounds,edge_count=after)
        end
    end
end

function moment_coordinates(H,rows,groups,gs; budget=Budget())
    words=Set{eltype(rows)}()
    for C in groups
        check(budget,length(words))
        for j in eachindex(C), i in 1:j
            _,w=product(rows[C[i]],rows[C[j]])
            push!(words,w)
        end
        check(budget,length(words))
    end
    issubset(Set(monomials(H)),words) || error("objective is not covered by retained blocks")
    ordered=sort!(collect(words)); coordinates=Dict{eltype(rows),Tuple{Int,Int}}()
    nvars=0; zero_orbits=0
    for seed in ordered
        haskey(coordinates,seed) && continue
        orbit=[seed]; signs=Dict(seed=>1); i=1; forcedzero=false
        while i<=length(orbit)
            check(budget,length(words))
            for g in gs
                s,v=N._act_monomial(g,orbit[i]); s*=signs[orbit[i]]
                v in words || error("moment support not closed")
                if haskey(signs,v)
                    forcedzero |= signs[v]!=s
                else
                    signs[v]=s; push!(orbit,v)
                end
            end
            i+=1
        end
        forcedzero ? (zero_orbits+=1) : (nvars+=1)
        for w in orbit
            coordinates[w]=(forcedzero ? 0 : nvars, signs[w])
        end
    end
    (;coordinates,nvars,zero_orbits)
end

function projected_forms(rows,C,transforms,coords,y; atol=1e-12)
    # M = sum_k K[k] y[k]. Project numerical sparse coefficients, not a
    # dense matrix of AffExpr objects. Product phases remain exact here.
    n=length(C)
    ii=[Int[] for _ in y]; jj=[Int[] for _ in y]
    vv=[ComplexF64[] for _ in y]
    for j in 1:n, i in 1:j
        phase,w=product(rows[C[i]],rows[C[j]])
        k,s=coords.coordinates[w]
        k==0 && continue
        push!(ii[k],i);push!(jj[k],j);push!(vv[k],phase*s)
        if i!=j
            push!(ii[k],j);push!(jj[k],i);push!(vv[k],conj(phase)*s)
        end
    end
    coefficients=[sparse(ii[k],jj[k],vv[k],n,n) for k in eachindex(y)]
    forms=[]
    for U in transforms
        m=size(U,1)
        A=[AffExpr(0.0) for _ in 1:m, _ in 1:m]
        B=[AffExpr(0.0) for _ in 1:m, _ in 1:m]
        work=zeros(ComplexF64,m,n); projected=zeros(ComplexF64,m,m)
        discarded=zeros(m,m)
        for (k,K) in enumerate(coefficients)
            mul!(work,U,K);mul!(projected,work,U')
            @assert norm(projected-projected')<1e-8
            for j in 1:m, i in 1:m
                a=real(projected[i,j]); b=imag(projected[i,j])
                abs(a)>atol ? add_to_expression!(A[i,j],a,y[k]) : (discarded[i,j]+=abs(a))
                abs(b)>atol ? add_to_expression!(B[i,j],b,y[k]) : (discarded[i,j]+=abs(b))
            end
        end
        @assert maximum(discarded)<1e-8
        push!(forms,A+im*B)
    end
    forms
end

function assemble(H,rows,groups,gs,optimizer; reduced=true,dual_form=false,budget=Budget(),
    coordinates=nothing,block_plans=nothing,numeric_projection=false)
    coords=isnothing(coordinates) ? moment_coordinates(H,rows,groups,gs;budget) : coordinates
    model=Model(optimizer); set_silent(model)
    @variable(model,y[1:coords.nvars])
    valueword(w)=begin
        i,s=coords.coordinates[w]
        i==0 ? AffExpr(0.0) : s*y[i]+0.0
    end
    @constraint(model,valueword(one(first(rows)))==1)
    @objective(model,Min,sum(real(c)*valueword(w) for (c,w) in terms(H)))
    @assert all(c->iszero(imag(c)),coefficients(H))
    plans=if !isnothing(block_plans)
        block_plans
    elseif reduced
        actions=actions_on(rows,gs)
        [(o.representative,decompose_signed(o.stabilizers,length(o.representative)))
            for o in block_orbits(groups,actions,length(rows);budget)]
    else
        [(C,[Matrix{ComplexF64}(I,length(C),length(C))]) for C in groups]
    end
    sizes=Int[]
    realblocks=[]
    for (C,transforms) in plans
        check(budget)
        projected=if numeric_projection && reduced
            projected_forms(rows,C,transforms,coords,y)
        else
            M=[begin
                phase,w=product(rows[i],rows[j]); phase*valueword(w)
            end for i in C,j in C]
            [reduced ? U*M*U' : M for U in transforms]
        end
        for (U,R) in zip(transforms,projected)
            # Realification of a Hermitian PSD matrix preserves positivity.
            A=real.(R); B=imag.(R)
            L=Symmetric([A -B; B A])
            if dual_form
                push!(realblocks,L)
            else
                @constraint(model,L in PSDCone())
            end
            push!(sizes,size(U,1))
        end
    end
    if !dual_form
        return (;model,sizes,coordinates=coords,moment_values=()->value.(y))
    end
    # Explicit conic dual: few coefficient equations instead of one affine
    # equation per entry of every large unreduced PSD block.
    dualmodel=Model(optimizer);set_silent(dualmodel)
    @variable(dualmodel,t)
    @objective(dualmodel,Max,t)
    accum=[AffExpr(0.0) for _ in y]
    variable_index=Dict(v=>k for (k,v) in enumerate(y))
    for L in realblocks
        check(budget)
        Z=@variable(dualmodel,[1:size(L,1),1:size(L,1)],PSD)
        for j in axes(L,2), i in 1:j
            expression=L[i,j]
            @assert iszero(constant(expression))
            for (coefficient,variable) in linear_terms(expression)
                add_to_expression!(accum[variable_index[variable]],
                    (i==j ? 1 : 2)*coefficient,Z[i,j])
            end
        end
    end
    identity_index,identity_sign=coords.coordinates[one(first(rows))]
    equations=[@constraint(dualmodel,
        accum[k]+(k==identity_index ? identity_sign*t : 0)==
        coefficient(objective_function(model),y[k])) for k in eachindex(y)]
    (;model=dualmodel,sizes,coordinates=coords,moment_values=()->-dual.(equations))
end
end

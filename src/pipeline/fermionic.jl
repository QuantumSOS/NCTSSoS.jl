"""
    FermionicConditions(modes; conditions=(:P, :Q, :G))

Add fermionic positivity conditions on a nonempty, sorted set of modes.
Supported conditions are `:P`, `:Q`, `:G`, `:T1`, `:T2`, and `:T2prime`.
P and Q use annihilator and creator pairs; G uses particle-hole rows and
the identity row. T1 uses annihilator triples. T2 uses creator-annihilator-
annihilator triples, including repeated creator/annihilator modes. T2prime
augments T2 with annihilator rows; only its cubic/cubic entries use the
anticommutator. It is not the same condition as T2.

All these conditions are valid for each state and take no particle-number
value. Use `Sector(N, q)` for fixed-particle-number identities. Block
generators close these mode families; moment-only maps add relations without
adding blocks. The selected modes need not be the whole problem registry.
"""
struct FermionicConditions <: StrengtheningKind
    modes::Vector{Int}
    conditions::Vector{Symbol}
    function FermionicConditions(modes;conditions=(:P,:Q,:G))
        ss=_strengthening_sites(modes)
        cs=collect(conditions)
        !isempty(cs) && allunique(cs) && all(c->c in (:P,:Q,:G,:T1,:T2,:T2prime),cs) ||
            throw(ArgumentError("FermionicConditions needs distinct conditions from :P, :Q, :G, :T1, :T2, :T2prime."))
        new(ss,cs)
    end
end
_validity(::FermionicConditions) = "valid for each state"
# Cubic anticommutators temporarily contain six letters before cancellation.
_strengthening_width(k::FermionicConditions,problem,order) = any(c->c in (:T1,:T2,:T2prime),k.conditions) ? 6 : 4

function _fermionic_rows(W,modes,condition)
    buffer=PolynomialBuffer{W,Int}(1)
    function word(ls)
        result=one(W)
        for l in ls
            empty!(buffer)
            # These tuples are already normal ordered, so each product has
            # exactly one term with coefficient one, including repeated modes.
            mul!(buffer,result,fermi_letter(W,l))
            result=buffer.words[1]
        end
        result
    end
    rows=if condition==:P
        W[word((p,q)) for p in modes for q in modes if p<q]
    elseif condition==:Q
        W[word((-q,-p)) for p in modes for q in modes if p<q]
    elseif condition==:G
        W[one(W);[word((-p,q)) for p in modes for q in modes]]
    elseif condition==:T1
        W[word((p,q,r)) for p in modes for q in modes for r in modes if p<q<r]
    else
        triples=W[word((-r,p,q)) for p in modes for q in modes for r in modes if p<q]
        condition==:T2prime && append!(triples,[word((p,)) for p in modes])
        triples
    end
    sort!(rows)
end

# Accumulate exact normal words BEFORE registering orbit moments. Otherwise
# the cancelling degree-six terms would leave unused 3RDM variables behind.
struct FermionicProduct{W}
    condition::Symbol
    scratch::SectorTerms{W,Complex{Int}}
end
function (product::FermionicProduct)(sink,a,b,polynomial,coefficient=1)
    product.condition in (:P,:Q,:G) && return _block_product!(sink,a,b,polynomial,coefficient)
    empty!(product.scratch.terms)
    phase,left=star(a)
    c=scale(phase,one(Complex{Int}))
    mul!(product.scratch,left,b,c)
    if product.condition!=:T2prime || degree(a)==degree(b)==3
        mul!(product.scratch,b,left,c)
    end
    for (w,d) in product.scratch.terms
        degree(w)<=4 || error("Fermionic cubic anticommutator did not cancel degree-six terms.")
        push_term!(sink,coefficient*d,w)
    end
    sink
end

function _fermionic_strengthen!(program,table,problem,kinds,closure)
    any(k->k isa FermionicConditions,kinds) || return nothing
    registry=problem.objective.registry
    registry.algebra isa Fermionic || throw(ArgumentError("FermionicConditions requires fermionic variables."))
    W=eltype(table.representatives.keys)
    representatives=RowBlock{W,ComplexF64}[]
    initial_count=closed_count=0
    for (source,kind) in enumerate(kinds)
        kind isa FermionicConditions || continue
        maximum(kind.modes)<=registry.nsites || throw(ArgumentError("A strengthening mode is outside the registry; select existing modes."))
        for condition in kind.conditions
            rows=_fermionic_rows(W,kind.modes,condition)
            isempty(rows) && continue
            block=RowBlock(rows,[(1.0+0im,one(W))],HPSD,source)
            family,index=close_family([block],closure.block_generators)
            initial_count+=1; closed_count+=length(family)
            product=FermionicProduct(condition,SectorTerms(Dict{W,Complex{Int}}()))
            for (root,_,stabilizers,transport) in _family_orbits(family,index,closure.block_generators)
                representative=family[root]
                sectors,real_blocks=_family_sectors(representative,stabilizers,transport,index,closure.moment_only)
                _append_family_sectors!(program,table,representative,sectors;real_blocks,product! = product)
                push!(representatives,representative)
            end
        end
    end
    (initial_cones=initial_count,closed_cones=closed_count,representatives=representatives)
end

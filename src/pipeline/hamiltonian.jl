# For a ground-space density operator, this form is
# tr(rho * a† * (H-E0) * b). It remains positive for every mixture in
# that space. No full row matrix is constructed before projection.
function _optimality_product!(sink,a,b,hamiltonian,coefficient=1)
    p,left=star(a)
    q,product=mul(left,b)
    for (c,h) in hamiltonian
        r,v=mul(left,h); s,u=mul(v,b)
        push_term!(sink,scale(p*r*s,coefficient*c),u)
        r,u=mul(h,product)
        push_term!(sink,scale(p*q*r,-coefficient*c/2),u)
        r,u=mul(product,h)
        push_term!(sink,scale(p*q*r,-coefficient*c/2),u)
    end
end

function _check_hamiltonian_strengthening(problem,kinds)
    any(k->k isa Union{Stationarity,StateOptimality},kinds) || return
    problem.objective.registry.algebra isa Pauli ||
        throw(ArgumentError("Stationarity and StateOptimality require Pauli variables in this delivery."))
    isempty(problem.constraints) || throw(ArgumentError(
        "Stationarity and StateOptimality require an unconstrained problem: a constraint restricts the state set, so its minimizer need not be an eigenstate of the objective. Remove these strengthening kinds."))
end

function _stationarity!(program,table,hamiltonian,closure,source)
    W=eltype(table.representatives.keys)
    words=Set{W}()
    for block in closure.family,a in block.rows,b in block.rows
        push!(words,last(mul(a,b)))
    end
    # An antiunitary changes the sign of i[H,u]. A zero moment of u
    # therefore does NOT imply a zero commutator (H=Z, u=Y, conjugation).
    unitary=filter(g->!is_antiunitary(g),(closure.block_generators...,closure.moment_only...))
    orbits=OrbitMoments(W,unitary)
    seen=Set{Int32}(); operators=W[]
    sink=OrbitFormSink(table,Dict{Int32,ComplexF64}())
    for u in sort!(collect(words))
        coordinate=orbit_coordinates!(orbits,u)
        coordinate.re==0 && continue
        coordinate.re in seen && continue
        push!(seen,coordinate.re)
        empty!(sink.terms)
        for (c,h) in hamiltonian
            p,w=mul(h,u); push_term!(sink,scale(p,-im*c),w)
            p,w=mul(u,h); push_term!(sink,scale(p,im*c),w)
        end
        all(iszero,values(sink.terms)) && continue
        _close_orbit_form!(program.forms,sink.terms;diagonal=true)
        push!(operators,u)
    end
    isempty(operators) || push_block!(program,ZERO,length(operators),source)
    operators
end

function _hamiltonian_strengthen!(program,table,problem,kinds,closure)
    any(k->k isa Union{Stationarity,StateOptimality},kinds) || return nothing
    W=eltype(table.representatives.keys)
    hamiltonian=[(ComplexF64(c),_pack_word(W,w)) for (c,w) in terms(problem.objective)]
    operators=W[]; optimality_orbits=Vector{Int}[]
    representatives=RowBlock{W,ComplexF64}[]
    index=Dict(_block_key(b)=>i for (i,b) in enumerate(closure.family))
    for (source,kind) in enumerate(kinds)
        if kind isa Stationarity
            append!(operators,_stationarity!(program,table,hamiltonian,closure,source))
        elseif kind isa StateOptimality
            for (root,orbit,stabilizers,transport) in _family_orbits(closure.family,index,closure.block_generators)
                block=closure.family[root]
                sectors,real_blocks=_family_sectors(block,stabilizers,transport,index,closure.moment_only)
                extra=RowBlock(block.rows,hamiltonian,HPSD,source)
                _append_family_sectors!(program,table,extra,sectors;real_blocks,product! = _optimality_product!)
                push!(optimality_orbits,orbit)
                push!(representatives,extra)
            end
        end
    end
    (stationarity_operators=operators,optimality_orbits=optimality_orbits,representatives=representatives)
end

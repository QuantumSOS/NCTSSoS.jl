import Graphs

"Correlative graph on Pauli sites, fermionic modes, or individual site-word letters; constraints form cliques."
struct Correlative{A}
    algorithm::A
end
Correlative(;algorithm=MF()) = Correlative(algorithm)

"""
Support-graph splitting, with an explicit number of support-extension steps.

`level` counts support-extension steps. A chordal algorithm such as `MMD()` or
`MF()` depends on the row order of `NCTSSoS.Next`. The blocks and the bound at a
given level can differ from the old package at the same level. Each level is a
valid relaxation, and the dense program implies it.
"""
struct TermSplit{A}
    algorithm::A
    level::Int
    function TermSplit(;algorithm=MaximalElimination(),level=1)
        level isa Integer && level>=1 || throw(ArgumentError("TermSplit level must be a positive integer."))
        new{typeof(algorithm)}(algorithm,Int(level))
    end
end

"Nested geometric regions or nested capped row blocks; high levels are dense."
struct PhysicsInformed
    level::Int
    function PhysicsInformed(;level=1)
        level isa Integer && level>=0 || throw(ArgumentError("PhysicsInformed level must be a nonnegative integer."))
        new(Int(level))
    end
end
next_level(p::PhysicsInformed) = PhysicsInformed(level=p.level+1)
next_level(p::TermSplit) = TermSplit(algorithm=p.algorithm,level=p.level+1)

"Rows with contiguous occupied sites in each group's site order."
struct Contiguous
    periodic::Bool
end
Contiguous(;periodic=false) = Contiguous(periodic)

function _interaction_graph(problem,n=problem.objective.registry.nsites,support_of=w->collect(sites(w)))
    graph=Graphs.SimpleGraph(n)
    supports=[sort!(unique!(support_of(w))) for (_,w) in terms(problem.objective)]
    for c in problem.constraints
        if c isa ExpectationConstraint
            append!(supports,[sort!(unique!(support_of(w))) for (_,w) in terms(c.polynomial)])
        else
            push!(supports,sort!(unique!(reduce(vcat,[support_of(w) for (_,w) in terms(c.polynomial)];init=Int[]))))
        end
    end
    for support in supports,i in eachindex(support),j in (i+1):length(support)
        Graphs.add_edge!(graph,support[i],support[j])
    end
    graph,supports
end

_site_groups(problem,::Nothing) = [collect(1:problem.objective.registry.nsites)]
_site_groups(problem,method,::Type{W}) where {W} = _site_groups(problem,method)
function _site_groups(problem,method::Correlative)
    graph,_=_interaction_graph(problem)
    sort!.(clique_decomp(graph,method.algorithm))
end
function _site_groups(problem,method::Correlative,::Type{W}) where {W<:SiteWord}
    registry=problem.objective.registry
    # Site-major numbering matches the current package's letter graph.
    letters=W[w for site in 1:registry.nsites for w in basis_rows(W,registry.algebra,[site],1) if !isone(w)]
    index=Dict(w=>i for (i,w) in enumerate(letters))
    support(w)=[index[v] for v in word_letters(_pack_word(W,w))]
    graph,_=_interaction_graph(problem,length(letters),support)
    [letters[sort!(clique)] for clique in clique_decomp(graph,method.algorithm)]
end
function _site_groups(problem,method::PhysicsInformed)
    graph,supports=_interaction_graph(problem)
    n=Graphs.nv(graph)
    method.level>=n && return [collect(1:n)]
    regions=[sort!(Graphs.neighborhood(graph,i,method.level)) for i in 1:n]
    append!(regions,filter(!isempty,supports))
    unique(regions)
end

_group_rows(::Nothing,W,registry,group,order) = basis_rows(W,registry.algebra,group,order)
_group_sites(group) = group
_group_sites(group::AbstractVector{<:SiteWord}) = sort!(unique!(reduce(vcat,collect.(sites.(group));init=Int[])))
_group_region(group,region) = region
_group_region(group::AbstractVector{<:SiteWord},region) = filter(w->only(sites(w)) in region,group)
_in_group(w,group) = issubset(sites(w),group)
_in_group(w::SiteWord,group::AbstractVector{<:SiteWord}) = all(in(group),word_letters(w))
function _group_rows(rule::Contiguous,W,registry,group,order)
    rows=Set([one(W)])
    physical_sites=_group_sites(group)
    n=length(physical_sites)
    for width in 1:min(order,n),start in 1:(rule.periodic ? n : n-width+1)
        region=[physical_sites[mod1(start+i,n)] for i in 0:width-1]
        for w in basis_rows(W,registry.algebra,_group_region(group,region),order)
            length(collect(sites(w)))==width && push!(rows,w)
        end
    end
    sort!(collect(rows))
end
function _group_rows(rule::Function,W,registry,group,order)
    rows=collect(rule(W,registry,copy(group),order))
    all(w->w isa W && degree(w)<=order && _in_group(w,group),rows) ||
        throw(ArgumentError("The row rule must return words of the requested type, on its group, and of degree at most the order."))
    allunique(rows) || throw(ArgumentError("The row rule returned duplicate rows."))
    sort!(rows)
end

struct WordSupportSink{W}
    support::Set{W}
end
function push_term!(sink::WordSupportSink,c::Number,w)
    iszero(c) || push!(sink.support,first(moment_key(w)))
    sink
end

function _close_starting_rows(rows::Vector{W},symmetry) where {W}
    isempty(symmetry.generators) && isempty(symmetry.continuous) && return rows
    generators=map(g->_group_action(W,g),symmetry.generators)
    seen=Set(rows)
    queue=copy(rows)
    cursor=1
    while cursor<=length(queue)
        w=queue[cursor]
        for g in generators
            image=act(g,w)[2]
            image in seen && continue
            push!(seen,image); push!(queue,image)
        end
        for d in symmetry.continuous,image in keys(_derivative(d,w))
            image in seen && continue
            push!(seen,image); push!(queue,image)
        end
        cursor+=1
    end
    sort!(queue)
end

struct MissingSupportSink{W}
    words::Set{W}
end
function push_term!(sink::MissingSupportSink,c::Number,w)
    iszero(c) || delete!(sink.words,first(moment_key(w)))
    sink
end
function _check_starting_coverage(block,objective,registry)
    W=eltype(block.rows)
    sink=MissingSupportSink(Set{W}(first(moment_key(w)) for (c,w) in objective if !iszero(c)))
    # Only track required words; never construct the unsplit moment matrix.
    for j in eachindex(block.rows),i in j:length(block.rows)
        isempty(sink.words) && return nothing
        _block_product!(sink,block.rows[i],block.rows[j],block.polynomial)
    end
    isempty(sink.words) || throw(ArgumentError(
        "Objective words $(join((_word_string(registry,w) for w in sort!(collect(sink.words))), ", ")) are not covered by the starting rows; add rows whose adjoint products cover these words before using TermSplit."))
    nothing
end

function _initial_support(blocks,objective)
    W=eltype(first(blocks).rows)
    sink=WordSupportSink(Set{W}())
    for (c,w) in objective; push_term!(sink,c,w); end
    for block in blocks
        for (c,w) in block.polynomial; push_term!(sink,c,w); end
        block.source==0 || continue
        for w in block.rows; _block_product!(sink,w,w,block.polynomial); end
    end
    sink.support
end

function _support_graph(block,active)
    graph=Graphs.SimpleGraph(length(block.rows))
    sink=WordSupportSink(empty(active))
    for j in eachindex(block.rows),i in (j+1):length(block.rows)
        empty!(sink.support)
        _block_product!(sink,block.rows[i],block.rows[j],block.polynomial)
        !isdisjoint(sink.support,active) && Graphs.add_edge!(graph,i,j)
    end
    graph
end

_split_family(blocks,objective,::Nothing) = blocks
function _split_family(blocks,objective,method::TermSplit)
    active=_initial_support(blocks,objective)
    selected=empty(blocks)
    for _ in 1:method.level
        empty!(selected)
        next=WordSupportSink(copy(active))
        for block in blocks
            graph=_support_graph(block,active)
            for clique in clique_decomp(graph,method.algorithm)
                rows=sort!(block.rows[clique])
                push!(selected,RowBlock(rows,block.polynomial,block.kind,block.source))
                for j in eachindex(rows),i in j:length(rows)
                    _block_product!(next,rows[i],rows[j],block.polynomial)
                end
            end
        end
        active=next.support
    end
    unique(_block_key,selected)
end

function _split_family(blocks,objective,method::PhysicsInformed)
    active=_initial_support(blocks,objective)
    selected=empty(blocks)
    for block in blocks
        n=length(block.rows)
        cap=method.level>=ceil(Int,log2(max(1,n))) ? n : min(n,2^method.level)
        if cap==n
            push!(selected,block)
            continue
        end
        for start in 1:cap:n
            push!(selected,RowBlock(block.rows[start:min(n,start+cap-1)],block.polynomial,block.kind,block.source))
        end
        # Keep each input-supported pair at every level. The capped chunks
        # nest under doubling, and the final level contains the dense block.
        graph=_support_graph(block,active)
        for edge in Graphs.edges(graph)
            rows=sort!(block.rows[[Graphs.src(edge),Graphs.dst(edge)]])
            push!(selected,RowBlock(rows,block.polynomial,block.kind,block.source))
        end
    end
    unique(_block_key,selected)
end

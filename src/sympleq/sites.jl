# Find one verified witness to each new orbit point under a stabilizer of
# the preceding sites. Their union generates the verified subgroup if the
# search completes: the pointwise stabilizer of the entire order is trivial.
function _pair_generators(evidence, pair_labels, verified; max_trials, max_generators=128)
    graph = evidence.graph
    n = Graphs.nv(graph)
    order = _automorphism_order(evidence)
    candidates = _automorphism_candidates(evidence)
    distances = reduce(hcat,(Graphs.gdistances(graph,v) for v in 1:n))
    images = zeros(Int,n)
    used = falses(n)
    generators = Vector{Int}[]
    trials = 0
    consistent(v,w,d) = all(1:d-1) do k
        u = order[k]; t = images[u]
        distances[u,v] == distances[t,w] && pair_labels[u,v] == pair_labels[t,w]
    end

    function witness(d, neighbor_move, stop)
        d > n && return verified(images)
        v = order[d]
        for w in candidates[v]
            used[w] && continue
            neighbor_move && distances[v,w] != 1 && continue
            trials >= stop && return false
            trials += 1
            consistent(v,w,d) || continue
            images[v] = w; used[w] = true
            witness(d+1,neighbor_move,stop) && return true
            used[w] = false
        end
        false
    end

    for (d,v) in enumerate(order)
        first_generator = length(generators)+1
        orbit = Set([v])
        for w in sort(candidates[v];by=w -> distances[v,w],alg=Base.Sort.MergeSort)
            (trials >= max_trials || length(generators) >= max_generators) && break
            (w in orbit || used[w]) && continue
            trials += 1
            consistent(v,w,d) || continue
            images[v] = w; used[w] = true
            # Give a seed edge a short search for a witness that moves every
            # site along an edge; reserve the remaining budget for any map.
            if (d == 1 && distances[v,w] == 1 && witness(d+1,true,min(max_trials,trials+8n))) ||
                witness(d+1,false,max_trials)
                push!(generators,copy(images))
                for k in d+1:n
                    used[images[order[k]]] = false
                end
                # Close this orbit under the witnesses that fix this prefix.
                queue = collect(orbit)
                for u in queue, k in first_generator:length(generators)
                    t = generators[k][u]
                    t in orbit && continue
                    push!(orbit,t); push!(queue,t)
                end
            end
            used[w] = false
        end
        images[v] = v; used[v] = true
        (trials >= max_trials || length(generators) >= max_generators) && break
    end
    @debug "Labelled pair automorphism witnesses" sites=n trials generators=length(generators)
    (trials >= max_trials || length(generators) >= max_generators) &&
        @warn "Labelled pair automorphism generator search reached its budget; returning verified partial results." trials generators=length(generators)
    generators
end

# Coefficients (including sign), axes and term incidences are graph colors.
function _site_symmetry_candidates(::Type{W},n,polynomials,verified) where {W}
    n<2 && return CliffordMap{W}[]
    # Retain the legacy probes outside the bounded pair-graph path only.
    maps = [mod1.(2:n+1,n),collect(n:-1:1)]
    terms = [(j,w,c) for (j,p) in enumerate(polynomials) for (w,c) in p if !isone(w)]
    axis(u) = all(iszero,u.x) ? 2 : all(iszero,u.z) ? 0 : 1
    spec = [(j,c,[(only(sites(u)),axis(u)) for u in word_letters(w)]) for (j,w,c) in terms]
    # Avoid factorial automorphism searches on Hamiltonians with unused sites.
    active = Set(i for (_,_,xs) in spec for (i,_) in xs)
    vertices = n + length(spec) + sum(t->length(t[3]),spec;init=0)
    trials = 2_000_000
    if length(active)==n && n<=512 && all(t->length(t[3])<=2,spec)
        # Color both endpoints and directed bonds. The latter encode the
        # actual coefficient and ordered axes, not just incident multisets.
        graph = Graphs.SimpleGraph(n)
        labels = [Tuple[] for _ in 1:n]
        bonds = Dict{Tuple{Int,Int},Vector{Tuple}}()
        for (j,c,xs) in spec
            i,a = first(xs)
            if length(xs)==1
                push!(labels[i], (j,a,-1,real(c),imag(c)))
            else
                k,b = last(xs)
                Graphs.add_edge!(graph,i,k)
                for (u,v,s,t) in ((i,k,a,b),(k,i,b,a))
                    color = (j,s,t,real(c),imag(c))
                    push!(labels[u],color)
                    push!(get!(() -> Tuple[],bonds,(u,v)),color)
                end
            end
        end
        evidence = SympleQGraph(graph,Any[Tuple(sort!(l)) for l in labels],collect(1:n),Int[],Vector{Int}[])
        colors = Dict{Tuple,Int}()
        pair_labels = zeros(Int,n,n)
        for ((i,k),terms) in bonds
            pair_labels[i,k] = get!(() -> length(colors)+1,colors,Tuple(sort!(terms)))
        end
        maps = _pair_generators(evidence,pair_labels,u -> verified(CliffordMap(W,n,SitePermutation(u)));max_trials=trials)
    elseif length(active)==n && vertices<=512
        graph = Graphs.SimpleGraph(vertices)
        colors = Any[:site for _ in 1:n]
        append!(colors,[(:term,j,c) for (j,c,_) in spec])
        vertex = length(colors)
        for (j,(_,_,xs)) in enumerate(spec), (i,a) in xs
            vertex += 1
            push!(colors,(:axis,a))
            Graphs.add_edge!(graph,n+j,vertex)
            Graphs.add_edge!(graph,i,vertex)
        end
        evidence = SympleQGraph(graph,colors,collect(1:n),Int[],Vector{Int}[])
        append!(maps,[u[1:n] for u in _enumerate_color_preserving_automorphisms(evidence;max_automorphisms=128,max_trials=trials)])
    end
    filter!(!isone,unique([CliffordMap(W,n,SitePermutation(u)) for u in maps]))
end

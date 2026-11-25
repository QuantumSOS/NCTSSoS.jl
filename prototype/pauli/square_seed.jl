# Each incident Hamiltonian word is visited once, including when both its
# sites occur in u. This also works when a bond spans packed-word limbs.
_square_incident(incidence, u) =
    (h for i in NC.sites(u) for (h,j) in incidence[i] if i < j || j ∉ NC.sites(u))

function _square_add_spin_rows!(rows::Set{W}, sites) where {W}
    for axes in Iterators.product(ntuple(_ -> 0:2, length(sites))...)
        w = one(W)
        for (i,a) in zip(sites,axes)
            w = last(NC.mul(w, NC.pauli_letter(W,i,a)))
        end
        push!(rows,w)
    end
end

function _square_motif_regions(motifs, L, site, neighbors)
    regions = Set{Tuple{Vararg{Int}}}()
    for motif in motifs
        offsets = [(x,y) for (x,y) in motif]
        isempty(offsets) && throw(ArgumentError("A square motif must be nonempty."))
        support = Set(site(x,y) for (x,y) in offsets)
        length(support) == length(offsets) || throw(ArgumentError("Motif sites must be distinct on the torus."))
        # Check radius in the induced interaction subgraph, not in the torus
        # with intermediate sites outside the motif.
        any(support) do root
            reached = Set([root])
            for _ in 1:2
                union!(reached, [j for i in reached for j in neighbors[i] if j in support])
            end
            return reached == support
        end || throw(ArgumentError("A square motif must be connected with radius at most two."))
        for swap in (false,true), sx in (-1,1), sy in (-1,1), dy in 0:L-1, dx in 0:L-1
            sites = [swap ? site(dx+sx*y,dy+sy*x) : site(dx+sx*x,dy+sy*y) for (x,y) in offsets]
            push!(regions,Tuple(sort!(sites)))
        end
    end
    return regions
end

struct _SquareMomentOrbits{W}
    representatives::Set{W}
    shifts::Matrix{Int}
    to_origin::Vector{Int}
end

Base.eltype(::Type{_SquareMomentOrbits{W}}) where {W} = W
Base.length(domain::_SquareMomentOrbits) = length(domain.representatives)
Base.iterate(domain::_SquareMomentOrbits, state...) = iterate(domain.representatives,state...)
Base.in(w::W, domain::_SquareMomentOrbits{W}) where {W} =
    all(<=(length(domain.to_origin)),NC.sites(w)) && _square_representative(w,domain) in domain.representatives
Base.in(::Any, ::_SquareMomentOrbits) = false

function _square_translate(w::NC.PauliWord{K}, shifts, t) where {K}
    x = z = ntuple(_ -> UInt64(0),Val(K))
    for i in NC.sites(w)
        j = @inbounds shifts[i,t]
        a,b = ((i-1) >>> 6)+1, ((j-1) >>> 6)+1
        inbit, outbit = UInt64(1) << ((i-1)&63), UInt64(1) << ((j-1)&63)
        @inbounds begin
            iszero(w.x[a] & inbit) || (x = Base.setindex(x,x[b] | outbit,b))
            iszero(w.z[a] & inbit) || (z = Base.setindex(z,z[b] | outbit,b))
        end
    end
    return NC.PauliWord(x,z)
end

function _square_representative(w, domain::_SquareMomentOrbits)
    # Every nonidentity orbit has an image containing site 1. It suffices to
    # compare the shifts bringing one occupied site to 1, not all L² shifts.
    best = w
    first = true
    for i in NC.sites(w)
        candidate = _square_translate(w,domain.shifts,domain.to_origin[i])
        if first || (candidate.x,candidate.z) < (best.x,best.z)
            best = candidate
            first = false
        end
    end
    return best
end

function _square_product_orbits(rows::Vector{W}, L) where {W}
    xy = [(x,y) for y in 0:L-1 for x in 0:L-1]
    shifts = [mod(x+dx,L)+1 + L*mod(y+dy,L) for (x,y) in xy, (dx,dy) in xy]
    to_origin = [mod(-x,L)+1 + L*mod(-y,L) for (x,y) in xy]
    domain = _SquareMomentOrbits(Set{W}(),shifts,to_origin)
    reps = sort!(collect(Set(_square_representative(w,domain) for w in rows)))
    shifted = [_square_translate(w,shifts,t) for t in eachindex(xy), w in reps]
    # Translation removes the common displacement in every pair. Swapping
    # its factors changes only the discarded phase and inverts displacement.
    for j in eachindex(reps), i in 1:j, t in eachindex(xy)
        w = last(NC.mul(reps[i],shifted[t,j]))
        push!(domain.representatives,_square_representative(w,domain))
    end
    return domain
end

function _square_stationarity(H::NC.Polynomial{W}, J, domain, incidence) where {W}
    sources = W[]
    equations = typeof(H)[]
    terms = Tuple{NC.Phase,W}[]
    # If v occurs in [H,u], then u is the phase-free product of v and some
    # anticommuting Hamiltonian word h. Distinct h give distinct outputs for
    # fixed u, so a rejected term cannot cancel. Enumerate each accepted u
    # from its least output v only; this avoids storing rejected candidates.
    for v in domain, h in _square_incident(incidence,v)
        NC._symplectic(h,v) || continue
        u = last(NC.mul(h,v))
        empty!(terms)
        for g in _square_incident(incidence,u)
            NC._symplectic(g,u) || continue
            phase,w = NC.mul(g,u)
            if w < v || w ∉ domain
                empty!(terms)
                break
            end
            push!(terms,(phase,w))
        end
        isempty(terms) && continue
        push!(sources,u)
        push!(equations,_wf_polynomial(H.registry,W,((NC.scale(p,2im*J),w) for (p,w) in terms)))
    end
    perm = sortperm(sources)
    return sources[perm], equations[perm]
end

function _square_other_axes(w::NC.PauliWord)
    x,z = w.x,w.z
    y = map(xor,x,z)
    return (NC.PauliWord(z,x),NC.PauliWord(x,y),NC.PauliWord(y,x),
        NC.PauliWord(y,z),NC.PauliWord(z,y))
end

function _square_stationarity(H::NC.Polynomial{W}, J, domain::_SquareMomentOrbits{W}, incidence) where {W}
    sources = Set{W}()
    seen = Set{W}()
    # Global axis permutations commute with translations, preserve the
    # anticommutation test and phase-free products, and permute the XXX
    # bonds and the domain (every row support has all axis assignments).
    # Every valid source orbit can therefore invert its least output over
    # both actions. Restore all axis images before emitting exact phases.
    for v in sort!(collect(domain.representatives))
        v in seen && continue
        for w in _square_other_axes(v)
            push!(seen,_square_representative(w,domain))
        end
        for h in _square_incident(incidence,v)
            NC._symplectic(h,v) || continue
            u = last(NC.mul(h,v))
            covered = true
            for g in _square_incident(incidence,u)
                NC._symplectic(g,u) || continue
                rep = _square_representative(last(NC.mul(g,u)),domain)
                if rep < v || rep ∉ domain.representatives
                    covered = false
                    break
                end
            end
            covered || continue
            rep = _square_representative(u,domain)
            rep in sources && continue
            push!(sources,rep)
            for w in _square_other_axes(rep)
                push!(sources,_square_representative(w,domain))
            end
        end
    end
    words = sort!(collect(sources))
    equations = [_wf_polynomial(H.registry,W,
        ((NC.scale(p,2im*J),w) for g in _square_incident(incidence,u) if NC._symplectic(g,u)
            for (p,w) in (NC.mul(g,u),))) for u in words]
    return words, equations
end

"""
    build_square_xxx_seed(H, L; motifs=(), rdm=true, rdm_shape=(2,2), stationarity=true,
                         stability=false, translation_quotient=false)

Build phases 1–2 for a periodic `L × L` square with site order `x + L*y + 1`
for zero-based coordinates. `H` must have uniform nonzero XXX nearest-neighbor
coupling (either sign, with an optional scalar shift). For `L=2` each distinct
bond is counted once in the geometry and each distinct face once in the RDMs.

Return sorted packed `rows` for one ordinary Gram block, the shared moment
`domain`, `plaquettes`, and `stationarity_words` in equality order. Append
`strengthening` to the **same** moment table as the Gram block. `plaquettes`
lists the actual sites in each RDM block's tensor order. Use `rdm_shape=(2,3)`
for 64-dimensional rectangles instead of the default 16-dimensional squares;
the rectangles imply the embedded square marginals by partial trace. Neither
adds Gram rows. Optional scalar `stability` is valid for unrestricted global
ground states.

Use `translation_quotient=true` **only** with a shared moment table that
identifies both primitive translations of this torus. `rows` still contains
the full Gram basis. `domain` iterates translation representatives (and its
length counts orbits), but `w in domain` accepts any translated member. The
builder returns one RDM cone per orientation (one square or two rectangles)
and one scalar constraint per translation orbit. The domain covers the full
D4 orbit without requiring the table to identify point-group or spin orbits.

Each additional `motifs` entry is a connected site support written as integer
`(dx,dy)` offsets, e.g. `(((0,0),(1,0),(0,1)),)` for one corner motif. The
builder adds its full translation and D4 orbit and all spin assignments,
then recomputes the moment domain and stationarity. It does not apply a
second Liouvillian shell to the added motifs.
"""
function build_square_xxx_seed(H::NC.Polynomial{W}, L::Integer;
        motifs=(), rdm::Bool=true, rdm_shape::Tuple{Integer,Integer}=(2,2),
        stationarity::Bool=true, stability::Bool=false,
        translation_quotient::Bool=false) where {W<:NC.PauliWord}
    n = H.registry.nsites
    L >= 2 && n == L^2 || throw(ArgumentError("Expected L ≥ 2 and exactly L² Pauli sites."))
    rdm_shape in ((2,2),(2,3)) || throw(ArgumentError("Expected rdm_shape=(2,2) or (2,3)."))
    (!rdm || L >= maximum(rdm_shape)) || throw(ArgumentError("The RDM region must fit on the torus without repeated sites."))
    NC._checked_hermitian(H,"Hamiltonian")
    HQ = _wf_exact(H)
    site(x,y) = mod(x,L) + L*mod(y,L) + 1
    bonds = Set(minmax(site(x,y),site(x+dx,y+dy))
        for x in 0:L-1 for y in 0:L-1 for (dx,dy) in ((1,0),(0,1)))
    letters = [NC.pauli_letter(W,i,a) for i in 1:n, a in 0:2]
    neighbors = [Int[] for _ in 1:n]
    incidence = [Tuple{W,Int}[] for _ in 1:n]
    expected = Set{W}()
    for (i,j) in bonds
        push!(neighbors[i],j)
        push!(neighbors[j],i)
        for a in axes(letters,2)
            h = last(NC.mul(letters[i,a],letters[j,a]))
            push!(expected,h)
            push!(incidence[i],(h,j))
            push!(incidence[j],(h,i))
        end
    end
    actual = Dict(w => real(c) for (c,w) in NC.terms(HQ) if !isone(w))
    keys(actual) == expected || throw(ArgumentError("Hamiltonian terms must be exactly the periodic square XXX bonds."))
    J = first(values(actual))
    all(==(J),values(actual)) || throw(ArgumentError("Periodic square XXX coupling must be uniform."))

    # Minimal factors I,X_i,Y_i,Z_i already cover every two-body term. Take
    # exactly one shell and close each of its supports under spin rotations:
    # individual Pauli words, unlike the chiral polynomial combinations,
    # need all nine axis pairs on a bond for the monomial span to be invariant.
    support = Set{Tuple{Vararg{Int}}}(Tuple(NC.sites(last(NC.mul(h,u))))
        for u in letters for h in _square_incident(incidence,u) if NC._symplectic(h,u))
    union!(support,_square_motif_regions(motifs,L,site,neighbors))
    selected = Set{W}(letters)
    push!(selected,one(W))
    for sites in support
        _square_add_spin_rows!(selected,sites)
    end
    rows = sort!(collect(selected))
    plaquettes = if !rdm
        Vector{Int}[]
    else
        orientations = rdm_shape == (2,2) ? (rdm_shape,) : (rdm_shape,reverse(rdm_shape))
        origins = 0:(translation_quotient ? 0 : L-1)
        sort!(unique!([sort!([site(x+dx,y+dy) for dy in 0:ny-1 for dx in 0:nx-1])
            for (nx,ny) in orientations for y in origins for x in origins]))
    end
    rdm_words = (w for sites in plaquettes for w in NC.basis_rows(W,NC.Pauli(),sites,length(sites)))
    domain = if translation_quotient
        orbits = _square_product_orbits(rows,L)
        union!(orbits.representatives,(_square_representative(w,orbits) for w in rdm_words))
        orbits
    else
        words = Set{W}(last(NC.mul(rows[i],rows[j])) for j in eachindex(rows) for i in 1:j)
        union!(words,rdm_words)
    end
    strengthening = workflow_strengthening(HQ,2; rdm_regions=plaquettes, energy_target=:global_ground)
    stationarity_words = W[]
    if stationarity
        stationarity_words, equations = _square_stationarity(HQ,J,domain,incidence)
        append!(strengthening.equalities,equations)
        push!(strengthening.rules,(; kind=:stationarity, validity=:target_eigenstates))
    end
    if stability
        for u in (translation_quotient ? view(letters,1,:) : letters)
            # u[H,u] = uHu-H; only anticommuting terms survive and the
            # result uses Hamiltonian words already covered by Gram products.
            push!(strengthening.inequalities,_wf_polynomial(H.registry,W,
                ((-2J,h) for h in _square_incident(incidence,u) if NC._symplectic(h,u))))
        end
        push!(strengthening.rules,(; kind=:one_site_stability, validity=:global_ground_states))
    end
    strengthening = merge(strengthening,(; hamiltonian_dependent=stationarity || stability))
    return (; rows,domain,plaquettes,stationarity_words,strengthening)
end

# Included by the Pauli workflow module after `import NCTSSoS as NC`.

include("square_moments.jl")

_symmetry_action(::Type{W}, n, u) where {W} = NC.CliffordMap(W,n,NC.SitePermutation(u))

function _symmetry_preserves(g, terms; sector=false)
    image = empty(terms)
    for (w,c) in terms
        p,v = NC.act(g,w)
        image[v] = NC.scale(p,NC.is_antiunitary(g) ? conj(c) : c)
    end
    image == terms || (sector && keys(image)==keys(terms) && all(image[w]==-c for (w,c) in terms))
end

function _symmetry_family_preserves(g, family)
    isnothing(family) && return true
    all(family) do block
        rows = block isa NC.RowBlock ? block.rows : block
        image = Set(last(NC.act(g,w)) for w in rows)
        any(family) do candidate
            candidate isa NC.RowBlock ? candidate.rows⊆image && length(candidate.rows)==length(image) &&
                candidate.source==block.source && candidate.kind==block.kind :
                candidate⊆image && length(candidate)==length(image)
        end
    end
end

# The image is stored on *all* commutator outputs, not projected onto the
# one-site generator space (which would miss the entire two-body obstruction).
struct WorkflowLocalRotation{W} <: NC.AbstractDerivation
    columns::Vector{Tuple{W,Rational{BigInt}}}
end
function NC.derive!(sink, d::WorkflowLocalRotation{W}, w::W, c::Number=1) where {W}
    for (u,a) in d.columns
        NC._symplectic(u,w) || continue
        p,v = NC.mul(u,w)
        NC.push_term!(sink,NC.scale(NC.Phase(3)*p,c*a),v)
    end
    sink
end

function _symmetry_derivative(d, terms::Dict{W}) where {W}
    sink = NC.SectorTerms(W)
    for (w,c) in terms
        NC.derive!(sink,d,w,c)
    end
    sink.terms
end

function _workflow_lie(registry, letters::Vector{W}, polynomials, family) where {W}
    closed(d) = isnothing(family) || all(family) do block
        rows = Set(block isa NC.RowBlock ? block.rows : block)
        all(w->all(in(rows),keys(NC._derivative(d,w))),rows)
    end
    valid(d) = all(p->isempty(_symmetry_derivative(d,p)),polynomials) && closed(d)
    collective = filter(valid,collect(NC.spin_rotations(registry)))
    # Keep a collective factor when available; its known representation is
    # what enables the analytical path. Two generators imply the third.
    length(collective) in (1,3) && return collective

    R = Rational{BigInt}
    pivots = Dict{Int,Dict{Int,R}}()
    units = [WorkflowLocalRotation([(w,one(R))]) for w in letters]
    for p in polynomials
        equations = Dict{W,Dict{Int,R}}()
        for (j,d) in enumerate(units), (w,c) in _symmetry_derivative(d,p)
            get!(() -> Dict{Int,R}(),equations,w)[j] = real(c)
        end
        for w in sort!(collect(keys(equations)))
            NC._insert_generator_equation!(pivots,equations[w];constant=false)
        end
    end
    if !isnothing(family)
        for block in family
            rows = Set(block isa NC.RowBlock ? block.rows : block)
            for w in sort!(collect(rows))
                equations = Dict{W,Dict{Int,R}}()
                for (j,d) in enumerate(units), (v,c) in NC._derivative(d,w)
                    v in rows && continue
                    get!(() -> Dict{Int,R}(),equations,v)[j] = real(c)
                end
                for v in sort!(collect(keys(equations)))
                    NC._insert_generator_equation!(pivots,equations[v];constant=false)
                end
            end
        end
    end
    result = WorkflowLocalRotation{W}[]
    for j in eachindex(letters)
        haskey(pivots,j) && continue
        x = zeros(R,length(letters)); x[j] = 1
        for i in sort!(collect(keys(pivots));rev=true)
            x[i] = -sum((c*x[k] for (k,c) in pivots[i] if k!=i);init=zero(R))
        end
        d = WorkflowLocalRotation([(w,c) for (w,c) in zip(letters,x) if !iszero(c)])
        valid(d) && push!(result,d)  # independent full-output recheck
    end
    result
end

function _workflow_parities(::Type{W}, n, polynomials, rotations) where {W}
    words = sort!(unique([w for p in polynomials for w in keys(p)]))
    # A parity must give a common sign to all sites of a collective axis.
    for d in rotations
        d isa NC.SpinRotation || continue
        v = NC.pauli_letter(W,1,d.axis)
        append!(words,[last(NC.mul(v,NC.pauli_letter(W,i,d.axis))) for i in 2:n])
    end
    # Pair columns with (X_i,Z_i) conjugating letters, i.e. [z | x].
    letters = [NC.pauli_letter(W,i,a) for a in (0,2) for i in 1:n]
    A = UInt8[NC._symplectic(w,u) for w in words,u in letters]
    as_word(bits) = foldl((w,u)->last(NC.mul(w,u)),letters[findall(!iszero,bits)];init=one(W))
    parity = [NC.CliffordMap(NC.PauliConjugation(as_word(v),false)) for v in NC._gf2_nullspace_basis(A)]
    # All other diagonal antisymmetries, if any, differ by a unitary parity.
    conjugation = NC.AntiunitaryPauli(one(W))
    reversal = NC.AntiunitaryPauli(as_word(ones(UInt8,2n)))
    rhs = UInt8[first(NC.act(conjugation,w)).k==2 for w in words]
    solution = NC._gf2_affine_solution_space(A,rhs)
    anti = NC.CliffordMap{W}[NC.CliffordMap(conjugation),NC.CliffordMap(reversal)]
    isnothing(solution) || push!(anti,NC.CliffordMap(NC.AntiunitaryPauli(as_word(first(solution)))))
    parity,unique(anti)
end

function _workflow_translations(::Type{W}, n, lattice_shape, spatial_mode) where {W}
    lattice_shape isa Tuple{Integer,Integer} && all(>=(2),lattice_shape) &&
        prod(BigInt.(lattice_shape))==n || throw(ArgumentError("lattice_shape must be (Lx,Ly), both at least 2, with Lx*Ly=$n."))
    spatial_mode in (:auto,:lattice,:circulant_dft) || throw(ArgumentError("lattice_shape requires a lattice spatial mode."))
    nx,ny = Int.(lattice_shape)
    tx = [mod1(x+1,nx)+(y-1)*nx for y in 1:ny for x in 1:nx]
    ty = [x+(mod1(y+1,ny)-1)*nx for y in 1:ny for x in 1:nx]
    [_symmetry_action(W,n,tx),_symmetry_action(W,n,ty)]
end

function _workflow_axes(::Type{W}, n) where {W}
    # Odd axis permutations need a compensating sign on the remaining axis.
    # (X,Y,Z) -> (Y,Z,X) and (Y,X,-Z) are proper rotations.
    [NC.CliffordMap([NC.pauli_letter(W,i,a) for a in axes for i in 1:n],
        BitVector([b for b in negative for _ in 1:n]))
        for (axes,negative) in (((1,0),(false,false)),((1,2),(false,true)),
            ((2,0),(false,false)),((0,1),(true,false)),((2,1),(false,false)))]
end

function _workflow_normalizes(g, rotations, registry, letters::Vector{W}) where {W}
    isempty(rotations) && return true
    first(rotations) isa NC.SpinRotation || return true # full joint centralizer of the checked operators
    allowed = Set(d.axis for d in rotations)
    n = registry.nsites
    all(rotations) do d
        p,w = NC.act(g,NC.pauli_letter(W,1,d.axis))
        j = only(NC.sites(w))
        a = only(a for a in 0:2 if w==NC.pauli_letter(W,j,a))
        a in allowed && all(i->begin
            q,u = NC.act(g,NC.pauli_letter(W,i,d.axis))
            q==p && u==NC.pauli_letter(W,only(NC.sites(u)),a)
        end,1:n)
    end
end

"""
    workflow_symmetry(H; observable=nothing, target=nothing, spatial_mode=:auto,
                      family=nothing, strengthen=(), parity=true, antiunitary=true,
                      spin=true, axis_permutations=true, lattice_shape=nothing)

Discover a verified subgroup for a Pauli relaxation. `family`, when supplied,
is a collection of packed word blocks (or `NC.RowBlock`s); no new blocks are
added. `target` is `nothing` or an `NC.Sector`. `lattice_shape=(Lx,Ly)` names
the two periodic shifts directly, with site `x+(y-1)*Lx`. Both must preserve
the problem and family. Without a shape, spatial discovery is bounded.
"""
function workflow_symmetry(H; observable=nothing, target=nothing, spatial_mode=:auto,
    family=nothing, strengthen=(), parity=true, antiunitary=true, spin=true, axis_permutations=true,
    lattice_shape=nothing, uniform_parity_blocks=false)
    spatial_mode in (:auto,:none,:lattice,:circulant_dft,:clique_orbits) || throw(ArgumentError("Unknown spatial mode $spatial_mode."))
    registry = H.registry
    registry.algebra isa NC.Pauli || throw(ArgumentError("Expected a Pauli Hamiltonian."))
    NC.Problem(H)
    n = registry.nsites
    W = eltype(H.data.words)
    translations = isnothing(lattice_shape) ? nothing : _workflow_translations(W,n,lattice_shape,spatial_mode)
    (isnothing(target) || target isa NC.Sector) || throw(ArgumentError("Use NC.Sector for a state target."))
    checked = [H]
    if !isnothing(observable)
        NC._checked_hermitian(observable,"Observable")
        # Preserve the original coefficient types for exact verification.
        checked = [H,observable]
    end
    targets = NC.Sector[k for k in strengthen if k isa NC.Sector]
    isnothing(target) || push!(targets,target)
    all(p->p.registry===registry,Iterators.flatten((checked,(k.operator for k in targets)))) ||
        throw(ArgumentError("The Hamiltonian, observable and target must share a registry."))
    polynomials = [NC._sector_terms(p,W).terms for p in checked]
    operators = [NC._sector_terms(k.operator,W).terms for k in targets]
    shifted = map(targets,operators) do k,p
        s = NC.SectorTerms(copy(p))
        NC.push_term!(s,-Rational{BigInt}(k.eigenvalue),one(W))
        s.terms
    end
    if !isnothing(family)
        localizers = [NC.SectorTerms(Dict(w=>Complex{Rational{BigInt}}(c) for (c,w) in b.polynomial)).terms
            for b in family if b isa NC.RowBlock && b.source!=0]
        append!(polynomials,localizers)
    end
    allops = [polynomials;operators]
    letters = [NC.pauli_letter(W,i,a) for i in 1:n for a in 0:2]
    # Unknown strengthening subspaces cannot be assumed invariant. Invariant
    # cones select a fixed labelled subspace, not an entire local algebra.
    if !all(k->k isa Union{NC.DensityMatrix,NC.Sector,NC.Stationarity,NC.StateOptimality},strengthen)
        isnothing(translations) || throw(ArgumentError("Cannot verify the named translations for this strengthening."))
        return NC.FamilySymmetry(NC.CliffordMap{W}[];detected=true)
    end
    rotations = !spin || n==0 ? NC.AbstractDerivation[] : _workflow_lie(registry,letters,allops,family)
    function verified(g)
        all(p->_symmetry_preserves(g,p),polynomials) && all(p->_symmetry_preserves(g,p;sector=true),shifted) &&
            _symmetry_family_preserves(g,family) && _workflow_normalizes(g,rotations,registry,letters) &&
            all(k->begin
                k isa NC.DensityMatrix || return true
                image = sort!([only(NC.sites(last(NC.act(g,NC.pauli_letter(W,i,0))))) for i in k.sites])
                any(l->l isa NC.DensityMatrix && l.sites==image,strengthen)
            end,strengthen)
    end
    n==0 && return NC.FamilySymmetry(NC.CliffordMap{W}[];detected=true)
    signs,anti = parity || antiunitary ? _workflow_parities(W,n,allops,rotations) :
        (NC.CliffordMap{W}[],NC.CliffordMap{W}[])
    parity ? filter!(verified,signs) : empty!(signs)
    antiunitary ? filter!(verified,anti) : empty!(anti)
    sites = if !isnothing(translations)
        all(verified,translations) || throw(ArgumentError("A coordinate translation of lattice_shape=$lattice_shape changes the problem or selected family."))
        # The abelian Fourier path accepts antiunitaries commuting with each
        # named shift. An incompatible antiunitary is not needed for coverage.
        filter!(g->all(t->NC._conjugate_action(g,t)==t,translations),anti)
        translations
    elseif spatial_mode==:none
        NC.CliffordMap{W}[]
    else
        filter(verified,NC._site_symmetry_candidates(W,n,allops,verified))
    end
    axes = axis_permutations ? filter(verified,_workflow_axes(W,n)) : NC.CliffordMap{W}[]
    uniform = parity && uniform_parity_blocks ? filter(verified,_wang_uniform_parities(H,allops)) : NC.CliffordMap{W}[]
    filter!(g -> g ∉ uniform,signs)
    # Keep diagonal parities out of the enumerated spatial point group. They
    # remain exact relations, except for uniform Wang signs (at most Klein 4).
    NC.FamilySymmetry(unique([sites;axes;uniform]);moment_only=Tuple([signs;first(anti,min(1,length(anti)))]),
        continuous=rotations,detected=true)
end

# A local spin-1 vector is represented by the orthonormal operator letters
# Z, (X+iY)/√2, (X-iY)/√2. The charge label is [Σ Z/2, ·], not a state sector.
function _workflow_axial_rows(rows::Vector{W}, a) where {W}
    columns = Dict{Int,Vector{Vector{Tuple{ComplexF64,W}}}}()
    b,c = mod(a+1,3),mod(a+2,3)
    for support in unique(Tuple(NC.sites(w)) for w in rows)
        # Tensor products with disjoint supports have phase +1.
        local_columns(i) = ([(1.0+0im,NC.pauli_letter(W,i,a))],
            [(complex(inv(sqrt(2.0))),NC.pauli_letter(W,i,b)),(im/sqrt(2.0),NC.pauli_letter(W,i,c))],
            [(complex(inv(sqrt(2.0))),NC.pauli_letter(W,i,b)),(-im/sqrt(2.0),NC.pauli_letter(W,i,c))])
        for labels in Iterators.product((1:3 for _ in support)...)
            m = sum((label==2 ? 1 : label==3 ? -1 : 0 for label in labels);init=0)
            col = [(1.0+0im,one(W))]
            for (i,label) in zip(support,labels)
                col = [(x*y,last(NC.mul(u,v))) for (x,u) in col for (y,v) in local_columns(i)[label]]
            end
            push!(get!(() -> Vector{Tuple{ComplexF64,W}}[],columns,m),col)
        end
    end
    columns
end

function _workflow_isotropic_rows(rows::Vector{W}) where {W}
    columns = Dict(J=>Vector{Tuple{ComplexF64,W}}[] for J in 0:2)
    for support in unique(Tuple(NC.sites(w)) for w in rows)
        if isempty(support)
            push!(columns[0],[(1+0im,one(W))])
        elseif length(support)==1
            push!(columns[1],[(1+0im,NC.pauli_letter(W,only(support),2))])
        else
            i,j = support
            w(a,b) = last(NC.mul(NC.pauli_letter(W,i,a),NC.pauli_letter(W,j,b)))
            push!(columns[0],[(inv(sqrt(3.0)),w(a,a)) for a in 0:2])
            push!(columns[1],[(im/sqrt(2.0),w(0,1)),(-im/sqrt(2.0),w(1,0))])
            push!(columns[2],[(a==2 ? 2/sqrt(6.0) : -inv(sqrt(6.0)),w(a,a)) for a in 0:2])
        end
    end
    columns
end

function _workflow_full_axes(rows)
    counts = Dict{Tuple,Int}()
    for w in rows
        s = Tuple(NC.sites(w))
        counts[s] = get(counts,s,0)+1
    end
    all(kv->big(3)^length(first(kv))==last(kv),counts)
end

function _workflow_grade_rows(rows::Vector{W}, parity, generators) where {W}
    signature(w) = Tuple(first(NC.act(g,w)).k for g in parity)
    sectors = Dict{Tuple,Vector{W}}()
    for w in rows
        push!(get!(() -> W[],sectors,signature(w)),w)
    end
    result = Vector{Vector{Tuple{ComplexF64,W}}}[]
    seen = Set{Tuple}()
    for key in sort!(collect(keys(sectors)))
        key in seen && continue
        # A verified normalizer transports the *whole* character block; no
        # dimension heuristic can identify redundant anisotropic sectors.
        orbit = [key]; push!(seen,key)
        i=1
        while i<=length(orbit)
            for g in generators
                image = signature(last(NC.act(g,first(sectors[orbit[i]]))))
                image in seen && continue
                push!(seen,image); push!(orbit,image)
            end
            i+=1
        end
        push!(result,[[(1+0im,w)] for w in sectors[key]])
    end
    result
end

# A discovered diagonal antiunitary either fixes each column up to sign or
# requires the generic complex lift. Retain sparse columns in the fixed case.
function _workflow_real_rows(sectors,anti)
    fixed = BitVector()
    for sector in sectors
        signs = map(sector) do col
            a,w = first(col)
            p,_ = NC.act(anti,w)
            sign = NC.scale(p,conj(a))/a
            all(col) do (b,v)
                q,u = NC.act(anti,v)
                u==v && NC.scale(q,conj(b))==sign*b
            end ? sign : 0im
        end
        ok = all(s->s==1 || s==-1,signs)
        push!(fixed,ok)
        ok || continue
        for (i,s) in enumerate(signs)
            s == -1 && (sector[i] = [(im*c,w) for (c,w) in sector[i]])
        end
    end
    fixed
end

function _workflow_reverses_charge(g, axis_words)
    !isempty(axis_words) && all(axis_words) do q
        p,u = NC.act(g,q)
        p.k==2 && u in axis_words
    end
end

"""
    workflow_symmetry_program(family, objective, symmetry; spatial_mode=:auto,
                              realify=true, lattice_shape=nothing)

Assemble an *unstrengthened* selected family with a single orbit table. Use
the redesign builder when adding constraints or strengthenings before exact
continuous-coordinate elimination. `symmetry` must be verified against the
same problem and family; this entry point rejects an unclosed family.
"""
function workflow_symmetry_program(family::AbstractVector{<:NC.RowBlock{W}}, objective,
    symmetry; spatial_mode=:auto, realify=true, lattice_shape=nothing) where {W<:NC.PauliWord}
    spatial_mode in (:auto,:none,:lattice,:circulant_dft,:clique_orbits) || throw(ArgumentError("Unknown spatial mode $spatial_mode."))
    isempty(family) && throw(ArgumentError("Select at least one moment block."))
    all(b->b.source==0 && b.kind==NC.HPSD,family) || throw(ArgumentError("Use moment blocks only; use BlockRelaxation for localizers."))
    gs,ms,ds = symmetry.generators,symmetry.moment_only,symmetry.continuous
    allgs = (gs...,ms...)
    NC._check_objective(objective,allgs)
    all(g->_symmetry_family_preserves(g,family),allgs) || throw(ArgumentError("A finite action moves the selected family."))
    NC._continuous_family(family,ds;custom=true)
    n = maximum((i for b in family for w in b.rows for i in NC.sites(w));init=0)
    lattice = if !isnothing(lattice_shape)
        length(family)==1 || throw(ArgumentError("Named lattice translations require one closed row block."))
        ts = _workflow_translations(W,n,lattice_shape,spatial_mode)
        all(in(gs),ts) || throw(ArgumentError("Discover symmetry with the same lattice_shape before assembly."))
        translations = NC.TranslationGroup(ts,collect(lattice_shape))
        point = filter(g->!haskey(translations.index,g),gs)
        # Discovery with a shape supplies only coordinate translations and
        # global axis actions, which normalize those translations.
        isempty(point) && push!(point,one(first(ts)))
        NC.abelian_lattice_program(only(family).rows,objective,translations,point;moment_generators=ms)
    elseif spatial_mode in (:none,:clique_orbits) || length(family)!=1
        nothing
    else
        letters = [[NC.pauli_letter(W,i,a) for a in 0:2] for i in 1:n]
        NC._whole_lattice_family(only(family),objective,gs,ms,letters)
    end
    continuous_blocks = false
    if !isnothing(lattice)
        program,table,reduction = lattice
        report = (;route=:lattice,reduction,lattice_shape)
        continuous_blocks = true
    else
        spatial_mode in (:lattice,:circulant_dft) && throw(ArgumentError("No verified regular translations on the selected single block."))
        row_sites = Set(i for b in family for w in b.rows for i in NC.sites(w))
        collective = !isempty(ds) && all(d->d isa NC.SpinRotation && row_sites⊆Set(NC.MaskSites(d.support)),ds)
        # Support closure is needed even for the longitudinal components
        # disconnected by one axial generator.
        whole = all(b->_workflow_full_axes(b.rows),family)
        axial = collective && whole && length(ds)==1
        isotropic = collective && whole && length(ds)==3 &&
            all(b->maximum(NC.degree,b.rows)<=2,family)
        parity = filter(g->!NC.is_antiunitary(g),ms)
        grade = isempty(ds) && !isempty(parity) && length(family)==1 &&
            all(g->all(w->last(NC.act(g,w))==w,only(family).rows),parity)
        if axial || isotropic || grade
            blocks,index = NC.close_family(family,gs;custom=true)
            table = NC.OrbitMoments(W,allgs;annihilators=isempty(ds) ? ms : ())
            program = NC.MomentProgram{ComplexF64}()
            labels = Vector{Int}[]
            anti = findfirst(NC.is_antiunitary,ms)
            for (root,_,stabilizers,_) in NC._family_orbits(blocks,index,gs)
                b = blocks[root]
                if grade
                    sectors = _workflow_grade_rows(b.rows,parity,stabilizers)
                    push!(labels,collect(eachindex(sectors)))
                else
                    columns = axial ? _workflow_axial_rows(b.rows,only(ds).axis) : _workflow_isotropic_rows(b.rows)
                    ks = sort!(collect(keys(columns)))
                    if axial
                        # Charge is a real eigenvalue of [Σσ/2, ·]. K fixes
                        # Z charge. Require the entire axial sum on this
                        # block to change sign, including its site support.
                        axis_words = Set(NC.pauli_letter(W,i,only(ds).axis) for w in b.rows for i in NC.sites(w))
                        reverses(g) = _workflow_reverses_charge(g,axis_words)
                        pairing = any(reverses,stabilizers) || any(g->reverses(g) && _symmetry_family_preserves(g,[b]),ms)
                        pairing && filter!(>=(0),ks)
                    end
                    sectors = [columns[k] for k in ks if !isempty(columns[k])]
                    push!(labels,[k for k in ks if !isempty(columns[k])])
                end
                real_blocks = false
                if !isnothing(anti)
                    g = ms[anti]
                    real_blocks = _workflow_real_rows(sectors,g)
                end
                NC._append_family_sectors!(program,table,b,sectors;real_blocks)
            end
            program,table = NC._orbit_objective!(program,table,objective)
            report = (;route=axial ? :axial_charge : isotropic ? :spin_j_m0 : :pauli_grading,labels)
        else
            program,table,report = NC.family_program(family,objective;
                generators=gs,moment_generators=ms,custom=true,check_coverage=false)
            continuous_blocks = true
        end
    end
    if !isempty(ds)
        program,table,continuous = NC._continuous_program(program,table,ds)
        if continuous_blocks
            program,blocks = NC._continuous_blocks(program)
            continuous = (;continuous...,blocks)
        end
        report = (;report...,continuous)
    end
    used = Set(program.forms.ids)
    all(id->id==1 || id in used,program.objective_ids) || throw(ArgumentError("The invariant cones do not cover the objective."))
    if realify
        sides = [b.side for b in program.blocks]
        program = NC.real_program(program)
        report = (;report...,doubled=[i for (i,b) in enumerate(program.blocks) if b.side!=sides[i]])
    end
    NC.refresh_coverage!(table,program)
    return program,table,report
end

# Square points are based at site (0,0). On L=2 some abstract D4
# elements have the same action; use the faithful permutation image.
function _square_point_generators(::Type{W}, L) where {W}
    site(x,y) = 1 + mod(x,L) + L*mod(y,L)
    [_symmetry_action(W,L^2,[site(f(x,y)...) for y in 0:L-1 for x in 0:L-1])
        for f in ((x,y)->(-y,x),(x,y)->(-x,y))]
end

function _square_pack(w::NC.PauliWord, H, ::Type{W}) where {W}
    row = NC._pack_word(W,w)
    all(<=(H.registry.nsites),NC.sites(row)) || throw(ArgumentError("Square row exceeds the registry."))
    row
end

_square_row(w::NC.PauliWord, H, ::Type{W}) where {W} =
    NC.SectorTerms(Dict(_square_pack(w,H,W)=>one(Complex{Rational{BigInt}})))

function _square_row(p::NC.Polynomial, H, ::Type{W}) where {W}
    p.registry === H.registry || throw(ArgumentError("Square rows must share the Hamiltonian registry."))
    NC._sector_terms(p,W)
end

function _square_action(g, row::NC.SectorTerms{W,C}) where {W,C}
    out = NC.SectorTerms(Dict{W,C}())
    for (w,c) in row.terms
        NC.act!(out,g,w,NC.is_antiunitary(g) ? conj(c) : c)
    end
    out
end

function _square_derivative(d, row::NC.SectorTerms{W}) where {W}
    out = NC.SectorTerms(W)
    for (w,c) in row.terms
        NC.derive!(out,d,w,c)
    end
    out
end

function _square_dependent!(pivots, row::NC.SectorTerms{W,C}, index) where {W,C}
    NC._insert_generator_equation!(pivots,Dict{Int,C}(index[w]=>c for (w,c) in row.terms);constant=false)
end

function _square_fixed_rows(rows, anti, index)
    pivots = Dict{Int,Dict{Int,Complex{Rational{BigInt}}}}()
    fixed = eltype(rows)[]
    for row in rows
        image = _square_action(anti,row)
        for odd in (false,true)
            out = NC.SectorTerms(keytype(row.terms))
            for (w,c) in row.terms
                NC.push_term!(out,odd ? im*c/2 : c/2,w)
            end
            for (w,c) in image.terms
                NC.push_term!(out,odd ? -im*c/2 : c/2,w)
            end
            _square_dependent!(pivots,out,index) || push!(fixed,out)
        end
    end
    length(fixed)==length(rows) || error("Incomplete K-fixed row decomposition.")
    fixed
end

# Apply either m² = -D_z² or the spin Casimir = -Σ D_a² without
# expanding the input span to all axis words on any of its supports.
function _square_spectral_factor(row::NC.SectorTerms{W}, ds, eigenvalue, denominator) where {W}
    out = NC.SectorTerms(W)
    for (w,c) in row.terms
        NC.push_term!(out,-c*eigenvalue/denominator,w)
    end
    for d in ds, (w,c) in _square_derivative(d,_square_derivative(d,row)).terms
        NC.push_term!(out,-c/denominator,w)
    end
    out
end

function _square_spin_rows(rows, ds, degree, index; check_dimension=true)
    isempty(ds) && return [(nothing,rows)]
    # Every site contributes spin 0 or 1, so 0:D contains the entire
    # integer-spin spectrum and each J copy contributes exactly one m=0.
    mz0 = map(rows) do row
        for m in 1:degree
            isempty(row.terms) && break
            row = _square_spectral_factor(row,(last(ds),),m^2,-m^2)
        end
        row
    end
    R = Complex{Rational{BigInt}}
    blocks = map(0:degree) do J
        pivots = Dict{Int,Dict{Int,R}}()
        basis = eltype(rows)[]
        λ = J*(J+1)
        for row in mz0
            for j in 0:degree
                (j==J || isempty(row.terms)) && continue
                μ = j*(j+1)
                row = _square_spectral_factor(row,ds,μ,λ-μ)
            end
            _square_dependent!(pivots,row,index) || push!(basis,row)
        end
        J,basis
    end
    !check_dimension || sum((2J+1)*length(b) for (J,b) in blocks) == length(rows) || error("Incomplete spin decomposition.")
    filter(p->!isempty(last(p)),blocks)
end

function _square_numeric_basis(A)
    F = LinearAlgebra.svd(A)
    # Reject an ambiguous numerical rank; check total multiplicity again
    # against the exact dimension after summing the momentum-star blocks.
    any(s->1e-12<s<1e-8,F.S) && error("Ambiguous square symmetry projection rank.")
    F.U[:,findall(>=(1e-8),F.S)]
end

_square_columns(A) = [[(A[i,j],i) for i in axes(A,1) if abs(A[i,j])>1e-12] for j in axes(A,2)]

function _square_fourier_lift(registry, orbits::NC.AbelianRows{W}, characters, cols) where {W}
    [_wf_polynomial(registry,W,((NC.scale(orbits.phases[a][t],c*characters[orbits.positions[a][t]]/sqrt(length(ws))),w)
        for (c,a) in col for ws in (orbits.words[a],) for (t,w) in enumerate(ws))) for col in cols]
end

function _square_lattice_sectors!(program, table, spin_blocks, translations, points, anti, recovery, folded, registry, packed_orbits)
    point = NC.OrbitGroup(points)
    permutations = NC._character_permutations(translations,point)
    labels = Tuple{Int,Int}[]
    report = NamedTuple[]
    sink = NC.OrbitFormSink(table,Dict{Int32,ComplexF64}())
    projector_cache = Dict{Tuple{Vararg{Int}},Tuple{typeof(point),Vector{Vector{ComplexF64}},Vector{Int}}}()
    for (J,rows) in spin_blocks
        orbits = if isnothing(packed_orbits)
            words = sort!(unique!([w for row in rows for w in keys(row.terms)]))
            NC.AbelianRows(words,translations)
        else
            packed_orbits
        end
        index = Dict(w=>(a,t) for (a,ws) in enumerate(orbits.words) for (t,w) in enumerate(ws))
        cache = Dict{Tuple{Int,Int},NC.LinearForms{ComplexF64}}()
        seen = Set{Int}()
        dimension = 0
        for k in eachindex(translations.coordinates)
            k in seen && continue
            star = Set(p[k] for p in permutations)
            union!(seen,star)
            characters = cispi.(-2 .* translations.characters[k,:] ./ translations.modulus)
            M = zeros(ComplexF64,length(orbits.words),length(rows))
            for (j,row) in enumerate(rows), (w,c) in row.terms
                a,t = index[w]
                orbits.admissible[k,a] || continue
                phase = conj(orbits.phases[a][t])
                M[a,j] += NC.scale(phase,ComplexF64(c)*conj(characters[orbits.positions[a][t]]))/sqrt(length(orbits.words[a]))
            end
            Q = _square_numeric_basis(M)
            isempty(Q) && continue
            dimension += length(star)*size(Q,2)
            stabilizer = Tuple(findall(p->p[k]==k,permutations))
            little,projections,copies = get!(projector_cache,stabilizer) do
                group = NC._little_group(point.elements[collect(stabilizer)])
                projections,ranks = NC.shared_projections(group)
                # A minimal projector of rank r in irrep dimension d
                # represents d/r copies of its retained range.
                sw = NC.SymbolicWedderburn
                irr = sw.irreducible_characters(sw.CharacterTable(Rational{Int},group))
                real_irr,_ = sw.affordable_real(irr,ones(Int,length(irr)))
                copies = [Int(real(χ(one(group)))) ÷ r for (χ,r) in zip(real_irr,ranks)]
                group,projections,copies
            end
            actions = [NC._momentum_action(orbits,index,h,k) for h in little.elements]
            sectors = map(projections) do coefficients
                P = zeros(ComplexF64,size(Q))
                for (action,c) in zip(actions,coefficients), a in axes(Q,1)
                    phase,b = action[a]
                    @views P[b,:] .+= (c*phase).*Q[a,:]
                end
                _square_columns(_square_numeric_basis(P))
            end
            sum(c*length(s) for (c,s) in zip(copies,sectors)) == size(Q,2) || error("Incomplete little-group decomposition.")
            original_sectors = sectors
            fixed = false
            if !isnothing(anti)
                negative = findfirst(==(mod.(-translations.coordinates[k],translations.orders)),translations.coordinates)
                candidates = point.elements[findall(p->p[negative]==k,permutations)]
                sectors,fixed = NC._real_character_sectors(sectors,orbits,index,k,anti,candidates)
            end
            NC._assemble_character!(program,sink,orbits,characters,k,sectors,labels,cache;real_blocks=fixed,original_sectors)
            for (s,cols) in enumerate(sectors)
                isempty(cols) && continue
                isnothing(recovery) || push!(recovery,()->_square_fourier_lift(registry,orbits,characters,cols))
                isnothing(folded) || push!(folded,(;orbits, characters, momentum=Tuple(translations.coordinates[k]), cols))
                push!(report,(; J,momentum=Tuple(translations.coordinates[k]),star_size=length(star),
                    little_order=length(little),irrep=s,irrep_copies=copies[s],side=length(cols)))
            end
        end
        !isnothing(packed_orbits) || dimension == length(rows) || error("Incomplete translation-star decomposition.")
    end
    for (i,b) in enumerate(program.blocks)
        program.blocks[i] = NC.ConeBlock(b.kind,b.side,b.first_form,0)
    end
    report
end

"""
    workflow_square_symmetry_program(H, rows, L; strengthening=nothing,
        coordinate_words=(), spin=true, spatial=true, real_gauge=true,
        recover_gram_rows=false, recover_folded_gram=false, compact_moments=true)

Compile a closed linear row span on the square torus. Rows may be packed
Pauli words or polynomials on `H`'s registry. Keep all spin copies and one
representative of every momentum star and little-group irrep; these are
operator sectors, not restrictions on the ground state's quantum numbers.
`real_gauge` averages complex conjugation K, then uses K-fixed rows (iP
for an odd-Y Pauli P); it does not split by Y parity. If a further block
has no proven fixed basis it retains its Hermitian cone. Strengthening
forms are appended before exact spin moment elimination; the caller must
supply a symmetry-compatible collection of strengthening constraints.
With `recover_gram_rows=true`, `report.gram_rows()` returns the polynomial
rows of each original moment cone, in conic order, before any real lifting
or pruning. The lift does not contain strengthening recipes or constitute
a rigorous certificate; audit its residual in ambient Pauli coordinates.
`recover_folded_gram=true` instead retains compact Fourier orbit recipes
for G5 (requires spatial symmetry and the compact square moment table).
`compact_moments=false` selects the generic alias table as a reference.
"""
function workflow_square_symmetry_program(H, rows, L::Int; strengthening=nothing,
    coordinate_words=(), spin=true, spatial=true, real_gauge=true, recover_gram_rows=false,
    recover_folded_gram=false, compact_moments=true)
    _check_hamiltonian(H)
    L>=2 && H.registry.nsites==L^2 || throw(ArgumentError("Expected an L×L square registry with L ≥ 2."))
    recover_folded_gram && !(spatial && compact_moments) &&
        throw(ArgumentError("Folded Gram recovery requires compact spatial square moments."))
    W = eltype(H.data.words)
    R = Complex{Rational{BigInt}}
    objective = [(ComplexF64(c),w) for (c,w) in NC.terms(H)]
    ds = spin ? NC.spin_rotations(H.registry) : ()
    shifts = spatial ? _workflow_translations(W,L^2,(L,L),:auto) : NC.CliffordMap{W}[]
    points = spatial ? _square_point_generators(W,L) : NC.CliffordMap{W}[]
    anti = real_gauge ? NC.CliffordMap(NC.AntiunitaryPauli(one(W))) : nothing
    generators = (shifts...,points...,(isnothing(anti) ? () : (anti,))...)
    NC._check_objective(objective,generators)
    NC._check_derivation_polynomial(H,W,ds,"square Hamiltonian")
    translations = spatial ? NC.TranslationGroup(shifts,[L,L]) : nothing
    monomials = rows isa AbstractVector{<:NC.PauliWord}
    packed_orbits = nothing
    closed_error() = throw(ArgumentError("Square row span is not closed under the requested spatial, spin or K generators."))
    if monomials
        words = sort!(unique!([_square_pack(w,H,W) for w in rows]))
        isempty(words) && throw(ArgumentError("Select at least one square row."))
        index = Dict(w=>i for (i,w) in enumerate(words))
        if spatial
            packed_orbits = NC.AbelianRows(words,translations)
            representatives = first.(packed_orbits.words)
        else
            representatives = words
        end
        # Translations are checked by AbelianRows. Each point map normalizes
        # those translations, and every spin derivation commutes with them.
        for w in representatives
            all(g->haskey(index,last(NC.act(g,w))),points) || closed_error()
            all(d->all(in(keys(index)),keys(NC._derivative(d,w))),ds) || closed_error()
        end
        basis = [NC.SectorTerms(Dict(w=>one(R))) for w in representatives]
        row_dimension = length(words)
    else
        input = [_square_row(row,H,W) for row in rows]
        words = sort!(unique!([w for row in input for w in keys(row.terms)]))
        index = Dict(w=>i for (i,w) in enumerate(words))
        pivots = Dict{Int,Dict{Int,R}}()
        basis = filter(row->!_square_dependent!(pivots,row,index),input)
        isempty(basis) && throw(ArgumentError("Square rows span the zero space."))
        for row in basis, image in Iterators.flatten(((_square_action(g,row) for g in generators),(_square_derivative(d,row) for d in ds)))
            all(in(keys(index)),keys(image.terms)) && _square_dependent!(copy(pivots),image,index) || closed_error()
        end
        row_dimension = length(basis)
    end
    spin_blocks = _square_spin_rows(basis,ds,maximum(NC.degree,words),index;check_dimension=isnothing(packed_orbits))
    table = if spatial && compact_moments
        SquareOrbitMoments(W,L,real_gauge)
    else
        cyclic = if spatial
            step = first(shifts)
            other = (shifts[2:end]...,points...)
            normalizes(g) = NC._conjugate_action(g,step) in (step,inv(step))
            (;generator=step,normalizers=Tuple(filter(normalizes,other)),
                remaining=(filter(!normalizes,other)...,(isnothing(anti) ? () : (anti,))...))
        end
        NC.OrbitMoments(W,generators;cyclic)
    end
    program = NC.MomentProgram{ComplexF64}()
    recovery = recover_gram_rows ? Function[] : nothing
    folded_grams = recover_folded_gram ? NamedTuple[] : nothing
    if spatial
        sectors = _square_lattice_sectors!(program,table,spin_blocks,translations,points,anti,recovery,folded_grams,H.registry,packed_orbits)
        sum(s.star_size*s.irrep_copies*s.side*(isnothing(s.J) ? 1 : 2s.J+1) for s in sectors) == row_dimension ||
            error("Square sectors do not recover the input row dimension.")
    else
        families = map(spin_blocks) do (_,b)
            fixed = isnothing(anti) ? b : _square_fixed_rows(b,anti,index)
            [[(ComplexF64(c),w) for (w,c) in sort!(collect(row.terms);by=first)] for row in fixed]
        end
        block = NC.RowBlock(words,[(1.0+0im,one(W))],NC.HPSD,0)
        NC._append_family_sectors!(program,table,block,families;real_blocks=!isnothing(anti))
        if !isnothing(recovery)
            for cols in families
                push!(recovery,()->[_wf_polynomial(H.registry,W,col) for col in cols])
            end
        end
        sectors = [(;J,momentum=nothing,star_size=1,little_order=1,irrep=1,irrep_copies=1,side=length(b)) for (J,b) in spin_blocks]
    end
    for w in coordinate_words
        NC.orbit_coordinates!(table,NC._pack_word(W,w))
    end
    program,table = NC._orbit_objective!(program,table,objective)
    isnothing(strengthening) || (program = _append_strengthening!(program,table,strengthening))
    continuous = nothing
    if spin
        program,table,continuous = NC._continuous_program(program,table,ds)
    end
    sym = NC.FamilySymmetry([shifts;points];moment_only=isnothing(anti) ? () : (anti,),continuous=ds,detected=true)
    gram_rows = isnothing(recovery) ? nothing : ()->[expand() for expand in recovery]
    return program,table,(;route=:square_spin_little_group,sectors,row_dimension,
        continuous,symmetry=sym,symmetry_generators=generators,gram_rows,folded_grams)
end

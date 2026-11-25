# Lift interface: witnesses are expressed in the original Pauli algebra, not
# in realified/quotiented moment or solver coordinates. Verification does not
# trust a solver objective, PSD tolerance, or any claimed lift map.

function _wf_gram_psd!(Q)
    size(Q,1) == size(Q,2) && LinearAlgebra.ishermitian(Q) ||
        throw(ArgumentError("A rational Gram witness must be square and Hermitian."))
    for i in axes(Q,1)
        pivot = real(Q[i,i])
        pivot >= 0 || throw(ArgumentError("Gram witness is not positive semidefinite."))
        if iszero(pivot)
            all(iszero, view(Q,i+1:size(Q,1),i)) ||
                throw(ArgumentError("A zero Gram pivot has a nonzero column."))
        else
            for j in i+1:size(Q,1), k in i+1:size(Q,1)
                Q[j,k] -= Q[j,i]*conj(Q[k,i])/pivot
            end
        end
    end
    return nothing
end

function _wf_witness_polynomial(reference, p)
    p isa NC.Polynomial && p.registry === reference.registry ||
        throw(ArgumentError("A Pauli witness polynomial uses the wrong registry."))
    return _wf_exact(p)
end

function _wf_certificate_bound(p, level, sectors; squares=(), grams=(), sector_multipliers=(),
                               inequalities=(), window_multipliers=())
    λ = _PauliQ(level)
    residual = p - λ
    for f in squares
        a = _wf_witness_polynomial(p,f)
        residual -= a'*a
    end
    for witness in grams
        rows = [_wf_witness_polynomial(p,v) for v in witness.rows]
        Q = Matrix{_PauliCQ}(witness.gram)
        size(Q) == (length(rows),length(rows)) || throw(ArgumentError("Gram dimensions do not match its rows."))
        _wf_gram_psd!(copy(Q))
        for j in eachindex(rows), i in eachindex(rows)
            iszero(Q[i,j]) || (residual -= Q[i,j]*(rows[i]'*rows[j]))
        end
    end
    if !isempty(sector_multipliers)
        length(sector_multipliers) == length(sectors) || throw(ArgumentError("Provide one multiplier per sector."))
        for (s,u) in zip(sectors, sector_multipliers)
            a = s*_wf_witness_polynomial(p,u)
            residual -= (a+a')*(1//2)
        end
    end
    if !isempty(window_multipliers)
        length(window_multipliers) == length(inequalities) ||
            throw(ArgumentError("Provide one multiplier per energy-window inequality."))
        for (g, μ) in zip(inequalities, window_multipliers)
            q = _PauliQ(μ)
            q >= 0 || throw(ArgumentError("Energy-window multipliers must be nonnegative."))
            residual -= q*g
        end
    end
    NC._checked_hermitian(residual, "Certificate residual")
    # Every Hermitian Pauli is unitary; all residual coefficients are REAL
    # rationals. The norm enclosure and the returned bound are exact.
    enclosure = sum((abs(real(c)) for (c,_) in NC.terms(residual)); init=zero(_PauliQ))
    return (; bound=λ-enclosure, level=λ, residual_l1=enclosure,
            residual=Tuple(NC.terms(residual)), status=:verified_rational)
end

_wf_signature(p) = (p.registry, Tuple(NC.terms(p)))

struct PauliSOSCertificate{S,T,R}
    hamiltonian::S
    target::T
    bound::_PauliQ
    level::_PauliQ
    residual_l1::_PauliQ
    residual::R
    status::Symbol
    function PauliSOSCertificate(H, level; target_mz_sector=nothing, target_spin_j=nothing,
                                squares=(), grams=(), sector_multipliers=())
        HQ = _wf_exact(NC._checked_hermitian(H,"Hamiltonian"))
        target = _wf_target(HQ; mz=target_mz_sector, spin_j=target_spin_j)
        sig = _wf_signature(HQ)
        label = (; mz=target.mz, spin_j=target.spin_j)
        result = _wf_certificate_bound(HQ,level,target.sectors; squares,grams,sector_multipliers)
        return new{typeof(sig),typeof(label),typeof(result.residual)}(
            sig,label,result.bound,result.level,result.residual_l1,result.residual,result.status)
    end
end

"""
    certify_rational_sos(H, level; squares=(), grams=(), target_mz_sector=nothing,
                         target_spin_j=nothing, sector_multipliers=())

Verify a lifted SOS witness in the original Pauli algebra. `squares` contains
polynomials f for f†f; each `grams` entry is `(;rows,gram)` with a Hermitian
PSD rational matrix tested by exact Schur elimination. Binary floats are
interpreted as their exact stored values, not as nearby rational numbers. Use
rational factors to avoid introducing indefinite rounded Gram matrices. The
optional one-sided sector terms are Hermitian parts of (Q-q)u. The corrected
bound subtracts the exact ℓ₁ residual; it does not certify the uncorrected
level. An empty witness gives a norm bound, not a lifted solver dual.
"""
certify_rational_sos(H,level;kwargs...) = PauliSOSCertificate(H,level;kwargs...)

"""
    pauli_dual_blocks(form, model; sos_variables=nothing)

Read actual JuMP multipliers in `MomentStandardForm.cones` order, without
requiring a solver status to serve as a proof. For an explicit moment model
read the constraint duals. For the SOS maximization pass its primal variables
in `svec` order. Return `nothing` when the relevant point is unavailable.
The model must have been built from this form, with no extra constraints.
"""
function pauli_dual_blocks(form, model; sos_variables=nothing)
    all(c -> c.kind in (:psd,:zero), form.cones) ||
        throw(ArgumentError("G5 supports real PSD cones and scalar equalities only."))
    if isnothing(sos_variables)
        JuMP.has_duals(model) || return nothing
        psd = JuMP.all_constraints(model, Vector{JuMP.AffExpr}, JuMP.MOI.PositiveSemidefiniteConeTriangle)
        eq = JuMP.all_constraints(model, JuMP.AffExpr, JuMP.MOI.EqualTo{Float64})
        length(psd) == count(c -> c.kind === :psd, form.cones) &&
            length(eq) == sum((c.length for c in form.cones if c.kind === :zero); init=0) &&
            length(JuMP.all_constraints(model; include_variable_in_set_constraints=false)) == length(psd)+length(eq) ||
            throw(ArgumentError("JuMP constraints do not match the original moment form."))
        p = e = 0
        return map(form.cones) do c
            if c.kind === :psd
                Matrix(JuMP.dual(psd[p+=1]))
            else
                [JuMP.dual(eq[e+=1]) for _ in 1:c.length]
            end
        end
    end
    JuMP.has_values(model) || return nothing
    z = JuMP.value.(sos_variables)
    length(z) == length(form.b) || throw(ArgumentError("SOS point does not match the original moment form."))
    all(isfinite,z) || throw(ArgumentError("The supplied SOS point must be finite."))
    return map(form.cones) do c
        c.kind === :zero && return z[c.first:c.first+c.length-1]
        [_PauliQ(z[c.first+NC.tri_offset(c.side,max(i,j),min(i,j))]) /
            (i==j ? 1 : _PauliQ(sqrt(2.0))) for i in 1:c.side, j in 1:c.side]
    end
end

# A floating factor only proposes a decomposition. The exact Gershgorin
# lower enclosure of D-F'F proves D+shift*I is PSD by Weyl. Replaying the
# shifted D (not F'F) keeps all original off-diagonal dual entries.
function _wf_dual_shift(D)
    function shift(A)
        return max(zero(_PauliQ), maximum((sum(abs,view(A,:,i))-abs(A[i,i])-A[i,i] for i in axes(A,1)); init=zero(_PauliQ)))
    end
    direct = shift(D)
    iszero(direct) && return direct
    candidate = Float64.(D)
    all(isfinite,candidate) || return direct
    eig = LinearAlgebra.eigen(LinearAlgebra.Symmetric(candidate))
    F = _PauliQ.(sqrt.(max.(eig.values,0.0)) .* eig.vectors')
    return min(direct,shift(D-F'*F))
end

function _wf_dual_entries(H, lift, side)
    layout = get(lift,:layout,:real)
    layout in (:real,:realified) || throw(ArgumentError("Unsupported real PSD lift layout: $layout."))
    if lift.kind === :gram
        rows = [_wf_witness_polynomial(H,p) for p in lift.rows]
        n = length(rows)
        entry = (i,j) -> rows[i]'*rows[j]
    elseif lift.kind === :rdm
        sites = Int[lift.sites...]
        1 <= length(sites) <= 6 && allunique(sites) && all(s -> 1 <= s <= H.registry.nsites,sites) ||
            throw(ArgumentError("An RDM lift needs 1–6 distinct existing sites in matrix order."))
        D = _wf_density(H,sites)
        V = get(lift,:columns,nothing)
        if !isnothing(V)
            V = Matrix{_PauliCQ}(V)
            size(V,1) == size(D,1) || throw(ArgumentError("RDM compression has the wrong input size."))
            D = [sum((conj(V[a,i])*V[b,j]*D[a,b]
                for a in axes(V,1) for b in axes(V,1) if !iszero(V[a,i]*V[b,j])); init=zero(H))
                for i in axes(V,2), j in axes(V,2)]
        end
        n = size(D,1)
        entry = (i,j) -> D[i,j]
    else
        throw(ArgumentError("Unsupported PSD lift $(lift.kind); provide a physical Gram or RDM map."))
    end
    full_side = layout === :real ? n : 2n
    indices = get(lift,:indices,1:full_side)
    length(indices) == side && allunique(indices) && all(i -> 1 <= i <= full_side,indices) ||
        throw(ArgumentError("PSD lift indices do not match the cone side."))
    # Re M and [Re M -Im M; Im M Re M] are PSD on real vectors whenever M is
    # Hermitian PSD. A principal restriction, including a permutation, is PSD.
    return function (i,j)
        a,b = indices[i],indices[j]
        p = entry(mod1(a,n),mod1(b,n))
        return (a>n)==(b>n) ? (p+p')*(1//2) : (p-p')*((a>n ? -im : im)*(1//2))
    end
end

# Reynolds-average a Hermitian residual over signed Pauli orbits. The group
# need not be enumerated: every image word has the same number of preimages,
# or cancels if both signs occur. Antiunitaries act on real coefficients here.
function _wf_average_dual_residual(H::NC.Polynomial{W}, residual, generators) where {W}
    isempty(generators) && return residual
    terms = Dict(w => c for (c,w) in NC.terms(H))
    checked = map(generators) do g
        g isa NC.CliffordMap{W} && length(g.columns) <= 2H.registry.nsites ||
            throw(ArgumentError("G5 averaging needs signed Clifford symmetries on H's registry."))
        # Revalidate the symplectic columns, which are mutable vectors.
        a = NC.CliffordMap(g.columns,g.negative; antiunitary=NC.is_antiunitary(g))
        _symmetry_preserves(a,terms) || throw(ArgumentError("An averaging symmetry does not exactly preserve H."))
        return a
    end
    remaining = Dict(w => real(c) for (c,w) in NC.terms(residual))
    out = NC.PolynomialBuffer{W,_PauliCQ}()
    while !isempty(remaining)
        w = first(keys(remaining))
        orbit = Dict(w => 1)
        queue = W[w]
        cancels = false
        for u in queue, g in checked
            p,v = NC.act(g,u)
            s = NC.scale(p,orbit[u])
            if haskey(orbit,v)
                cancels |= orbit[v] != s
            else
                orbit[v] = s
                push!(queue,v)
            end
        end
        c = sum(s*pop!(remaining,v,zero(_PauliQ)) for (v,s) in orbit)/length(orbit)
        cancels && continue
        for (v,s) in orbit
            NC.push_term!(out,c*s,v)
        end
    end
    return NC._polynomial(H.registry,out)
end

# On degree-d Pauli words collective rotations act as d tensor factors of
# spin 1. Their Casimir C=-Σ D_axis^2 has eigenvalues j(j+1), 0≤j≤d.
# The product of (I-C/[j(j+1)]) is the exact Haar projector, with no
# numerical nullspace or irrational rotation matrices in the replay.
function _wf_spin_average_dual_residual(H::NC.Polynomial{W}, residual, max_weight) where {W}
    max_weight >= 0 || throw(ArgumentError("The exact spin-projection weight limit must be nonnegative."))
    rotations = NC.spin_rotations(H.registry)
    function derivative(d, terms)
        out = NC.SectorTerms(W)
        for (w,c) in terms
            NC.derive!(out,d,w,c)
        end
        return out.terms
    end
    objective = NC._sector_terms(H,W).terms
    all(d -> isempty(derivative(d,objective)),rotations) ||
        throw(ArgumentError("Collective spin averaging requires exact spin invariance of H."))
    weight = NC.degree(residual)
    if weight > max_weight
        # Collective quarter-turns around x and z generate the 24 proper
        # signed permutations of the three Pauli axes. This is a rigorous
        # finite average, not a claim of full spin invariance of the result.
        n = H.registry.nsites
        x = [NC.pauli_letter(W,i,0) for i in 1:n]
        y = [NC.pauli_letter(W,i,1) for i in 1:n]
        z = [NC.pauli_letter(W,i,2) for i in 1:n]
        rx = NC.CliffordMap([x;y],[falses(n);trues(n)]) # X→X, Z→−Y
        rz = NC.CliffordMap([y;z])                      # X→Y, Z→Z
        return _wf_average_dual_residual(H,residual,(rx,rz)), :clifford24
    end
    values = NC._sector_terms(residual,W).terms
    for j in weight:-1:1
        next = copy(values)
        for d in rotations
            for (w,c) in derivative(d,derivative(d,values))
                NC._sector_add!(next,w,c/(j*(j+1)))
            end
        end
        values = next
    end
    return _wf_polynomial(H.registry,W,((c,w) for (w,c) in values)), :su2_exact
end

# Audit the map in the actual quotient coordinates. This threshold only
# detects incorrect provenance; it is NOT used in the bound or PSD checks.
function _wf_audit_dual_row!(expected, actual, table, p, scale)
    empty!(expected)
    for (c,w) in NC.terms(p)
        terms = NC.moment_terms(table,w)
        isnothing(terms) && throw(ArgumentError("Physical lift uses a word missing from the solved moment table."))
        for (id,a) in terms
            NC._sector_add!(expected, Int(id), scale*real(c*_PauliCQ(a)))
        end
    end
    error = max(maximum((abs(v-get(actual,i,0)) for (i,v) in expected); init=zero(_PauliQ)),
        maximum((abs(v-get(expected,i,0)) for (i,v) in actual); init=zero(_PauliQ)))
    magnitude = maximum(abs,values(actual); init=one(_PauliQ))
    error <= _PauliQ(1e-10)*magnitude ||
        throw(ArgumentError("Physical lift does not match a row of the solved conic form."))
    return error
end

"""
    verify_pauli_dual(H, form, table, dual, lifts;
                      symmetry_generators=(), spin_average=false, spin_max_weight=6)

Replay a Pauli SDP dual in original Hamiltonian units. Supply one lift per
cone, in actual conic order: `(;kind=:gram, rows)` (polynomial Gram rows),
`(;kind=:rdm, sites)` (physical reduced density matrix, optionally with
compression `columns`), or `(;kind=:stationarity, rows)` (Hermitian test
operators for the scalar equations `i[H, rows[k]]`). PSD lifts accept
`layout=:real` (default) or `:realified`, and optional principal `indices`
in that real layout. Give an explicit lift even for a zero multiplier.
Unsupported strengthenings and unmapped affine rows fail closed.

The proof uses the ambient Gram and density matrices and the exact
stationarity of a ground state. Optionally pass exact signed Clifford
symmetries of H (including antiunitaries) to average the residual and
recover cancellations of the moment quotient. `spin_average=true` also
applies the exact collective spin Haar projector after checking that H is
spin invariant. If the residual exceeds `spin_max_weight`, it uses the
exact 24-Clifford average instead (`spin_projection=:clifford24`). Otherwise
the result reports `:su2_exact`. Each average preserves the ground-space trace. The proof
does not assume that a ground state is symmetric or trust a rounded quotient.
Binary inputs are converted exactly. For repaired dual
blocks let `level = offset - b'z`. The **only** certified level is
`bound = level + residual_identity - residual_nonidentity_l1`, where the
residual is `H - level*I - dual_contribution` in the full Pauli algebra,
averaged if requested.
`raw_dual_objective`, both balance norms and the lift matching error are
diagnostics, not bounds. This proof targets the unrestricted ground energy.
"""
function verify_pauli_dual(H, form, table, dual, lifts;
                           symmetry_generators=(), spin_average=false, spin_max_weight::Int=6)
    HQ = _wf_exact(NC._checked_hermitian(H,"Hamiltonian"))
    length(form.cones) == length(dual) == length(lifts) ||
        throw(ArgumentError("Provide a dual block and physical lift for every cone."))
    # Traverse the original affine form by row without retaining an exact
    # dictionary for every matrix entry.
    transposed = copy(transpose(form.A))
    actual = Dict{Int,_PauliQ}()
    expected = empty(actual)
    function audit_row(k,p,scale)
        empty!(actual)
        iszero(form.b[k]) || (actual[1] = _PauliQ(form.b[k]))
        for t in transposed.colptr[k]:transposed.colptr[k+1]-1
            actual[transposed.rowval[t]+1] = -_PauliQ(transposed.nzval[t])
        end
        return _wf_audit_dual_row!(expected,actual,table,p,scale)
    end
    sqrt2 = _PauliQ(sqrt(2.0))
    z = zeros(_PauliQ,length(form.b))
    raw = similar(z)
    shifts = zeros(_PauliQ,length(dual))
    max_lift_error = zero(_PauliQ)
    W = eltype(HQ.data.words)
    contribution = NC.SectorTerms(W)
    function accumulate!(q,p)
        iszero(q) && return nothing
        for (c,w) in NC.terms(p)
            NC.push_term!(contribution,q*c,w)
        end
    end
    for (index,(cone,D,lift)) in enumerate(zip(form.cones,dual,lifts))
        all(isfinite,D) || throw(ArgumentError("The supplied dual point must be finite."))
        k = cone.first
        if cone.kind === :psd
            D isa AbstractMatrix && size(D) == (cone.side,cone.side) &&
                eltype(D) <: Real && LinearAlgebra.issymmetric(D) ||
                throw(ArgumentError("PSD dual must be real symmetric and match the cone side."))
            Q = Matrix{_PauliQ}(D)
            δ = shifts[index] = _wf_dual_shift(Q)
            entry = _wf_dual_entries(HQ,lift,cone.side)
            for j in 1:cone.side, i in j:cone.side
                p = entry(i,j)
                scale = i==j ? one(_PauliQ) : sqrt2
                max_lift_error = max(max_lift_error,audit_row(k,p,scale))
                raw[k] = scale*Q[i,j]
                q = i==j ? Q[i,j]+δ : Q[i,j]
                z[k] = scale*q
                accumulate!(i==j ? q : 2q,p)
                k += 1
            end
        elseif cone.kind === :zero && lift.kind === :stationarity
            D isa AbstractVector{<:Real} && length(D) == cone.length == length(lift.rows) ||
                throw(ArgumentError("Stationarity rows and duals must match the cone length."))
            for (d,u) in zip(D,lift.rows)
                U = NC._checked_hermitian(_wf_witness_polynomial(HQ,u),"Stationarity test operator")
                p = im*(HQ*U-U*HQ)
                max_lift_error = max(max_lift_error,audit_row(k,p,one(_PauliQ)))
                z[k] = raw[k] = _PauliQ(d)
                accumulate!(z[k],p)
                k += 1
            end
        else
            throw(ArgumentError("Unsupported cone/lift pair $(cone.kind)/$(lift.kind)."))
        end
        k == cone.first+cone.length || throw(ArgumentError("Malformed conic block length."))
    end
    modeled_objective = Dict(i+1 => _PauliQ(c) for (i,c) in pairs(form.c) if !iszero(c))
    iszero(form.offset) || (modeled_objective[1] = _PauliQ(form.offset))
    max_lift_error = max(max_lift_error,_wf_audit_dual_row!(expected,modeled_objective,table,HQ,one(_PauliQ)))
    objective(v) = _PauliQ(form.offset)-sum((_PauliQ(b)*q for (b,q) in zip(form.b,v)); init=zero(_PauliQ))
    function balance_inf(v)
        mapreduce(max,axes(form.A,2); init=zero(_PauliQ)) do j
            terms = form.A.colptr[j]:form.A.colptr[j+1]-1
            abs(_PauliQ(form.c[j])+sum((v[form.A.rowval[t]]*_PauliQ(form.A.nzval[t]) for t in terms); init=zero(_PauliQ)))
        end
    end
    level = objective(z)
    residual = HQ-level-_wf_polynomial(HQ.registry,W,((c,w) for (w,c) in contribution.terms))
    NC._checked_hermitian(residual,"Dual certificate residual")
    residual = _wf_average_dual_residual(HQ,residual,symmetry_generators)
    spin_projection = :none
    if spin_average
        residual,spin_projection = _wf_spin_average_dual_residual(HQ,residual,spin_max_weight)
    end
    residual_identity = sum((real(c) for (c,w) in NC.terms(residual) if isone(w)); init=zero(_PauliQ))
    residual_nonidentity_l1 = sum((abs(real(c)) for (c,w) in NC.terms(residual) if !isone(w)); init=zero(_PauliQ))
    return (; bound=level+residual_identity-residual_nonidentity_l1, level,
        residual_identity, residual_nonidentity_l1, residual=Tuple(NC.terms(residual)),
        shifts, spin_projection, raw_dual_objective=objective(raw), raw_balance_inf=balance_inf(raw),
        repaired_balance_inf=balance_inf(z), max_lift_error, status=:verified_rational)
end

"""
    square_lifts(H, seed, before, real, pruned, pruning, symmetry_report; form)

Map a recovered square Gram basis and the seed's physical strengthening
sources to the final real moment cones. `before` includes the appended
strengthening, `real` is before term sparsity, and `pruned` is what the solver
uses. Pass the already-built `form` to avoid building it again. Principal
cliques and recorded constant-singleton deletions are supported; aggregate
chordal auxiliaries and non-stationarity/non-RDM strengthening are not.
"""
function square_lifts(H, seed, before, real, pruned, pruning, symmetry_report;
                      form=NC.moment_standard_form(pruned))
    strengthen = seed.strengthening
    isempty(strengthen.inequalities) || throw(ArgumentError("G5 does not support square inequality strengthening."))
    pruning.algorithm in (:None,:Maximal,:MMD) && isempty(pruning.auxiliary_ids) ||
        throw(ArgumentError("G5 square lifts require principal cliques without chordal auxiliaries."))
    isnothing(symmetry_report.gram_rows) && throw(ArgumentError("Build square symmetry with recover_gram_rows=true."))
    grams = symmetry_report.gram_rows()
    words, sites = seed.stationarity_words, seed.plaquettes
    ngram, neq, nrdm = length(grams),length(words),length(sites)
    neq == length(strengthen.equalities) && nrdm == length(strengthen.psd_blocks) &&
        all(b -> b.kind === :rdm_positivity,strengthen.psd_blocks) ||
        throw(ArgumentError("Square strengthening sources do not match the seed's stationarity/RDM labels."))
    length(before.blocks) == length(real.blocks) == ngram+neq+nrdm ||
        throw(ArgumentError("Realification changed the square cone order or count."))
    removed = Set(get(pruning,:pruned_singletons,Tuple{Int,Int}[]))
    lifts = NamedTuple[]
    p = c = 1
    for (i,(source,block)) in enumerate(zip(before.blocks,real.blocks))
        label = i <= ngram ? 0 : i-ngram
        source.label == block.label == label || throw(ArgumentError("Unexpected square strengthening label/order."))
        if ngram < i <= ngram+neq
            source.kind == block.kind == NC.ZERO && source.side == block.side == 1 ||
                throw(ArgumentError("A square stationarity source must be one scalar equality."))
            p <= length(pruned.blocks) &&
                (pruned.blocks[p].kind,pruned.blocks[p].label,pruned.blocks[p].side) == (NC.ZERO,label,1) ||
                throw(ArgumentError("Pruning changed a square stationarity block."))
            p += 1
            # Standard form drops projected zero equalities, not their labels.
            if c <= length(form.cones) && form.cones[c].kind === :zero && form.cones[c].label == label
                form.cones[c].length == 1 || throw(ArgumentError("Expected a scalar conic stationarity row."))
                push!(lifts,(;kind=:stationarity,rows=[_wf_word(H,words[label])]))
                c += 1
            end
            continue
        end
        source.kind in (NC.PSD,NC.HPSD) && block.kind == NC.PSD && block.side in (source.side,2source.side) ||
            throw(ArgumentError("Unexpected square PSD realification."))
        layout = block.side == source.side ? :real : :realified
        lift = if i <= ngram
            length(grams[i]) == source.side || throw(ArgumentError("Recovered square Gram rows do not match the cone side."))
            (;kind=:gram,rows=grams[i],layout)
        else
            region = sites[label-neq]
            source.side == 2^length(region) || throw(ArgumentError("Square RDM cone is not the seed's full density matrix."))
            (;kind=:rdm,sites=region,layout)
        end
        cliques = get(pruning.cliques,i,(1:block.side,))
        for indices in cliques
            if length(indices) == 1 && (i,only(indices)) in removed
                delete!(removed,(i,only(indices)))
                continue
            end
            side = length(indices)
            p <= length(pruned.blocks) && c <= length(form.cones) &&
                (pruned.blocks[p].kind,pruned.blocks[p].label,pruned.blocks[p].side) == (NC.PSD,label,side) &&
                (form.cones[c].kind,form.cones[c].label,form.cones[c].side) == (:psd,label,side) ||
                throw(ArgumentError("Principal cliques do not match the final square cone order/sides."))
            push!(lifts,(;lift...,indices))
            p += 1
            c += 1
        end
    end
    isempty(removed) && p == length(pruned.blocks)+1 && c == length(form.cones)+1 ||
        throw(ArgumentError("Some pruned square blocks or singleton removals have no matching lift."))
    return lifts
end

# Apply a Pauli polynomial exactly, without constructing a 2^N × 2^N matrix.
function _wf_apply(p, ket)
    out = zeros(_PauliCQ,length(ket))
    region = 1:p.registry.nsites
    for (c,w) in NC.terms(p)
        x,z = _wf_local_masks(w,region)
        for j in 0:length(ket)-1
            out[(j ⊻ x)+1] += c*_wf_pauli_phase(x,z,j)*ket[j+1]
        end
    end
    return out
end

struct PauliEnergyWindow{S,T}
    hamiltonian::S
    target::T
    lower::_PauliQ
    upper::_PauliQ
    function PauliEnergyWindow(H, cert::PauliSOSCertificate, ket::AbstractVector;
                               target_mz_sector=cert.target.mz, target_spin_j=cert.target.spin_j)
        HQ = _wf_exact(NC._checked_hermitian(H,"Hamiltonian"))
        signature = _wf_signature(HQ)
        signature == cert.hamiltonian || throw(ArgumentError("Energy certificate belongs to another Hamiltonian."))
        target = _wf_target(HQ; mz=target_mz_sector, spin_j=target_spin_j)
        label = (; mz=target.mz,spin_j=target.spin_j)
        cert.target in ((;mz=nothing,spin_j=nothing),label) ||
            throw(ArgumentError("Energy lower certificate does not cover this target."))
        big(2)^HQ.registry.nsites == length(ket) ||
            throw(ArgumentError("A variational ket needs one entry per computational basis vector."))
        a = _PauliCQ[x for x in ket]
        normsq = sum(abs2,a)
        normsq > 0 || throw(ArgumentError("The variational ket must be nonzero."))
        all(s -> all(iszero,_wf_apply(s,a)), target.sectors) ||
            throw(ArgumentError("The variational ket is outside the requested sector."))
        ub = real(LinearAlgebra.dot(a,_wf_apply(HQ,a)))/normsq
        return new{typeof(signature),typeof(label)}(signature,label,cert.bound,ub)
    end
end

"Construct a ground-energy enclosure from an exact lower certificate and an exact target-compatible Rayleigh quotient."
certify_energy_window(H, cert, ket; kwargs...) = PauliEnergyWindow(H,cert,ket;kwargs...)

_wf_window_bounds(window, H, target) = throw(ArgumentError("Pass a PauliEnergyWindow from certify_energy_window."))
function _wf_window_bounds(window::PauliEnergyWindow, H, target)
    _wf_signature(H) == window.hamiltonian || throw(ArgumentError("Energy window belongs to another Hamiltonian."))
    (;mz=target.mz,spin_j=target.spin_j) == window.target ||
        throw(ArgumentError("Energy window belongs to another sector."))
    return window.lower, window.upper
end

"""
    certify_observable_interval(H, O, window; lower, upper)

Verify endpoint witnesses for O and -O, respectively. Supply NamedTuples
with `level` and optional `squares`, `grams`, `sector_multipliers` and
`window_multipliers` (in order H-lo, hi-H). Returns exact rationals valid for
all sector-compatible states in the energy window, including non-invariant
states. The certificate never uses symmetry equalities, stationarity, or KMS.
No solver estimate is promoted to a verified endpoint without this check.
"""
function certify_observable_interval(H, O, window::PauliEnergyWindow; lower, upper)
    HQ = _wf_exact(NC._checked_hermitian(H,"Hamiltonian"))
    OQ = _wf_witness_polynomial(HQ, NC._checked_hermitian(O,"Observable"))
    target = _wf_target(HQ; window.target...)
    lo, hi = _wf_window_bounds(window,HQ,target)
    inequalities = (HQ-lo,hi-HQ)
    function endpoint(p; level, squares=(), grams=(), sector_multipliers=(), window_multipliers=())
        return _wf_certificate_bound(p,level,target.sectors;
            inequalities,squares,grams,sector_multipliers,window_multipliers)
    end
    lb, negub = endpoint(OQ;lower...),endpoint(-OQ;upper...)
    return (; lower=lb.bound, upper=-negub.bound, status=:verified_rational,
            energy_target=:energy_window, lower_certificate=lb, negated_upper_certificate=negub)
end

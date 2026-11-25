const _PauliQ = Rational{BigInt}
const _PauliCQ = Complex{_PauliQ}

function _wf_polynomial(registry, ::Type{W}, terms) where {W<:NC.PauliWord}
    data = NC.PolynomialBuffer{W,_PauliCQ}()
    for (c, w) in terms
        NC.push_term!(data, c, w)
    end
    return NC._polynomial(registry, data)
end

_wf_exact(p::NC.Polynomial{W}) where {W<:NC.PauliWord} =
    _wf_polynomial(p.registry, W, NC.terms(p))

_wf_word(p::NC.Polynomial{W}, w::W) where {W<:NC.PauliWord} =
    _wf_polynomial(p.registry, W, ((1, w),))

function _wf_target(H::NC.Polynomial{W}; mz=nothing, spin_j=nothing,
                    energy_target::Symbol=:target_ground) where {W<:NC.PauliWord}
    energy_target in (:target_ground, :global_ground, :energy_window) ||
        throw(ArgumentError("Unknown energy target $energy_target."))
    n = H.registry.nsites
    sectors = typeof(_wf_exact(H))[]
    if !isnothing(mz)
        mz isa Integer && -n <= mz <= n && iseven(n - mz) ||
            throw(ArgumentError("Magnetization is the eigenvalue of ∑Z, with values -N:2:N."))
        operator = _wf_polynomial(H.registry, W, ((1, NC.pauli_letter(W, i, 2)) for i in 1:n))
        push!(sectors, operator - mz)
    end
    if !isnothing(spin_j)
        spin_j isa Real && isfinite(spin_j) || throw(ArgumentError("Total spin must be finite."))
        j = _PauliQ(spin_j)
        0 <= 2j <= n && isinteger(2j) && iseven(n - numerator(2j)) ||
            throw(ArgumentError("Total spin must have the parity of N/2 and lie in [0,N/2]."))
        isnothing(mz) || abs(mz) <= 2j || throw(ArgumentError("These spin and magnetization sectors do not intersect."))
        generator = ((1//2, last(NC.mul(NC.pauli_letter(W, i, a), NC.pauli_letter(W, k, a))))
                     for i in 1:n for k in i+1:n for a in 0:2)
        push!(sectors, _wf_polynomial(H.registry, W, generator) + n//4*3 - j*(j+1))
        spin_j = j
    end
    energy_target == :global_ground && !isempty(sectors) &&
        throw(ArgumentError("A sector target does not certify unrestricted global ground states."))
    # Without commutation, compression and ground-state stationarity of H
    # differ. This prototype supports conserved sectors only.
    HQ = _wf_exact(H)
    all(s -> HQ*s == s*HQ, sectors) ||
        throw(ArgumentError("The Hamiltonian must commute exactly with each requested sector operator."))
    return (; mz, spin_j, sectors=Tuple(sectors), energy_target)
end

_wf_nonzero!(out, p) = iszero(p) ? out : push!(out, p)
function _wf_split_zero!(out, p)
    _wf_nonzero!(out, (p+p')*(1//2))
    _wf_nonzero!(out, (p-p')*(-im*(1//2)))
end

function _wf_regions(n, k)
    k == 0 && return Vector{Int}[]
    regions = Vector{Int}[]
    function choose!(sites, start)
        length(sites) == k && return push!(regions, copy(sites))
        for i in start:n-(k-length(sites))+1
            push!(sites, i)
            choose!(sites, i+1)
            pop!(sites)
        end
    end
    choose!(Int[], 1)
    return regions
end

# Kronecker order: the first listed site is the most significant ket bit.
function _wf_local_masks(w, region)
    x = z = 0
    for (site, axis) in NC._locations(w)
        bit = 1 << (length(region) - findfirst(==(site), region))
        axis != 2 && (x |= bit)
        axis != 0 && (z |= bit)
    end
    return x, z
end

_wf_pauli_phase(x, z, j) = im^count_ones(x & z) * (-1)^count_ones(z & j)

function _wf_density(H::NC.Polynomial{W}, region) where {W<:NC.PauliWord}
    d = 1 << length(region)
    entries = [NC.PolynomialBuffer{W,_PauliCQ}(d) for _ in 1:d, _ in 1:d]
    # Every nonzero Pauli matrix entry has one of these four exact coefficients.
    unit = _PauliCQ(1//d)
    coefficients = ntuple(k -> NC.scale(NC.Phase(k-1),unit), Val(4))
    for w in NC.basis_rows(W, NC.Pauli(), region, length(region))
        x, z = _wf_local_masks(w, region)
        ny = count_ones(x & z)
        for j in 0:d-1
            k = mod(ny + 2count_ones(z & j), 4)
            NC.push_term!(entries[(j ⊻ x)+1,j+1], coefficients[k+1], w)
        end
    end
    return map(b -> NC._polynomial(H.registry, b), entries)
end

"""
Replace a single *ordinary* moment block by a local density block, only if
its packed rows span the complete specified local Pauli algebra. Do not use
after projecting or dropping rows, or to replace a localizer.
"""
function workflow_density_rewrite(H::NC.Polynomial{W}, rows, region) where {W<:NC.PauliWord}
    sites = _wf_checked_region(H.registry.nsites, region)
    Set(rows) == Set(NC.basis_rows(W, NC.Pauli(), sites, length(sites))) ||
        throw(ArgumentError("Exact density replacement requires a complete ordinary local moment block."))
    return (;kind=:rdm_positivity, matrix=_wf_density(H,sites))
end

"""
    workflow_strengthening(H, order; kwargs...)

Register supplemental packed rows, Hermitian scalar equality/inequality
polynomials, and Hermitian polynomial PSD blocks before symmetry reduction.
The compiler must map all entries to the same global moment table. Supply
`moment_rows` as one packed row vector or a vector of such vectors to keep
sector relations and default KMS perturbations local to the selected family.
Stationarity also uses the degree-limited words in each requested RDM region,
without adding them as moment rows. Omit `moment_rows` to enumerate dense bases.
`rdm_regions` selects ordered subsets of up to six sites;
`rdm_positivity_k=k` alone requests every subset of size k.

`:energy_window` covers arbitrary states whose energy expectation is in the
window, not just eigenstates. Stationarity and KMS are rejected on that domain.
An explicit centered second-moment cut further narrows the domain, unless
verified for every target ground state by a certified energy window.
"""
function workflow_strengthening(H::NC.Polynomial{W}, order::Integer;
        commutator_linear_order::Integer=0, kms_psd_order::Integer=0,
        kms_frustration_guard::Bool=true, kms_null_tol::Real=1e-7, kms_rows=nothing,
        rdm_positivity_k::Integer=0, rdm_regions=nothing, lightcone_depth::Integer=0,
        rg_coarse_grain_depth::Integer=0, rg_maps=(), moment_rows=nothing,
        target_mz_sector=nothing, target_spin_j=nothing, energy_target::Symbol=:target_ground,
        energy_window=nothing, energy_upper_bound=nothing,
        variance_center=nothing, variance_upper_bound=nothing,
        rational_sos_certify::Bool=false) where {W<:NC.PauliWord}
    NC._checked_hermitian(H, "Hamiltonian")
    all(>=(0), (order, commutator_linear_order, kms_psd_order, rdm_positivity_k,
                lightcone_depth, rg_coarse_grain_depth)) || throw(ArgumentError("Orders must be nonnegative."))
    isfinite(kms_null_tol) && kms_null_tol >= 0 || throw(ArgumentError("Invalid KMS null tolerance."))
    isnothing(energy_upper_bound) || throw(ArgumentError("Use a checked energy_window instead of a bare energy_upper_bound."))
    target = _wf_target(H; mz=target_mz_sector, spin_j=target_spin_j, energy_target)
    energy_target == :energy_window && (commutator_linear_order > 0 || kms_psd_order > 0) &&
        throw(ArgumentError("Stationarity and KMS do not hold for arbitrary energy-window states."))
    HQ = _wf_exact(H)
    P = typeof(HQ)
    equalities = P[]
    inequalities = P[]
    psd_blocks = NamedTuple{(:kind,:matrix),Tuple{Symbol,Matrix{P}}}[]
    rules = NamedTuple{(:kind,:validity),Tuple{Symbol,Symbol}}[]
    rows = Set{W}()
    n = H.registry.nsites
    basis(d) = NC.basis_rows(W, NC.Pauli(), 1:n, d)
    families = isnothing(moment_rows) ? nothing :
        eltype(moment_rows) <: W ? (moment_rows,) : moment_rows

    if lightcone_depth > 0
        power = one(HQ)
        for _ in 1:lightcone_depth
            power *= HQ
            union!(rows, (w for (_, w) in NC.terms(power)))
        end
        push!(rules, (; kind=:lightcone, validity=:all_physical_states))
    end
    if !isempty(target.sectors)
        # All products of a single selected family are supported. Never put Q-q
        # between two unrestricted words: the right word can leave the sector.
        selected = isnothing(families) ? (basis(order),) : families
        domain = Set{W}((one(W),))
        for block in Iterators.flatten((selected, (rows,)))
            union!(domain, block, (last(NC.mul(u,v)) for u in block for v in block))
        end
        ordered = sort!(collect(domain))
        for s in target.sectors
            push!(equalities, s*s)
            for w in ordered
                _wf_split_zero!(equalities, s*_wf_word(HQ, w))
            end
        end
        push!(rules, (; kind=:sector_support, validity=:selected_sector))
    end

    rdm_positivity_k <= min(n,6) || throw(ArgumentError("rdm_positivity_k must not exceed the site count or six."))
    regions = isnothing(rdm_regions) ? _wf_regions(n, rdm_positivity_k) :
        [_wf_checked_region(n, region) for region in rdm_regions]
    all(sites -> rdm_positivity_k in (0,length(sites)), regions) ||
        throw(ArgumentError("RDM region does not match rdm_positivity_k."))
    if commutator_linear_order > 0
        d = commutator_linear_order
        words = if isnothing(families)
            basis(d)
        else
            selected = Set{W}(w for block in families for w in block if NC.degree(w) <= d)
            # Stationarity tests an operator, not a Gram row. For example,
            # a four-site RDM region contains punctured triples that are absent
            # from a contiguous degree-three moment basis.
            for sites in regions
                union!(selected, NC.basis_rows(W, NC.Pauli(), sites, d))
            end
            sort!(collect(selected))
        end
        for w in words
            u = _wf_word(HQ, w)
            _wf_nonzero!(equalities, im*(HQ*u-u*HQ))
        end
        push!(rules, (; kind=:stationarity, validity=:target_eigenstates))
    end
    if kms_psd_order > 0
        isnothing(kms_rows) && !isempty(target.sectors) &&
            throw(ArgumentError("Sector KMS requires explicit sector-preserving kms_rows."))
        perturbbases = if isnothing(kms_rows)
            selected = isnothing(families) ? (basis(kms_psd_order),) : families
            ([_wf_word(HQ,w) for w in block if NC.degree(w) <= kms_psd_order] for block in selected)
        else
            (P[_wf_exact(u) for u in kms_rows],)
        end
        for perturbations in perturbbases
            all(u -> u.registry === H.registry && NC.degree(u) <= kms_psd_order, perturbations) ||
                throw(ArgumentError("KMS perturbations must use the Hamiltonian registry and fit kms_psd_order."))
            all(u*s == s*u for u in perturbations for s in target.sectors) ||
                throw(ArgumentError("KMS perturbations must commute exactly with each selected sector operator."))
            if kms_frustration_guard
                perturbations = _wf_kms_rows(HQ, perturbations, _PauliQ(kms_null_tol))
            end
            isempty(perturbations) && continue
            M = Matrix{P}(undef, length(perturbations), length(perturbations))
            for j in eachindex(perturbations), i in 1:j
                u, v = perturbations[i], perturbations[j]
                uv = u'*v
                M[i,j] = u'*HQ*v - (HQ*uv + uv*HQ)*(1//2)
                M[j,i] = M[i,j]'
            end
            push!(psd_blocks, (;kind=:kms_psd, matrix=M))
        end
        kms_frustration_guard && push!(rules, (; kind=:kms_guard, validity=:principal_subblock))
        push!(rules, (; kind=:kms_psd, validity=isempty(target.sectors) ? :global_ground_states : :sector_ground_states))
    elseif !isnothing(kms_rows)
        throw(ArgumentError("kms_rows requires a positive kms_psd_order."))
    end

    for sites in regions
        D = _wf_density(HQ, sites)
        if isnothing(target.mz)
            push!(psd_blocks, (;kind=:rdm_positivity, matrix=D))
        else
            # Fixed global Mz implies the partial trace commutes with local Mz.
            # Keep BOTH signs of each nonzero charge and enforce zero cross blocks.
            charge = count_ones.(0:size(D,1)-1)
            for q in 0:length(sites)
                ids = findall(==(q), charge)
                push!(psd_blocks, (;kind=:rdm_positivity, matrix=D[ids,ids]))
            end
            for j in axes(D,2), i in 1:j-1
                charge[i] == charge[j] || _wf_split_zero!(equalities,D[i,j])
            end
        end
    end
    isempty(regions) || push!(rules, (;kind=:rdm_positivity, validity=:all_physical_states))
    isnothing(target.mz) || isempty(regions) || push!(rules, (;kind=:rdm_charge_blocks, validity=:selected_sector))

    rg_coarse_grain_depth > 0 && isempty(rg_maps) && throw(ArgumentError("RG needs explicit local isometries in rg_maps."))
    rg_coarse_grain_depth == 0 && !isempty(rg_maps) && throw(ArgumentError("rg_maps requires a positive RG depth."))
    for map in rg_maps
        sites = _wf_checked_region(n, map.sites)
        position = map.depth
        V = _PauliCQ.(map.isometry)
        1 <= position <= rg_coarse_grain_depth && size(V,1) == 1 << length(sites) &&
            0 < size(V,2) <= size(V,1) && V'*V == LinearAlgebra.I ||
            throw(ArgumentError("RG maps must provide a depth and an exact isometry on the designated sites."))
        D = _wf_density(HQ, sites)
        out = [_wf_exact(NC.polynomial(H.registry, (conj(V[a,i])*V[b,j]*D[a,b]
                    for a in axes(V,1) for b in axes(V,1) if !iszero(V[a,i]*V[b,j]))))
               for i in axes(V,2), j in axes(V,2)]
        push!(psd_blocks, (;kind=:rg_compression, matrix=out))
    end
    isempty(rg_maps) || push!(rules, (;kind=:rg_compression, validity=:all_physical_states))

    lo, hi = isnothing(energy_window) ? (nothing,nothing) : _wf_window_bounds(energy_window,HQ,target)
    energy_target == :energy_window && isnothing(lo) && throw(ArgumentError("An explicit energy-window target needs a checked energy window."))
    if !isnothing(lo)
        push!(inequalities, HQ-lo, hi-HQ)
        push!(rules, (;kind=:energy_window, validity=:energy_window_states))
    end
    isnothing(variance_center) == isnothing(variance_upper_bound) ||
        throw(ArgumentError("A second-moment cut needs both a center and an upper bound."))
    if !isnothing(variance_center)
        c, v = _PauliQ(variance_center), _PauliQ(variance_upper_bound)
        v >= 0 || throw(ArgumentError("Second-moment upper bound must be nonnegative."))
        ground_cut = energy_target != :energy_window
        ground_cut && (isnothing(lo) || max((lo-c)^2,(hi-c)^2) > v) &&
            throw(ArgumentError("The supplied window does not prove this variance cut for all target ground states."))
        push!(inequalities, v-(HQ-c)^2)
        push!(rules, (;kind=:centered_second_moment, validity=ground_cut ? :target_ground_states : :second_moment_window_states))
    end
    hamiltonian_dependent = !isempty(target.sectors) || commutator_linear_order > 0 || kms_psd_order > 0 ||
        lightcone_depth > 0 || !isnothing(lo) || !isnothing(variance_center)
    return (; rows=sort!(collect(rows)), equalities=unique(equalities), inequalities,
             psd_blocks, target, rules, hamiltonian_dependent, rational_sos_certify)
end

function _wf_checked_region(n, region)
    sites = NC._strengthening_sites(region)
    length(sites) <= 6 && last(sites) <= n || throw(ArgumentError("Local maps require 1–6 existing sites."))
    return sites
end

# Greedy, deterministic rank guard on the full ad_H image, not on its
# projection back to the perturbation basis. Remove exact nulls and small
# entries in the reduced image. Keeping actual source rows makes this a
# principal subblock even when a tolerance drops a nonzero direction.
function _wf_kms_rows(H::NC.Polynomial{W}, candidates, tol) where {W<:NC.PauliWord}
    kept = empty(candidates)
    pivots = Pair{W,typeof(H)}[]
    scale = tol*sum((abs(real(c)) + abs(imag(c)) for (c,_) in NC.terms(H));init=zero(_PauliQ))
    for u in candidates
        r = H*u-u*H
        for (w,p) in pivots
            for (c,v) in NC.terms(r)
                v == w || continue
                r -= c*p
                break
            end
        end
        best = scale
        pivot = nothing
        for (c,w) in NC.terms(r)
            mag = abs(real(c)) + abs(imag(c))
            mag > best || continue
            best, pivot = mag, (c,w)
        end
        isnothing(pivot) && continue
        c,w = pivot
        push!(pivots,w => inv(c)*r)
        push!(kept,u)
    end
    return kept
end

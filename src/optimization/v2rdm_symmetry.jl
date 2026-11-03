# Finite mode symmetries of real structured RDM data. The original maps and
# equality sector must be invariant; every discarded off-block coefficient is
# checked, rather than inferring correctness from an unchanged objective value.
function reduce_v2rdm_symmetry(
    linear::MomentLinearData{K,Float64,M},
    generators::AbstractVector{<:FermionicModePermutation};
    moment_generators::AbstractVector{<:FermionicModePermutation}=generators,
    atol::Float64=1e-9,
) where {T<:Signed,K<:Vector{T},M<:NormalMonomial{FermionicAlgebra,T}}
    domain = sort!(unique(Int(abs(i)) for k in linear.moments for i in k))
    spec = SymmetrySpec(; fermionic_generators=generators)
    spatial = _enumerate_symmetry_group(spec, domain)
    full_spec = SymmetrySpec(; fermionic_generators=[generators; moment_generators])
    group = _enumerate_symmetry_group(full_spec, domain)
    sw_group = _sw_fermionic_mode_permutation_group(spec, spatial, domain)
    F = LinearMomentForm{K,Float64}

    function image(g, key)
        sign, mono = _act_monomial(g, linear.key_to_monomial[key])
        adj = _v2rdm_adjoint_key(mono.word)
        word = key_lt(adj, mono.word) ? adj : mono.word
        haskey(linear.moment_index, word) || throw(ArgumentError("moment basis is not symmetry closed"))
        return sign, word
    end
    actions = [Dict(k => image(g, k) for k in linear.moments) for g in group]
    act(form, action) = F(action[k][2] => c * action[k][1] for (k,c) in form)
    checks = [_convert_fermionic_mode_permutation(Int, g) for g in [generators; moment_generators]]
    check_actions = [Dict(k => image(g, k) for k in linear.moments) for g in checks]
    for action in check_actions
        act(linear.objective_lin, action).terms == linear.objective_lin.terms ||
            throw(ArgumentError("objective is not invariant under the supplied mode permutation"))
    end

    # Equality rows need preserve their span, not each row separately.
    E = zeros(length(linear.zero_constraints), length(linear.moments))
    for (i,row) in enumerate(linear.zero_constraints), (k,c) in row.form
        E[i,linear.moment_index[k]] = c
    end
    decomposition = svd(E; full=false)
    rank = count(s -> s > atol, decomposition.S)
    rowspace = decomposition.V[:,1:rank]
    for action in check_actions, row in linear.zero_constraints
        v = zeros(length(linear.moments))
        for (k,c) in act(row.form, action)
            v[linear.moment_index[k]] = c
        end
        norm(v - rowspace * (rowspace' * v), Inf) <= atol * max(1, norm(v,Inf)) ||
            throw(ArgumentError("equality sector is not invariant under the supplied mode permutation"))
    end

    # Verify the complete cone set, including generators that exchange spin
    # blocks. Only the spatial subgroup must preserve each individual block.
    lookups = [Dict(row => i for (i,row) in enumerate(b.meta.row_labels)) for b in linear.psd_blocks_lin]
    for (g,action) in zip(checks,check_actions), block in linear.psd_blocks_lin
        images = [_act_monomial(g, row) for row in block.meta.row_labels]
        target = findfirst(eachindex(lookups)) do i
            candidate = linear.psd_blocks_lin[i]
            candidate.meta.origin.family == block.meta.origin.family &&
                candidate.size == block.size && all(haskey(lookups[i], m) for (_,m) in images)
        end
        isnothing(target) && throw(ArgumentError("positive block bases are not symmetry closed"))
        other = linear.psd_blocks_lin[target]
        indices = [lookups[target][m] for (_,m) in images]
        signs = first.(images)
        for j in 1:block.size, i in 1:block.size
            expected = other.entries[indices[i],indices[j]].terms
            actual = act(block.entries[i,j],action).terms
            sign = signs[i]*signs[j]
            # Multiplication by a fermionic permutation sign preserves the
            # canonical term order and cannot introduce zeros or duplicates.
            length(actual) == length(expected) && all(
                a.first == e.first && a.second == sign*e.second
                for (a,e) in zip(actual,expected)) ||
                throw(ArgumentError("positive matrix map is not covariant under the supplied mode permutation"))
        end
    end

    reconstruction = Dict{K,Tuple{Int,K}}()
    for key in linear.moments
        orbit = [a[key] for a in actions]
        representative = first(sort!(last.(orbit); lt=key_lt))
        signs = [s for (s,k) in orbit if k == representative]
        reconstruction[key] = (all(==(first(signs)), signs) ? first(signs) : 0, representative)
    end
    project(form) = F(reconstruction[k][2] => c*reconstruction[k][1] for (k,c) in form)
    moments = sort!(unique([k for (s,k) in values(reconstruction) if s != 0]); lt=key_lt)
    blocks = PSDBlockLin{K,Float64,M}[]
    transforms = Vector{Matrix{Float64}}[]
    max_offblock = 0.0
    for block in linear.psd_blocks_lin
        _check_basis_closure("structured RDM rows", block.meta.row_labels, spatial)
        # semisimple=true keeps a complete orthogonal basis, including all
        # copies of multidimensional irreducible representations.
        bases = if block.size == 1 || _basis_action_is_trivial(block.meta.row_labels, sw_group)
            [Matrix{Float64}(I,block.size,block.size)]
        else
            summands = SymbolicWedderburn.symmetry_adapted_basis(Float64, sw_group,
                NCFermionicModePermutationAction{T}(), block.meta.row_labels; semisimple=true)
            [Matrix(_wedderburn_row_basis(s)) for s in summands]
        end
        U = reduce(vcat, bases)
        size(U) == (block.size,block.size) && norm(U*U' - I,Inf) <= atol ||
            error("symmetry basis is not a complete orthogonal transformation")
        push!(transforms,bases)
        ranges = UnitRange{Int}[]
        start = 1
        for basis in bases
            push!(ranges, start:(start+size(basis,1)-1))
            start += size(basis,1)
        end
        coefficients = Dict{K,Matrix{Float64}}()
        for j in 1:block.size, i in 1:block.size, (key,c) in project(block.entries[i,j])
            A = get!(coefficients,key) do
                zeros(block.size,block.size)
            end
            A[i,j] += c
        end
        pairs = [[Pair{K,Float64}[] for _ in r, _ in r] for r in ranges]
        for (key,A) in coefficients
            B = U*A*U'
            for i in eachindex(ranges), j in (i+1):length(ranges)
                violation = maximum(abs, B[ranges[i],ranges[j]]; init=0.0)
                max_offblock = max(max_offblock,violation)
                violation <= atol || error("nonzero symmetry off-block coefficient")
            end
            for (p,r) in zip(pairs,ranges), j in eachindex(r), i in eachindex(r)
                c = (B[r[i],r[j]] + B[r[j],r[i]])/2
                abs(c) > atol/100 && push!(p[i,j],key=>c)
            end
        end
        for p in pairs
            n = size(p,1)
            meta = _default_block_meta(M,:PSD,length(blocks)+1,n)
            push!(blocks,PSDBlockLin{K,Float64,M}(n,map(F,p),meta))
        end
    end
    equalities = ScalarLinearConstraint{K,Float64}[]
    seen = Set{Any}()
    for row in linear.zero_constraints
        form = project(row.form)
        isempty(form) && continue
        normalized = F(k => c/first(form).second for (k,c) in form)
        signature = Tuple((Tuple(k),c) for (k,c) in normalized)
        signature in seen && continue
        push!(seen,signature)
        push!(equalities,ScalarLinearConstraint(normalized,row.kind,row.origin))
    end
    adjoints = Dict(k=>k for k in moments)
    pivots = _discover_linear_pivots(blocks,adjoints)
    reduced = MomentLinearData{K,Float64,M}(moments,Dict(k=>i for (i,k) in enumerate(moments)),
        linear.identity,Dict(k=>linear.key_to_monomial[k] for k in moments),adjoints,
        blocks,collect(1:length(blocks)),equalities,project(linear.objective_lin),
        pivots,_build_pivot_at(pivots),K[k for k in moments if !haskey(pivots,k)])
    return (; linear=reduced, reconstruction, transforms, max_offblock,
        spatial_group_order=length(spatial), moment_group_order=length(group))
end

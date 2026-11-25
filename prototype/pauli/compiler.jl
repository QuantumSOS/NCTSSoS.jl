"""Read-only compiler output, except for the private per-optimizer cache."""
struct CompiledSDPPlan{H,W,S,T,P,F,D,O,R}
    hamiltonian::H
    coefficients::Vector{Tuple{W,Float64}}
    symmetry::S
    table::T
    program::P
    standard::F
    sos::Union{Nothing,NC.SOSStandardForm}
    objective_map::D
    objective_data::Dict{Int,Float64}
    constraint_ids::BitSet
    moi_form::Symbol
    options::O
    reports::R
    provenance::BoundProvenance
    h_dependent::Bool
    certify_default::Bool
    cache::IdDict{Any,Any}
end

function _check_hamiltonian(H)
    H isa NC.Polynomial && H.registry.algebra isa NC.Pauli ||
        throw(ArgumentError("Supply a Hermitian NCTSSoS Pauli polynomial."))
    NC._checked_hermitian(H, "Hamiltonian")
    # The promoted block compiler stores Float64 coefficients. Reject input
    # that would change the polynomial before symmetry checks and assembly.
    all(isfinite(Float64(real(c))) && Float64(real(c)) == real(c) for (c,_) in NC.terms(H)) ||
        throw(ArgumentError("Hamiltonian coefficients must be exactly representable as finite Float64 values."))
    return H
end

function _symmetry_options(spec)
    isnothing(spec.drop_irreps) || throw(ArgumentError("Explicit irrep dropping has no supported lift witness."))
    spec.time_reversal_gauge == spec.anti_sign_realify || throw(ArgumentError(
        "The prototype selects one verified antiunitary route; set time_reversal_gauge and anti_sign_realify together."))
    spec.rotated_lie_detection == spec.cartan_u1_weyl == spec.su2_mz0_tree || throw(ArgumentError(
        "The prototype selects a spin route from the detected Lie algebra; set rotated_lie_detection, cartan_u1_weyl and su2_mz0_tree together."))
    return (; spec.spatial_mode, parity=spec.klein4_sign, antiunitary=spec.time_reversal_gauge,
        spin=spec.rotated_lie_detection, axis_permutations=spec.xyz_permutation)
end

_pauli_coefficients(H) = [(w, Float64(real(c))) for (c,w) in NC.terms(H)]

# Round the exact binary-input Pauli norm estimate down, not to nearest.
function _pauli_floor(H)
    bound = sum(NC.terms(H); init=zero(Rational{BigInt})) do (c,w)
        a = Rational{BigInt}(real(c))
        isone(w) ? a : -abs(a)
    end
    return Float64(bound, RoundDown)
end

function _append_strengthening!(program, table, strengthening)
    W = eltype(table.representatives.keys)
    sink = NC.OrbitFormSink(table, Dict{Int32,ComplexF64}())
    function append_form(p, coordinates=nothing; diagonal=false, round_constant_up=false)
        empty!(sink.terms)
        for (c,w) in NC.terms(p)
            a = round_constant_up && isone(w) ? Float64(real(c),RoundUp) : c
            u = NC._pack_word(W,w)
            if isnothing(coordinates)
                NC.push_term!(sink,a,u)
            else
                m = get!(coordinates,u) do
                    NC.orbit_coordinates!(table,u)
                end
                NC.push_term!(sink,a,m)
            end
        end
        for (id,c) in sort!(collect(sink.terms); by=first)
            diagonal && !isreal(c) && error("Nonreal scalar or diagonal strengthening entry.")
            iszero(c) || NC.push_term!(program.forms, c, id)
        end
        NC.close_form!(program.forms)
    end
    label = 0
    for (polynomials,kind) in ((strengthening.equalities,NC.ZERO),(strengthening.inequalities,NC.NONNEGATIVE))
        for p in polynomials
            append_form(p; diagonal=true, round_constant_up=kind==NC.NONNEGATIVE)
            NC.push_block!(program, kind, 1, label+=1)
        end
    end
    for (;kind,matrix) in strengthening.psd_blocks
        n = size(matrix,1)
        # Six-site RDMs repeatedly use the same 4^6 local Pauli words. Keep
        # a per-block cache, including zeros after their Ward heads register.
        coordinates = kind==:rdm_positivity && n==64 && table isa SquareOrbitMoments ?
            Dict{W,NC.OrbitCoordinates}() : nothing
        for j in 1:n, i in j:n
            append_form(matrix[i,j],coordinates; diagonal=i==j)
        end
        NC.push_block!(program, NC.HPSD, n, label+=1)
    end
    program.nvars = length(table.labels)
    NC.refresh_coverage!(table, program)
    return NC.validate(program)
end

function _append_strengthening!(program::NC.MomentProgram{Float64}, table, strengthening)
    # A lattice can already have real momentum blocks, but a local density
    # or KMS entry can be complex on the same real moment coordinates.
    complex_program = NC.MomentProgram{ComplexF64}(program.nvars)
    append!(complex_program.objective_ids,program.objective_ids)
    append!(complex_program.objective_coefs,program.objective_coefs)
    for b in program.blocks
        NC._copy_cone!(complex_program,program,b)
    end
    return _append_strengthening!(complex_program,table,strengthening)
end

function _real_workflow_program(program)
    all(b -> b.kind in (NC.PSD,NC.HPSD),program.blocks) && return NC.real_program(program)
    real_program = NC.MomentProgram{Float64}(program.nvars)
    append!(real_program.objective_ids,program.objective_ids)
    append!(real_program.objective_coefs,program.objective_coefs)
    for b in program.blocks
        if b.kind == NC.HPSD
            block_program = NC.MomentProgram{ComplexF64}(program.nvars)
            NC._copy_cone!(block_program,program,b)
            block_program_real = NC.real_program(block_program)
            NC._copy_cone!(real_program,block_program_real,only(block_program_real.blocks))
        else
            NC._copy_cone!(real_program,program,b)
        end
    end
    return NC.validate(real_program)
end

function _local_rdm_regions(support, k, shape)
    if support.mode === :contiguous_lattice
        if length(shape) == 1
            n = only(shape)
            return unique([sort!([mod1(i+j,n) for j in 0:k-1]) for i in 1:n])
        end
        Lx,Ly = shape
        site(x,y) = mod1(x,Lx) + (mod1(y,Ly)-1)*Lx
        return unique([sort!(unique!([site(x+dx,y+dy) for (dx,dy) in ((0,0),(1,0),(0,1),(1,1))[1:min(k,4)]]))
            for x in 1:Lx for y in 1:Ly])
    end
    # At most one bounded neighborhood per selected support region. Never
    # enumerate binomial(n,k) subsets to instantiate a scalable preset.
    return unique([collect(region[1:min(k,length(region))]) for region in support.regions if !isempty(region)])
end

function _coordinate_data(p, table)
    data = Dict{Int,Float64}()
    for (c,w) in NC.terms(p)
        pairs = NC.moment_terms(table,w)
        isnothing(pairs) && throw(ArgumentError("A requested polynomial uses unrepresented word $w; increase order or change supports."))
        for (id,factor) in pairs
            data[id] = get(data,id,0.0) + real(c*factor)
        end
    end
    filter!(pair -> !iszero(last(pair)), data)
    return data
end

"""
    compile_relaxation(H; order=2, preset, observable, symmetry, sparsity, strengthening, moi_form=:auto)

Compile Pauli supports and invariant moments, block reductions, post-build
term sparsity and sparse conic data. Use a preset or individual specs, not both.
The plan reports reductions that were actually applied; it is not a certificate.
"""
function compile_relaxation(H; order::Int=2, preset::Union{Nothing,Symbol}=nothing,
    observable=nothing, symmetry::Union{Nothing,PauliSymmetrySpec}=nothing,
    sparsity::Union{Nothing,PauliSparsitySpec}=nothing,
    strengthening::Union{Nothing,PauliStrengtheningSpec}=nothing,
    moi_form::Symbol=:auto, lattice_shape=nothing, max_rows::Int=1_000_000,
    rdm_regions=nothing, kms_rows=nothing, rg_maps=(), energy_target::Symbol=:target_ground,
    energy_window=nothing, variance_center=nothing, variance_upper_bound=nothing,
    wang_matched::Bool=false, wang_lso_span::Int=0)

    options = (; order, preset, observable, symmetry, sparsity, strengthening, moi_form,
        lattice_shape, max_rows, rdm_regions, kms_rows, rg_maps, energy_target,
        energy_window, variance_center, variance_upper_bound, wang_matched, wang_lso_span)
    if !isnothing(preset)
        wang_matched && throw(ArgumentError("Use individual specs with wang_matched."))
        all(isnothing, (symmetry,sparsity,strengthening)) || throw(ArgumentError("Use either a preset or individual specs, not both."))
        symmetry, sparsity, strengthening = resolve_workflow_preset(preset)
    else
        symmetry = something(symmetry, wang_matched ? PauliSymmetrySpec(rotated_lie_detection=false,
            cartan_u1_weyl=false,su2_mz0_tree=false) : PauliSymmetrySpec())
        sparsity = something(sparsity, wang_matched ? PauliSparsitySpec(ts_algo=:None,
            rewrite_complete_local_rdm=false) : PauliSparsitySpec())
        strengthening = something(strengthening, PauliStrengtheningSpec())
    end

    # Stage 0: validate inputs and choose the target/analysis policy. Symmetry
    # discovery consumes the actual selected family and target below.
    _check_hamiltonian(H)
    order >= 0 || throw(ArgumentError("order must be nonnegative."))
    moi_form in (:auto,:dualized,:primalized) || throw(ArgumentError("moi_form must be :auto, :dualized or :primalized."))
    if !isnothing(observable)
        _check_hamiltonian(observable)
        observable.registry === H.registry || throw(ArgumentError("Observable and Hamiltonian must use the same registry."))
    end
    sym_options = _symmetry_options(symmetry)
    n = H.registry.nsites
    shape = isnothing(lattice_shape) ? (n,) : Tuple(lattice_shape)
    named_shape = length(shape)==2 ? shape : nothing
    mode = sparsity.support_mode === :auto && sparsity.cs_algo === :MF ? :correlative_cliques : sparsity.support_mode
    direct = mode in (:graph_distance,:local_products)
    mode !== :correlative_cliques && sparsity.cs_algo !== :NoElimination &&
        throw(ArgumentError("cs_algo applies only to :correlative_cliques."))
    mode ∉ (:correlative_cliques,:graph_distance) && sparsity.power_graph_radius != 2 &&
        throw(ArgumentError("power_graph_radius applies only to :correlative_cliques or :graph_distance."))
    mode !== :correlative_cliques && sparsity.min_cluster_size != 4 &&
        throw(ArgumentError("min_cluster_size applies only to :correlative_cliques."))
    mode !== :local_products && sparsity.local_product_depth != 2 &&
        throw(ArgumentError("local_product_depth applies only to :local_products."))
    direct && sparsity.ts_algo ∉ (:None,:Maximal) &&
        throw(ArgumentError("Direct supports do not accept chordal term methods."))
    mode !== :contiguous_lattice && sparsity.two_body_range_r != 4 &&
        throw(ArgumentError("two_body_range_r applies only to :contiguous_lattice."))
    mode !== :correlative_cliques && sparsity.ts_algo !== :AggregateDualChordal && sparsity.merge_threshold != 8 &&
        throw(ArgumentError("merge_threshold applies only to correlative supports or :AggregateDualChordal."))
    !isnothing(lattice_shape) && mode !== :contiguous_lattice && throw(ArgumentError("lattice_shape requires contiguous lattice support."))
    if wang_matched
        mode in (:auto,:full,:graph_distance,:local_products) && sparsity.ts_algo === :None &&
            sparsity.ts_step == 0 && isnothing(lattice_shape) && !sym_options.spin ||
            throw(ArgumentError("Wang matching requires hint-free supports, no term decomposition and no continuous spin projection."))
        all(isnothing,(rdm_regions,kms_rows,energy_window,strengthening.target_mz_sector,
            strengthening.target_spin_j,variance_center,variance_upper_bound)) &&
            energy_target === :target_ground && strengthening.lightcone_depth == strengthening.rg_coarse_grain_depth == 0 ||
            throw(ArgumentError("Wang matching supports global ground states with discovered LSO, PSO and RDM pools only."))
        strengthening.commutator_linear_order in (0,3,5,7) && strengthening.kms_psd_order in 0:order &&
            strengthening.rdm_positivity_k in (0,8,9,10) && 0 <= wang_lso_span <= n ||
            throw(ArgumentError("Use Wang LSO degree 0/3/5/7, PSO up to order, RDM 0/8/9/10 and a valid span."))
    end
    wang_lso_span != 0 && !(wang_matched && strengthening.commutator_linear_order > 0) &&
        throw(ArgumentError("wang_lso_span requires Wang LSO."))
    strengthening.kms_psd_order == 0 && (!strengthening.kms_frustration_guard || strengthening.kms_null_tol != 1e-7) &&
        throw(ArgumentError("KMS guard options require a positive kms_psd_order."))
    !strengthening.kms_frustration_guard && strengthening.kms_null_tol != 1e-7 &&
        throw(ArgumentError("kms_null_tol requires the KMS frustration guard."))

    # Stage 1: select supports first so discovery is restricted to actions
    # that actually preserve these blocks. Light-cone rows are an augmentation.
    row_options = direct ? (;observable,sparsity.local_product_depth) : (;)
    support = workflow_support(H,order; support_mode=mode, sparsity.power_graph_radius,
        sparsity.cs_algo, sparsity.min_cluster_size, sparsity.merge_threshold, sparsity.two_body_range_r,
        lattice_shape=shape,max_rows,row_options...)
    W = eltype(H.data.words)
    rows = [W[only(NC.terms(p))[2] for p in block] for block in support.family]
    requested_k = strengthening.rdm_positivity_k
    k = isnothing(preset) ? requested_k : min(requested_k,n)
    local_rdm = !isnothing(preset) && k > 0 && isnothing(rdm_regions)
    rdm_regions = local_rdm ? _local_rdm_regions(support,k,shape) : rdm_regions
    !isnothing(rdm_regions) && k == 0 && throw(ArgumentError("rdm_regions requires rdm_positivity_k > 0."))
    extra = workflow_strengthening(H,order; moment_rows=rows,
        commutator_linear_order=wang_matched ? 0 : strengthening.commutator_linear_order,
        kms_psd_order=wang_matched ? 0 : strengthening.kms_psd_order,
        strengthening.kms_frustration_guard, strengthening.kms_null_tol, kms_rows,
        rdm_positivity_k=local_rdm || wang_matched ? 0 : k, rdm_regions, strengthening.lightcone_depth,
        strengthening.rg_coarse_grain_depth, rg_maps, strengthening.energy_upper_bound,
        strengthening.target_mz_sector, strengthening.target_spin_j, energy_target,
        energy_window, variance_center, variance_upper_bound,
        strengthening.rational_sos_certify)
    if !isempty(extra.rows)
        length(extra.rows) + sum(length,rows) <= max_rows || throw(ArgumentError("Light-cone augmentation exceeds max_rows."))
        # Union each supplemental row with the identity so its moment is
        # present, without replacing the original sparse family by a full one.
        push!(rows, sort!(unique!([one(W); extra.rows])))
        unique!(rows)
    end
    sectors = Tuple(NC.Sector(s,0) for s in extra.target.sectors)
    sym = workflow_symmetry(H; observable, family=rows, strengthen=sectors, lattice_shape=named_shape,
        uniform_parity_blocks=wang_matched,sym_options...)
    layout = nothing
    if wang_matched
        sym,layout = _wang_spatial_symmetry(H,observable,rows,sym)
    end
    blocks = [NC.RowBlock(block,[(1.0+0im,one(W))],NC.HPSD,0) for block in rows]
    objective = [(ComplexF64(c),w) for (c,w) in NC.terms(H)]

    # Stages 2–3 share NC's global invariant table with the block assembler.
    # For added forms, defer continuous elimination until *all* cones use it.
    augmented = !isempty(extra.equalities) || !isempty(extra.inequalities) || !isempty(extra.psd_blocks)
    if augmented
        finite = NC.Symmetry(sym.generators...; moment_only=sym.moment_only)
        program,table,sym_report = workflow_symmetry_program(blocks,objective,finite;
            sym_options.spatial_mode, lattice_shape=named_shape, realify=false)
        program = _append_strengthening!(program,table,extra)
        if !isempty(sym.continuous)
            program,table,continuous = NC._continuous_program(program,table,sym.continuous; algebra=H.registry.algebra)
            sym_report = merge(sym_report,(; continuous))
        end
    else
        program,table,sym_report = workflow_symmetry_program(blocks,objective,sym;
            sym_options.spatial_mode, lattice_shape=named_shape, realify=false)
    end
    wang = nothing
    if wang_matched
        program,wang = _append_wang!(program,table,H,objective,rows,sym,sym_report,strengthening,wang_lso_span,layout)
        !isempty(wang.pso.sides) && push!(extra.rules,(;kind=:kms_psd,validity=:global_ground_states))
        !isempty(wang.rdm.sides) && push!(extra.rules,(;kind=:rdm_positivity,validity=:all_physical_states))
        !isnothing(wang.lso) && wang.lso.rows > 0 && push!(extra.rules,(;kind=:stationarity,validity=:target_eigenstates))
    end
    density_report = nothing
    if sparsity.rewrite_complete_local_rdm
        candidates = workflow_local_density_candidates(H,support.family)
        program,density_report = workflow_local_density((;program,table),candidates)
    end
    before = [b.side for b in program.blocks]
    program = _real_workflow_program(program)
    sym_report = merge(sym_report,(;doubled=findall(i -> program.blocks[i].side != before[i],eachindex(before))))
    observable_data = isnothing(observable) ? Dict{Int,Float64}() : _coordinate_data(observable,table)

    # Stage 4: protect the observable and every non-moment cone.
    split_blocks = findall(b -> b.label == 0 && b.kind in (NC.PSD,NC.HPSD),program.blocks)
    program,ts_report = workflow_term_sparsity(program; sparsity.ts_algo, sparsity.ts_step,
        sparsity.merge_threshold, split_blocks, protected_ids=keys(observable_data))
    NC.refresh_coverage!(table,program)

    # Stage 5: lower shared affine terms once and retain a sparse objective map.
    standard = NC.moment_standard_form(program)
    all(isfinite,standard.c) && isfinite(standard.offset) || throw(ArgumentError("The projected objective overflows Float64."))
    form = moi_form === :auto ? (2length(standard.c) <= length(standard.b) ? :dualized : :primalized) : moi_form
    sos = form === :primalized ? NC.sos_standard_form(standard) : nothing
    coefficients = _pauli_coefficients(H)
    objective_map = Dict(w => Tuple(NC.moment_terms(table,w)) for (w,_) in coefficients)
    data = Dict(id => c for (id,c) in zip(Int.(program.objective_ids),program.objective_coefs))
    A = standard.A
    constraint_ids = BitSet((1,))
    for j in axes(A,2)
        any(!iszero(A.nzval[i]) for i in A.colptr[j]:A.colptr[j+1]-1) && push!(constraint_ids,j+1)
    end
    relaxations = (support.relation === :exact ? () : (support.mode,),
        ts_report.relation === :exact ? () : (ts_report.algorithm,))
    active_relaxations = (relaxations[1]...,relaxations[2]...)
    active_strengthenings = Tuple(unique(rule.kind for rule in extra.rules if rule.validity !== :selected_sector))
    sector_targets = (isnothing(extra.target.mz) ? () : (:mz,),isnothing(extra.target.spin_j) ? () : (:spin_j,))
    targets = (sector_targets[1]...,sector_targets[2]...)
    exact = isempty(active_relaxations) && isempty(active_strengthenings) && isempty(targets)
    witnesses = ((NC._has_symmetry(sym) ? (:averaging,) : ())...,
        (ts_report.zero_lift_audit === true ? (:zero_extension,) : ())...,
        (isempty(ts_report.auxiliary_ids) ? () : (:aggregate_chordal_lift,))...,
        (isnothing(density_report) || isempty(density_report.applied) ? () : (:local_density_equivalence,))...)
    provenance = BoundProvenance(:unavailable,false,exact,witnesses,
        active_relaxations,active_strengthenings,targets,_pauli_floor(H))
    reports = (; analysis=(; order, target=extra.target, rules=Tuple(extra.rules),wang,
            scalar_equalities=length(extra.equalities) + (isnothing(wang) || isnothing(wang.lso) ? 0 : wang.lso.rows),
            active_equality_rows=sum(c.length for c in standard.cones if c.kind === :zero && c.label!=0;init=0),
            rdm=(; requested_k,k,selection=k == 0 ? :disabled : wang_matched ? :wang_region : local_rdm ? :local_regions : isnothing(rdm_regions) ? :all_subsets : :explicit_regions,
                regions=k == 0 ? 0 : wang_matched ? 1 : isnothing(rdm_regions) ? binomial(big(n),k) : length(rdm_regions),
                count=wang_matched ? length(wang.rdm.sides) : count(b -> b.kind === :rdm_positivity,extra.psd_blocks))),
        support=(; support.mode,support.relation,support.regions,row_counts=length.(rows),
            density_rewrite=density_report,generation=get(support,:generation,nothing)),
        invariants=(; coordinates=program.nvars-1), symmetry=sym_report, term_sparsity=ts_report,
        conic=(; form, rows=size(standard.A,1), coordinates=size(standard.A,2), cones=length(standard.cones)),
        original_rows=rows)
    h_dependent = extra.hamiltonian_dependent || (wang_matched &&
        (strengthening.kms_psd_order > 0 || strengthening.commutator_linear_order > 0 || k > 0))
    return CompiledSDPPlan(H,coefficients,sym,table,program,standard,sos,objective_map,data,constraint_ids,
        form,options,reports,provenance,h_dependent,strengthening.rational_sos_certify,IdDict{Any,Any}())
end

function _guarded_objective(plan,H)
    incoming = _pauli_coefficients(H)
    incoming == plan.coefficients && return (:reused,plan.objective_data)
    first.(incoming) == first.(plan.coefficients) || return (:support_changed,nothing)
    plan.h_dependent && return (:hamiltonian_constraints_changed,nothing)
    sym = plan.symmetry
    original = Dict(incoming)
    for g in Iterators.flatten((sym.generators,sym.moment_only)), (w,c) in incoming
        phase,mapped = NC.act(g,w)
        get(original,mapped,0.0) == NC.scale(phase,c) || return (:symmetry_changed,nothing)
    end
    if !isempty(sym.continuous)
        try
            NC._check_derivation_polynomial(H,eltype(H.data.words),sym.continuous,"new Hamiltonian")
        catch error
            error isa ArgumentError || rethrow()
            return (:symmetry_changed,nothing)
        end
    end
    data = Dict{Int,Float64}()
    for (w,c) in incoming, (id,factor) in plan.objective_map[w]
        data[id] = get(data,id,0.0) + c * real(factor)
    end
    filter!(pair -> !iszero(last(pair)),data)
    all(isfinite,values(data)) || throw(ArgumentError("The projected objective overflows Float64."))
    all(in(plan.constraint_ids),keys(data)) || return (:term_sparsity_changed,nothing)
    return (:reused,data)
end

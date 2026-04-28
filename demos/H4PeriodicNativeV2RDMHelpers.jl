module H4PeriodicNativeV2RDMHelpers

using JuMP
using LinearAlgebra
using NCTSSoS

const MOI = JuMP.MOI

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: load_nk2_asset
include(joinpath(@__DIR__, "H4PeriodicPQGHelpers.jl"))
const PQGHelpers = H4PeriodicPQGHelpers

export load_native_problem_data,
       native_size_summary,
       build_native_jump_model,
       jump_model_summary,
       refinement_label

struct SpinOrbitalLabel
    k::Int
    spin::Symbol
    orb::Int
end

spin_weight(spin::Symbol) = spin === :up ? 1 :
                            spin === :dn ? -1 :
                            throw(ArgumentError("Unknown spin label $(spin)."))

pair_two_sz(m1, m2) = spin_weight(m1.spin) + spin_weight(m2.spin)
ph_transfer_two_sz(m1, m2) = spin_weight(m1.spin) - spin_weight(m2.spin)

function refinement_label(refinement::Symbol)
    refinement === :k_only && return "K-only"
    refinement === :spin_resolved && return "spin-resolved"
    throw(ArgumentError("Unknown refinement $(refinement). Use :k_only or :spin_resolved."))
end

function build_spin_orbitals(; nk::Int, norb::Int)
    modes = Vector{SpinOrbitalLabel}(undef, 2 * nk * norb)
    idx = 0
    for family in 0:2nk-1
        k = family ÷ 2
        spin = iseven(family) ? :up : :dn
        for orb in 1:norb
            idx += 1
            modes[idx] = SpinOrbitalLabel(k, spin, orb)
        end
    end
    return modes
end

function pair_sector_keys(nk::Int; refinement::Symbol)
    if refinement === :k_only
        return collect(0:nk-1)
    elseif refinement === :spin_resolved
        return [(K, two_Sz) for K in 0:nk-1 for two_Sz in (-2, 0, 2)]
    end
    throw(ArgumentError("Unknown refinement $(refinement). Use :k_only or :spin_resolved."))
end

function ph_sector_keys(nk::Int; refinement::Symbol)
    if refinement === :k_only
        return collect(0:nk-1)
    elseif refinement === :spin_resolved
        return [(Δk, two_Sz_transfer) for Δk in 0:nk-1 for two_Sz_transfer in (-2, 0, 2)]
    end
    throw(ArgumentError("Unknown refinement $(refinement). Use :k_only or :spin_resolved."))
end

function pair_sector_label(modes, p::Int, q::Int, nk::Int; refinement::Symbol)
    K = mod(modes[p].k + modes[q].k, nk)
    if refinement === :k_only
        return K
    elseif refinement === :spin_resolved
        return (K, pair_two_sz(modes[p], modes[q]))
    end
    throw(ArgumentError("Unknown refinement $(refinement). Use :k_only or :spin_resolved."))
end

function ph_sector_label(modes, p::Int, q::Int, nk::Int; refinement::Symbol)
    Δk = mod(modes[p].k - modes[q].k, nk)
    if refinement === :k_only
        return Δk
    elseif refinement === :spin_resolved
        return (Δk, ph_transfer_two_sz(modes[p], modes[q]))
    end
    throw(ArgumentError("Unknown refinement $(refinement). Use :k_only or :spin_resolved."))
end

function signed_token(x::Int)
    x < 0 && return "m$(abs(x))"
    x > 0 && return "p$(x)"
    return "0"
end

pair_sector_slug(sector::Int) = "K$(sector)"
pair_sector_slug(sector::Tuple{Int, Int}) = "K$(sector[1])_2Sz$(signed_token(sector[2]))"
ph_sector_slug(sector::Int) = "dK$(sector)"
ph_sector_slug(sector::Tuple{Int, Int}) = "dK$(sector[1])_d2Sz$(signed_token(sector[2]))"

pair_sector_display(sector::Int) = "K=$(sector)"
pair_sector_display(sector::Tuple{Int, Int}) = "K=$(sector[1]), 2S_z=$(sector[2])"
ph_sector_display(sector::Int) = "Δk=$(sector)"
ph_sector_display(sector::Tuple{Int, Int}) = "Δk=$(sector[1]), Δ(2S_z)=$(sector[2])"

function build_pair_partition(modes, nk::Int; refinement::Symbol)
    sector_keys = pair_sector_keys(nk; refinement)
    sector_type = typeof(first(sector_keys))
    by_sector = Dict(key => NTuple{2, Int}[] for key in sector_keys)
    index = Dict(key => Dict{NTuple{2, Int}, Int}() for key in sector_keys)
    sector_of = Dict{NTuple{2, Int}, sector_type}()

    for p in 1:length(modes)-1, q in p+1:length(modes)
        sector = pair_sector_label(modes, p, q, nk; refinement)
        pair = (p, q)
        push!(by_sector[sector], pair)
        index[sector][pair] = length(by_sector[sector])
        sector_of[pair] = sector
    end

    rows = refinement === :k_only ?
        [(; K = key, dim = length(by_sector[key])) for key in sector_keys] :
        [(; K = key[1], two_Sz = key[2], dim = length(by_sector[key])) for key in sector_keys]

    labels = [pair_sector_display(key) for key in sector_keys]

    return (; refinement,
             sector_keys,
             basis = by_sector,
             index,
             sector_of,
             block_sizes = [length(by_sector[key]) for key in sector_keys],
             rows,
             labels)
end

function build_ph_partition(modes, nk::Int; refinement::Symbol)
    sector_keys = ph_sector_keys(nk; refinement)
    sector_type = typeof(first(sector_keys))
    by_sector = Dict(key => Tuple{Int, Int}[] for key in sector_keys)
    index = Dict(key => Dict{Tuple{Int, Int}, Int}() for key in sector_keys)
    sector_of = Dict{Tuple{Int, Int}, sector_type}()

    for p in 1:length(modes), q in 1:length(modes)
        sector = ph_sector_label(modes, p, q, nk; refinement)
        pair = (p, q)
        push!(by_sector[sector], pair)
        index[sector][pair] = length(by_sector[sector])
        sector_of[pair] = sector
    end

    rows = refinement === :k_only ?
        [(; delta_k = key, dim = length(by_sector[key])) for key in sector_keys] :
        [(; delta_k = key[1], two_Sz_transfer = key[2], dim = length(by_sector[key])) for key in sector_keys]

    labels = [ph_sector_display(key) for key in sector_keys]

    return (; refinement,
             sector_keys,
             basis = by_sector,
             index,
             sector_of,
             block_sizes = [length(by_sector[key]) for key in sector_keys],
             rows,
             labels)
end

function build_nct_pair_partition(spin_orbitals, nk::Int; refinement::Symbol)
    NM = typeof(monomials(spin_orbitals[1].a * spin_orbitals[2].a)[1])
    sector_keys = pair_sector_keys(nk; refinement)
    blocks = Dict(key => NM[] for key in sector_keys)

    for p in 1:length(spin_orbitals)-1, q in p+1:length(spin_orbitals)
        sector = pair_sector_label(spin_orbitals, p, q, nk; refinement)
        push!(blocks[sector], monomials(spin_orbitals[p].a * spin_orbitals[q].a)[1])
    end

    return (; sector_keys,
             block_bases = [blocks[key] for key in sector_keys],
             block_sizes = [length(blocks[key]) for key in sector_keys])
end

canonical_pair(p::Int, q::Int) = p == q ? (0, (0, 0)) : p < q ? (1, (p, q)) : (-1, (q, p))

real_lift_psd_triangle_rows(n::Int) = (2n) * (2n + 1) ÷ 2

function native_partitions(data, refinement::Symbol)
    if refinement === :k_only
        return (; pair = data.pair_partition_k, ph = data.ph_partition_k)
    elseif refinement === :spin_resolved
        return (; pair = data.pair_partition_spin, ph = data.ph_partition_spin)
    end
    throw(ArgumentError("Unknown refinement $(refinement). Use :k_only or :spin_resolved."))
end

function load_native_problem_data()
    asset = load_nk2_asset()
    nk = asset.nk
    norb = asset.n_active_orb
    M = asset.total_spin_orbital_modes
    N = asset.total_active_electrons

    registry, vars, ham = PQGHelpers.build_h4_nk2_hamiltonian(asset.h1e, asset.eri; nk, norb)
    nct_spin_orbitals = PQGHelpers.build_spin_orbitals(vars; nk, norb)
    pqg_basis_data = PQGHelpers.build_pqg_basis(nct_spin_orbitals; nk, norb)

    d_nct_partition_k = build_nct_pair_partition(nct_spin_orbitals, nk; refinement = :k_only)
    D_block_bases = d_nct_partition_k.block_bases
    D_moment_basis = vcat(D_block_bases...)
    D_term_sparsity = NCTSSoS.TermSparsity{pqg_basis_data.NM}(copy(D_moment_basis), D_block_bases)

    d_nct_partition_spin = build_nct_pair_partition(nct_spin_orbitals, nk; refinement = :spin_resolved)
    D_block_bases_spin = d_nct_partition_spin.block_bases
    D_term_sparsity_spin = NCTSSoS.TermSparsity{pqg_basis_data.NM}(copy(D_moment_basis), D_block_bases_spin)

    modes = build_spin_orbitals(nk = nk, norb = norb)
    mode_id = Dict((m.k, m.spin, m.orb) => i for (i, m) in enumerate(modes))

    pair_partition_k = build_pair_partition(modes, nk; refinement = :k_only)
    pair_partition_spin = build_pair_partition(modes, nk; refinement = :spin_resolved)
    ph_partition_k = build_ph_partition(modes, nk; refinement = :k_only)
    ph_partition_spin = build_ph_partition(modes, nk; refinement = :spin_resolved)

    trace_targets = (
        D = N * (N - 1) ÷ 2,
        Q = (M - N) * (M - N - 1) ÷ 2,
        G = N * (M - N + 1),
    )

    return (; asset, nk, norb, M, N, registry, vars, ham,
             pqg_basis_data,
             D_block_bases,
             D_block_bases_spin,
             D_moment_basis,
             D_term_sparsity,
             D_term_sparsity_spin,
             modes,
             mode_id,
             pair_partition_k,
             pair_partition_spin,
             ph_partition_k,
             ph_partition_spin,
             pair_basis = pair_partition_k.basis,
             pair_index = pair_partition_k.index,
             ph_basis = ph_partition_k.basis,
             ph_index = ph_partition_k.index,
             D_block_sizes = pair_partition_k.block_sizes,
             Q_block_sizes = pair_partition_k.block_sizes,
             G_block_sizes = ph_partition_k.block_sizes,
             D_block_sizes_spin = pair_partition_spin.block_sizes,
             Q_block_sizes_spin = pair_partition_spin.block_sizes,
             G_block_sizes_spin = ph_partition_spin.block_sizes,
             trace_targets)
end

function native_size_summary(data; refinement::Symbol = :k_only,
                             include_trace_constraints::Bool = true)
    partitions = native_partitions(data, refinement)
    D = partitions.pair.block_sizes
    Q = copy(D)
    G = partitions.ph.block_sizes
    hermitian_block_sizes = vcat(D, Q, G)
    scalar_equalities = include_trace_constraints ? 5 : 2
    solver_psd_rows = [real_lift_psd_triangle_rows(n) for n in hermitian_block_sizes]
    d_complex_slots = sum(n^2 for n in D)

    return (
        refinement = refinement,
        pair_sector_rows = partitions.pair.rows,
        ph_sector_rows = partitions.ph.rows,
        pair_sector_labels = partitions.pair.labels,
        ph_sector_labels = partitions.ph.labels,
        d_only_block_sizes = D,
        q_block_sizes = Q,
        g_block_sizes = G,
        d_only_moment_vector_complex_variables = d_complex_slots,
        d_only_moment_vector_real_variables = 2 * d_complex_slots,
        native_hermitian_d_real_variables = d_complex_slots,
        jump_variables = d_complex_slots,
        jump_psd_blocks = length(hermitian_block_sizes),
        jump_scalar_equalities = scalar_equalities,
        jump_constraints_total = length(hermitian_block_sizes) + scalar_equalities,
        hermitian_block_sizes = hermitian_block_sizes,
        solver_real_lift_block_sizes = [2n for n in hermitian_block_sizes],
        solver_psd_rows = solver_psd_rows,
        solver_psd_total_rows = sum(solver_psd_rows),
        solver_total_rows = sum(solver_psd_rows) + scalar_equalities,
    )
end

function _make_model(optimizer_factory)
    if optimizer_factory === nothing
        return Model()
    end
    opt = MOI.instantiate(optimizer_factory; with_bridge_type = Float64, with_cache_type = Float64)
    return direct_model(opt)
end

function alloc_expr_matrix(n::Int, prototype)
    mat = Matrix{typeof(prototype)}(undef, n, n)
    for j in 1:n, i in 1:n
        mat[i, j] = zero(prototype)
    end
    return mat
end

function build_native_jump_model(; optimizer_factory = nothing,
                                 include_trace_constraints::Bool = true,
                                 refinement::Symbol = :k_only)
    data = load_native_problem_data()
    partitions = native_partitions(data, refinement)
    pair_partition = partitions.pair
    ph_partition = partitions.ph
    size_summary = native_size_summary(data;
        refinement = refinement,
        include_trace_constraints = include_trace_constraints)

    model = _make_model(optimizer_factory)

    total_t0 = time_ns()

    D_blocks = Dict{typeof(first(pair_partition.sector_keys)), Any}()
    for sector in pair_partition.sector_keys
        n = length(pair_partition.basis[sector])
        D_blocks[sector] = @variable(model,
            [1:n, 1:n] in HermitianPSDCone(),
            base_name = "D_$(pair_sector_slug(sector))")
    end
    prototype = D_blocks[pair_partition.sector_keys[1]][1, 1]

    function D_entry(p::Int, q::Int, r::Int, s::Int)
        sign1, pair1 = canonical_pair(p, q)
        sign2, pair2 = canonical_pair(r, s)
        if sign1 == 0 || sign2 == 0
            return zero(prototype)
        end
        sector1 = pair_partition.sector_of[pair1]
        sector2 = pair_partition.sector_of[pair2]
        sector1 == sector2 || return zero(prototype)
        i = pair_partition.index[sector1][pair1]
        j = pair_partition.index[sector2][pair2]
        return sign1 * sign2 * D_blocks[sector1][i, j]
    end

    t_one = time_ns()
    one_rdm = Matrix{typeof(prototype)}(undef, data.M, data.M)
    for p in 1:data.M, s in 1:data.M
        acc = zero(prototype)
        for q in 1:data.M
            acc += D_entry(p, q, s, q)
        end
        one_rdm[p, s] = acc / (data.N - 1)
    end
    one_rdm_time = (time_ns() - t_one) / 1e9

    t_q = time_ns()
    Q_blocks = Dict{typeof(first(pair_partition.sector_keys)), Any}()
    for sector in pair_partition.sector_keys
        basis = pair_partition.basis[sector]
        n = length(basis)
        mat = alloc_expr_matrix(n, prototype)
        for j in 1:n, i in 1:j
            p, q = basis[i]
            s, t = basis[j]
            expr = zero(prototype)
            expr += ((p == s && q == t) ? 1.0 + 0.0im : 0.0 + 0.0im)
            expr -= ((p == t && q == s) ? 1.0 + 0.0im : 0.0 + 0.0im)
            if q == t
                expr -= one_rdm[s, p]
            end
            if p == t
                expr += one_rdm[s, q]
            end
            if q == s
                expr += one_rdm[t, p]
            end
            if p == s
                expr -= one_rdm[t, q]
            end
            expr += D_entry(s, t, p, q)
            mat[i, j] = expr
        end
        Q_blocks[sector] = Hermitian(mat, :U)
        @constraint(model, Q_blocks[sector] in HermitianPSDCone())
    end
    q_time = (time_ns() - t_q) / 1e9

    t_g = time_ns()
    G_blocks = Dict{typeof(first(ph_partition.sector_keys)), Any}()
    for sector in ph_partition.sector_keys
        basis = ph_partition.basis[sector]
        n = length(basis)
        mat = alloc_expr_matrix(n, prototype)
        for j in 1:n, i in 1:j
            p, q = basis[i]
            s, t = basis[j]
            expr = zero(prototype)
            if q == t
                expr += one_rdm[p, s]
            end
            expr -= D_entry(p, t, s, q)
            mat[i, j] = expr
        end
        G_blocks[sector] = Hermitian(mat, :U)
        @constraint(model, G_blocks[sector] in HermitianPSDCone())
    end
    g_time = (time_ns() - t_g) / 1e9

    t_scalar = time_ns()
    up_modes = [i for (i, m) in enumerate(data.modes) if m.spin === :up]
    dn_modes = [i for (i, m) in enumerate(data.modes) if m.spin === :dn]
    @constraint(model, real(sum(one_rdm[p, p] for p in up_modes)) == data.N ÷ 2)
    @constraint(model, real(sum(one_rdm[p, p] for p in dn_modes)) == data.N ÷ 2)
    if include_trace_constraints
        @constraint(model,
            real(sum(D_blocks[sector][i, i]
                     for sector in pair_partition.sector_keys
                     for i in eachindex(pair_partition.basis[sector]))) == data.trace_targets.D)
        @constraint(model,
            real(sum(Q_blocks[sector][i, i]
                     for sector in pair_partition.sector_keys
                     for i in eachindex(pair_partition.basis[sector]))) == data.trace_targets.Q)
        @constraint(model,
            real(sum(G_blocks[sector][i, i]
                     for sector in ph_partition.sector_keys
                     for i in eachindex(ph_partition.basis[sector]))) == data.trace_targets.G)
    end
    scalar_time = (time_ns() - t_scalar) / 1e9

    t_obj = time_ns()
    obj_expr = let acc = zero(prototype)
        for k in 0:data.nk-1
            h_block = data.asset.h1e[k] / data.nk
            for p in 1:data.norb, s in 1:data.norb, spin in (:up, :dn)
                coeff = h_block[p, s]
                iszero(coeff) && continue
                ip = data.mode_id[(k, spin, p)]
                is = data.mode_id[(k, spin, s)]
                acc += coeff * one_rdm[ip, is]
            end
        end
        for ((k1, k2, k3, k4), Vblock) in data.asset.eri
            Vnorm = Vblock / data.nk^2
            for p in 1:data.norb, r in 1:data.norb, q in 1:data.norb, s in 1:data.norb
                coeff = Vnorm[p, r, q, s]
                iszero(coeff) && continue
                for (σ, τ) in ((:up, :up), (:up, :dn), (:dn, :up), (:dn, :dn))
                    i1 = data.mode_id[(k1, σ, p)]
                    i2 = data.mode_id[(k2, τ, q)]
                    i3 = data.mode_id[(k3, σ, r)]
                    i4 = data.mode_id[(k4, τ, s)]
                    acc += 0.5 * coeff * D_entry(i1, i2, i3, i4)
                end
            end
        end
        acc
    end
    @objective(model, Min, real(obj_expr))
    objective_time = (time_ns() - t_obj) / 1e9

    total_time = (time_ns() - total_t0) / 1e9

    build_times = (
        one_rdm = one_rdm_time,
        q_blocks = q_time,
        g_blocks = g_time,
        scalar_constraints = scalar_time,
        objective = objective_time,
        total = total_time,
    )

    return (; model, data, refinement, size_summary, build_times)
end

function jump_model_summary(model)
    constraint_counts = Dict(
        string(F) * " in " * string(S) => num_constraints(model, F, S)
        for (F, S) in list_of_constraint_types(model)
    )
    return (
        n_variables = num_variables(model),
        n_constraints = sum(values(constraint_counts)),
        constraint_counts = constraint_counts,
    )
end

end # module

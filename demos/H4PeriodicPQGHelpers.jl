module H4PeriodicPQGHelpers

using NCTSSoS

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: load_nk2_asset

export build_h4_pqg_problem_data,
       build_h4_nk2_hamiltonian,
       build_spin_orbitals,
       build_pqg_basis,
       build_route_moment_problem,
       evaluate_polynomial,
       moment_problem_stats,
       route_block_bases,
       route_label

function build_h4_nk2_hamiltonian(h1e, eri; nk::Int, norb::Int)
    registry, ((c_up_k0, c_up_k0_dag),
               (c_dn_k0, c_dn_k0_dag),
               (c_up_k1, c_up_k1_dag),
               (c_dn_k1, c_dn_k1_dag)) = create_fermionic_variables([
        ("c_up_k0", 1:norb),
        ("c_dn_k0", 1:norb),
        ("c_up_k1", 1:norb),
        ("c_dn_k1", 1:norb),
    ])

    ann = Dict(
        (0, :up) => c_up_k0,
        (0, :dn) => c_dn_k0,
        (1, :up) => c_up_k1,
        (1, :dn) => c_dn_k1,
    )
    dag = Dict(
        (0, :up) => c_up_k0_dag,
        (0, :dn) => c_dn_k0_dag,
        (1, :up) => c_up_k1_dag,
        (1, :dn) => c_dn_k1_dag,
    )
    spin_channels = ((:up, :up), (:up, :dn), (:dn, :up), (:dn, :dn))

    ham = (0.0 + 0.0im) * (c_up_k0_dag[1] * c_up_k0[1])

    for k in 0:nk-1
        h_block = h1e[k] / nk
        for p in 1:norb, s in 1:norb, sigma in (:up, :dn)
            coeff = h_block[p, s]
            iszero(coeff) && continue
            ham += coeff * dag[(k, sigma)][p] * ann[(k, sigma)][s]
        end
    end

    for (k1, k2, k3, k4) in sort!(collect(keys(eri)))
        block = eri[(k1, k2, k3, k4)] / nk^2
        for p in 1:norb, r in 1:norb, q in 1:norb, s in 1:norb
            coeff = block[p, r, q, s]
            iszero(coeff) && continue
            for (sigma, tau) in spin_channels
                ham += 0.5 * coeff *
                    dag[(k1, sigma)][p] * dag[(k2, tau)][q] *
                    ann[(k4, tau)][s] * ann[(k3, sigma)][r]
            end
        end
    end

    return registry,
           ((c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
            (c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag)),
           0.5 * (ham + adjoint(ham))
end

function build_spin_orbitals(vars; nk::Int, norb::Int)
    (c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
    (c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag) = vars

    ann = Dict(
        (0, :up) => c_up_k0,
        (0, :dn) => c_dn_k0,
        (1, :up) => c_up_k1,
        (1, :dn) => c_dn_k1,
    )
    dag = Dict(
        (0, :up) => c_up_k0_dag,
        (0, :dn) => c_dn_k0_dag,
        (1, :up) => c_up_k1_dag,
        (1, :dn) => c_dn_k1_dag,
    )

    spin_orbitals = NamedTuple[]
    for family in 0:2nk-1
        k = family ÷ 2
        spin = iseven(family) ? :up : :dn
        for orb in 1:norb
            push!(spin_orbitals,
                  (; a = ann[(k, spin)][orb],
                     adag = dag[(k, spin)][orb],
                     k,
                     spin,
                     orb,
                     family))
        end
    end
    return spin_orbitals
end

function mode_info(mode_id::Integer; norb::Int)
    family, orb0 = divrem(Int(mode_id) - 1, norb)
    k = family ÷ 2
    spin = iseven(family) ? :up : :dn
    return (; k, spin, orb = orb0 + 1, family)
end

function basis_label(m::NormalMonomial; nk::Int, norb::Int)
    ΔN = 0
    K = 0
    for idx in m.word
        info = mode_info(abs(Int(idx)); norb)
        sign = idx < 0 ? 1 : -1
        ΔN += sign
        K = mod(K + sign * info.k, nk)
    end
    return (; ΔN, K, deg = length(m.word))
end

function build_pqg_basis(spin_orbitals; nk::Int, norb::Int)
    NM = typeof(monomials(spin_orbitals[1].adag * spin_orbitals[1].a)[1])
    identity_mono = one(NM)

    D_basis = NM[]
    Q_basis = NM[]
    G_basis = NM[]

    M = length(spin_orbitals)
    for i in 1:M-1, j in i+1:M
        push!(D_basis, monomials(spin_orbitals[i].a * spin_orbitals[j].a)[1])
        push!(Q_basis, monomials(spin_orbitals[i].adag * spin_orbitals[j].adag)[1])
    end
    for i in 1:M, j in 1:M
        push!(G_basis, monomials(spin_orbitals[i].adag * spin_orbitals[j].a)[1])
    end

    pqg_basis = NM[identity_mono; D_basis; Q_basis; G_basis]
    unique!(pqg_basis)

    labels = [basis_label(m; nk, norb) for m in pqg_basis]
    sort_perm = sortperm(eachindex(pqg_basis); by = i -> (
        labels[i].ΔN,
        labels[i].K,
        labels[i].deg,
        string(pqg_basis[i]),
    ))
    pqg_basis = pqg_basis[sort_perm]
    labels = labels[sort_perm]

    sectors = sort!(unique([(l.ΔN, l.K) for l in labels]); by = s -> (s[1], s[2]))
    paper_blocks = Dict(s => NM[] for s in sectors)
    for (basis_elem, sector) in zip(pqg_basis, [(l.ΔN, l.K) for l in labels])
        push!(paper_blocks[sector], basis_elem)
    end

    sector_rows = [(; ΔN = s[1], K = s[2], dim = length(paper_blocks[s])) for s in sectors]
    paper_block_bases = [paper_blocks[s] for s in sectors]

    return (; NM,
             identity_mono,
             pqg_basis,
             D_basis,
             Q_basis,
             G_basis,
             labels,
             sectors,
             paper_blocks,
             paper_block_bases,
             sector_rows)
end

function build_trace_operators(spin_orbitals, ham)
    D_trace = zero(ham)
    Q_trace = zero(ham)
    M = length(spin_orbitals)

    for i in 1:M-1, j in i+1:M
        D_trace += spin_orbitals[j].adag * spin_orbitals[i].adag * spin_orbitals[i].a * spin_orbitals[j].a
        Q_trace += spin_orbitals[j].a * spin_orbitals[i].a * spin_orbitals[i].adag * spin_orbitals[j].adag
    end

    G_trace = zero(ham)
    for i in 1:M, j in 1:M
        G_trace += spin_orbitals[j].adag * spin_orbitals[i].a * spin_orbitals[i].adag * spin_orbitals[j].a
    end

    return (; D_trace, Q_trace, G_trace)
end

function build_particle_number_operators(vars; norb::Int)
    (c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
    (c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag) = vars

    n_up = 1.0 * sum(c_up_k0_dag[i] * c_up_k0[i] + c_up_k1_dag[i] * c_up_k1[i] for i in 1:norb)
    n_dn = 1.0 * sum(c_dn_k0_dag[i] * c_dn_k0[i] + c_dn_k1_dag[i] * c_dn_k1[i] for i in 1:norb)
    return (; n_up, n_dn, n_total = n_up + n_dn)
end

function build_trace_targets(n_modes::Int, n_electrons::Int)
    n_holes = n_modes - n_electrons
    return (
        n_modes = n_modes,
        n_electrons = n_electrons,
        n_holes = n_holes,
        D = n_electrons * (n_electrons - 1) ÷ 2,
        Q = n_holes * (n_holes - 1) ÷ 2,
        G = n_electrons * (n_holes + 1),
    )
end

function build_h4_pqg_problem_data()
    asset = load_nk2_asset()
    registry, vars, ham = build_h4_nk2_hamiltonian(asset.h1e, asset.eri;
        nk = asset.nk, norb = asset.n_active_orb)

    spin_orbitals = build_spin_orbitals(vars; nk = asset.nk, norb = asset.n_active_orb)
    basis_data = build_pqg_basis(spin_orbitals; nk = asset.nk, norb = asset.n_active_orb)
    trace_ops = build_trace_operators(spin_orbitals, ham)
    number_ops = build_particle_number_operators(vars; norb = asset.n_active_orb)
    trace_targets = build_trace_targets(asset.total_spin_orbital_modes, asset.total_active_electrons)

    I = one(ham)
    moment_eq_constraints = [
        number_ops.n_up - 4.0 * I,
        number_ops.n_dn - 4.0 * I,
        trace_ops.D_trace - trace_targets.D * I,
        trace_ops.Q_trace - trace_targets.Q * I,
        trace_ops.G_trace - trace_targets.G * I,
    ]

    pop = polyopt(ham, registry; moment_eq_constraints = moment_eq_constraints)
    config = SolverConfig(
        optimizer = nothing,
        moment_basis = basis_data.pqg_basis,
        cs_algo = NoElimination(),
        ts_algo = NoElimination(),
    )
    sparsity = compute_sparsity(pop, config)

    @assert length(sparsity.corr_sparsity.cliques) == 1
    @assert length(sparsity.cliques_term_sparsities[1]) == 1

    default_ts = sparsity.cliques_term_sparsities[1][1]
    paper_ts = NCTSSoS.TermSparsity{basis_data.NM}(
        default_ts.term_sparse_graph_supp,
        basis_data.paper_block_bases,
    )

    return (; asset,
             registry,
             vars,
             spin_orbitals,
             ham,
             basis_data,
             trace_ops,
             number_ops,
             trace_targets,
             pop,
             sparsity,
             fat_cliques_term_sparsities = sparsity.cliques_term_sparsities,
             paper_cliques_term_sparsities = [[paper_ts]])
end

route_label(route::Symbol) =
    route === :fat ? "single PQG block" :
    route === :paper ? "paper-style (ΔN, K) blocks" :
    throw(ArgumentError("Unknown route $(route)"))

function route_block_bases(data, route::Symbol)
    if route === :fat
        return data.fat_cliques_term_sparsities[1][1].block_bases
    elseif route === :paper
        return data.paper_cliques_term_sparsities[1][1].block_bases
    end
    throw(ArgumentError("Unknown route $(route)"))
end

function build_route_moment_problem(data, route::Symbol)
    cliques_term_sparsities = route === :fat ?
        data.fat_cliques_term_sparsities :
        route === :paper ?
            data.paper_cliques_term_sparsities :
            throw(ArgumentError("Unknown route $(route)"))
    return NCTSSoS.moment_relax(data.pop, data.sparsity.corr_sparsity, cliques_term_sparsities)
end

function moment_problem_stats(moment_problem)
    psd_block_sizes = sort(
        [size(mat, 1) for (cone, mat) in moment_problem.constraints if cone === :HPSD];
        rev = true,
    )
    n_zero_scalar_eqs = sum(length(mat) for (cone, mat) in moment_problem.constraints if cone === :Zero)
    return (
        n_psd_blocks = length(psd_block_sizes),
        max_psd_block = isempty(psd_block_sizes) ? 0 : maximum(psd_block_sizes),
        sum_sq = sum(d^2 for d in psd_block_sizes),
        sum_tri_slots = sum(d * (d + 1) ÷ 2 for d in psd_block_sizes),
        psd_block_sizes = psd_block_sizes,
        n_unique_moments = moment_problem.n_unique_moment_matrix_elements,
        total_basis_len = length(moment_problem.total_basis),
        n_zero_scalar_eqs = n_zero_scalar_eqs,
    )
end

function evaluate_polynomial(poly, monomap)
    value = zero(ComplexF64)
    for (coef, mono) in zip(coefficients(poly), monomials(poly))
        canon = NCTSSoS.symmetric_canon(NCTSSoS.expval(mono))
        value += coef * get(monomap, canon, 0.0 + 0.0im)
    end
    return value
end

end # module H4PeriodicPQGHelpers

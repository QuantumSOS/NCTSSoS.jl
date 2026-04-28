using Printf
using LinearAlgebra
using SparseArrays
using NCTSSoS
using JuMP
using COSMO
using Clarabel

const MOI = JuMP.MOI
BLAS.set_num_threads(1)

function flushln(args...)
    println(args...)
    flush(stdout)
    Base.Libc.flush_cstdio()
    return nothing
end
flush(stdout)
Base.Libc.flush_cstdio()

function format_bytes(b)
    return b < 1024        ? @sprintf("%d B",   b) :
           b < 1024^2      ? @sprintf("%.2f KiB", b / 1024) :
           b < 1024^3      ? @sprintf("%.2f MiB", b / 1024^2) :
                             @sprintf("%.2f GiB", b / 1024^3)
end

summarize_blocks(blocks) = length.(blocks)

function block_stats(sizes)
    isempty(sizes) && return (n_blocks = 0, max_size = 0, total_scalar_entries = 0, sum_sq = 0)
    return (
        n_blocks = length(sizes),
        max_size = maximum(sizes),
        total_scalar_entries = sum(s * (s + 1) ÷ 2 for s in sizes),
        sum_sq = sum(s^2 for s in sizes),
    )
end

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: fermionic_order2_basis_size,
                         load_nk2_asset

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
        (0, :up) => c_up_k0, (0, :dn) => c_dn_k0,
        (1, :up) => c_up_k1, (1, :dn) => c_dn_k1,
    )
    dag = Dict(
        (0, :up) => c_up_k0_dag, (0, :dn) => c_dn_k0_dag,
        (1, :up) => c_up_k1_dag, (1, :dn) => c_dn_k1_dag,
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

function build_h4_nk2_moment_problem(; report_prefix = nothing)
    prefix = isnothing(report_prefix) ? "" : string(report_prefix)

    flushln(prefix, "="^78)
    flushln(prefix, "H4 periodic (Nk=2) — order-2 V2RDM — term sparsity probe")
    flushln(prefix, "="^78)

    asset = load_nk2_asset()
    flushln(prefix,
            "Asset preflight: Nk = ", asset.nk,
            ", norb_per_k = ", asset.n_active_orb,
            ", total spin-orbital modes = ", asset.total_spin_orbital_modes,
            ", HF constant shift = ", round(asset.preflight["hf_constant_shift"]; digits = 10),
            " Ha")

    build_time = @timed begin
        registry, vars, ham = build_h4_nk2_hamiltonian(
            asset.h1e, asset.eri; nk = asset.nk, norb = asset.n_active_orb)
        (c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
        (c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag) = vars

        n_up = 1.0 * sum(
            c_up_k0_dag[i] * c_up_k0[i] + c_up_k1_dag[i] * c_up_k1[i]
            for i in 1:asset.n_active_orb)
        n_dn = 1.0 * sum(
            c_dn_k0_dag[i] * c_dn_k0[i] + c_dn_k1_dag[i] * c_dn_k1[i]
            for i in 1:asset.n_active_orb)

        polyopt(ham, registry;
                moment_eq_constraints = [n_up - 4.0 * one(ham),
                                         n_dn - 4.0 * one(ham)])
    end
    pop = build_time.value
    flushln(prefix, @sprintf("[stage] build polyopt        : %.2f s  |  alloc %s",
                             build_time.time, format_bytes(build_time.bytes)))

    flushln(prefix, "\n[stage] compute_sparsity     : running order-2 TS with MMD() ...")
    sparsity_time = @timed compute_sparsity(pop, SolverConfig(optimizer = nothing, order = 2, ts_algo = MMD()))
    sparsity = sparsity_time.value
    flushln(prefix, @sprintf("[stage] compute_sparsity     : %.2f s  |  alloc %s",
                             sparsity_time.time, format_bytes(sparsity_time.bytes)))

    blocks = sparsity.cliques_term_sparsities[1][1].block_bases
    sizes = summarize_blocks(blocks)
    stats = block_stats(sizes)
    dense_block_dim = Int(fermionic_order2_basis_size(asset.total_spin_orbital_modes))
    flushln(prefix, "\n-- block footprint at order 2 --")
    flushln(prefix, @sprintf("  dense baseline PSD block          : %d × %d", dense_block_dim, dense_block_dim))
    flushln(prefix, @sprintf("  TS PSD blocks                     : %d", stats.n_blocks))
    flushln(prefix, @sprintf("  TS largest block                  : %d × %d", stats.max_size, stats.max_size))
    flushln(prefix, @sprintf("  sum(block_size^2) proxy           : %d", stats.sum_sq))
    flushln(prefix, @sprintf("  sum of upper-tri PSD scalar slots : %d", stats.total_scalar_entries))

    flushln(prefix, "\n[stage] moment_relax         : assembling symbolic moment problem ...")
    mp_time = @timed NCTSSoS.moment_relax(pop, sparsity.corr_sparsity,
                                          sparsity.cliques_term_sparsities)
    moment_problem = mp_time.value
    flushln(prefix, @sprintf("[stage] moment_relax         : %.2f s  |  alloc %s",
                             mp_time.time, format_bytes(mp_time.bytes)))
    flushln(prefix, @sprintf("  unique moment variables           : %d",
                             moment_problem.n_unique_moment_matrix_elements))

    return (; asset, pop, sparsity, stats, dense_block_dim, moment_problem,
              build_time, sparsity_time, mp_time)
end

function build_basis_data(moment_problem)
    Mty = eltype(moment_problem.total_basis)
    basis = [NCTSSoS.symmetric_canon(NCTSSoS.expval(m)) for m in moment_problem.total_basis]
    NCTSSoS.sorted_unique!(basis)
    basis_to_idx = Dict(m => i for (i, m) in enumerate(basis))
    one_sym = NCTSSoS.symmetric_canon(NCTSSoS.expval(one(Mty)))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("Expected identity moment to be present in basis")
    return (; basis, basis_to_idx, idx_one)
end

@inline re_var_idx(idx::Int, n_basis::Int) = idx
@inline im_var_idx(idx::Int, n_basis::Int) = n_basis + idx

function add_poly_triplets!(I, J, V, row::Int, poly, basis_to_idx, n_basis::Int;
                            part::Symbol = :re, scale::Float64 = 1.0)
    for (coef, mono) in zip(coefficients(poly), monomials(poly))
        canon_mono = NCTSSoS.symmetric_canon(NCTSSoS.expval(mono))
        idx = get(basis_to_idx, canon_mono, 0)
        idx == 0 && continue

        c_re = real(coef)
        c_im = imag(coef)
        if part == :re
            iszero(c_re) || begin
                push!(I, row); push!(J, re_var_idx(idx, n_basis)); push!(V, scale * c_re)
            end
            iszero(c_im) || begin
                push!(I, row); push!(J, im_var_idx(idx, n_basis)); push!(V, -scale * c_im)
            end
        elseif part == :im
            iszero(c_im) || begin
                push!(I, row); push!(J, re_var_idx(idx, n_basis)); push!(V, scale * c_im)
            end
            iszero(c_re) || begin
                push!(I, row); push!(J, im_var_idx(idx, n_basis)); push!(V, scale * c_re)
            end
        else
            error("Unexpected polynomial part $(part)")
        end
    end
    return nothing
end

function direct_objective(moment_problem, basis_to_idx, n_basis::Int)
    I = Int[]
    J = Int[]
    V = Float64[]
    add_poly_triplets!(I, J, V, 1, moment_problem.objective, basis_to_idx, n_basis; part = :re, scale = 1.0)
    A = sparse(I, J, V, 1, 2 * n_basis)
    return vec(Matrix(A))
end

function build_direct_cosmo_constraints(moment_problem; report_every_fraction = 0.25)
    (; basis, basis_to_idx, idx_one) = build_basis_data(moment_problem)
    n_basis = length(basis)
    n_vars = 2 * n_basis
    constraints = Vector{COSMO.Constraint{Float64}}()

    total_rows = 0
    total_nnz = 0
    n_constraint_blocks = length(moment_problem.constraints)
    report_every = max(1, round(Int, n_constraint_blocks * report_every_fraction))
    flushln(@sprintf("  basis length / unique moments     : %d", n_basis))
    flushln(@sprintf("  direct COSMO variable count       : %d", n_vars))
    flushln(@sprintf("  constraint blocks to assemble     : %d", n_constraint_blocks))

    norm_A = spzeros(Float64, 2, n_vars)
    norm_A[1, re_var_idx(idx_one, n_basis)] = 1.0
    norm_A[2, im_var_idx(idx_one, n_basis)] = 1.0
    norm_b = [-1.0, 0.0]
    push!(constraints, COSMO.Constraint(norm_A, norm_b, COSMO.ZeroSet))
    total_rows += 2
    total_nnz += nnz(norm_A)

    for (k, (cone, mat)) in enumerate(moment_problem.constraints)
        dim = size(mat, 1)
        if cone == :Zero
            n_rows = 2 * dim * dim
            I = Int[]
            J = Int[]
            V = Float64[]
            b = zeros(Float64, n_rows)
            row = 0
            for j in 1:dim, i in 1:dim
                row += 1
                add_poly_triplets!(I, J, V, row, mat[i, j], basis_to_idx, n_basis; part = :re, scale = 1.0)
            end
            for j in 1:dim, i in 1:dim
                row += 1
                add_poly_triplets!(I, J, V, row, mat[i, j], basis_to_idx, n_basis; part = :im, scale = 1.0)
            end
            A = sparse(I, J, V, n_rows, n_vars)
            push!(constraints, COSMO.Constraint(A, b, COSMO.ZeroSet))
            total_rows += n_rows
            total_nnz += nnz(A)
        elseif cone == :HPSD
            emb_dim = 2 * dim
            n_rows = emb_dim * (emb_dim + 1) ÷ 2
            I = Int[]
            J = Int[]
            V = Float64[]
            b = zeros(Float64, n_rows)
            row = 0
            for col in 1:emb_dim
                for r in 1:col
                    row += 1
                    scale = r == col ? 1.0 : sqrt(2.0)
                    if col <= dim
                        add_poly_triplets!(I, J, V, row, mat[r, col], basis_to_idx, n_basis; part = :re, scale = scale)
                    elseif r <= dim
                        add_poly_triplets!(I, J, V, row, mat[r, col - dim], basis_to_idx, n_basis; part = :im, scale = -scale)
                    else
                        add_poly_triplets!(I, J, V, row, mat[r - dim, col - dim], basis_to_idx, n_basis; part = :re, scale = scale)
                    end
                end
            end
            A = sparse(I, J, V, n_rows, n_vars)
            push!(constraints, COSMO.Constraint(A, b, COSMO.PsdConeTriangle(n_rows)))
            total_rows += n_rows
            total_nnz += nnz(A)
        else
            error("Unexpected cone type $(cone) for complex moment problem")
        end
        if k % report_every == 0 || k == n_constraint_blocks
            flushln(@sprintf("    ... direct constraint block %d / %d", k, n_constraint_blocks))
        end
    end

    q = direct_objective(moment_problem, basis_to_idx, n_basis)
    return (; constraints, q, basis, basis_to_idx, idx_one, n_basis, n_vars, total_rows, total_nnz)
end

function build_primal_jump_model(moment_problem, optimizer;
                                 typed_exprs::Bool,
                                 report_every_fraction = 0.25)
    Cr = real(eltype(coefficients(moment_problem.objective)))
    model = GenericModel{Cr}()
    (; basis, basis_to_idx, idx_one) = build_basis_data(moment_problem)
    n_basis = length(basis)
    @variable(model, y_re[1:n_basis], set_string_name = false)
    @variable(model, y_im[1:n_basis], set_string_name = false)
    @constraint(model, y_re[idx_one] == 1)
    @constraint(model, y_im[idx_one] == 0)

    flushln(@sprintf("  basis length / unique moments     : %d", n_basis))
    n_constraint_blocks = length(moment_problem.constraints)
    flushln(@sprintf("  constraint blocks to assemble     : %d", n_constraint_blocks))
    report_every = max(1, round(Int, n_constraint_blocks * report_every_fraction))

    sample_expr = zero(Cr) * y_re[1]
    ExprT = typeof(sample_expr)

    for (k, (cone, mat)) in enumerate(moment_problem.constraints)
        dim = size(mat, 1)
        if typed_exprs
            mat_re = Matrix{ExprT}(undef, dim, dim)
            mat_im = Matrix{ExprT}(undef, dim, dim)
        else
            mat_re = Matrix{Any}(undef, dim, dim)
            mat_im = Matrix{Any}(undef, dim, dim)
        end
        for i in 1:dim, j in 1:dim
            re_expr, im_expr = NCTSSoS._substitute_complex_poly(mat[i, j], basis_to_idx, y_re, y_im)
            mat_re[i, j] = re_expr
            mat_im[i, j] = im_expr
        end
        if cone == :Zero
            @constraint(model, [mat_re[i, j] for i in 1:dim, j in 1:dim] .== 0)
            @constraint(model, [mat_im[i, j] for i in 1:dim, j in 1:dim] .== 0)
        elseif cone == :HPSD
            embedded = [
                [mat_re[i, j] for i in 1:dim, j in 1:dim]   [-mat_im[i, j] for i in 1:dim, j in 1:dim]
                [mat_im[i, j] for i in 1:dim, j in 1:dim]   [mat_re[i, j] for i in 1:dim, j in 1:dim]
            ]
            @constraint(model, embedded in PSDCone())
        else
            error("Unexpected cone type $(cone) for complex problem")
        end
        if k % report_every == 0 || k == n_constraint_blocks
            flushln(@sprintf("    ... constraint block %d / %d", k, n_constraint_blocks))
        end
    end

    obj_re, _ = NCTSSoS._substitute_complex_poly(moment_problem.objective, basis_to_idx, y_re, y_im)
    @objective(model, Min, obj_re)
    set_optimizer(model, optimizer)
    return (; model, n_basis)
end

function jump_model_summary(model)
    return (
        n_variables = num_variables(model),
        n_constraints = sum(num_constraints(model, F, S)
                            for (F, S) in list_of_constraint_types(model); init = 0),
    )
end

function energy_summary(asset, obj_active)
    hf_shift = asset.preflight["hf_constant_shift"]
    e_total = obj_active + hf_shift
    return (; hf_shift, e_total,
              figure_v2rdm = asset.preflight["figure_v2rdm_nk2"],
              figure_hf = asset.preflight["figure_hf_nk2"],
              figure_unc = asset.preflight["figure_uncertainty_ha"])
end

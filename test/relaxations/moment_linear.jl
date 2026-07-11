using Test, NCTSSoS, JuMP
import Clarabel

@testset "Moment linear form cache primitives" begin
    keys = [Int[2], Int[], Int[1, 2], Int[1]]
    @test sort(keys; lt=NCTSSoS.key_lt) == [Int[], Int[1], Int[1, 2], Int[2]]
    @test NCTSSoS.key_lt(Int[1], Int[1, 0])
    @test !NCTSSoS.key_lt(Int[1, 0], Int[1])
    @test NCTSSoS.key_isequal(Int[3, 4], Int[3, 4])

    form = NCTSSoS.LinearMomentForm{Vector{Int},Float64}([
        Int[2] => 1,
        Int[1] => 2,
        Int[2] => -0.25,
        Int[1] => -2,
        Int[3] => 0,
    ])

    @test collect(form) == [Int[2] => 0.75]
    @test length(form) == 1
    @test !isempty(form)
    @test NCTSSoS._assert_linear_moment_form_invariants(form) === nothing

    complex_form = NCTSSoS.LinearMomentForm{Vector{Int},ComplexF64}([
        Int[2] => 1,
        Int[1] => im,
        Int[1] => -im,
    ])
    @test collect(complex_form) == [Int[2] => 1.0 + 0.0im]
    @test NCTSSoS._within_abs_atol(3.0 + 4.0im, 5.0)
    @test !NCTSSoS._within_abs_atol(3.0 + 4.0im, 4.999)
    @test NCTSSoS._within_abs_atol(-1.0e-13, 1.0e-12)

    two_term_form = NCTSSoS.LinearMomentForm{Vector{Int},Float64}([
        Int[2] => 1.0,
        Int[1] => 2.0,
    ])
    @test collect(two_term_form) == [Int[1] => 2.0, Int[2] => 1.0]

    two_term_cancel = NCTSSoS.LinearMomentForm{Vector{Int},Float64}([
        Int[1] => 2.0,
        Int[1] => -2.0,
    ])
    @test isempty(two_term_cancel)

    acc = Dict{Vector{Int},Float64}(
        Int[3] => 0.0,
        Int[2] => 1.5,
        Int[1] => -2.0,
    )
    acc_form = NCTSSoS._linear_moment_form_from_accumulator!(acc; atol=1e-12)
    @test collect(acc_form) == [Int[1] => -2.0, Int[2] => 1.5]
    @test NCTSSoS._assert_linear_moment_form_invariants(acc_form) === nothing

    duplicate_unique_pairs = Pair{Vector{Int},Float64}[Int[2] => 1.0, Int[2] => 2.0]
    duplicate_form = NCTSSoS._linear_moment_form_from_unique_owned_pairs!(
        duplicate_unique_pairs,
    )
    @test collect(duplicate_form) == [Int[2] => 3.0]

    distinct_pairs = Pair{Vector{Int},Float64}[
        Int[3] => 3.0,
        Int[1] => 1.0,
        Int[2] => 0.0,
    ]
    distinct_form = NCTSSoS._linear_moment_form_from_distinct_owned_pairs!(
        distinct_pairs,
    )
    @test collect(distinct_form) == [Int[1] => 1.0, Int[3] => 3.0]

    left_subtract = NCTSSoS.LinearMomentForm{Vector{Int},Float64}([
        Int[1] => 1.0,
        Int[3] => 4.0,
        Int[5] => 1e-13,
    ])
    right_subtract = NCTSSoS.LinearMomentForm{Vector{Int},Float64}([
        Int[2] => 2.0,
        Int[3] => 1.0,
        Int[4] => 1.0,
    ])
    subtract_form = NCTSSoS._subtract_linear_forms(
        left_subtract,
        right_subtract;
        atol=1e-12,
    )
    @test collect(subtract_form) == [
        Int[1] => 1.0,
        Int[2] => -2.0,
        Int[3] => 3.0,
        Int[4] => -1.0,
    ]
    subtract_cancel = NCTSSoS._subtract_linear_forms(
        NCTSSoS.LinearMomentForm{Vector{Int},Float64}([Int[1] => 2.0]),
        NCTSSoS.LinearMomentForm{Vector{Int},Float64}([Int[1] => 2.0]);
        atol=1e-12,
    )
    @test isempty(subtract_cancel)

    cancel_acc = Dict{Vector{Int},Float64}()
    cancel_form = NCTSSoS.LinearMomentForm{Vector{Int},Float64}([Int[1] => 2.0])
    NCTSSoS._add_scaled_linear_form_terms!(
        cancel_acc,
        1.0,
        cancel_form,
        Float64;
        atol=1e-12,
    )
    NCTSSoS._add_scaled_linear_form_terms!(
        cancel_acc,
        -1.0,
        cancel_form,
        Float64;
        atol=1e-12,
    )
    @test isempty(cancel_acc)

    near_cancel_acc = Dict{Vector{Int},ComplexF64}()
    near_cancel_form = NCTSSoS.LinearMomentForm{Vector{Int},ComplexF64}([
        Int[1] => 1.0 + 0.0im,
    ])
    NCTSSoS._add_scaled_linear_form_terms!(
        near_cancel_acc,
        1.0 + 0.0im,
        near_cancel_form,
        ComplexF64;
        atol=1e-12,
    )
    NCTSSoS._add_scaled_linear_form_terms!(
        near_cancel_acc,
        -1.0 + 5.0e-13im,
        near_cancel_form,
        ComplexF64;
        atol=1e-12,
    )
    @test isempty(near_cancel_acc)

    entries = Matrix{NCTSSoS.LinearMomentForm{Vector{Int},Float64}}(undef, 3, 3)
    for j in 1:3, i in 1:3
        entries[i, j] = NCTSSoS.LinearMomentForm{Vector{Int},Float64}([
            Int[i] => i + j,
            Int[j + 3] => i - j,
        ])
    end
    entries[2, 3] = NCTSSoS.LinearMomentForm{Vector{Int},Float64}(Pair{Vector{Int},Float64}[])
    transform = [
        1.0 0.0 0.0
        0.0 1.0 1.0
        0.0 1.0 -1.0
    ]
    sparse_transform = NCTSSoS._pauli_sparse_transform_rows(transform; atol=1e-12)
    dense = NCTSSoS._transform_linear_block(entries, transform, transform; atol=1e-12)
    sparse = NCTSSoS._transform_linear_block(entries, sparse_transform, sparse_transform; atol=1e-12)
    @test [collect(dense[i, j]) for i in axes(dense, 1), j in axes(dense, 2)] ==
        [collect(sparse[i, j]) for i in axes(sparse, 1), j in axes(sparse, 2)]
    selected = NCTSSoS._sparse_transform_linear_block_entries(
        entries,
        sparse_transform,
        [1, 3],
        [2, 3];
        atol=1e-12,
    )
    @test [collect(selected[i, j]) for i in axes(selected, 1), j in axes(selected, 2)] ==
        [collect(dense[row, col]) for row in [1, 3], col in [2, 3]]
    empty_transformed = NCTSSoS._sparse_transformed_linear_form(
        entries,
        [(2, 1.0)],
        [(3, 1.0)],
        Float64;
        atol=1e-12,
    )
    @test isempty(empty_transformed)
    single_transformed = NCTSSoS._sparse_transformed_linear_form(
        entries,
        [(1, 2.0)],
        [(1, 0.5)],
        Float64;
        atol=1e-12,
    )
    @test collect(single_transformed) == [Int[1] => 2.0]

    indexed_entries, index_to_key = NCTSSoS._indexed_linear_block_entries(entries)
    @test index_to_key == sort(unique(vcat([
        first.(collect(form)) for form in entries
    ]...)); lt=NCTSSoS.key_lt)
    indexed_selected = NCTSSoS._sparse_transform_linear_block_entries(
        indexed_entries,
        sparse_transform,
        [1, 3],
        [2, 3];
        atol=1e-12,
    )
    rekeyed_selected = NCTSSoS._rekey_indexed_linear_block_entries(
        indexed_selected,
        index_to_key,
    )
    @test [collect(rekeyed_selected[i, j]) for i in axes(rekeyed_selected, 1), j in axes(rekeyed_selected, 2)] ==
        [collect(selected[i, j]) for i in axes(selected, 1), j in axes(selected, 2)]

    symmetric_entries = NCTSSoS._symmetrize_real_linear_block(entries)
    symmetric_sparse = NCTSSoS._transform_real_symmetric_linear_block(
        symmetric_entries,
        sparse_transform;
        atol=1e-12,
    )
    symmetric_reference = NCTSSoS._symmetrize_real_linear_block(
        NCTSSoS._transform_linear_block(symmetric_entries, sparse_transform, sparse_transform; atol=1e-12),
    )
    @test [collect(symmetric_sparse[i, j]) for i in axes(symmetric_sparse, 1), j in axes(symmetric_sparse, 2)] ==
        [collect(symmetric_reference[i, j]) for i in axes(symmetric_reference, 1), j in axes(symmetric_reference, 2)]
end

@testset "No-symmetry relaxation can build MomentLinearData directly" begin
    reg, (x,) = create_noncommutative_variables([("x", 1:2)])
    objective = 1.0 * x[1]^2 + 0.5 * x[2]^2
    pop = polyopt(
        objective,
        reg;
        eq_constraints=[1.0 * x[1] - x[2]],
        ineq_constraints=[1.0 - x[1]^2],
    )
    sparsity = compute_sparsity(pop, SolverConfig(optimizer=nothing, order=1))
    symbolic = NCTSSoS.moment_relax(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )
    direct = NCTSSoS._moment_relax_linear(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )

    direct_psd_entries = [
        collect(direct.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct.psd_blocks_lin)
        for i in 1:direct.psd_blocks_lin[block_idx].size
        for j in 1:direct.psd_blocks_lin[block_idx].size
    ]
    symbolic_psd_entries = [
        collect(symbolic.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic.linear.psd_blocks_lin)
        for i in 1:symbolic.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic.linear.psd_blocks_lin[block_idx].size
    ]
    direct_zero_entries = [(collect(z.form), z.kind, repr(z.origin)) for z in direct.zero_constraints]
    symbolic_zero_entries =
        [(collect(z.form), z.kind, repr(z.origin)) for z in symbolic.linear.zero_constraints]

    @test NCTSSoS.assert_moment_linear_data_invariants(direct) === nothing
    @test direct.moments == symbolic.linear.moments
    @test collect(direct.objective_lin) == collect(symbolic.linear.objective_lin)
    @test direct_psd_entries == symbolic_psd_entries
    @test direct_zero_entries == symbolic_zero_entries
    @test direct.psd_block_constraint_idx == symbolic.linear.psd_block_constraint_idx
end

@testset "No-symmetry direct MomentLinearData includes moment equality rows" begin
    reg, (x,) = create_noncommutative_variables([("x", 1:2)])
    objective = 1.0 * x[1]^2 + 0.5 * x[2]^2
    pop = polyopt(
        objective,
        reg;
        moment_eq_constraints=[1.0 * x[1] - x[2]],
    )
    sparsity = compute_sparsity(pop, SolverConfig(optimizer=nothing, order=1))
    symbolic = NCTSSoS.moment_relax(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )
    direct = NCTSSoS._moment_relax_linear(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )

    direct_zero_entries = [(collect(z.form), z.kind, repr(z.origin)) for z in direct.zero_constraints]
    symbolic_zero_entries =
        [(collect(z.form), z.kind, repr(z.origin)) for z in symbolic.linear.zero_constraints]

    @test NCTSSoS.assert_moment_linear_data_invariants(direct) === nothing
    @test direct.moments == symbolic.linear.moments
    @test collect(direct.objective_lin) == collect(symbolic.linear.objective_lin)
    @test direct_zero_entries == symbolic_zero_entries
    @test direct.psd_block_constraint_idx == symbolic.linear.psd_block_constraint_idx
end

@testset "No-symmetry direct MomentLinearData includes fermionic parity rows" begin
    reg, (c, c_dag) = create_fermionic_variables(1:1)
    n = c_dag[1] * c[1]
    pop = polyopt(1.0 * n, reg)
    sparsity = compute_sparsity(pop, SolverConfig(optimizer=nothing, order=1))
    symbolic = NCTSSoS.moment_relax(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )
    direct = NCTSSoS._moment_relax_linear(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )

    direct_psd_entries = [
        collect(direct.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct.psd_blocks_lin)
        for i in 1:direct.psd_blocks_lin[block_idx].size
        for j in 1:direct.psd_blocks_lin[block_idx].size
    ]
    symbolic_psd_entries = [
        collect(symbolic.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic.linear.psd_blocks_lin)
        for i in 1:symbolic.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic.linear.psd_blocks_lin[block_idx].size
    ]
    direct_zero_entries = [(collect(z.form), z.kind, repr(z.origin)) for z in direct.zero_constraints]
    symbolic_zero_entries =
        [(collect(z.form), z.kind, repr(z.origin)) for z in symbolic.linear.zero_constraints]

    @test !isempty(direct.zero_constraints)
    @test NCTSSoS.assert_moment_linear_data_invariants(direct) === nothing
    @test direct.moments == symbolic.linear.moments
    @test collect(direct.objective_lin) == collect(symbolic.linear.objective_lin)
    @test direct_psd_entries == symbolic_psd_entries
    @test direct_zero_entries == symbolic_zero_entries
    @test direct.psd_block_constraint_idx == symbolic.linear.psd_block_constraint_idx
end

@testset "No-symmetry direct MomentLinearData splits complex equality rows" begin
    reg, (σx, σy, _) = create_pauli_variables(1:1)
    objective = 1.0 * σx[1]
    eq = 1.0 * σx[1] + 1.0im * σy[1]
    pop = polyopt(objective, reg; eq_constraints=[eq])
    sparsity = compute_sparsity(pop, SolverConfig(optimizer=nothing, order=1))
    symbolic = NCTSSoS.moment_relax(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )
    direct = NCTSSoS._moment_relax_linear(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )

    direct_zero_entries = [(collect(z.form), z.kind, repr(z.origin)) for z in direct.zero_constraints]
    symbolic_zero_entries =
        [(collect(z.form), z.kind, repr(z.origin)) for z in symbolic.linear.zero_constraints]

    @test NCTSSoS.assert_moment_linear_data_invariants(direct) === nothing
    @test direct.moments == symbolic.linear.moments
    @test collect(direct.objective_lin) == collect(symbolic.linear.objective_lin)
    @test direct_zero_entries == symbolic_zero_entries
    @test direct.psd_block_constraint_idx == symbolic.linear.psd_block_constraint_idx
end

@testset "No-symmetry direct MomentLinearData splits complex moment equality rows" begin
    reg, (σx, σy, _) = create_pauli_variables(1:1)
    objective = 0.0 * σx[1]
    meq = 1.0 * σx[1] + 1.0im * σy[1]
    pop = polyopt(objective, reg; moment_eq_constraints=[meq])
    sparsity = compute_sparsity(pop, SolverConfig(optimizer=nothing, order=1))
    symbolic = NCTSSoS.moment_relax(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )
    direct = NCTSSoS._moment_relax_linear(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )

    direct_zero_entries = [(collect(z.form), z.kind, repr(z.origin)) for z in direct.zero_constraints]
    symbolic_zero_entries =
        [(collect(z.form), z.kind, repr(z.origin)) for z in symbolic.linear.zero_constraints]

    @test NCTSSoS.assert_moment_linear_data_invariants(direct) === nothing
    @test direct.moments == symbolic.linear.moments
    @test collect(direct.objective_lin) == collect(symbolic.linear.objective_lin)
    @test direct_zero_entries == symbolic_zero_entries
    @test direct.psd_block_constraint_idx == symbolic.linear.psd_block_constraint_idx
end

@testset "No-symmetry direct MomentLinearData matches bosonic PBW path" begin
    reg, (b, b_dag) = create_bosonic_variables(1:1)
    n = b_dag[1] * b[1]
    pop = polyopt(1.0 * n, reg)
    sparsity = compute_sparsity(pop, SolverConfig(optimizer=nothing, order=1))
    symbolic = NCTSSoS.moment_relax(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )
    direct = NCTSSoS._moment_relax_linear(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )

    direct_psd_entries = [
        collect(direct.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct.psd_blocks_lin)
        for i in 1:direct.psd_blocks_lin[block_idx].size
        for j in 1:direct.psd_blocks_lin[block_idx].size
    ]
    symbolic_psd_entries = [
        collect(symbolic.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic.linear.psd_blocks_lin)
        for i in 1:symbolic.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic.linear.psd_blocks_lin[block_idx].size
    ]

    @test NCTSSoS.assert_moment_linear_data_invariants(direct) === nothing
    @test direct.moments == symbolic.linear.moments
    @test collect(direct.objective_lin) == collect(symbolic.linear.objective_lin)
    @test direct_psd_entries == symbolic_psd_entries
    @test direct.zero_constraints == symbolic.linear.zero_constraints
    @test direct.psd_block_constraint_idx == symbolic.linear.psd_block_constraint_idx
end

@testset "MomentLinearBuilder stages owned forms before finalization" begin
    _, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    M = typeof(one_σ)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))
    C = ComplexF64

    builder = NCTSSoS.MomentLinearBuilder(K, C, M)
    x_key = NCTSSoS.register_moment!(builder, σx[1])
    scratch_key = copy(x_key)

    NCTSSoS.add_objective_terms!(
        builder,
        [scratch_key => 2.0 + 0im, scratch_key => -0.5 + 0im],
    )
    NCTSSoS.add_psd_block!(
        builder,
        :HPSD,
        reshape([[builder.identity => 1.0 + 0im]], 1, 1),
        NCTSSoS.BlockMeta(:HPSD, NCTSSoS.GlobalOrigin(1), [one_σ]);
        constraint_idx=1,
    )
    NCTSSoS.add_zero_constraint!(
        builder,
        [scratch_key => 1.0 + 0im],
        NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar),
    )

    form_key = copy(x_key)
    staged_form = NCTSSoS.LinearMomentForm{K,C}([form_key => 1.0 + 0im])
    NCTSSoS.add_psd_block!(
        builder,
        :HPSD,
        reshape([staged_form], 1, 1),
        NCTSSoS.BlockMeta(:HPSD, NCTSSoS.GlobalOrigin(2), [σx[1]]);
        constraint_idx=2,
    )

    zero_form_key = copy(x_key)
    staged_zero_form = NCTSSoS.LinearMomentForm{K,C}([zero_form_key => 0.25 + 0im])
    NCTSSoS.add_zero_constraint!(
        builder,
        staged_zero_form,
        NCTSSoS.ZeroMatrixOrigin(3, 1, 1, :scalar),
    )

    scratch_key[1] = 0xff
    form_key[1] = 0xfe
    zero_form_key[1] = 0xfd
    linear = NCTSSoS.finalize!(builder)

    @test NCTSSoS.assert_moment_linear_data_invariants(linear) === nothing
    @test collect(linear.objective_lin) == [x_key => 1.5 + 0im]
    @test haskey(linear.pivots, builder.identity)
    @test linear.psd_block_constraint_idx == [1, 2]
    @test collect(only(linear.psd_blocks_lin[2].entries)) == [x_key => 1.0 + 0im]
    @test length(linear.zero_constraints) == 2
    @test collect(linear.zero_constraints[1].form) == [x_key => 1.0 + 0im]
    @test collect(linear.zero_constraints[2].form) == [x_key => 0.25 + 0im]
    @test_throws ArgumentError NCTSSoS.add_objective_terms!(builder, [x_key => 1.0 + 0im])
    @test_throws ArgumentError NCTSSoS.finalize!(builder)
end

@testset "MomentLinearBuilder validates unregistered public zero rows" begin
    _, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    M = typeof(one_σ)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))

    builder = NCTSSoS.MomentLinearBuilder(K, Float64, M)
    x_key = convert(K, symmetric_canon(NCTSSoS.expval(σx[1])))
    NCTSSoS.add_psd_block!(
        builder,
        :PSD,
        reshape([[builder.identity => 1.0]], 1, 1),
        NCTSSoS.BlockMeta(:PSD, NCTSSoS.GlobalOrigin(1), [one_σ]);
        constraint_idx=1,
    )
    NCTSSoS.add_zero_constraint!(
        builder,
        [x_key => 1.0],
        NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar),
    )

    @test_throws ArgumentError NCTSSoS.finalize!(builder)
end

@testset "MomentLinearBuilder closes complex adjoint keys" begin
    _, (b, b_dag) = create_bosonic_variables(1:1)
    one_b = one(b[1])
    M = typeof(one_b)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_b)))
    C = ComplexF64

    builder = NCTSSoS.MomentLinearBuilder(K, C, M)
    b_key = NCTSSoS.register_moment!(builder, b[1])
    b_dag_key = symmetric_canon(NCTSSoS.expval(b_dag[1]))

    entries = Matrix{Vector{Pair{K,C}}}(undef, 2, 2)
    entries[1, 1] = [builder.identity => 1.0 + 0im]
    entries[1, 2] = [b_key => 1.0 + 0im]
    entries[2, 1] = [b_dag_key => 1.0 + 0im]
    entries[2, 2] = [builder.identity => 1.0 + 0im]

    NCTSSoS.add_psd_block!(
        builder,
        :HPSD,
        entries,
        NCTSSoS.BlockMeta(:HPSD, NCTSSoS.GlobalOrigin(1), [one_b, b[1]]);
        constraint_idx=1,
    )
    linear = NCTSSoS.finalize!(builder)

    @test NCTSSoS.assert_moment_linear_data_invariants(linear) === nothing
    @test haskey(linear.moment_index, b_key)
    @test haskey(linear.moment_index, b_dag_key)
    @test haskey(linear.pivots, b_key)
    @test haskey(linear.pivots, b_dag_key)
    @test linear.pivots[b_dag_key].adjoint
    @test linear.adjoint_key[b_key] == b_dag_key
end

@testset "MomentProblem linear cache construction" begin
    reg, (b, b_dag) = create_bosonic_variables(1:1)
    one_b = one(b[1])
    objective = zero(typeof(1.0 * one_b))
    P = typeof(objective)

    block = Matrix{P}(undef, 2, 2)
    block[1, 1] = 1.0 * one_b
    block[1, 2] = 1.0 * b[1]
    block[2, 1] = 1.0 * b_dag[1]
    block[2, 2] = 1.0 * one_b

    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block)],
        [one_b, b[1], b_dag[1]],
        3,
    )

    key(m) = symmetric_canon(NCTSSoS.expval(m))
    @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
    @test length(mp.linear.psd_blocks_lin) == 1
    @test mp.linear.psd_block_constraint_idx == [1]
    @test isempty(mp.linear.free_keys)
    @test haskey(mp.linear.pivots, key(b[1]))
    @test haskey(mp.linear.pivots, key(b_dag[1]))
    @test mp.linear.pivots[key(b_dag[1])].adjoint
    @test length(mp.linear.pivot_at[(1, 1, 2)]) == 2
end

@testset "MomentProblem linear cache keeps zero-only keys free" begin
    reg, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    objective = zero(typeof(1.0 * one_σ))
    P = typeof(objective)

    block = reshape([1.0 * one_σ], 1, 1)
    zero_mat = reshape([1.0 * σx[1]], 1, 1)
    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block), (:Zero, zero_mat)],
        [one_σ, σx[1]],
        1,
    )

    σx_key = symmetric_canon(NCTSSoS.expval(σx[1]))
    @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
    @test σx_key in mp.linear.free_keys
    @test !haskey(mp.linear.pivots, σx_key)
    @test length(mp.linear.zero_constraints) == 1

    int_identity = Polynomial([(1, one_σ)])
    int_zero = zero(typeof(int_identity))
    int_block = reshape([int_identity], 1, 1)
    int_zero_mat = reshape([Polynomial([(1, σx[1])])], 1, 1)
    int_mp = NCTSSoS.MomentProblem(
        int_zero,
        [(:HPSD, int_block), (:Zero, int_zero_mat)],
        [one_σ, σx[1]],
        1,
    )
    @test NCTSSoS.assert_moment_linear_data_invariants(int_mp.linear, int_mp.constraints) === nothing
end

@testset "MomentProblem linear cache trusts total_basis support" begin
    reg, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    objective = zero(typeof(1.0 * one_σ))
    P = typeof(objective)

    block = reshape([1.0 * one_σ], 1, 1)
    zero_mat = reshape([1.0 * σx[1]], 1, 1)

    @test_throws ArgumentError NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block), (:Zero, zero_mat)],
        [one_σ],
        1,
    )
end

@testset "_add_parity_constraints! appends linear rows in place" begin
    reg, (c, c_dag) = create_fermionic_variables(1:1)
    one_c = one(c[1])
    objective = zero(typeof((1.0 + 0im) * one_c))
    P = typeof(objective)

    block = Matrix{P}(undef, 2, 2)
    block[1, 1] = (1.0 + 0im) * one_c
    block[1, 2] = (1.0 + 0im) * c[1]
    block[2, 1] = (1.0 + 0im) * c_dag[1]
    block[2, 2] = (1.0 + 0im) * one_c
    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block)],
        [one_c, c[1], c_dag[1]],
        3,
    )

    before_linear = mp.linear
    before_constraints = length(mp.constraints)
    before_zero_forms = length(mp.linear.zero_constraints)
    NCTSSoS._add_parity_constraints!(mp)

    @test mp.linear === before_linear
    @test length(mp.constraints) > before_constraints
    @test length(mp.linear.zero_constraints) > before_zero_forms
    @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
end

@testset "_add_moment_eq_constraints! appends linear rows in place" begin
    reg, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    objective = zero(typeof(1.0 * one_σ))
    P = typeof(objective)

    block = Matrix{P}(undef, 2, 2)
    block[1, 1] = 1.0 * one_σ
    block[1, 2] = 1.0 * σx[1]
    block[2, 1] = 1.0 * σx[1]
    block[2, 2] = 1.0 * one_σ
    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block)],
        [one_σ, σx[1]],
        2,
    )
    pop = polyopt(objective, reg; moment_eq_constraints=[1.0 * σx[1]])

    before_linear = mp.linear
    before_constraints = length(mp.constraints)
    before_zero_forms = length(mp.linear.zero_constraints)
    NCTSSoS._add_moment_eq_constraints!(mp, pop, [one_σ, σx[1]], [0, 1])

    @test mp.linear === before_linear
    @test length(mp.constraints) > before_constraints
    @test length(mp.linear.zero_constraints) > before_zero_forms
    @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
end

@testset "PSD-block lowering consumes cached linear zero rows" begin
    MOI = JuMP.MOI
    reg, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    objective = zero(typeof(1.0 * one_σ))
    P = typeof(objective)

    block = Matrix{P}(undef, 2, 2)
    block[1, 1] = 1.0 * one_σ
    block[1, 2] = 1.0 * σx[1]
    block[2, 1] = 1.0 * σx[1]
    block[2, 2] = 1.0 * one_σ
    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block)],
        [one_σ, σx[1]],
        2,
    )

    base_model, _ = build_jump_model(mp; formulation=:psd_blocks, representation=:complex)
    base_backend = JuMP.backend(base_model)
    base_real_equalities = MOI.get(
        base_backend,
        MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}(),
    )

    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))
    x_key = symmetric_canon(NCTSSoS.expval(σx[1]))
    form = NCTSSoS.LinearMomentForm{K,ComplexF64}([x_key => 1.0 + 0im])
    push!(
        mp.linear.zero_constraints,
        NCTSSoS.ScalarLinearConstraint(form, :zero, NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar)),
    )

    model, _ = build_jump_model(mp; formulation=:psd_blocks, representation=:complex)
    backend = JuMP.backend(model)
    @test MOI.get(
        backend,
        MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}(),
    ) == base_real_equalities + 1
end

@testset "MomentLinearData lowers directly to PSD-block JuMP model" begin
    MOI = JuMP.MOI
    _, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    M = typeof(one_σ)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))
    C = ComplexF64

    builder = NCTSSoS.MomentLinearBuilder(K, C, M)
    x_key = NCTSSoS.register_moment!(builder, σx[1])
    NCTSSoS.add_objective_terms!(builder, [x_key => -1.0 + 0im])
    NCTSSoS.add_psd_block!(
        builder,
        :HPSD,
        reshape(
            [
                [builder.identity => 1.0 + 0im],
                [x_key => 1.0 + 0im],
                [x_key => 1.0 + 0im],
                [builder.identity => 1.0 + 0im],
            ],
            2,
            2,
        ),
        NCTSSoS.BlockMeta(:HPSD, NCTSSoS.GlobalOrigin(1), [one_σ, σx[1]]);
        constraint_idx=1,
    )
    NCTSSoS.add_zero_constraint!(
        builder,
        [x_key => 1.0 + 0im],
        NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar),
    )

    linear = NCTSSoS.finalize!(builder)
    model, extract = build_jump_model(linear; formulation=:psd_blocks, representation=:complex)
    backend = JuMP.backend(model)

    @test extract isa Function
    @test MOI.get(
        backend,
        MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}(),
    ) == 3

    real_builder = NCTSSoS.MomentLinearBuilder(K, Float64, M)
    real_x_key = NCTSSoS.register_moment!(real_builder, σx[1])
    NCTSSoS.add_objective_terms!(real_builder, [real_x_key => -1.0])
    NCTSSoS.add_psd_block!(
        real_builder,
        :PSD,
        reshape(
            [
                [real_builder.identity => 1.0],
                [real_x_key => 1.0],
                [real_x_key => 1.0],
                [real_builder.identity => 1.0],
            ],
            2,
            2,
        ),
        NCTSSoS.BlockMeta(:PSD, NCTSSoS.GlobalOrigin(1), [one_σ, σx[1]]);
        constraint_idx=1,
    )
    NCTSSoS.add_zero_constraint!(
        real_builder,
        [real_x_key => 1.0],
        NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar),
    )

    real_linear = NCTSSoS.finalize!(real_builder)
    real_model, real_extract = build_jump_model(
        real_linear;
        formulation=:psd_blocks,
        representation=:real,
    )
    real_types = JuMP.list_of_constraint_types(real_model)
    @test real_extract isa Function
    @test any(
        pair -> last(pair) == MOI.PositiveSemidefiniteConeTriangle,
        real_types,
    )

    orphan_builder = NCTSSoS.MomentLinearBuilder(K, Float64, M)
    orphan_key = NCTSSoS.register_moment!(orphan_builder, σx[1])
    NCTSSoS.add_objective_terms!(orphan_builder, [orphan_key => 1.0])
    NCTSSoS.add_psd_block!(
        orphan_builder,
        :PSD,
        reshape([[orphan_builder.identity => 1.0]], 1, 1),
        NCTSSoS.BlockMeta(:PSD, NCTSSoS.GlobalOrigin(1), [one_σ]);
        constraint_idx=1,
    )
    orphan_linear = NCTSSoS.finalize!(orphan_builder)
    @test !isempty(orphan_linear.free_keys)
    @test_throws ArgumentError build_jump_model(
        orphan_linear;
        formulation=:psd_blocks,
        representation=:real,
    )
    orphan_model, _ = build_jump_model(
        orphan_linear;
        formulation=:psd_blocks,
        representation=:real,
        orphan_policy=:free_variables,
    )
    @test JuMP.num_variables(orphan_model) > 1
end

@testset "MomentLinearData lowers directly to moment-variable JuMP models" begin
    MOI = JuMP.MOI
    _, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    M = typeof(one_σ)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))

    real_builder = NCTSSoS.MomentLinearBuilder(K, Float64, M)
    real_x_key = NCTSSoS.register_moment!(real_builder, σx[1])
    NCTSSoS.add_objective_terms!(real_builder, [real_x_key => -1.0])
    NCTSSoS.add_psd_block!(
        real_builder,
        :PSD,
        reshape(
            [
                [real_builder.identity => 1.0],
                [real_x_key => 1.0],
                [real_x_key => 1.0],
                [real_builder.identity => 1.0],
            ],
            2,
            2,
        ),
        NCTSSoS.BlockMeta(:PSD, NCTSSoS.GlobalOrigin(1), [one_σ, σx[1]]);
        constraint_idx=1,
    )
    NCTSSoS.add_zero_constraint!(
        real_builder,
        [real_x_key => 1.0],
        NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar),
    )
    real_linear = NCTSSoS.finalize!(real_builder)
    real_model, real_extract = build_jump_model(
        real_linear;
        formulation=:moment_variables,
        representation=:real,
    )
    real_backend = JuMP.backend(real_model)

    @test real_extract isa Function
    @test JuMP.num_variables(real_model) == 2
    @test MOI.get(
        real_backend,
        MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}(),
    ) == 2

    complex_builder = NCTSSoS.MomentLinearBuilder(K, ComplexF64, M)
    complex_x_key = NCTSSoS.register_moment!(complex_builder, σx[1])
    NCTSSoS.add_objective_terms!(complex_builder, [complex_x_key => -1.0 + 0im])
    NCTSSoS.add_psd_block!(
        complex_builder,
        :HPSD,
        reshape(
            [
                [complex_builder.identity => 1.0 + 0im],
                [complex_x_key => 1.0 + 0im],
                [complex_x_key => 1.0 + 0im],
                [complex_builder.identity => 1.0 + 0im],
            ],
            2,
            2,
        ),
        NCTSSoS.BlockMeta(:HPSD, NCTSSoS.GlobalOrigin(1), [one_σ, σx[1]]);
        constraint_idx=1,
    )
    NCTSSoS.add_zero_constraint!(
        complex_builder,
        [complex_x_key => 1.0 + 0im],
        NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar),
    )
    complex_linear = NCTSSoS.finalize!(complex_builder)
    complex_model, complex_extract = build_jump_model(
        complex_linear;
        formulation=:moment_variables,
        representation=:real,
    )
    complex_backend = JuMP.backend(complex_model)

    @test complex_extract isa Function
    @test JuMP.num_variables(complex_model) == 4
    @test MOI.get(
        complex_backend,
        MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}(),
    ) == 3
end

@testset "MomentLinearData dualizes directly to SOS model" begin
    _, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    M = typeof(one_σ)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))
    C = ComplexF64

    builder = NCTSSoS.MomentLinearBuilder(K, C, M)
    x_key = NCTSSoS.register_moment!(builder, σx[1])
    NCTSSoS.add_objective_terms!(builder, [x_key => -1.0 + 0im])
    NCTSSoS.add_psd_block!(
        builder,
        :HPSD,
        reshape(
            [
                [builder.identity => 1.0 + 0im],
                [x_key => 1.0 + 0im],
                [x_key => 1.0 + 0im],
                [builder.identity => 1.0 + 0im],
            ],
            2,
            2,
        ),
        NCTSSoS.BlockMeta(:HPSD, NCTSSoS.GlobalOrigin(1), [one_σ, σx[1]]);
        constraint_idx=1,
    )
    NCTSSoS.add_zero_constraint!(
        builder,
        [x_key => 1.0 + 0im],
        NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar),
    )

    linear = NCTSSoS.finalize!(builder)
    sos = NCTSSoS.sos_dualize(linear)

    @test length(sos.psd_dual_blocks) == 1
    @test length(sos.zero_duals) == 1
    @test JuMP.num_variables(sos.model) == 11
end

@testset "Pauli translation block emits linear entries without polynomial matrix" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    one_σ = one(first(ops[1]))
    M = typeof(one_σ)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))
    PReal = Polynomial{PauliAlgebra,UInt8,Float64}
    PComplex = Polynomial{PauliAlgebra,UInt8,ComplexF64}

    orbit_reps = NCTSSoS._pauli_contiguous_chain_orbit_representatives(ops, 1; periodic=true)
    block_basis = [mono for mono in orbit_reps if !isone(mono)]
    translated = Dict{M,Vector{M}}(
        rep => [NCTSSoS._translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
        for rep in block_basis
    )
    translated[one_σ] = fill(one_σ, n)
    rep_cache = Dict{M,M}()
    product_cache = NCTSSoS._TranslationProductCache(
        Dict{Tuple{M,M},Vector{Tuple{Int,ComplexF64,M}}}(),
        0,
        0,
    )

    symbolic = NCTSSoS._realify_hermitian_block(
        NCTSSoS._translation_momentum_block(
            block_basis,
            1,
            n,
            translated,
            rep_cache,
            PComplex,
            product_cache,
        ),
        PReal;
        atol=1e-12,
    )
    direct = NCTSSoS._translation_momentum_block_linear_entries(
        block_basis,
        1,
        n,
        translated,
        rep_cache,
        K,
        Float64,
        product_cache;
        real_moment_matrix=true,
        atol=1e-12,
    )

    @test size(direct) == size(symbolic)
    @test [
        collect(NCTSSoS._linearize_moment_polynomial(K, Float64, symbolic[i, j]))
        for i in axes(symbolic, 1), j in axes(symbolic, 2)
    ] == [collect(direct[i, j]) for i in axes(direct, 1), j in axes(direct, 2)]
end

@testset "Pauli translation base relaxation can build MomentLinearData directly" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
    )

    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.n_unique_moment_matrix_elements ==
        length(direct_linear.moments) ==
        length(symbolic_mp.linear.moments)
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test collect(direct_linear.objective_lin) == collect(symbolic_mp.linear.objective_lin)
    @test isempty(direct_linear.zero_constraints)
    @test length(direct_linear.psd_blocks_lin) == length(symbolic_mp.linear.psd_blocks_lin)
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
end

translation_zero_signature(zero) = (
    collect(zero.form),
    zero.kind,
    repr(zero.origin.label),
    zero.origin.row,
    zero.origin.col,
    zero.origin.part,
)

psd_entry_signature(linear) = [
    collect(linear.psd_blocks_lin[block_idx].entries[i, j])
    for block_idx in eachindex(linear.psd_blocks_lin)
    for i in 1:linear.psd_blocks_lin[block_idx].size
    for j in 1:linear.psd_blocks_lin[block_idx].size
]

@testset "MomentLinearBuilder supports finalized fermionic HPSD data" begin
    reg, (c, c_dag) = create_fermionic_variables(1:1)
    one_c = one(c[1])
    M = typeof(one_c)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_c)))
    P = typeof((1.0 + 0im) * one_c)

    symbolic_block = Matrix{P}(undef, 2, 2)
    symbolic_block[1, 1] = (1.0 + 0im) * one_c
    symbolic_block[1, 2] = (1.0 + 0im) * c[1]
    symbolic_block[2, 1] = (1.0 + 0im) * c_dag[1]
    symbolic_block[2, 2] = (1.0 + 0im) * one_c
    symbolic_mp = NCTSSoS.MomentProblem(
        zero(P),
        [(:HPSD, symbolic_block)],
        [one_c, c[1], c_dag[1]],
        3,
    )

    builder = NCTSSoS.MomentLinearBuilder(K, ComplexF64, M)
    c_key = NCTSSoS.register_moment!(builder, c[1])
    c_dag_key = symmetric_canon(NCTSSoS.expval(c_dag[1]))
    direct_entries = Matrix{Vector{Pair{K,ComplexF64}}}(undef, 2, 2)
    direct_entries[1, 1] = [builder.identity => 1.0 + 0im]
    direct_entries[1, 2] = [c_key => 1.0 + 0im]
    direct_entries[2, 1] = [c_dag_key => 1.0 + 0im]
    direct_entries[2, 2] = [builder.identity => 1.0 + 0im]
    NCTSSoS.add_psd_block!(
        builder,
        :HPSD,
        direct_entries,
        NCTSSoS.BlockMeta(:HPSD, NCTSSoS.GlobalOrigin(1), [one_c, c[1]]);
        constraint_idx=1,
    )
    direct_linear = NCTSSoS.finalize!(builder)

    @test NCTSSoS.assert_moment_linear_data_invariants(direct_linear) === nothing
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test psd_entry_signature(direct_linear) == psd_entry_signature(symbolic_mp.linear)

    model, extract = build_jump_model(
        direct_linear;
        formulation=:psd_blocks,
        representation=:complex,
    )
    @test extract isa Function
    @test JuMP.num_variables(model) > 0

    sos = NCTSSoS.sos_dualize(direct_linear)
    @test length(sos.psd_dual_blocks) == 1
    @test isempty(sos.zero_duals)
    @test NCTSSoS.sos_zero_duals(direct_linear, sos) == []
    @test NCTSSoS.sos_zero_dual_values(direct_linear, sos) == []
    @test NCTSSoS.sos_zero_dual_diagnostics(direct_linear, sos) == []

    set_optimizer(sos.model, Clarabel.Optimizer)
    set_silent(sos.model)
    optimize!(sos.model)
    @test termination_status(sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)

    residual = NCTSSoS.sos_dual_certificate_residual(direct_linear, sos)
    @test residual.max_abs_residual <= 1e-7
    @test abs(residual.identity_residual) <= 1e-7
    @test residual.moment_count == length(direct_linear.moments)
    @test residual.identity_moment == direct_linear.identity
    @test haskey(residual.residual_by_moment, residual.max_residual_moment)
    @test residual.residual_by_moment[residual.max_residual_moment] ==
        residual.max_residual_value
    @test abs(residual.max_residual_value) == residual.max_abs_residual

    block_values = NCTSSoS.sos_dual_block_values(direct_linear, sos)
    @test length(block_values) == 1
    @test only(block_values).representation == :hermitian_lift
    @test only(block_values).transform_family === nothing
    @test size(only(block_values).native_value) == (2, 2)

    diagnostics = NCTSSoS.sos_dual_block_diagnostics(direct_linear, sos)
    @test length(diagnostics) == 1
    @test only(diagnostics).representation == :hermitian_lift
    @test only(diagnostics).transform_family === nothing
    @test only(diagnostics).native_hermitian_residual <= 1e-7
    @test only(diagnostics).native_min_eigenvalue >= -1e-7

    native_sos = NCTSSoS.sos_dualize(direct_linear; hermitian_representation=:native)
    @test length(native_sos.psd_dual_blocks) == 1
    @test isempty(native_sos.zero_duals)
    @test only(NCTSSoS.sos_dual_blocks(native_sos)).representation == :hermitian
    @test JuMP.num_variables(native_sos.model) < JuMP.num_variables(sos.model)

    set_optimizer(native_sos.model, Clarabel.Optimizer)
    set_silent(native_sos.model)
    optimize!(native_sos.model)
    @test termination_status(native_sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
    @test objective_value(native_sos.model) ≈ objective_value(sos.model) atol = 1e-7

    native_residual = NCTSSoS.sos_dual_certificate_residual(direct_linear, native_sos)
    @test native_residual.max_abs_residual <= 1e-7

    native_values = NCTSSoS.sos_dual_block_values(direct_linear, native_sos)
    @test only(native_values).representation == :hermitian
    @test size(only(native_values).native_value) == (2, 2)
end

@testset "MomentLinearData SOS helpers accept direct real data" begin
    _, (σx, _, _) = create_pauli_variables(1:1)
    one_σ = one(σx[1])
    M = typeof(one_σ)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))

    builder = NCTSSoS.MomentLinearBuilder(K, Float64, M)
    x_key = NCTSSoS.register_moment!(builder, σx[1])
    NCTSSoS.add_psd_block!(
        builder,
        :PSD,
        reshape([[builder.identity => 1.0]], 1, 1),
        NCTSSoS.BlockMeta(:PSD, NCTSSoS.GlobalOrigin(1), [one_σ]);
        constraint_idx=1,
    )
    NCTSSoS.add_zero_constraint!(
        builder,
        [x_key => 1.0],
        NCTSSoS.ZeroMatrixOrigin(2, 1, 1, :scalar),
    )
    direct_linear = NCTSSoS.finalize!(builder)
    sos = NCTSSoS.sos_dualize(direct_linear)

    set_optimizer(sos.model, Clarabel.Optimizer)
    set_silent(sos.model)
    optimize!(sos.model)
    @test termination_status(sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)

    residual = NCTSSoS.sos_dual_certificate_residual(direct_linear, sos)
    @test residual.max_abs_residual <= 1e-7
    @test abs(residual.identity_residual) <= 1e-7
    @test residual.moment_count == length(direct_linear.moments)
    @test residual.identity_moment == direct_linear.identity
    @test haskey(residual.residual_by_moment, residual.max_residual_moment)
    @test residual.residual_by_moment[residual.max_residual_moment] ==
        residual.max_residual_value
    @test abs(residual.max_residual_value) == residual.max_abs_residual

    block_values = NCTSSoS.sos_dual_block_values(direct_linear, sos)
    @test length(block_values) == 1
    @test only(block_values).representation == :psd
    @test size(only(block_values).native_value) == (1, 1)

    @test isdefined(NCTSSoS, :sos_zero_dual_diagnostics)
    zero_duals = NCTSSoS.sos_zero_duals(direct_linear, sos)
    zero_values = NCTSSoS.sos_zero_dual_values(direct_linear, sos)
    zero_diagnostics = NCTSSoS.sos_zero_dual_diagnostics(direct_linear, sos)
    @test length(zero_duals) == 1
    @test length(zero_values) == 1
    @test length(zero_diagnostics) == 1
    @test only(zero_duals).label === nothing
    @test only(zero_duals).feature === nothing
    @test only(zero_duals).decomposition === nothing
    @test only(zero_duals).reason === nothing
    @test only(zero_duals).coefficient_domain === nothing
    @test only(zero_duals).exact_coefficient_domain === nothing
    @test only(zero_duals).term_count == 1
    @test only(zero_values).label === nothing
    @test only(zero_values).feature === nothing
    @test only(zero_values).decomposition === nothing
    @test only(zero_values).reason === nothing
    @test only(zero_values).coefficient_domain === nothing
    @test only(zero_values).exact_coefficient_domain === nothing
    @test only(zero_values).term_count == 1
    @test only(zero_diagnostics).origin == only(zero_values).origin
    @test only(zero_diagnostics).label === nothing
    @test only(zero_diagnostics).feature === nothing
    @test only(zero_diagnostics).decomposition === nothing
    @test only(zero_diagnostics).reason === nothing
    @test only(zero_diagnostics).coefficient_domain === nothing
    @test only(zero_diagnostics).exact_coefficient_domain === nothing
    @test only(zero_diagnostics).kind == :zero
    @test only(zero_diagnostics).term_count == 1
    @test only(zero_diagnostics).max_abs_coefficient == 1.0
    @test only(zero_diagnostics).value == only(zero_values).value
    @test isfinite(only(zero_diagnostics).value)

    diagnostics = NCTSSoS.sos_dual_block_diagnostics(direct_linear, sos)
    @test length(diagnostics) == 1
    @test only(diagnostics).representation == :psd
    @test only(diagnostics).native_min_eigenvalue >= -1e-7

    @test isdefined(NCTSSoS, :sos_dual_certificate_diagnostics)
    certificate = NCTSSoS.sos_dual_certificate_diagnostics(direct_linear, sos)
    @test certificate.residual == residual
    @test certificate.psd_blocks == diagnostics
    @test certificate.zero_duals == zero_diagnostics
    @test certificate.moment_count == residual.moment_count
    @test certificate.psd_block_count == length(diagnostics)
    @test certificate.zero_dual_count == length(zero_diagnostics)
    @test certificate.max_abs_residual == residual.max_abs_residual
    @test certificate.identity_residual == residual.identity_residual
    @test certificate.max_residual_moment == residual.max_residual_moment
    @test certificate.max_residual_value == residual.max_residual_value
end

@testset "Pauli translation direct builder supports scalar constraints" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    _, _, σz = ops
    objective = heisenberg_chain_hamiltonian(ops)
    pop = polyopt(
        objective,
        registry;
        eq_constraints=[sum(σz)],
        ineq_constraints=[one(objective) + objective],
    )

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count == 1
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test last(direct_report.block_labels) == (feature=:scalar_inequality, index=1)
    @test collect(only(last(direct_linear.psd_blocks_lin).entries)) ==
        collect(only(last(symbolic_mp.linear.psd_blocks_lin).entries))
    @test [
        (collect(zero.form), zero.kind, zero.origin)
        for zero in direct_linear.zero_constraints
    ] == [
        (collect(zero.form), zero.kind, zero.origin)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

@testset "Pauli translation direct builder supports complex scalar constraints" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    _, _, σz = ops
    objective = heisenberg_chain_hamiltonian(ops)
    pop = polyopt(
        objective,
        registry;
        eq_constraints=[sum(σz)],
        ineq_constraints=[one(objective) + objective],
    )

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        real_moment_matrix=false,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        real_moment_matrix=false,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test collect(only(last(direct_linear.psd_blocks_lin).entries)) ==
        collect(only(last(symbolic_mp.linear.psd_blocks_lin).entries))
    @test [
        translation_zero_signature(zero)
        for zero in direct_linear.zero_constraints
    ] == [
        translation_zero_signature(zero)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

moment_eq_zero_signature(zero) = (
    collect(zero.form),
    zero.kind,
    zero.origin.label.feature,
    zero.origin.label.index,
    zero.origin.label.row,
    repr(zero.origin.label.row_monomial),
    zero.origin.row,
    zero.origin.col,
    zero.origin.part,
)

@testset "Pauli translation direct builder supports moment equality rows" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    _, _, σz = ops
    objective = heisenberg_chain_hamiltonian(ops)
    meq = sum(one(objective) - σz[i] * σz[mod1(i + 1, n)] for i in 1:n)
    pop = polyopt(objective, registry; moment_eq_constraints=[meq])

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        moment_eq_zero_signature(zero)
        for zero in direct_linear.zero_constraints
    ] == [
        moment_eq_zero_signature(zero)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

@testset "Pauli translation direct builder supports complex moment equality rows" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    _, _, σz = ops
    objective = heisenberg_chain_hamiltonian(ops)
    meq = sum(one(objective) - σz[i] * σz[mod1(i + 1, n)] for i in 1:n)
    pop = polyopt(objective, registry; moment_eq_constraints=[meq])

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        real_moment_matrix=false,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        real_moment_matrix=false,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        moment_eq_zero_signature(zero)
        for zero in direct_linear.zero_constraints
    ] == [
        moment_eq_zero_signature(zero)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

@testset "Pauli fixed-momentum reflection transform emits linear entries" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    one_σ = one(first(ops[1]))
    M = typeof(one_σ)
    K = typeof(symmetric_canon(NCTSSoS.expval(one_σ)))
    PReal = Polynomial{PauliAlgebra,UInt8,Float64}
    PComplex = Polynomial{PauliAlgebra,UInt8,ComplexF64}

    orbit_reps = NCTSSoS._pauli_contiguous_chain_orbit_representatives(ops, 1; periodic=true)
    block_basis = [mono for mono in orbit_reps if !isone(mono)]
    translated = Dict{M,Vector{M}}(
        rep => [NCTSSoS._translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
        for rep in block_basis
    )
    translated[one_σ] = fill(one_σ, n)
    rep_cache = Dict{M,M}()
    product_cache = NCTSSoS._TranslationProductCache(
        Dict{Tuple{M,M},Vector{Tuple{Int,ComplexF64,M}}}(),
        0,
        0,
    )
    complex_symbolic = NCTSSoS._translation_momentum_block(
        block_basis,
        0,
        n,
        translated,
        rep_cache,
        PComplex,
        product_cache,
    )
    complex_linear = NCTSSoS._translation_momentum_block_linear_entries(
        block_basis,
        0,
        n,
        translated,
        rep_cache,
        K,
        ComplexF64,
        product_cache;
        real_moment_matrix=false,
    )

    for adapted in NCTSSoS._translation_reflection_adapted_blocks(block_basis, 0, n; atol=1e-12)
        symbolic = NCTSSoS._realify_hermitian_block(
            NCTSSoS._hermitianize_pauli_polynomial_block(
                NCTSSoS._transform_polynomial_block(
                    complex_symbolic,
                    adapted.row_basis,
                    adapted.row_basis;
                    atol=1e-12,
                ),
                PComplex,
            ),
            PReal;
            atol=1e-12,
        )
        direct = NCTSSoS._realify_hermitian_linear_block(
            NCTSSoS._hermitianize_pauli_linear_block(
                NCTSSoS._transform_linear_block(
                    complex_linear,
                    adapted.row_basis,
                    adapted.row_basis;
                    atol=1e-12,
                ),
            ),
            Float64;
            atol=1e-12,
        )

        @test size(direct) == size(symbolic)
        @test [
            collect(NCTSSoS._linearize_moment_polynomial(K, Float64, symbolic[i, j]))
            for i in axes(symbolic, 1), j in axes(symbolic, 2)
        ] == [collect(direct[i, j]) for i in axes(direct, 1), j in axes(direct, 2)]
    end
end

@testset "Pauli translation direct builder supports fixed reflection sectors" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        reflection_symmetry=true,
        momenta=[0, 2],
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        reflection_symmetry=true,
        momenta=[0, 2],
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test length(direct_linear.psd_blocks_lin) == length(symbolic_mp.linear.psd_blocks_lin)
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
end

@testset "Pauli translation direct builder supports real reflection sectors" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        reflection_symmetry=true,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        reflection_symmetry=true,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
end

@testset "Pauli translation direct builder supports full contiguous RDM blocks" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        contiguous_rdm_k=2,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        contiguous_rdm_k=2,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
end

@testset "Pauli translation direct builder supports complex full contiguous RDM blocks" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        real_moment_matrix=false,
        contiguous_rdm_k=2,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        real_moment_matrix=false,
        contiguous_rdm_k=2,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.real_moment_matrix == symbolic_report.real_moment_matrix == false
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test psd_entry_signature(direct_linear) == psd_entry_signature(symbolic_mp.linear)
end

@testset "Pauli translation direct builder supports RDM support extension" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        contiguous_rdm_k=2,
        contiguous_rdm_support=:extend,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        contiguous_rdm_k=2,
        contiguous_rdm_support=:extend,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.n_unique_moment_matrix_elements ==
        symbolic_report.n_unique_moment_matrix_elements
    @test direct_report.linear_moment_count == symbolic_report.linear_moment_count
    @test direct_report.n_unique_moment_matrix_elements < direct_report.linear_moment_count
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
end

@testset "Pauli translation direct builder supports U1-decomposed RDM blocks" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        u1_symmetry=true,
        contiguous_rdm_k=3,
        contiguous_rdm_decomposition=:u1,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        u1_symmetry=true,
        contiguous_rdm_k=3,
        contiguous_rdm_decomposition=:u1,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
    @test [
        (collect(zero.form), zero.kind, zero.origin)
        for zero in direct_linear.zero_constraints
    ] == [
        (collect(zero.form), zero.kind, zero.origin)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

@testset "Pauli translation direct builder supports complex U1-decomposed RDM blocks" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        real_moment_matrix=false,
        u1_symmetry=true,
        contiguous_rdm_k=3,
        contiguous_rdm_decomposition=:u1,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        real_moment_matrix=false,
        u1_symmetry=true,
        contiguous_rdm_k=3,
        contiguous_rdm_decomposition=:u1,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_report.real_moment_matrix == symbolic_report.real_moment_matrix == false
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test psd_entry_signature(direct_linear) == psd_entry_signature(symbolic_mp.linear)
    @test [
        translation_zero_signature(zero)
        for zero in direct_linear.zero_constraints
    ] == [
        translation_zero_signature(zero)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

@testset "Pauli translation direct builder supports SU2-decomposed RDM blocks" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        su2_symmetry=true,
        contiguous_rdm_k=2,
        contiguous_rdm_decomposition=:su2,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        su2_symmetry=true,
        contiguous_rdm_k=2,
        contiguous_rdm_decomposition=:su2,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.block_coefficient_domains == symbolic_report.block_coefficient_domains
    @test direct_report.block_exact_coefficient_domains == symbolic_report.block_exact_coefficient_domains
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
    @test [
        (collect(zero.form), zero.kind, zero.origin)
        for zero in direct_linear.zero_constraints
    ] == [
        (collect(zero.form), zero.kind, zero.origin)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

@testset "Pauli SU2 RDM dense accumulator matches Dict fallback" begin
    n = 4
    k = 3
    T = UInt16
    M = NormalMonomial{PauliAlgebra,T}
    K = Vector{T}
    C = ComplexF64
    orbit_monos, orbit_indices = NCTSSoS._contiguous_rdm_reduced_orbit_data(n, k, M)
    code_actions = NCTSSoS._pauli_rdm_code_actions(C, k)
    schur_transform, states = NCTSSoS._pauli_su2_schur_matrix(k)
    schur_rows = NCTSSoS._pauli_sparse_transform_rows(transpose(schur_transform); atol=1e-12)
    columns = NCTSSoS._pauli_su2_state_columns(states)
    block = first(pauli_su2_rdm_blocks(k))
    reference_m2 = -block.spin2
    reference_rows = [
        columns[(block.spin2, reference_m2, mult)]
        for mult in 1:block.multiplicity
    ]
    transform = NCTSSoS._select_sparse_transform_rows(schur_rows, reference_rows)
    rows_by_state = NCTSSoS._sparse_transform_columns_by_state(transform)
    state_sign_parities = NCTSSoS._pauli_rdm_state_sign_parities(k)
    for sign_mask in 0:((1 << k) - 1), state0 in 0:((1 << k) - 1)
        @test state_sign_parities[state0 + 1, sign_mask + 1] ==
              isodd(count_ones(sign_mask & state0))
    end

    dense_entries = NCTSSoS._translation_contiguous_rdm_su2_reduced_linear_block_dense(
        orbit_monos,
        orbit_indices,
        k,
        transform,
        rows_by_state,
        code_actions,
        K,
        C;
        atol=1e-12,
        state_sign_parities,
    )
    dense_stage_times = Dict{Symbol,Int}()
    dense_timed_entries = NCTSSoS._translation_contiguous_rdm_su2_reduced_linear_block_dense(
        orbit_monos,
        orbit_indices,
        k,
        transform,
        rows_by_state,
        code_actions,
        K,
        C;
        atol=1e-12,
        stage_times_ns=dense_stage_times,
        state_sign_parities,
    )
    dict_entries = NCTSSoS._translation_contiguous_rdm_su2_reduced_linear_block_dict(
        orbit_monos,
        orbit_indices,
        k,
        transform,
        rows_by_state,
        code_actions,
        K,
        C;
        atol=1e-12,
        state_sign_parities,
    )
    dict_stage_times = Dict{Symbol,Int}()
    dict_timed_entries = NCTSSoS._translation_contiguous_rdm_su2_reduced_linear_block_dict(
        orbit_monos,
        orbit_indices,
        k,
        transform,
        rows_by_state,
        code_actions,
        K,
        C;
        atol=1e-12,
        stage_times_ns=dict_stage_times,
        state_sign_parities,
    )

    @test size(dense_entries) == size(dict_entries)
    @test [
        collect(dense_entries[row, col])
        for col in axes(dense_entries, 2)
        for row in axes(dense_entries, 1)
    ] == [
        collect(dict_entries[row, col])
        for col in axes(dict_entries, 2)
        for row in axes(dict_entries, 1)
    ]
    @test [
        collect(dense_timed_entries[row, col])
        for col in axes(dense_timed_entries, 2)
        for row in axes(dense_timed_entries, 1)
    ] == [
        collect(dense_entries[row, col])
        for col in axes(dense_entries, 2)
        for row in axes(dense_entries, 1)
    ]
    @test [
        collect(dict_timed_entries[row, col])
        for col in axes(dict_timed_entries, 2)
        for row in axes(dict_timed_entries, 1)
    ] == [
        collect(dict_entries[row, col])
        for col in axes(dict_entries, 2)
        for row in axes(dict_entries, 1)
    ]
    @test dense_stage_times[:su2_extend_rdm_block_accumulate] > 0
    @test dense_stage_times[:su2_extend_rdm_block_finalize] > 0
    @test dict_stage_times[:su2_extend_rdm_block_accumulate] > 0
    @test dict_stage_times[:su2_extend_rdm_block_finalize] > 0
end

@testset "Pauli translation direct builder supports complex SU2-decomposed RDM blocks" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        real_moment_matrix=false,
        su2_symmetry=true,
        contiguous_rdm_k=2,
        contiguous_rdm_decomposition=:su2,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        1;
        sign_symmetry=false,
        real_moment_matrix=false,
        su2_symmetry=true,
        contiguous_rdm_k=2,
        contiguous_rdm_decomposition=:su2,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.block_coefficient_domains == symbolic_report.block_coefficient_domains
    @test direct_report.block_exact_coefficient_domains ==
        symbolic_report.block_exact_coefficient_domains
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_report.real_moment_matrix == symbolic_report.real_moment_matrix == false
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test psd_entry_signature(direct_linear) == psd_entry_signature(symbolic_mp.linear)
    @test [
        translation_zero_signature(zero)
        for zero in direct_linear.zero_constraints
    ] == [
        translation_zero_signature(zero)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

linear_state_opt_zero_signature(zero) = (
    collect(zero.form),
    zero.kind,
    zero.origin.label.feature,
    zero.origin.label.width,
    repr(zero.origin.label.test_monomial),
    zero.origin.row,
    zero.origin.col,
    zero.origin.part,
)

@testset "Pauli translation direct builder supports linear state-opt rows" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        linear_state_opt_width=2,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        linear_state_opt_width=2,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
    @test [
        linear_state_opt_zero_signature(zero)
        for zero in direct_linear.zero_constraints
    ] == [
        linear_state_opt_zero_signature(zero)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

@testset "Pauli translation direct builder supports complex linear state-opt rows" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        real_moment_matrix=false,
        linear_state_opt_width=2,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        real_moment_matrix=false,
        linear_state_opt_width=2,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.zero_constraint_count == symbolic_report.zero_constraint_count
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
    @test [
        linear_state_opt_zero_signature(zero)
        for zero in direct_linear.zero_constraints
    ] == [
        linear_state_opt_zero_signature(zero)
        for zero in symbolic_mp.linear.zero_constraints
    ]
end

@testset "Pauli translation direct builder supports PSD state-opt blocks" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        psd_state_opt_width=1,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        psd_state_opt_width=1,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
end

@testset "Pauli translation direct builder supports complex PSD state-opt blocks" begin
    n = 4
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

    symbolic_mp, symbolic_report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        real_moment_matrix=false,
        psd_state_opt_width=1,
    )
    direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        2;
        sign_symmetry=false,
        real_moment_matrix=false,
        psd_state_opt_width=1,
    )

    @test direct_report.block_labels == symbolic_report.block_labels
    @test direct_report.psd_block_sizes == symbolic_report.psd_block_sizes
    @test direct_report.logical_block_sizes == symbolic_report.logical_block_sizes
    @test direct_report.psd_block_sizes[end - 2:end] == fill(3, 3)
    @test direct_linear.moments == symbolic_mp.linear.moments
    @test all(block.meta.cone == :HPSD for block in direct_linear.psd_blocks_lin[end - 2:end])
    @test [
        collect(direct_linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(direct_linear.psd_blocks_lin)
        for i in 1:direct_linear.psd_blocks_lin[block_idx].size
        for j in 1:direct_linear.psd_blocks_lin[block_idx].size
    ] == [
        collect(symbolic_mp.linear.psd_blocks_lin[block_idx].entries[i, j])
        for block_idx in eachindex(symbolic_mp.linear.psd_blocks_lin)
        for i in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
        for j in 1:symbolic_mp.linear.psd_blocks_lin[block_idx].size
    ]
end

@testset "moment_relax attaches cache after symbolic mutations" begin
    reg, (c, c_dag) = create_fermionic_variables(1:1)
    n = c_dag[1] * c[1]
    objective = 1.0 * n
    pop = polyopt(objective, reg; moment_eq_constraints=[1.0 * (n - one(n))])
    sparsity = compute_sparsity(pop, SolverConfig(optimizer=nothing, order=1))
    mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)

    @test any(c -> c[1] == :Zero, mp.constraints)
    @test !isempty(mp.linear.zero_constraints)
    @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing

    P = typeof(mp.objective)
    odd_poly = convert(P, 1.0 * c[1])
    push!(mp.constraints, (:HPSD, reshape([odd_poly], 1, 1)))
    before_zero_forms = length(mp.linear.zero_constraints)
    NCTSSoS._add_parity_constraints!(mp)
    @test length(mp.linear.zero_constraints) > before_zero_forms
    @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
end

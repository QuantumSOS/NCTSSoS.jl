using Test, NCTSSoS

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

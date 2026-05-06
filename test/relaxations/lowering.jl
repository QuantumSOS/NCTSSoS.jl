using Test, NCTSSoS, JuMP

@testset "Moment lowering cached pivots" begin
    reg, (σx, σy, σz) = create_pauli_variables(1:2)
    m_pos = σx[1]
    m_neg = σy[1]
    m_ipos = σz[1]
    m_ineg = σx[2]

    P = typeof(1.0 * m_pos)
    z = zero(P)
    mat = fill(z, 4, 4)
    mat[1, 1] = 1.0 * m_pos
    mat[1, 2] = -1.0 * m_neg
    mat[1, 3] = 1.0im * m_ipos
    mat[1, 4] = -1.0im * m_ineg
    mat[2, 1] = 1.0 * m_pos  # duplicate pivot candidate loses to (1, 1)

    mp = NCTSSoS.MomentProblem(
        zero(P),
        [(:HPSD, mat)],
        [one(m_pos), m_pos, m_neg, m_ipos, m_ineg],
        4,
    )

    pivots = mp.linear.pivots

    key(m) = symmetric_canon(NCTSSoS.expval(m))
    @test pivots[key(m_pos)].phase == 1 + 0im
    @test pivots[key(m_neg)].phase == -1 + 0im
    @test pivots[key(m_ipos)].phase == 0 + 1im
    @test pivots[key(m_ineg)].phase == 0 - 1im
    @test (pivots[key(m_pos)].block, pivots[key(m_pos)].row, pivots[key(m_pos)].col) == (1, 1, 1)
    @test mp.linear.psd_block_constraint_idx[pivots[key(m_pos)].block] == 1
end

@testset "Moment lowering PSD-block complex formulation solves tiny Hermitian model" begin
    reg, (b, b_dag) = create_bosonic_variables(1:1)
    objective = -(1.0 * b[1] + 1.0 * b_dag[1])
    P = typeof(objective)

    block = Matrix{P}(undef, 2, 2)
    block[1, 1] = 1.0 * one(b[1])
    block[1, 2] = 1.0 * b[1]
    block[2, 1] = 1.0 * b_dag[1]
    block[2, 2] = 1.0 * one(b[1])

    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block)],
        [one(b[1]), b[1], b_dag[1]],
        3,
    )

    model, extract = build_jump_model(mp; formulation=:psd_blocks, representation=:complex)
    set_optimizer(model, SOLVER)
    set_silent(model)
    optimize!(model)
    NCTSSoS._check_solver_status(model)

    @test objective_value(model) ≈ -2.0 atol = 1e-6
    monomap = extract()
    @test monomap[symmetric_canon(NCTSSoS.expval(one(b[1])))] ≈ 1.0 atol = 1e-7
end

@testset "Moment lowering PSD-block resolver uses Hermitian adjoint pivots" begin
    reg, (b, b_dag) = create_bosonic_variables(1:1)
    objective = 1.0im * b_dag[1] - 1.0im * b[1]
    P = typeof(objective)

    block = Matrix{P}(undef, 2, 2)
    block[1, 1] = 1.0 * one(b[1])
    block[1, 2] = 1.0im * b[1]
    block[2, 1] = -1.0im * b_dag[1]
    block[2, 2] = 1.0 * one(b[1])

    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block)],
        [one(b[1]), b[1], b_dag[1]],
        3,
    )

    key(m) = symmetric_canon(NCTSSoS.expval(m))
    @test mp.linear.pivots[key(b_dag[1])].adjoint

    model, extract = build_jump_model(mp; formulation=:psd_blocks, representation=:complex)
    set_optimizer(model, SOLVER)
    set_silent(model)
    optimize!(model)
    NCTSSoS._check_solver_status(model)

    @test objective_value(model) ≈ -2.0 atol = 1e-6
    monomap = extract()
    @test monomap[key(b[1])] ≈ -1.0im atol = 1e-6
    @test monomap[key(b_dag[1])] ≈ 1.0im atol = 1e-6
end

@testset "Moment lowering PSD-block binds full HPSD and raw complex Zero rows" begin
    MOI = JuMP.MOI
    reg, (b, b_dag) = create_bosonic_variables(1:1)
    P = typeof(0.0 * one(b[1]))

    block = Matrix{P}(undef, 2, 2)
    block[1, 1] = 1.0 * one(b[1])
    block[1, 2] = 1.0 * b[1]
    block[2, 1] = zero(P)          # Deliberately non-Hermitian: lower entry must not be ignored.
    block[2, 2] = 1.0 * one(b[1])
    zero_mat = reshape([1.0 * b[1]], 1, 1)

    mp = NCTSSoS.MomentProblem(
        zero(P),
        [(:HPSD, block), (:Zero, zero_mat)],
        [one(b[1]), b[1], b_dag[1]],
        3,
    )

    model, _ = build_jump_model(mp; formulation=:psd_blocks, representation=:complex)
    backend = JuMP.backend(model)

    @test MOI.get(backend, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}()) == 2
    @test MOI.get(backend, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{ComplexF64},MOI.EqualTo{ComplexF64}}()) == 3
end

@testset "Moment lowering moment-variable formulation uses cached identity" begin
    reg, (b, b_dag) = create_bosonic_variables(1:1)
    objective = 0.0 * one(b[1])
    P = typeof(objective)

    block = Matrix{P}(undef, 1, 1)
    block[1, 1] = 1.0 * one(b[1])

    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block)],
        [b[1], b_dag[1]],
        1,
    )

    model, _ = build_jump_model(mp; formulation=:moment_variables, representation=:real)
    @test JuMP.num_variables(model) == 2 * length(mp.linear.moments)
end

@testset "Moment lowering cached free keys report orphans" begin
    reg, (σx, _, _) = create_pauli_variables(1:1)
    one_m = one(σx[1])
    orphan_m = σx[1]
    P = typeof(1.0 * orphan_m)

    mat = reshape([1.0 * one_m], 1, 1)
    mp = NCTSSoS.MomentProblem(
        1.0 * orphan_m,
        [(:HPSD, mat)],
        [one_m, orphan_m],
        1,
    )

    orphan_key = symmetric_canon(NCTSSoS.expval(orphan_m))
    @test orphan_key in NCTSSoS.orphan_keys(mp)
    @test orphan_key in mp.linear.free_keys
    @test !haskey(mp.linear.pivots, orphan_key)
end

@testset "Moment lowering PSD-block formulation represents orphans as free variables" begin
    reg, (b, b_dag) = create_bosonic_variables(1:1)
    linear = 1.0 * b[1] + 1.0 * b_dag[1]
    objective = -linear
    P = typeof(objective)

    block = Matrix{P}(undef, 2, 2)
    block[1, 1] = 1.0 * one(b[1])
    block[1, 2] = linear
    block[2, 1] = linear
    block[2, 2] = 1.0 * one(b[1])

    mp = NCTSSoS.MomentProblem(
        objective,
        [(:HPSD, block)],
        [one(b[1]), b[1], b_dag[1]],
        3,
    )

    @test length(NCTSSoS.orphan_keys(mp)) == 2
    @test_throws ArgumentError build_jump_model(mp; formulation=:psd_blocks, representation=:complex)

    aux_model, _ = build_jump_model(mp;
        formulation=:psd_blocks,
        representation=:complex,
        orphan_policy=:aux_psd_free,
    )
    @test JuMP.num_variables(aux_model) > 0

    model, extract = build_jump_model(mp;
        formulation=:psd_blocks,
        representation=:complex,
        orphan_policy=:free_variables,
    )
    set_optimizer(model, SOLVER)
    set_silent(model)
    optimize!(model)
    NCTSSoS._check_solver_status(model)

    @test objective_value(model) ≈ -1.0 atol = 1e-6
    monomap = extract()
    moment_sum = monomap[symmetric_canon(NCTSSoS.expval(b[1]))] +
                 monomap[symmetric_canon(NCTSSoS.expval(b_dag[1]))]
    @test real(moment_sum) ≈ 1.0 atol = 1e-6
end

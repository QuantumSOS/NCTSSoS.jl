using Test, NCTSSoS, JuMP

@testset "Moment lowering pivot discovery" begin
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

    pivots = NCTSSoS.discover_pivots(mp)

    key(m) = symmetric_canon(NCTSSoS.expval(m))
    @test pivots[key(m_pos)].phase == 1 + 0im
    @test pivots[key(m_neg)].phase == -1 + 0im
    @test pivots[key(m_ipos)].phase == 0 + 1im
    @test pivots[key(m_ineg)].phase == 0 - 1im
    @test (pivots[key(m_pos)].constraint_idx, pivots[key(m_pos)].i, pivots[key(m_pos)].j) == (1, 1, 1)
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

@testset "Moment lowering pivot discovery reports orphans" begin
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

    @test_throws ArgumentError NCTSSoS.discover_pivots(mp)
    @test symmetric_canon(NCTSSoS.expval(orphan_m)) in NCTSSoS.orphan_keys(mp)
end

using Test, NCTSSoS, LinearAlgebra

using NCTSSoS: get_ncbasis, monomials

if !@isdefined(MOI)
    const MOI = NCTSSoS.MOI
end

function _pipeline_moment_key(mono::NormalMonomial)
    return symmetric_canon(NCTSSoS.expval(mono))
end

function _pipeline_moment_lookup(monomap, mono)
    return get(monomap, _pipeline_moment_key(mono), 0.0)
end

function _pipeline_evaluate_monomial(matrices::Dict, mono::NormalMonomial)
    sample_matrix = first(values(matrices))
    T = eltype(sample_matrix)
    dim = size(sample_matrix, 1)
    result = Matrix{T}(I, dim, dim)
    for idx in mono.word
        result *= matrices[idx]
    end
    return result
end

function _pipeline_expectation(matrices::Dict, xi::AbstractVector, mono::NormalMonomial)
    return dot(xi, _pipeline_evaluate_monomial(matrices, mono) * xi)
end

function _solve_dense_moment_problem(pop, cfg)
    sparsity = compute_sparsity(pop, cfg)
    mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    result = NCTSSoS.solve_moment_problem(mp, SOLVER)
    NCTSSoS._check_solver_status(result.model)
    return result
end

function _dense_hankel_data(monomap, reg, H_deg, hankel_deg; atol=1e-6)
    full_basis = get_ncbasis(reg, H_deg)
    basis = get_ncbasis(reg, hankel_deg)
    hankel = NCTSSoS.hankel_matrix(monomap, full_basis)
    flatness = test_flatness(hankel, full_basis, basis; atol=atol)
    hankel_flat = flat_extend(hankel, full_basis, basis; atol=atol)
    flatness_flat = test_flatness(hankel_flat, full_basis, basis; atol=atol)

    return (
        full_basis=full_basis,
        basis=basis,
        hankel=hankel,
        flatness=flatness,
        hankel_flat=hankel_flat,
        flatness_flat=flatness_flat,
    )
end

function _solve_dense_example(pop, reg, order; H_deg=order, hankel_deg=order - 1, atol=1e-6)
    cfg = SolverConfig(
        optimizer=SOLVER,
        order=order,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
    )
    result = _solve_dense_moment_problem(pop, cfg)
    hankel_data = _dense_hankel_data(result.monomap, reg, H_deg, hankel_deg; atol=atol)
    return merge((result=result,), hankel_data)
end

function _solve_ball_example()
    reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
    X, Y = vars

    f = 2.0 - X^2 + X * Y^2 * X - Y^2
    g = 1.0 - X^2 - Y^2
    pop = polyopt(f, reg; ineq_constraints=[g])

    return merge(
        (
            reg=reg,
            X=X,
            Y=Y,
            f=f,
            g=g,
        ),
        _solve_dense_example(pop, reg, 3; H_deg=3, hankel_deg=2, atol=1e-6),
    )
end

@testset "Dense GNS Pipeline" begin
    @testset "Flat extension correctness" begin
        reg, (x,) = create_noncommutative_variables([("X", 1:1)])
        full_basis = get_ncbasis(reg, 2)
        basis = get_ncbasis(reg, 1)
        H = [1.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 1.0]
        H_rank_deficient = ones(3, 3)

        flatness = test_flatness(H, full_basis, basis; atol=1e-12)
        Hflat = flat_extend(H, full_basis, basis; atol=1e-12)
        flatness_rank_deficient = test_flatness(H_rank_deficient, full_basis, basis; atol=1e-12)
        Hflat_rank_deficient = flat_extend(H_rank_deficient, full_basis, basis; atol=1e-12)

        @test flatness.is_flat
        @test flatness.rank_principal == 2
        @test flatness.rank_full == 2
        @test flatness.err_flat ≤ 1e-12
        @test rank(Hflat) == 2
        @test Hflat ≈ H atol = 1e-12
        @test flatness_rank_deficient.is_flat
        @test flatness_rank_deficient.rank_principal == 1
        @test flatness_rank_deficient.rank_full == 1
        @test flatness_rank_deficient.err_flat ≤ 1e-12
        @test Hflat_rank_deficient ≈ H_rank_deficient atol = 1e-12

        basis_deg1 = get_ncbasis(reg, 1)
        Hdeg1 = Matrix{Float64}(I, length(basis_deg1), length(basis_deg1))
        flatness_deg1 = test_flatness(Hdeg1, basis_deg1, basis_deg1; atol=1e-12)

        @test flatness_deg1.is_flat
        @test flatness_deg1.rank_principal == 2
        @test flatness_deg1.rank_full == 2
        @test flatness_deg1.err_flat == 0.0
    end

    @testset "Cholesky GNS toy" begin
        reg, (x,) = create_noncommutative_variables([("X", 1:1)])
        H = [1.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 1.0]
        gns = gns_reconstruct(H, reg, 2; method=:cholesky, hankel_deg=1, atol=1e-12)
        Xmat = gns.matrices[reg[:X₁]]

        @test gns.rank == 2
        @test Xmat ≈ [0.0 1.0; 1.0 0.0] atol = 1e-12
        @test gns.xi ≈ [1.0, 0.0] atol = 1e-12

        for k in 0:4
            mono = k == 0 ? one(x[1]) : only(monomials(x[1]^k))
            expected = iseven(k) ? 1.0 : 0.0
            @test real(_pipeline_expectation(gns.matrices, gns.xi, mono)) ≈ expected atol = 1e-12
        end
    end

    @testset "Cholesky basis repair keeps identity" begin
        witness = [
            0.06193274031408013 -0.5958244153640522 1.0857940215432762 0.1759399913010747
            0.2784058141640002 0.04665938957338174 -1.5765649225859841 0.8653808054093252
        ]
        hankel_block = witness' * witness
        selected, svals = NCTSSoS._gns_cholesky_basis_indices(hankel_block; atol=1e-8)

        @test length(selected) == 2
        @test selected[1] == 1
        @test svals[2] > 1e-8
    end

    @testset "SVD vs Cholesky agreement" begin
        reg, (x,) = create_noncommutative_variables([("X", 1:1)])
        H = [1.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 1.0]

        gns_svd = gns_reconstruct(H, reg, 2; method=:svd, hankel_deg=1, atol=1e-12)
        gns_chol = gns_reconstruct(H, reg, 2; method=:cholesky, hankel_deg=1, atol=1e-12)

        for k in 0:4
            mono = k == 0 ? one(x[1]) : only(monomials(x[1]^k))
            moment_svd = _pipeline_expectation(gns_svd.matrices, gns_svd.xi, mono)
            moment_chol = _pipeline_expectation(gns_chol.matrices, gns_chol.xi, mono)
            @test abs(moment_svd - moment_chol) ≤ 1e-8
        end
    end

    @testset "Dense wrappers and method dispatch" begin
        reg, (x,) = create_noncommutative_variables([("X", 1:1)])
        H = [1.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 1.0]
        monomap = Dict(
            _pipeline_moment_key(one(x[1])) => 1.0,
            _pipeline_moment_key(only(monomials(x[1]))) => 0.0,
            _pipeline_moment_key(only(monomials(x[1]^2))) => 1.0,
            _pipeline_moment_key(only(monomials(x[1]^3))) => 0.0,
            _pipeline_moment_key(only(monomials(x[1]^4))) => 1.0,
        )

        gns_from_monomap = gns_reconstruct(monomap, reg, 2; method=:cholesky, hankel_deg=1, atol=1e-12)
        matrices_only = reconstruct(H, reg, 2; hankel_deg=1, atol=1e-12)

        @test gns_from_monomap.matrices[reg[:X₁]] ≈ [0.0 1.0; 1.0 0.0] atol = 1e-12
        @test matrices_only[reg[:X₁]] ≈ [0.0 1.0; 1.0 0.0] atol = 1e-12
        @test_throws ArgumentError gns_reconstruct(H, reg, 2; method=:bogus, hankel_deg=1)
        @test_throws ArgumentError gns_reconstruct(monomap, reg, 2; method=:bogus, hankel_deg=1)
    end

    @testset "Sparse GNS rejects inactive registry variables" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        pop = polyopt(1.0 + x[1]^2, reg)
        cfg = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        sparsity = compute_sparsity(pop, cfg)

        @test length(sparsity.corr_sparsity.cliques) == 1
        @test sparsity.corr_sparsity.cliques[1] == [reg[:x₁]]

        err = try
            gns_reconstruct(Dict{Any,Float64}(), sparsity, 2; hankel_deg=1)
            nothing
        catch caught
            caught
        end

        @test err isa ArgumentError
        @test occursin("every registry variable", sprint(showerror, err))
        @test occursin("9", sprint(showerror, err))
        @test occursin("13", sprint(showerror, err))
    end

    ball = _solve_ball_example()

    @testset "Ball extraction" begin
        @test ball.result.objective ≈ 1.0 atol = 1e-5

        flatness = test_flatness(ball.hankel_flat, ball.full_basis, ball.basis; atol=1e-6)
        gns = gns_reconstruct(
            ball.hankel_flat,
            ball.reg,
            3;
            method=:cholesky,
            hankel_deg=2,
            atol=1e-6,
        )

        A1 = gns.matrices[ball.reg[:X₁]]
        A2 = gns.matrices[ball.reg[:X₂]]
        F = Matrix(2I - A1^2 + A1 * A2^2 * A1 - A2^2)
        F = Symmetric((F + F') / 2)
        λmin = eigmin(F)
        identity_matrix = Matrix{eltype(A1)}(I, size(A1, 1), size(A1, 1))
        ball_matrix = Symmetric((identity_matrix - A1^2 - A2^2 + (identity_matrix - A1^2 - A2^2)') / 2)

        @test flatness.is_flat
        @test gns.rank == 5
        @test size(A1) == (5, 5)
        @test size(A2) == (5, 5)
        @test λmin ≈ 1.0 atol = 1e-5
        @test real(dot(gns.xi, F * gns.xi)) ≈ 1.0 atol = 1e-5
        @test eigmin(ball_matrix) ≥ -1e-6
    end

    @testset "Robustness diagnostics" begin
        gns = gns_reconstruct(
            ball.hankel_flat,
            ball.reg,
            3;
            method=:cholesky,
            hankel_deg=2,
            atol=1e-6,
        )
        report = robustness_report(gns, ball.hankel, ball.full_basis, ball.basis)
        report_convenience = robustness_report(gns, ball.hankel)

        @test report.sigma_min > 0
        @test report.sigma_max ≥ report.sigma_min
        @test isfinite(report.condition_number)
        @test report.dist_to_flat ≤ 0.1
        @test report.operator_error_bound ≤ 5.0
        @test report_convenience.dist_to_flat ≈ report.dist_to_flat atol = 1e-12
        @test report_convenience.operator_error_bound ≈ report.operator_error_bound atol = 1e-12
    end

    @testset "Verification suite" begin
        gns = gns_reconstruct(
            ball.hankel_flat,
            ball.reg,
            3;
            method=:cholesky,
            hankel_deg=2,
            atol=1e-6,
        )
        report = verify_gns(
            gns,
            ball.result.monomap,
            ball.reg;
            poly=ball.f,
            f_star=ball.result.objective,
            constraints=[ball.g],
            ball=true,
            atol=5e-5,
        )

        @test report.is_symmetric
        @test report.moment_max_error ≤ 5e-5
        @test report.objective_error ≤ 5e-5
        @test all(v -> v >= -5e-5, report.constraint_min_eigenvalues)
        @test report.ball_contained
    end

    # Example 4.22 on the nc ball is already covered above. The following
    # regression cases exercise the remaining dense GNS geometries from the
    # design notes without duplicating that baseline test.

    @testset "Polydisc extraction" begin
        reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
        X, Y = vars

        f = X * Y * X
        g1 = 1.0 - X^2
        g2 = 1.0 - Y^2
        pop = polyopt(f, reg; ineq_constraints=[g1, g2])
        problem = _solve_dense_example(pop, reg, 3; H_deg=3, hankel_deg=2, atol=1e-6)
        gns = gns_reconstruct(
            problem.hankel_flat,
            reg,
            3;
            method=:svd,
            hankel_deg=2,
            atol=1e-6,
        )

        A1 = gns.matrices[reg[:X₁]]
        A2 = gns.matrices[reg[:X₂]]
        F = Symmetric((A1 * A2 * A1 + (A1 * A2 * A1)') / 2)
        identity_matrix = Matrix{eltype(A1)}(I, size(A1, 1), size(A1, 1))
        g1_matrix = Symmetric((identity_matrix - A1^2 + (identity_matrix - A1^2)') / 2)
        g2_matrix = Symmetric((identity_matrix - A2^2 + (identity_matrix - A2^2)') / 2)
        report = verify_gns(
            gns,
            problem.result.monomap,
            reg;
            poly=f,
            f_star=problem.result.objective,
            constraints=[g1, g2],
            atol=5e-5,
        )

        @test problem.result.objective ≈ -1.0 atol = 1e-5
        @test !problem.flatness.is_flat
        @test problem.flatness_flat.is_flat
        @test gns.rank == 5
        @test size(A1) == (5, 5)
        @test size(A2) == (5, 5)
        @test eigmin(F) ≈ -1.0 atol = 1e-5
        @test real(dot(gns.xi, F * gns.xi)) ≈ -1.0 atol = 1e-5
        @test eigmin(g1_matrix) ≥ -5e-5
        @test eigmin(g2_matrix) ≥ -5e-5
        @test report.is_symmetric
        @test report.moment_max_error ≤ 5e-5
        @test report.objective_error ≤ 5e-5
        @test all(v -> v >= -5e-5, report.constraint_min_eigenvalues)
    end

    @testset "Mixed equality and inequality extraction" begin
        reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
        X, Y = vars

        f = 2.0 - X^2 + X * Y^2 * X - Y^2
        g = 4.0 - X^2 - Y^2
        h = X * Y + Y * X - 2.0
        pop = polyopt(f, reg; eq_constraints=[h], ineq_constraints=[g])
        # COSMO reaches a flat dense extraction for this mixed-constraint example
        # already at order 3, so keep the test at that smaller stable instance.
        problem = _solve_dense_example(pop, reg, 3; H_deg=3, hankel_deg=2, atol=1e-6)
        gns = gns_reconstruct(
            problem.hankel_flat,
            reg,
            3;
            method=:svd,
            hankel_deg=2,
            atol=1e-6,
        )

        A1 = gns.matrices[reg[:X₁]]
        A2 = gns.matrices[reg[:X₂]]
        identity_matrix = Matrix{eltype(A1)}(I, size(A1, 1), size(A1, 1))
        h_matrix = A1 * A2 + A2 * A1 - 2 * identity_matrix
        report = verify_gns(
            gns,
            problem.result.monomap,
            reg;
            poly=f,
            f_star=problem.result.objective,
            constraints=[g],
            atol=5e-5,
        )

        @test problem.result.objective ≈ -1.0 atol = 1e-5
        @test problem.flatness.is_flat
        @test problem.flatness_flat.is_flat
        @test report.is_symmetric
        @test report.moment_max_error ≤ 5e-5
        @test report.objective_error ≤ 5e-5
        @test all(v -> v >= -5e-5, report.constraint_min_eigenvalues)
        @test norm(h_matrix) ≤ 1e-5
    end

    @testset "Motzkin quartic ball extraction" begin
        reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
        X, Y = vars

        f = X * Y^4 * X + Y * X^4 * Y - 3.0 * X * Y^2 * X + 1.0
        g = 1.0 - X^4 - Y^4
        pop = polyopt(f, reg; ineq_constraints=[g])
        problem = _solve_dense_example(pop, reg, 3; H_deg=3, hankel_deg=2, atol=1e-6)
        gns = gns_reconstruct(
            problem.hankel_flat,
            reg,
            3;
            method=:svd,
            hankel_deg=2,
            atol=1e-6,
        )

        A1 = gns.matrices[reg[:X₁]]
        A2 = gns.matrices[reg[:X₂]]
        identity_matrix = Matrix{eltype(A1)}(I, size(A1, 1), size(A1, 1))
        f_operator = A1 * A2^4 * A1 + A2 * A1^4 * A2 - 3 * A1 * A2^2 * A1 + identity_matrix
        g_operator = identity_matrix - A1^4 - A2^4
        f_matrix = Symmetric((f_operator + f_operator') / 2)
        g_matrix = Symmetric((g_operator + g_operator') / 2)
        report = verify_gns(
            gns,
            problem.result.monomap,
            reg;
            poly=f,
            f_star=problem.result.objective,
            constraints=[g],
            atol=5e-5,
        )

        @test problem.result.objective ≈ -0.160813 atol = 1e-3
        @test !problem.flatness.is_flat
        @test problem.flatness_flat.is_flat
        @test real(dot(gns.xi, f_matrix * gns.xi)) ≈ problem.result.objective atol = 1e-5
        @test eigmin(g_matrix) ≥ -5e-5
        @test report.is_symmetric
        @test report.moment_max_error ≤ 5e-5
        @test report.objective_error ≤ 5e-5
        @test all(v -> v >= -5e-5, report.constraint_min_eigenvalues)
    end

    @testset "Unbounded dense problem does not yield a fake optimizer" begin
        reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
        X, Y = vars

        f = 2.0 - X^2 + X * Y^2 * X - Y^2
        pop = polyopt(f, reg)
        cfg = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        sparsity = compute_sparsity(pop, cfg)
        mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)

        outcome = try
            result = NCTSSoS.solve_moment_problem(mp, SOLVER)
            NCTSSoS._check_solver_status(result.model)
            result
        catch caught
            caught
        end

        if outcome isa NCTSSoS.SolverStatusError
            @test outcome.termination in (MOI.DUAL_INFEASIBLE, MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED) ||
                  outcome.primal == MOI.INFEASIBILITY_CERTIFICATE ||
                  outcome.dual == MOI.INFEASIBILITY_CERTIFICATE
        else
            problem = _dense_hankel_data(outcome.monomap, reg, 2, 1; atol=1e-6)
            @test !problem.flatness.is_flat
        end
    end
end

using Test, NCTSSoS, LinearAlgebra

using NCTSSoS: get_ncbasis, monomials

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

function _solve_ball_example()
    reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
    X, Y = vars

    f = 2.0 - X^2 + X * Y^2 * X - Y^2
    g = 1.0 - X^2 - Y^2
    pop = polyopt(f, reg; ineq_constraints=[g])
    cfg = SolverConfig(
        optimizer=SOLVER,
        order=3,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
    )

    sparsity = compute_sparsity(pop, cfg)
    mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    result = NCTSSoS.solve_moment_problem(mp, SOLVER)

    full_basis = get_ncbasis(reg, 3)
    basis = get_ncbasis(reg, 2)
    hankel = NCTSSoS.hankel_matrix(result.monomap, full_basis)
    hankel_flat = flat_extend(hankel, full_basis, basis; atol=1e-6)

    return (
        reg=reg,
        X=X,
        Y=Y,
        f=f,
        g=g,
        result=result,
        full_basis=full_basis,
        basis=basis,
        hankel=hankel,
        hankel_flat=hankel_flat,
    )
end

@testset "Dense GNS Pipeline" begin
    @testset "Flat extension correctness" begin
        reg, (x,) = create_noncommutative_variables([("X", 1:1)])
        full_basis = get_ncbasis(reg, 2)
        basis = get_ncbasis(reg, 1)
        H = [1.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 1.0]

        flatness = test_flatness(H, full_basis, basis; atol=1e-12)
        Hflat = flat_extend(H, full_basis, basis; atol=1e-12)

        @test flatness.is_flat
        @test flatness.rank_principal == 2
        @test flatness.rank_full == 2
        @test flatness.err_flat ≤ 1e-12
        @test rank(Hflat) == 2
        @test Hflat ≈ H atol = 1e-12
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

        @test report.sigma_min > 0
        @test report.sigma_max ≥ report.sigma_min
        @test isfinite(report.condition_number)
        @test report.dist_to_flat ≤ 0.1
        @test report.operator_error_bound ≤ 5.0
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
end

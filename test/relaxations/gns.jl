using Test, NCTSSoS, LinearAlgebra, COSMO

if !@isdefined(SOLVER)
    const SOLVER = optimizer_with_attributes(
        COSMO.Optimizer,
        "verbose" => false,
        "eps_abs" => 1e-7,
        "eps_rel" => 1e-7,
    )
end

function example_27_fixture()
    reg, (theta,) = create_noncommutative_variables([("theta", 1:2)])
    hankel = [
        1.00001 0.499907 0.500102 1.0483 -0.5483 -0.5483 1.0484;
        0.499907 1.0483 -0.548283 1.0627 -0.0144 -0.6090 0.0606;
        0.500102 -0.548283 1.04827 -0.0144 -0.5340 0.0606 0.9878;
        1.0483 1.0627 -0.0144 1.4622 -0.3995 -0.8006 0.7863;
        -0.5483 -0.0144 -0.5340 -0.3995 0.3852 0.1917 -0.7256;
        -0.5483 -0.6090 0.0606 -0.8006 0.1917 0.4411 -0.3804;
        1.0484 0.0606 0.9878 0.7863 -0.7256 -0.3804 1.3682;
    ]
    return reg, theta, hankel
end

@testset "GNS Construction" begin
    @testset "Dense Hankel Reconstruction" begin
        reg, theta, hankel = example_27_fixture()
        model = gns_reconstruct(hankel, reg, 2; atol=0.1, rtol=1e-6)

        @test model.report.flat
        @test model.report.rank_full == 2
        @test model.report.rank_quotient == 2
        @test model.report.max_moment_error < 1e-3
        @test model.quotient_basis == [one(theta[1]), theta[1]]
        @test model.cyclic_vector ≈ [1.0, 0.0] atol = 1e-12

        theta1_idx = theta[1].word[1]
        theta2_idx = theta[2].word[1]
        @test model.operators[theta1_idx] * model.cyclic_vector ≈ [0.0, 1.0] atol = 1e-9
        @test model.operators[theta2_idx] * model.cyclic_vector ≈ [1.0, -1.0] atol = 1e-3
        @test model.operators[theta1_idx]' * model.gram ≈ model.gram * model.operators[theta1_idx] atol = 1e-12
        @test model.operators[theta2_idx]' * model.gram ≈ model.gram * model.operators[theta2_idx] atol = 1e-12

        gram_factor = cholesky(Symmetric(model.gram)).U
        paper_sign = Diagonal([1.0, -1.0])
        theta1_paper = paper_sign * (gram_factor * model.operators[theta1_idx] / gram_factor) * paper_sign
        theta2_paper = paper_sign * (gram_factor * model.operators[theta2_idx] / gram_factor) * paper_sign

        @test theta1_paper ≈ [0.5019 -0.8931; -0.8931 0.1727] atol = 3e-3
        @test theta2_paper ≈ [0.4981 0.8939; 0.8939 0.0825] atol = 3e-3

        wrapped = reconstruct(hankel, reg, 2; atol=0.1, rtol=1e-6)
        @test wrapped[theta1_idx] ≈ model.operators[theta1_idx] atol = 1e-9
        @test wrapped[theta2_idx] ≈ model.operators[theta2_idx] atol = 1e-9
    end

    @testset "Non-Flat Report" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])
        hankel = Matrix{Float64}(I, 3, 3)
        model = gns_reconstruct(hankel, reg, 2; atol=1e-8, rtol=1e-8)

        @test !model.report.flat
        @test model.report.rank_full == 3
        @test model.report.rank_quotient == 2
        @test_throws ArgumentError gns_reconstruct(hankel, reg, 2; atol=1e-8, rtol=1e-8, require_flat=true)
    end

    @testset "Direct Solve Plumbing" begin
        reg, (u,) = create_unipotent_variables([("u", 1:1)])
        pop = polyopt(-1.0 * u[1], reg)

        direct = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=1); dualize=false)
        sol = dense_moment_solution(direct)
        model = gns_reconstruct(direct; atol=1e-7, rtol=1e-7)

        @test direct.objective ≈ -1.0 atol = 1e-6
        @test sol.order == 1
        @test sol.basis == [one(u[1]), u[1]]
        @test model.report.flat
        @test model.report.max_moment_error == 0.0
        @test model.quotient_basis == [one(u[1])]
        @test model.cyclic_vector == [1.0]
        @test model.operators[u[1].word[1]][1, 1] ≈ 1.0 atol = 1e-6

        dualized = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=1); dualize=true)
        @test_throws ArgumentError dense_moment_solution(dualized)
    end
end

using Test
using NCTSSoS
using LinearAlgebra

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

@testset "Moment Matrix Extraction" begin
    @testset "CHSH inequality" begin
        @ncpolyvar x[1:2]  # x = (A_1, A_2)
        @ncpolyvar y[1:2]  # y = (B_1, B_2)

        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(f, comm_gps=[x, y], is_unipotent=true)

        @testset "Primal moment formulation" begin
            solver_config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(pop, solver_config, dualize=false)

            moments = result.moment_matrices

            expected_moment = [1.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 -√2/2 -√2/2; 0.0 0.0 1.0 -√2/2 √2/2; 0.0 -√2/2 -√2/2 1.0 0.0; 0.0 -√2/2 √2/2 0.0 1.0]

            @test length(moments) == 1
            @test length(moments[1]) == 1
            @test moments[1][1] ≈ expected_moment atol = 1e-6
        end

        @testset "Primal moment formulation with Term Sparsity" begin
            solver_config = SolverConfig(optimizer=SOLVER, order=1, ts_algo=MMD())
            result = cs_nctssos(pop, solver_config, dualize=false)

            moments = result.moment_matrices

            @test length(moments) == 1
            @test length(moments[1]) == 3

            # TODO: NCTSSOS Gives different blocks than CliqueTrees.jl no way of comparing esaily
            expected_moment1 = [1.0 0.0 -√2/2; 0.0 1.0 -√2/2; -√2/2 -√2/2 1.0]
            expected_moment2 = [1.0 0.0 -√2/2; 0.0 1.0 √2/2; -√2/2 √2/2 1.0]

            @test moments[1][1] ≈ expected_moment1 atol = 1e-6
            @test moments[1][2] ≈ expected_moment2 atol = 1e-6
            @test moments[1][3] ≈ [1;;] atol = 1e-6
        end


        # Test with dual formulation (default)
        @testset "Dual SOS formulation" begin
            solver_config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(pop, solver_config, dualize=true)

            moments = result.moment_matrices

            expected_moment = [1.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 -√2/2 -√2/2; 0.0 0.0 1.0 -√2/2 √2/2; 0.0 -√2/2 -√2/2 1.0 0.0; 0.0 -√2/2 √2/2 0.0 1.0]

            @test moments[1][1] ≈ expected_moment atol = 1e-6
        end

        @testset "Dual SOS TS Formulation" begin
            solver_config = SolverConfig(optimizer=SOLVER, order=1, ts_algo=MMD())
            result = cs_nctssos(pop, solver_config, dualize=true)

            moments = result.moment_matrices

            expected_moment1 = [1.0 0.0 -√2/2; 0.0 1.0 -√2/2; -√2/2 -√2/2 1.0]
            expected_moment2 = [1.0 0.0 -√2/2; 0.0 1.0 √2/2; -√2/2 √2/2 1.0]

            @test moments[1][1] ≈ expected_moment1 atol = 1e-6
            @test moments[1][2] ≈ expected_moment2 atol = 1e-6
            @test moments[1][3] ≈ [1;;] atol = 1e-6
        end
    end

    @testset "Pauli algebra problem" begin
        @ncpolyvar x y z

        # Simple objective with Pauli constraints
        obj = x^2 + y^2 + z^2

        # Pauli algebra: σ_i^2 = 1, σ_i σ_j = -σ_j σ_i (i ≠ j)
        pop = polyopt(obj; eq_constraints=[x^2 - 1, y^2 - 1, z^2 - 1],
            comm_gps=[[x], [y], [z]])

        solver_config = SolverConfig(optimizer=Clarabel.Optimizer, order=1)
        result = cs_nctssos(pop, solver_config)

        # Extract moment matrices
        moments = get_moment_matrices(result)

        @test moments isa Vector
        @test length(moments) > 0
        @test moments[1] isa Vector
        @test moments[1][1] isa Matrix

        println("✓ Pauli algebra test passed")
        println("  Objective value: ", result.objective)
        println("  Number of cliques: ", length(moments))
    end
end

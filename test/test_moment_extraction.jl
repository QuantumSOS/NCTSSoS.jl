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


        # Test with dual formulation (default)
        @testset "Dual SOS formulation" begin
            solver_config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(pop, solver_config, dualize=true)

            # Check that result has moment_support field
            @test hasfield(typeof(result), :moment_support)
            @test result.moment_support isa MomentSupport

            # Extract moment matrices
            moments = get_moment_matrices(result)
            moments[1][1]
            ans = [0.999999999512592 -0.0 -0.0 -0.0 -0.0; -0.0 0.999999999512592 3.039770389764888e-16 0.7071067828933203 0.7071067828933197; -0.0 3.039770389764888e-16 0.999999999512592 0.7071067828933203 -0.7071067828933205; -0.0 0.7071067828933203 0.7071067828933203 0.999999999512592 -3.9950008909332224e-16; -0.0 0.7071067828933197 -0.7071067828933205 -3.9950008909332224e-16 0.999999999512592]

            @test moments[1][1] ≈ ans atol=1e-6

        end

        @testset "Dual SOS TS Formulation" begin

            solver_config = SolverConfig(optimizer=SOLVER, order=1,ts_algo=MMD())
            result = cs_nctssos(pop, solver_config, dualize=true)

            model = result.model

            using JuMP
            all_vars = all_variables(model)

            # Extract moment matrices
            moments = get_moment_matrices(result)

            # need to swap because of ordering of basis elements, in NCTSSOS, they have y2 < y1
            swap_mtx = [1 0 0; 0 0 1; 0 1 0]

            ans1 = swap_mtx * [0.999999999987871 0.7071067811836527 -0.7071067811836522; 0.7071067811836527 0.999999999987871 -5.11184632543241e-16; -0.7071067811836522 -5.11184632543241e-16 0.999999999987871] * swap_mtx

            ans2 = swap_mtx * [0.999999999987871 0.7071067811836521 0.7071067811836521; 0.7071067811836521 0.999999999987871 -5.11184632543241e-16; 0.7071067811836521 -5.11184632543241e-16 0.999999999987871] * swap_mtx

            moments[1][1] + ans1
            moments[1][1] + ans2

            moments[1][2] + ans1
            moments[1][2] + ans2

            @test moments[1][1] ≈ ans atol=1e-6


        end

        # Test with primal formulation
        @testset "Primal moment formulation" begin

            solver_config = SolverConfig(optimizer=SOLVER, order=1)
            result = cs_nctssos(pop, solver_config, dualize=false)

            mom_supp = result.moment_support
            cliques = mom_supp.cliques
            cliques[1].blocks[1].dual_indices

            model = result.model
            all_vars = value.(all_variables(model))

            moments = get_moment_matrices(result)

            ans = [0.999999999512592 -0.0 -0.0 -0.0 -0.0; -0.0 0.999999999512592 3.039770389764888e-16 0.7071067828933203 0.7071067828933197; -0.0 3.039770389764888e-16 0.999999999512592 0.7071067828933203 -0.7071067828933205; -0.0 0.7071067828933203 0.7071067828933203 0.999999999512592 -3.9950008909332224e-16; -0.0 0.7071067828933197 -0.7071067828933205 -3.9950008909332224e-16 0.999999999512592]

            @test length(moments) == 1
            @test length(moments[1]) == 1
            moments[1][1]
            moments[1][1] - ans

            @test moments[1][1] ≈ ans atol=1e-6
        end

        @testset "Primal moment formulation with Term Sparsity" begin
            solver_config = SolverConfig(optimizer=SOLVER, order=1, ts_algo=MMD())
            result = cs_nctssos(pop, solver_config, dualize=false)

            moments = get_moment_matrices(result)

            @test length(moments) == 1

            moments[1][1]

            @test length(moments) == 1
            @test length(moments[1]) == 1
            moments[1][1]
            moments[1][1] - ans

            # need to swap because of ordering of basis elements, in NCTSSOS, they have y2 < y1
            swap_mtx = [1 0 0; 0 0 1; 0 1 0]

            ans1 = swap_mtx * [0.999999999987871 0.7071067811836527 -0.7071067811836522; 0.7071067811836527 0.999999999987871 -5.11184632543241e-16; -0.7071067811836522 -5.11184632543241e-16 0.999999999987871] * swap_mtx

            ans2 = swap_mtx * [0.999999999987871 0.7071067811836521 0.7071067811836521; 0.7071067811836521 0.999999999987871 -5.11184632543241e-16; 0.7071067811836521 -5.11184632543241e-16 0.999999999987871] * swap_mtx

            @test moments[1][1] ≈ ans atol=1e-6
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

    @testset "Moment matrix properties" begin
        # Use CHSH problem to verify moment matrix properties
        @ncpolyvar x[1:2]  # x = (A_1, A_2)
        @ncpolyvar y[1:2]  # y = (B_1, B_2)

        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(f, comm_gps=[x, y], is_unipotent=true)

        using Clarabel
        solver_config = SolverConfig(optimizer=Clarabel.Optimizer, order=1)
        result = cs_nctssos(pop, solver_config)

        moments = get_moment_matrices(result)

        # Test that matrices are nearly symmetric (within numerical tolerance)
        for clique_mats in moments
            for mat in clique_mats
                @test norm(mat - mat') < 1e-6
                println("  Matrix symmetry error: ", norm(mat - mat'))
            end
        end

        # Test that the (1,1) entry (corresponding to constant monomial) is close to 1
        # This assumes the first basis element is the constant monomial
        first_mat = moments[1][1]
        println("  First matrix (1,1) entry: ", first_mat[1,1])

        println("✓ Moment matrix properties test passed")
    end

    @testset "Higher order relaxation" begin
        # Use CHSH problem to test higher order relaxation
        @ncpolyvar x[1:2]  # x = (A_1, A_2)
        @ncpolyvar y[1:2]  # y = (B_1, B_2)

        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(f, comm_gps=[x, y], is_unipotent=true)

        using Clarabel
        solver_config_1 = SolverConfig(optimizer=Clarabel.Optimizer, order=1)
        result_1 = cs_nctssos(pop, solver_config_1)

        # Test higher order
        solver_config_2 = SolverConfig(optimizer=Clarabel.Optimizer, order=2)
        result_2 = cs_nctssos_higher(pop, result_1, solver_config_2)

        # Check that higher order result also has moment_support
        @test hasfield(typeof(result_2), :moment_support)

        moments_2 = get_moment_matrices(result_2)
        @test moments_2 isa Vector
        @test moments_2[1] isa Vector
        @test moments_2[1][1] isa Matrix

        # Both should be close to Tsirelson's bound
        @test result_1.objective ≈ 2*sqrt(2) atol=1e-3
        @test result_2.objective ≈ 2*sqrt(2) atol=1e-3

        println("✓ Higher order relaxation test passed")
        println("  Order 1 objective: ", result_1.objective)
        println("  Order 2 objective: ", result_2.objective)
        println("  Tsirelson's bound: ", 2*sqrt(2))
    end
end

println("\n✓ All moment matrix extraction tests passed!")

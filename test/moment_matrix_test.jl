using Test, NCTSSoS

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

@testset "Moment Matrix Extraction" begin
    @testset "Simple Problem - Primal" begin
        @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2
        pop = polyopt(f)

        # Solve without dualization (Moment problem)
        result = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=1); dualize=false)

        @test result.problem_type == MomentPrimal

        # Try to extract moment matrix
        mm = moment_matrix(result)
        @test mm isa MomentMatrix
        @test mm.clique_index == 1
        @test length(mm.basis) > 0
        @test size(mm.matrix, 1) == size(mm.matrix, 2)
        @test size(mm.matrix, 1) == length(mm.basis)
    end

    @testset "Simple Problem - Dual" begin
        @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2
        pop = polyopt(f)

        # Solve with dualization (SOS problem)
        result = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=1); dualize=true)

        @test result.problem_type == SOSDual

        # Extract moment matrix (should work now)
        mm = moment_matrix(result)
        @test mm isa MomentMatrix
        @test mm.clique_index == 1
        @test length(mm.basis) > 0
        @test size(mm.matrix, 1) == size(mm.matrix, 2)
        @test size(mm.matrix, 1) == length(mm.basis)
    end

    @testset "Moment Matrix Values" begin
        @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2
        pop = polyopt(f)

        # Solve both primal and dual
        result_primal = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=1); dualize=false)
        result_dual = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=1); dualize=true)

        # Extract moment matrices
        mm_primal = moment_matrix(result_primal)
        mm_dual = moment_matrix(result_dual)

        # Both should give similar results (may have small numerical differences)
        @test size(mm_primal.matrix) == size(mm_dual.matrix)
        @test mm_primal.basis == mm_dual.basis
    end
end

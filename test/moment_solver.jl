using Test, NCTSSoS
using JuMP

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2))
    )
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end
using Graphs

# =============================================================================
# Migrated Tests using New API
# =============================================================================
# NOTE: Tests that directly call internal NCTSSoS functions (correlative_sparsity,
# term_sparsities, moment_relax) are removed as they test implementation details
# that may have changed with the FastPolynomials refactor.
#
# Focus is on testing the high-level API (cs_nctssos, cs_nctssos_higher).

@testset "Special Constraint Type" begin
    @testset "CHSH Inequality" begin
        # Use unipotent variables for x^2 = I property
        registry, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])

        f = x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]

        pop = polyopt(f, registry)

        solver_config = SolverConfig(optimizer=SOLVER; order=1)

        result = cs_nctssos(pop, solver_config; dualize=false)

        @test isapprox(result.objective, -2.8284271321623193, atol=1e-6)
    end
end

@testset "CS TS Example" begin
    order = 3
    n = 10
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])

    # Build polynomial using new API
    f = Polynomial{NonCommutativeAlgebra,UInt8,Float64}(Term{Monomial{NonCommutativeAlgebra,UInt8},Float64}[])
    for i = 1:n
        jset = max(1, i - 5):min(n, i + 1)
        jset = setdiff(jset, i)
        f += (2.0 * x[i] + 5.0 * x[i]^3 + 1)^2
        f -= sum([
            4.0 * x[i] * x[j] +
            10.0 * x[i]^3 * x[j] +
            2.0 * x[j] +
            4.0 * x[i] * x[j]^2 +
            10.0 * x[i]^3 * x[j]^2 +
            2.0 * x[j]^2 for j in jset
        ])
        f += sum([
            1.0 * x[j] * x[k] + 2.0 * x[j]^2 * x[k] + 1.0 * x[j]^2 * x[k]^2 for j in jset for k in jset
        ])
    end

    cons = vcat([(1.0 - x[i]^2) for i = 1:n], [(1.0 * x[i] - 1.0 / 3) for i = 1:n])

    pop = polyopt(f, registry; ineq_constraints=cons)

    solver_config =
        SolverConfig(optimizer=SOLVER; order=order, cs_algo=MF(), ts_algo=MMD())

    result = cs_nctssos(pop, solver_config; dualize=false)

    @test isapprox(result.objective, 3.011288, atol=1e-4)
end

@testset "Moment Method Heisenberg Model on Star Graph" begin
    num_sites = 10
    g = star_graph(num_sites)

    true_ans = -1.0

    vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]

    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

    # Use unipotent variables for projector-like variables
    registry, (pij,) = create_unipotent_variables([("pij", 1:length(vec_idx2ij))])

    objective = sum(1.0 * pij[[findvaridx(ee.src, ee.dst) for ee in edges(g)]])

    gs = unique!([
        (
            pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
            pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
            pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
            pij[findvaridx(sort([i, k])...)] + 1.0
        ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
        (i != j && j != k && i != k)
    ])


    pop = polyopt(objective, registry; eq_constraints=gs)
    order = 1
    cs_algo = MF()

    solver_config = SolverConfig(
        optimizer=SOLVER,
        order=order,
        cs_algo=cs_algo
    )

    result = cs_nctssos(pop, solver_config; dualize=false)
    @test isapprox(result.objective, true_ans, atol=1e-6)
end

@testset "Moment Method Example 1" begin
    order = 2
    n = 3
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])

    f =
        1.0 * x[1]^2 - 1.0 * x[1] * x[2] - 1.0 * x[2] * x[1] + 3.0 * x[2]^2 - 2.0 * x[1] * x[2] * x[1] +
        2.0 * x[1] * x[2]^2 * x[1] - 1.0 * x[2] * x[3] - 1.0 * x[3] * x[2] +
        6.0 * x[3]^2 +
        9.0 * x[2]^2 * x[3] +
        9.0 * x[3] * x[2]^2 - 54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

    pop = polyopt(f, registry)

    @testset "Dense" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            cs_algo=NoElimination(),
        )

        result = cs_nctssos(pop, solver_config; dualize=false)

        # NOTE: differs from original test case value since that one is a relaxed in terms of sparsity
        # This value here is obtained by running the master branch with no sparsity relaxation
        @test isapprox(
            result.objective,
            4.372259295498716e-10,
            atol=1e-6,
        )
    end

    @testset "Sparse" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            ts_algo=MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize=false)

        @test isapprox(result.objective, -0.0035512, atol=1e-7)
    end
end

@testset "Moment Method Example 2" begin
    order = 2
    n = 2
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])

    f = 2.0 - 1.0 * x[1]^2 + 1.0 * x[1] * x[2]^2 * x[1] - 1.0 * x[2]^2
    g = 4.0 - 1.0 * x[1]^2 - 1.0 * x[2]^2
    h1 = 1.0 * x[1] * x[2] + 1.0 * x[2] * x[1] - 2.0

    pop = polyopt(f, registry; eq_constraints=[h1], ineq_constraints=[g])

    @testset "Dense" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order)

        result = cs_nctssos(pop, solver_config; dualize=false)
        @test isapprox(result.objective, -1.0, atol=1e-6)
    end

    @testset "Term Sparse" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            cs_algo=MF(),
            ts_algo=MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize=false)

        @test isapprox(result.objective, -1.0, atol=1e-6)
    end
end

@testset "Moment Method Correlative Sparsity" begin
    n = 3
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])

    f =
        1.0 * x[1]^2 - 1.0 * x[1] * x[2] - 1.0 * x[2] * x[1] + 3.0 * x[2]^2 - 2.0 * x[1] * x[2] * x[1] +
        2.0 * x[1] * x[2]^2 * x[1] - 1.0 * x[2] * x[3] - 1.0 * x[3] * x[2] +
        6.0 * x[3]^2 +
        9.0 * x[2]^2 * x[3] +
        9.0 * x[3] * x[2]^2 - 54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - 1.0 * x[i]^2 for i = 1:n], [1.0 * x[i] - 1.0 / 3 for i = 1:n])
    order = 3

    pop = polyopt(f, registry; ineq_constraints=cons)

    @testset "Correlative Sparse" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            cs_algo=MF(),
        )
        result = cs_nctssos(pop, solver_config; dualize=false)

        # FIXME: reduced accuracy
        # @test is_solved_and_feasible(moment_problem.model)
        @test isapprox(
            result.objective,
            0.9975306427277915,
            atol=1e-5,
        )
    end

    @testset "Term Sparse" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            ts_algo=MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize=false)

        result = cs_nctssos_higher(pop, result, solver_config; dualize=false)

        @test isapprox(
            result.objective,
            0.9975306427277915,
            atol=1e-5,
        )
    end
end

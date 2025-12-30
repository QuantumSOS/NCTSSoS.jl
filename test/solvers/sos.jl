# SOS (Sum of Squares) Dualization Tests
# =======================================
# Tests the SOS dualization approach to polynomial optimization.

using Test, NCTSSoS
using SparseArrays, Graphs, CliqueTrees
using NCTSSoS: get_Cαj

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "setup.jl"))

if LOCAL_TESTING
    @testset "I_3322 inequality" begin
        # Use projector algebra for x and y (projective = P² = P)
        registry, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
        f =
            1.0 * x[1] * (y[1] + y[2] + y[3]) +
            x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]
        pop = polyopt(-f, registry)

        solver_config = SolverConfig(optimizer=SOLVER; order=3)

        result = cs_nctssos(pop, solver_config)

        @test isapprox(result.objective, -0.2508753049688358, atol=1e-6)
    end
end

# This test requires high solver precision - only run with LOCAL_TESTING
if LOCAL_TESTING
    @testset "CS TS Example" begin
        order = 3
        n = 10
        registry, (x,) = create_noncommutative_variables([("x", 1:n)])
        # Start with a zero polynomial (0.0 * x[1] creates polynomial, then - x[1] gives zero)
        f = 1.0 * x[1] - 1.0 * x[1]  # zero polynomial
        for i = 1:n
            jset = max(1, i - 5):min(n, i + 1)
            jset = setdiff(jset, i)
            f += (2x[i] + 5 * x[i]^3 + 1)^2
            f -= sum([
                4x[i] * x[j] +
                10x[i]^3 * x[j] +
                2x[j] +
                4x[i] * x[j]^2 +
                10x[i]^3 * x[j]^2 +
                2x[j]^2 for j in jset
            ])
            f += sum([
                x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
            ])
        end

        cons = vcat([(1 - x[i]^2) for i = 1:n], [(x[i] - 1 / 3) for i = 1:n])

        pop = polyopt(f, registry; ineq_constraints=cons)
        solver_config = SolverConfig(optimizer=SOLVER, order=order,
            cs_algo=MF(), ts_algo=MMD())

        result = cs_nctssos(pop, solver_config)

        @test isapprox(result.objective, 3.011288, atol=1e-3)
    end
end

@testset "Cαj" begin
    # Test get_Cαj with new API: Vector{Monomial}, Matrix{Polynomial}
    registry, (x,) = create_noncommutative_variables([("x", 1:3)])

    # Create a properly typed polynomial matrix
    poly_matrix = [
        1.0*x[1] 1.0*x[2]
        1.0*x[2] 1.0*x[1]+1.0*x[3]
    ]

    # Basis: the monomials themselves (sorted)
    basis_monomials = sort([one(x[1]), x[1], x[2], x[3]])

    C_α_js = get_Cαj(basis_monomials, poly_matrix)

    # Just verify structure - we have 5 non-zero entries
    @test length(C_α_js) == 5
    @test all(v -> v == 1.0, values(C_α_js))
end

@testset "Cαj complex" begin
    registry, (x,) = create_noncommutative_variables([("x", 1:2)])
    basis_polys = get_ncbasis(registry, 2)
    # For NC algebra, each basis poly is single-term; extract the monomial
    basis_monomials = [monomials(p)[1] for p in basis_polys]

    # Create properly typed polynomial matrix
    localizing_mtx = [
        1.0*x[1]-1.0*x[2]   1.0*x[2]^2-1.0*one(x[1])
        1.0*x[2]^2-1.0*one(x[1])   1.0*x[2]^3
    ]
    C_α_js = get_Cαj(basis_monomials, localizing_mtx)

    @test !isempty(C_α_js)
    @test all(v -> v == 1.0 || v == -1.0, values(C_α_js))
end

@testset "Dualization Trivial Example" begin
    n = 2
    true_min = 3.0
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])

    f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min

    pop = polyopt(f, registry)
    order = 2

    solver_config = SolverConfig(
        optimizer=SOLVER,
        order=order
    )

    result = cs_nctssos(pop, solver_config; dualize=true)

    @test isapprox(result.objective, true_min, atol=1e-6)
end

# This test requires high precision solver - COSMO gives Inf for one method
if LOCAL_TESTING
    @testset "Dualization Trivial Example 2" begin
        n = 2
        true_min = 3.0
        registry, (x,) = create_noncommutative_variables([("x", 1:n)])

        f = x[1]^2 + x[1] * x[2] + x[2] * x[1] + x[2]^2 + true_min
        r = -10.0
        g1 = r - x[1]
        g2 = r - x[2]
        g3 = x[1] - r
        g4 = x[2] - r

        pop = polyopt(f, registry; ineq_constraints=[g1, g2, g3, g4])
        order = 2

        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order
        )

        result_mom = cs_nctssos(pop, solver_config; dualize=false)
        result_sos = cs_nctssos(pop, solver_config; dualize=true)

        @test isapprox(result_mom.objective, result_sos.objective, atol=1e-3)
    end
end

@testset "Dualization Example 1" begin
    n = 3
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])

    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] +
        2.0 * x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    pop = polyopt(f, registry)
    order = 2

    solver_config = SolverConfig(
        optimizer=SOLVER,
        order=order,
    )

    result = cs_nctssos(pop, solver_config; dualize=true)

    @test isapprox(result.objective, 0.0, atol=1e-6)
end

@testset "Dualization Example 2" begin
    n = 2
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    pop = polyopt(f, registry; ineq_constraints=[g], eq_constraints=[h1])

    order = 2

    @testset "Dense" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
        )

        result = cs_nctssos(pop, solver_config; dualize=true)
        @test isapprox(result.objective, -1, atol=1e-6)
    end

    @testset "Term Sparse" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            ts_algo=MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize=true)

        @test isapprox(result.objective, -1.0, atol=1e-6)
    end
end

@testset "Dualization Heisenberg Model on Star Graph" begin
    num_sites = 8
    star = star_graph(num_sites)

    true_ans = -1.0

    vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]

    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

    # Use unipotent algebra (U² = I)
    registry, (pij,) = create_unipotent_variables([("p", 1:length(vec_idx2ij))])

    objective = sum(1.0 * pij[[findvaridx(ee.src, ee.dst) for ee in edges(star)]])

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

    solver_config = SolverConfig(
        optimizer=SOLVER,
        order=order,
    )

    result = cs_nctssos(pop, solver_config; dualize=true)

    @test isapprox(result.objective, true_ans, atol=1e-6)
end

@testset "SOS Method Correlative Sparsity" begin
    n = 3
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])
    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] +
        2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0 / 3 for i = 1:n])
    order = 3

    pop = polyopt(f, registry; ineq_constraints=cons)

    @testset "Correlative Sparsity" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            cs_algo=MF(),
        )

        result = cs_nctssos(pop, solver_config; dualize=true)

        @test isapprox(result.objective, 0.9975306427277915, atol=1e-5)
    end

    @testset "Term Sparsity" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            ts_algo=MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize=true)

        # Term sparsity gives same bound as correlative sparsity for this problem
        # Mosek achieves 1.42e-7 error, COSMO needs 1e-4 tolerance
        @test isapprox(result.objective, 0.9975306427277915, atol=1e-4)
    end
end

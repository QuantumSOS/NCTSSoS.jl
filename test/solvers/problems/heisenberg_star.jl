# =============================================================================
# Heisenberg Model on Star Graph Tests
# =============================================================================
# Consolidates all Heisenberg star graph tests:
#   - Moment method
#   - SOS dualization
#   - Various graph sizes
#
# The Heisenberg model on a star graph with singlet projection constraints.
# Expected optimal value: -1.0 (ground state energy per edge)
# Results verified against NCTSSOS oracles.
# =============================================================================

using Test, NCTSSoS
using Graphs

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__DIR__), "..", "setup.jl"))

# Load oracle values
include(joinpath(dirname(@__DIR__), "..", "oracles", "results", "heisenberg_star_oracles.jl"))

# Helper: flatten moment_matrix_sizes for comparison with oracle
flatten_sizes(sizes) = reduce(vcat, sizes)

const HEISENBERG_STAR_EXPECTED = -1.0

@testset "Heisenberg Star Graph" begin

    # =========================================================================
    # Problem Setup Helper
    # =========================================================================
    function create_heisenberg_star_problem(num_sites::Int)
        g = star_graph(num_sites)

        vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]
        findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

        # Use unipotent variables for projector-like variables
        reg, (pij,) = create_unipotent_variables([("pij", 1:length(vec_idx2ij))])

        objective = sum(1.0 * pij[[findvaridx(ee.src, ee.dst) for ee in edges(g)]])

        # Singlet projection constraints
        gs = unique!([
            (
                pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
                pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
                pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
                pij[findvaridx(sort([i, k])...)] + 1.0
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ])

        pop = polyopt(objective, reg; eq_constraints=gs)
        return pop, reg
    end

    # =========================================================================
    # Dense (No Sparsity, n=10) - HeisenbergStar_Dense_n10_d1
    # =========================================================================
    @testset "Dense (n=10, order=1)" begin
        pop, _ = create_heisenberg_star_problem(10)
        oracle = HEISENBERG_STAR_ORACLES["HeisenbergStar_Dense_n10_d1"]

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )

        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    # =========================================================================
    # Correlative Sparsity MF (n=10) - HeisenbergStar_CS_n10_d1
    # =========================================================================
    # NOTE: Block order may differ from NCTSSOS due to clique ordering.
    # We compare sorted block sizes instead of exact order.
    # =========================================================================
    @testset "Correlative Sparsity MF (n=10, order=1)" begin
        pop, _ = create_heisenberg_star_problem(10)
        oracle = HEISENBERG_STAR_ORACLES["HeisenbergStar_CS_n10_d1"]

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )

        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-6
        # Block sizes match (order may differ due to clique enumeration)
        @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    # =========================================================================
    # Dense (n=8) - HeisenbergStar_Dense_n8_n8_d1
    # =========================================================================
    @testset "Dense (n=8, order=1)" begin
        pop, _ = create_heisenberg_star_problem(8)
        oracle = HEISENBERG_STAR_ORACLES["HeisenbergStar_Dense_n8_n8_d1"]

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )

        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    # =========================================================================
    # SOS Dualization (n=8)
    # =========================================================================
    @testset "SOS Dualization (n=8)" begin
        pop, _ = create_heisenberg_star_problem(8)
        oracle = HEISENBERG_STAR_ORACLES["HeisenbergStar_Dense_n8_n8_d1"]

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )

        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ oracle.opt atol = 1e-6
    end

    # =========================================================================
    # Majumdar-Ghosh Model (ring with J1-J2 interactions)
    # =========================================================================
    @testset "Majumdar-Ghosh Model (n=6)" begin
        num_sites = 6
        J1_interactions =
            unique!([tuple(sort([i, mod1(i + 1, num_sites)])...) for i = 1:num_sites])
        J2_interactions =
            unique!([tuple(sort([i, mod1(i + 2, num_sites)])...) for i = 1:num_sites])

        J1 = 2.0
        J2 = 1.0
        true_ans = -num_sites / 4 * 6

        ij2idx_dict = Dict(
            zip(
                [(i, j) for i in 1:num_sites, j in 1:num_sites if j > i],
                1:(num_sites*(num_sites-1)÷2),
            ),
        )

        # Use projector algebra (P² = P)
        reg, (hij,) = create_projector_variables([("h", 1:(num_sites*(num_sites-1)÷2))])

        objective = (
            sum([J1 * hij[ij2idx_dict[(i, j)]] for (i, j) in J1_interactions]) +
            sum([J2 * hij[ij2idx_dict[(i, j)]] for (i, j) in J2_interactions])
        )

        gs = unique!([
            (
                hij[ij2idx_dict[tuple(sort([i, j])...)]] *
                hij[ij2idx_dict[tuple(sort([j, k])...)]] +
                hij[ij2idx_dict[tuple(sort([j, k])...)]] *
                hij[ij2idx_dict[tuple(sort([i, j])...)]] -
                0.5 * (
                    hij[ij2idx_dict[tuple(sort([i, j])...)]] +
                    hij[ij2idx_dict[tuple(sort([j, k])...)]] -
                    hij[ij2idx_dict[tuple(sort([i, k])...)]]
                )
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ])

        pop = polyopt(-objective, reg; eq_constraints=gs)

        config = SolverConfig(optimizer=SOLVER, order=1)
        result = cs_nctssos(pop, config)

        @test result.objective ≈ true_ans atol = 1e-4
    end
end

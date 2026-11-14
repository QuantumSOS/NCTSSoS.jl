# Heisenberg Model on Star Graph Tests
# Consolidates all Heisenberg star graph tests:
#   - Moment method
#   - SOS dualization
#   - Various graph sizes
#
# The Heisenberg model on a star graph with singlet projection constraints.
# The unipotent variables represent swaps, not singlet projectors.
# For n sites, sum_j SWAP(1,j) = (n-1)/2 + 2 S_1 ⋅ S_leaves.
# Coupling the maximal leaf spin (n-1)/2 to spin 1/2 gives minimum -1
# for the TOTAL objective, not the energy per edge.

using Test, NCTSSoS, JuMP
import Graphs

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

const HEISENBERG_STAR_EXPECTED = -1.0

@testset "Heisenberg Star Graph" begin

    # =========================================================================
    # Problem Setup Helper
    # =========================================================================
    function create_heisenberg_star_problem(num_sites::Int)
        g = Graphs.star_graph(num_sites)

        vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]
        findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

        # Swap operators satisfy U² = I.
        reg, (pij,) = create_unipotent_variables([("pij", 1:length(vec_idx2ij))])

        objective = sum(1.0 * pij[[findvaridx(ee.src, ee.dst) for ee in Graphs.edges(g)]])

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
    # Dense (No Sparsity, n=10)
    # =========================================================================
    @testset "Dense (n=10, order=1)" begin
        pop, _ = create_heisenberg_star_problem(10)

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )

        result = cs_nctssos(pop, config; dualize=false)
        @test is_solved_and_feasible(result.model)
        @test result.objective ≈ HEISENBERG_STAR_EXPECTED atol = 1e-6
    end

    # =========================================================================
    # Correlative Sparsity MF (n=10)
    # =========================================================================
    @testset "Correlative Sparsity MF (n=10, order=1)" begin
        pop, _ = create_heisenberg_star_problem(10)

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )

        result = cs_nctssos(pop, config; dualize=false)
        @test is_solved_and_feasible(result.model)
        @test result.objective ≈ HEISENBERG_STAR_EXPECTED atol = 1e-6
    end

    # =========================================================================
    # Dense (n=8)
    # =========================================================================
    @testset "Dense (n=8, order=1)" begin
        pop, _ = create_heisenberg_star_problem(8)

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )

        result = cs_nctssos(pop, config; dualize=false)
        @test is_solved_and_feasible(result.model)
        @test result.objective ≈ HEISENBERG_STAR_EXPECTED atol = 1e-6
    end

    # =========================================================================
    # SOS Dualization (n=8)
    # =========================================================================
    @testset "SOS Dualization (n=8)" begin
        pop, _ = create_heisenberg_star_problem(8)

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )

        result = cs_nctssos(pop, config; dualize=true)
        @test is_solved_and_feasible(result.model)
        @test result.objective ≈ HEISENBERG_STAR_EXPECTED atol = 1e-6
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
        # A product of nearest-neighbor singlets attains the MG minimum.
        # -sum J_ij P_singlet = sum J_ij S_i⋅S_j - sum J_ij/4 = -9.
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

        @test is_solved_and_feasible(result.model)
        @test result.objective ≈ true_ans atol = 1e-4
    end
end

# test/relaxations/sparsity.jl
# Non-correlated sparsity tests:
# - term sparsity initialization
# - PolyOptResult structural fields

using Test, NCTSSoS, JuMP

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

@testset "Term Sparsity" begin
    @testset "init_activated_supp" begin
        registry, (x,) = create_noncommutative_variables([("x", 1:2)])
        T = eltype(indices(registry))
        M = NormalMonomial{NonCommutativeAlgebra,T}
        P = Polynomial{NonCommutativeAlgebra,T,Float64}

        x_idx = [x[i].word[1] for i in 1:2]
        objective = x[1] + x[2]
        constraints = P[]

        one_mono = one(M)
        basis = M[one_mono, M([x_idx[1]]), M([x_idx[2]])]

        supp = NCTSSoS.init_activated_supp(objective, constraints, basis)
        @test !isempty(supp)
        @test length(supp) >= 2
    end
end

@testset "PolyOptResult Fields" begin
    @testset "moment_matrix_sizes and n_unique_moment_matrix_elements" begin
        oracle = expectations_oracle("expectations/relaxations_sparsity.json", "unipotent_n2_order1_dense")
        reg, (u,) = create_unipotent_variables([("u", 1:2)])
        objective = 1.0 * u[1] * u[2]
        pop = polyopt(objective, reg)

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "minimal single variable" begin
        oracle = expectations_oracle("expectations/relaxations_sparsity.json", "unipotent_n1_order1_dense")
        reg, (u,) = create_unipotent_variables([("u", 1:1)])
        objective = 1.0 * u[1]
        pop = polyopt(objective, reg)

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config)

        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end

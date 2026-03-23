# test/correlated_sparsity/runtests.jl
# Correlated sparsity test entrypoint and shared helpers.

using Test, NCTSSoS, JuMP, Graphs

using NCTSSoS:
    assign_constraint,
    clique_decomp,
    correlative_sparsity,
    get_correlative_graph

if !@isdefined(MOI)
    const MOI = JuMP.MOI
end

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

const CORRELATED_PIPELINE_ORACLES = (
    CS_d3 = expectations_oracle("expectations/correlated_pipeline.toml", "CS_d3"),
    TS_d3 = expectations_oracle("expectations/correlated_pipeline.toml", "TS_d3"),  # sides vary
    E8_sparse_d2 = expectations_oracle("expectations/correlated_pipeline.toml", "E8_sparse_d2"),
    E8_sparse_d3 = expectations_oracle("expectations/correlated_pipeline.toml", "E8_sparse_d3"),
    E8_dense_d2 = expectations_oracle("expectations/correlated_pipeline.toml", "E8_dense_d2"),
    E9_dense_n4_d2 = expectations_oracle("expectations/correlated_pipeline.toml", "E9_dense_n4_d2"),
    E9_sparse_n4_d2 = expectations_oracle("expectations/correlated_pipeline.toml", "E9_sparse_n4_d2"),
    E9_sparse_n8_d2 = expectations_oracle("expectations/correlated_pipeline.toml", "E9_sparse_n8_d2"),
    E9_sparse_n12_d2 = expectations_oracle("expectations/correlated_pipeline.toml", "E9_sparse_n12_d2"),
)

const CORRELATED_STRUCTURE_EXPECTATIONS =
    TestExpectations.expectations_load("expectations/correlated_structure.toml")

function correlated_structure_case(id::AbstractString)
    case = TestExpectations.expectations_case(CORRELATED_STRUCTURE_EXPECTATIONS, id)
    haskey(case, "expected") || error("Missing key `expected` for case $(repr(id)).")
    return case["expected"]
end

if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

var_indices(vars) = [v.word[1] for v in vars]
normalize_cliques(cliques) = sort(sort.(cliques))

function clique_symbol_names(registry, cliques)
    normalize_cliques([[string(registry[idx]) for idx in clique] for clique in cliques])
end

function clique_variable_positions(vars, cliques)
    position_by_index = Dict(var.word[1] => i for (i, var) in enumerate(vars))
    normalize_cliques([[position_by_index[idx] for idx in clique] for clique in cliques])
end

function graph_adjacency_by_index(graph::SimpleGraph, sorted_indices)
    Dict(
        sorted_indices[i] => sort([sorted_indices[j] for j in neighbors(graph, i)])
        for i in 1:nv(graph)
    )
end

function nc_chain_objective(x)
    1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
    2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
    x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
    9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
    54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]
end

function nc_bipartite_objective(x, y)
    1.0 * x[1] * (y[1] + y[2] + y[3]) +
    x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2.0 * y[1] - y[2]
end

function nc_large_scale_objective(x)
    n = length(x)
    poly_type = typeof(x[1] + x[2])
    objective = zero(poly_type)
    for i in 1:n
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)
        objective += (2.0 * x[i] + 5.0 * x[i]^3 + 1.0)^2
        objective -= sum([
            4.0 * x[i] * x[j] +
            10.0 * x[i]^3 * x[j] +
            2.0 * x[j] +
            4.0 * x[i] * x[j]^2 +
            10.0 * x[i]^3 * x[j]^2 +
            2.0 * x[j]^2 for j in jset
        ])
        objective += sum([
            x[j] * x[k] + 2.0 * x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
        ])
    end
    objective
end

function build_nc_correlative_problem()
    n = 3
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    objective = nc_chain_objective(x)
    constraints = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0 / 3 for i = 1:n])
    pop = polyopt(objective, reg; ineq_constraints=constraints)
    return pop
end

function build_global_constraint_fixture()
    reg, (x,) = create_noncommutative_variables([("x", 1:2)])
    objective = 1.0 * x[1]^2 + 1.0 * x[2]^2 + 0.5 * x[1] * x[2] + 0.5 * x[2] * x[1]
    constraint = 1.0 * x[1] + x[2]
    pop = polyopt(objective, reg; ineq_constraints=[constraint])
    config = SolverConfig(
        optimizer=SOLVER,
        order=1,
        cs_algo=NoElimination(),
        ts_algo=NoElimination()
    )
    sparsity = compute_sparsity(pop, config)
    corr = sparsity.corr_sparsity
    corr_with_global = typeof(corr)(
        corr.cliques,
        corr.registry,
        corr.cons,
        corr.clq_cons,
        [1],
        corr.clq_mom_mtx_bases,
        corr.clq_localizing_mtx_bases
    )
    return pop, sparsity, corr_with_global
end

function build_e8_polyball_problem()
    reg, (x,) = create_noncommutative_variables([("X", 1:4)])
    X1, X2, X3, X4 = x

    f1 =
        4.0 - X1 + 3.0 * X2 - 3.0 * X3 - 3.0 * X1^2 - 7.0 * X1 * X2 + 6.0 * X1 * X3 -
        X2 * X1 - 5.0 * X3 * X1 + 5.0 * X3 * X2 - 5.0 * X1^3 - 3.0 * X1^2 * X3 +
        4.0 * X1 * X2 * X1 - 6.0 * X1 * X2 * X3 + 7.0 * X1 * X3 * X1 +
        2.0 * X1 * X3 * X2 - X1 * X3^2 - X2 * X1^2 + 3.0 * X2 * X1 * X2 -
        X2 * X1 * X3 - 2.0 * X2^3 - 5.0 * X2^2 * X3 - 4.0 * X2 * X3^2 -
        5.0 * X3 * X1^2 + 7.0 * X3 * X1 * X2 + 6.0 * X3 * X2 * X1 -
        4.0 * X3 * X2^2 - X3^2 * X1 - 2.0 * X3^2 * X2 + 7.0 * X3^3

    f2 =
        -1.0 + 6.0 * X2 + 5.0 * X3 + 3.0 * X4 - 5.0 * X2^2 + 2.0 * X2 * X3 +
        4.0 * X2 * X4 - 4.0 * X3 * X2 + X3^2 - X3 * X4 + X4 * X2 - X4 * X3 +
        2.0 * X4^2 - 7.0 * X2^3 + 4.0 * X2 * X3^2 + 5.0 * X2 * X3 * X4 -
        7.0 * X2 * X4 * X3 - 7.0 * X2 * X4^2 + X3 * X2^2 + 6.0 * X3 * X2 * X3 -
        6.0 * X3 * X2 * X4 - 3.0 * X3^2 * X2 - 7.0 * X3^2 * X4 + 6.0 * X3 * X4 * X2 -
        3.0 * X3 * X4 * X3 - 7.0 * X3 * X4^2 + 3.0 * X4 * X2^2 - 7.0 * X4 * X2 * X3 -
        X4 * X2 * X4 - 5.0 * X4 * X3^2 + 7.0 * X4 * X3 * X4 + 6.0 * X4^2 * X2 -
        4.0 * X4^3

    # Examples 5.10/5.15 are stated as a self-adjoint eigenvalue benchmark.
    # The tabulated f₁ and f₂ list one orientation of the mixed words, so the
    # benchmark objective is reconstructed as f₁ + f₁† + f₂ + f₂†. Using only
    # f₁ + f₂ cuts the reported optimum in half and no longer matches the paper.
    objective = f1 + adjoint(f1) + f2 + adjoint(f2)
    constraints = [1.0 - X1^2 - X2^2 - X3^2, 1.0 - X2^2 - X3^2 - X4^2]
    return polyopt(objective, reg; ineq_constraints=constraints)
end

function build_e9_chained_singular_polydisc_problem(n::Int)
    n >= 4 || throw(ArgumentError("E9 requires n ≥ 4."))
    iseven(n) || throw(ArgumentError("E9 is defined on even n with 4-variable chained blocks."))

    reg, (x,) = create_noncommutative_variables([("X", 1:n)])
    poly_type = typeof(x[1] + x[2])
    objective = zero(poly_type)

    # Keep the literature SOHS formulation from the chained-singular benchmark.
    # Replacing these quartic pieces by naive NC fourth powers changes the
    # support and no longer matches the published correlative-sparsity instance.
    for i in 1:2:(n - 3)
        linear1 = x[i] + 10.0 * x[i + 1]
        linear2 = x[i + 2] - x[i + 3]
        quadratic1 = x[i + 1]^2 - 4.0 * x[i + 1] * x[i + 2] + 4.0 * x[i + 2]^2
        quadratic2 = x[i]^2 - 20.0 * x[i] * x[i + 3] + 100.0 * x[i + 3]^2

        objective += linear1' * linear1
        objective += 5.0 * (linear2' * linear2)
        objective += quadratic1' * quadratic1
        objective += 10.0 * (quadratic2' * quadratic2)
    end

    constraints = vcat([1.0 - x[i]^2 for i in 1:n], [x[i] - 1.0 / 3 for i in 1:n])
    return reg, x, polyopt(objective, reg; ineq_constraints=constraints)
end

@testset "Correlated Sparsity" begin
    include("graph_and_cliques.jl")
    include("core_pipeline_structure.jl")
    include("core_pipeline_numeric.jl")
    include("coverage_edges.jl")
    include("e9_chained_singular_polydisc.jl")
    include("literature_matrix.jl")
end

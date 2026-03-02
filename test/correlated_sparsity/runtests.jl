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
    CS_d3 = expectations_oracle("expectations/correlated_pipeline.json", "CS_d3"),
    TS_d3 = expectations_oracle("expectations/correlated_pipeline.json", "TS_d3"),  # sides vary
)

const CORRELATED_STRUCTURE_EXPECTATIONS =
    TestExpectations.expectations_load("expectations/correlated_structure.json")

function correlated_structure_case(id::AbstractString)
    case = TestExpectations.expectations_case(CORRELATED_STRUCTURE_EXPECTATIONS, id)
    haskey(case, "expected") || error("Missing key `expected` for case $(repr(id)).")
    return case["expected"]
end

json_int(v) = Int(v)
json_int_vec(values) = [Int(v) for v in values]
json_int_vec_vec(values) = [json_int_vec(row) for row in values]

if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

var_indices(vars) = [v.word[1] for v in vars]
normalize_cliques(cliques) = sort(sort.(cliques))

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

@testset "Correlated Sparsity" begin
    include("graph_and_cliques.jl")
    include("core_pipeline_structure.jl")
    include("core_pipeline_numeric.jl")
    include("coverage_edges.jl")
    include("literature_matrix.jl")
end

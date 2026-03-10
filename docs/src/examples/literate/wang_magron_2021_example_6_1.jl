# # [Wang-Magron 2021 Example 6.1: Fast Sparse REPL Checks](@id wang-magron-2021-example-6-1)
#
# This file turns the fast, fixture-backed part of Section 6.1 from
# Wang and Magron's NCTSSOS paper [wangExploitingTermSparsity2021](@cite)
# into manual REPL checks.
#
# Nothing executes when the file is included. That keeps it REPL-friendly:
#
# ```julia
# include("docs/src/examples/literate/wang_magron_2021_example_6_1.jl")
# using COSMO, JuMP
#
# cosmo = optimizer_with_attributes(
#     COSMO.Optimizer,
#     "verbose" => false,
#     "eps_abs" => 1e-7,
#     "eps_rel" => 1e-7,
# )
#
# test_table_3_sparse_n20(optimizer=cosmo)
# test_table_4_sparse_n20(optimizer=cosmo)
# test_table_5_sparse_n20(optimizer=cosmo)
# test_table_6_sparse_n20(optimizer=cosmo)
# test_all_fast_cases(optimizer=cosmo)
# ```
#
# The tracked JSON expectations reproduce the published sparse rows for
# Tables 3--6. Table 2 is omitted from the callable surface because its
# first published sparse row is too slow for the warm-REPL budget used by
# this example, and Table 7 is outside the fixture-backed scope.

# ## Setup
#
# The default backend is silent Mosek when `MosekTools` is available in the
# active environment. Tests and local smoke runs can pass a different
# `optimizer=` keyword, for example a silent COSMO optimizer.

using NCTSSoS, Printf, JSON3

const MOI = NCTSSoS.MOI
const EXAMPLE_6_1_DATA_DIR = joinpath(
    pkgdir(NCTSSoS), "docs", "src", "examples", "literate", "data"
)
const EXAMPLE_6_1_HAS_MOSEKTOOLS = Base.find_package("MosekTools") !== nothing

if EXAMPLE_6_1_HAS_MOSEKTOOLS
    import MosekTools
end

function example_6_1_flatten_sizes(sizes)
    return reduce(vcat, sizes)
end

function example_6_1_max_block_size(result)
    return maximum(example_6_1_flatten_sizes(result.moment_matrix_sizes))
end

function example_6_1_default_optimizer()
    EXAMPLE_6_1_HAS_MOSEKTOOLS || error(
        "MosekTools is not available in the active environment. " *
        "Pass `optimizer=` explicitly, for example a silent COSMO optimizer.",
    )
    return MOI.OptimizerWithAttributes(MosekTools.Optimizer, MOI.Silent() => true)
end

function example_6_1_resolve_optimizer(optimizer)
    return optimizer === nothing ? example_6_1_default_optimizer() : optimizer
end

function example_6_1_table_expectation_path(table::Int)
    3 <= table <= 6 || error("Only Tables 3--6 are fixture-backed in this example.")
    return joinpath(
        EXAMPLE_6_1_DATA_DIR, @sprintf("wang_magron_2021_example_6_1_table_%d.json", table)
    )
end

function example_6_1_load_table_expectation(table::Int)
    path = example_6_1_table_expectation_path(table)
    isfile(path) || error("Missing expectation JSON: $(repr(path))")
    return JSON3.read(read(path, String))
end

function example_6_1_sparse_reference_row(table::Int, n::Int=20)
    data = example_6_1_load_table_expectation(table)["expected_values"]
    for row in data["sparse_rows"]
        Int(row["n"]) == n || continue
        return (
            table=table,
            family=String(data["family"]),
            n=Int(row["n"]),
            paper=Float64(row["opt"]),
            paper_mb=Int(row["mb"]),
        )
    end
    error("No sparse expectation row for Table $(table) at n=$(n).")
end

function example_6_1_default_atol(table::Int)
    return table == 6 ? 1e-6 : 1e-3
end

function example_6_1_summarize_rows(rows)
    io = IOBuffer()
    println(
        io,
        @sprintf(
            "%-7s %-8s %-4s %12s %10s %10s %5s %8s %8s",
            "table",
            "family",
            "n",
            "objective",
            "paper",
            "abs err",
            "mb",
            "nuniq",
            "sec",
        ),
    )
    println(io, repeat("-", 84))
    for row in rows
        println(
            io,
            @sprintf(
                "%-7s %-8s %-4d %12.6f %10.4f %10.3e %5d %8d %8.2f",
                "T$(row.table)",
                row.family,
                row.n,
                row.objective,
                row.paper,
                row.abs_error,
                row.mb,
                row.n_unique,
                row.seconds,
            ),
        )
    end
    return String(take!(io))
end

function example_6_1_build_row(table::Int, result, expected, seconds::Real)
    return (
        table=table,
        family=expected.family,
        n=expected.n,
        objective=result.objective,
        paper=expected.paper,
        abs_error=abs(result.objective - expected.paper),
        mb=example_6_1_max_block_size(result),
        paper_mb=expected.paper_mb,
        n_unique=result.n_unique_moment_matrix_elements,
        seconds=round(seconds; digits=2),
    )
end

function example_6_1_verify_row(row, atol::Real)
    @assert isapprox(row.objective, row.paper; atol=atol) (
        "Objective mismatch for Table $(row.table): got $(row.objective), " *
        "expected $(row.paper) with atol=$(atol)."
    )
    @assert row.mb == row.paper_mb (
        "Moment block size mismatch for Table $(row.table): got $(row.mb), " *
        "expected $(row.paper_mb)."
    )
    return nothing
end

function example_6_1_run_case(table::Int, pop, config; n::Int=20, atol::Real=example_6_1_default_atol(table))
    expected = example_6_1_sparse_reference_row(table, n)
    seconds = @elapsed result = cs_nctssos(pop, config)
    row = example_6_1_build_row(table, result, expected, seconds)
    example_6_1_verify_row(row, atol)
    return (row=row, result=result)
end

# ## Problem Builders
#
# The builders below encode the mathematical formulations documented in the
# source notes for Tables 3--6. They are hand-transcribed into package-native
# Julia rather than parsed from Markdown snippets.

function build_table_3_problem(n::Int)
    n >= 4 || error("Table 3 requires n >= 4.")
    iseven(n) || error("Table 3 expects an even n.")
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    P = typeof(x[1] + x[2])
    f = zero(P)

    for i in 1:2:(n - 3)
        linear_1 = x[i] + 10.0 * x[i + 1]
        linear_2 = x[i + 2] - x[i + 3]

        left_quad_1 = x[i + 1]^2 - 4.0 * x[i + 2] * x[i + 1] + 4.0 * x[i + 2]^2
        right_quad_1 = x[i + 1]^2 - 4.0 * x[i + 1] * x[i + 2] + 4.0 * x[i + 2]^2
        left_quad_2 = x[i]^2 - 20.0 * x[i + 3] * x[i] + 100.0 * x[i + 3]^2
        right_quad_2 = x[i]^2 - 20.0 * x[i] * x[i + 3] + 100.0 * x[i + 3]^2

        f += linear_1 * linear_1
        f += 5.0 * linear_2 * linear_2
        f += left_quad_1 * right_quad_1
        f += 10.0 * left_quad_2 * right_quad_2
    end

    return polyopt(f, reg)
end

function build_table_4_problem(n::Int)
    n >= 2 || error("Table 4 requires n >= 2.")
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    f = Float64(n) * one(typeof(x[1]))

    for i in 2:n
        f += 100.0 * x[i - 1]^4 - 200.0 * x[i - 1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
    end

    return polyopt(f, reg)
end

function build_table_5_problem(n::Int)
    n >= 4 || error("Table 5 requires n >= 4.")
    n % 4 == 0 || error("Table 5 expects n divisible by 4.")
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    P = typeof(x[1] + x[2])
    f = one(P)

    for i in 1:2:(n - 3)
        a = x[i + 1] - x[i]^2
        b = 1.0 - x[i]
        c = x[i + 3] - x[i + 2]^2
        d = 1.0 - x[i + 2]
        e = x[i + 1] + x[i + 3] - 2.0
        g = x[i + 1] - x[i + 3]

        f += 100.0 * a * a
        f += b * b
        f += 90.0 * c * c
        f += d * d
        f += 10.0 * e * e
        f += 0.1 * g * g
    end

    return polyopt(f, reg)
end

function build_table_6_problem(n::Int)
    n >= 2 || error("Table 6 requires n >= 2.")
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    a = 3.0 * x[1] - 2.0 * x[1]^2 - 2.0 * x[2] + 1.0
    f = a * a

    for i in 2:(n - 1)
        b = 3.0 * x[i] - 2.0 * x[i]^2 - x[i - 1] - 2.0 * x[i + 1] + 1.0
        f += b * b
    end

    c = 3.0 * x[n] - 2.0 * x[n]^2 - x[n - 1] + 1.0
    f += c * c

    return polyopt(f, reg)
end

function table_3_sparse_config(pop; optimizer)
    basis = newton_chip_basis(pop, 2)
    return SolverConfig(
        optimizer=optimizer,
        moment_basis=basis,
        cs_algo=MF(),
        ts_algo=MMD(),
    )
end

function table_4_sparse_config(pop; optimizer)
    basis = newton_chip_basis(pop, 2)
    return SolverConfig(
        optimizer=optimizer,
        moment_basis=basis,
        cs_algo=MF(),
        ts_algo=MMD(),
    )
end

function table_5_sparse_config(pop; optimizer)
    return SolverConfig(
        optimizer=optimizer,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=MMD(),
    )
end

function table_6_sparse_config(pop; optimizer)
    basis = newton_chip_basis(pop, 2)
    return SolverConfig(
        optimizer=optimizer,
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=MMD(),
    )
end

function example_6_1_execute_table(
    table::Int,
    build_problem::Function,
    build_config::Function;
    optimizer=nothing,
    n::Int=20,
    atol::Real=example_6_1_default_atol(table),
    io::IO=stdout,
    print_summary::Bool=true,
)
    chosen_optimizer = example_6_1_resolve_optimizer(optimizer)
    pop = build_problem(n)
    report = example_6_1_run_case(
        table,
        pop,
        build_config(pop; optimizer=chosen_optimizer);
        n=n,
        atol=atol,
    )
    if print_summary
        print(io, example_6_1_summarize_rows([report.row]))
    end
    return report
end

# ## REPL Entry Points
#
# Each `test_table_*` function runs the first sparse published row (`n = 20`)
# for its table, checks the objective and maximal block size against the
# tracked JSON oracle, prints a one-row summary, and returns `(row, result)`.

function test_table_3_sparse_n20(; optimizer=nothing, atol::Real=example_6_1_default_atol(3), io::IO=stdout)
    return example_6_1_execute_table(
        3,
        build_table_3_problem,
        table_3_sparse_config;
        optimizer=optimizer,
        n=20,
        atol=atol,
        io=io,
    )
end

function test_table_4_sparse_n20(; optimizer=nothing, atol::Real=example_6_1_default_atol(4), io::IO=stdout)
    return example_6_1_execute_table(
        4,
        build_table_4_problem,
        table_4_sparse_config;
        optimizer=optimizer,
        n=20,
        atol=atol,
        io=io,
    )
end

function test_table_5_sparse_n20(; optimizer=nothing, atol::Real=example_6_1_default_atol(5), io::IO=stdout)
    return example_6_1_execute_table(
        5,
        build_table_5_problem,
        table_5_sparse_config;
        optimizer=optimizer,
        n=20,
        atol=atol,
        io=io,
    )
end

function test_table_6_sparse_n20(; optimizer=nothing, atol::Real=example_6_1_default_atol(6), io::IO=stdout)
    return example_6_1_execute_table(
        6,
        build_table_6_problem,
        table_6_sparse_config;
        optimizer=optimizer,
        n=20,
        atol=atol,
        io=io,
    )
end

function test_all_fast_cases(; optimizer=nothing, io::IO=stdout)
    chosen_optimizer = example_6_1_resolve_optimizer(optimizer)
    reports = [
        example_6_1_execute_table(
            3,
            build_table_3_problem,
            table_3_sparse_config;
            optimizer=chosen_optimizer,
            n=20,
            atol=example_6_1_default_atol(3),
            io=io,
            print_summary=false,
        ),
        example_6_1_execute_table(
            4,
            build_table_4_problem,
            table_4_sparse_config;
            optimizer=chosen_optimizer,
            n=20,
            atol=example_6_1_default_atol(4),
            io=io,
            print_summary=false,
        ),
        example_6_1_execute_table(
            5,
            build_table_5_problem,
            table_5_sparse_config;
            optimizer=chosen_optimizer,
            n=20,
            atol=example_6_1_default_atol(5),
            io=io,
            print_summary=false,
        ),
        example_6_1_execute_table(
            6,
            build_table_6_problem,
            table_6_sparse_config;
            optimizer=chosen_optimizer,
            n=20,
            atol=example_6_1_default_atol(6),
            io=io,
            print_summary=false,
        ),
    ]
    print(io, example_6_1_summarize_rows([report.row for report in reports]))
    return reports
end

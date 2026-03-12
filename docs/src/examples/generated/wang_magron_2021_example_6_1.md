```@meta
EditURL = "../literate/wang_magron_2021_example_6_1.jl"
```

# [Wang-Magron 2021 Example 6.1: Reproducing `n = 20` Results](@id wang-magron-2021-example-6-1)

This example reproduces the `n = 20` rows of Tables 2--6 from Section 6.1
of Wang and Magron's NCTSSOS paper [wangExploitingTermSparsity2021](@cite),
comparing both dense and sparse relaxations.

Including this file defines the functions below but does not run any solves.
To run individual cases or the full suite interactively:

```julia
include("docs/src/examples/literate/wang_magron_2021_example_6_1.jl")

test_table_4_dense_n20()              # single case (Mosek auto-detected)
test_table_4_sparse_n20()
test_all_n20_reference_cases()        # all 10 cases at once
```

Each function checks the computed objective and maximal moment-matrix block
size (`mb`) against tracked JSON expectations derived from the paper.
Table 7 is not covered.

## Recorded Mosek results (`n = 20`)

The table below summarises the silent-Mosek solves recorded in
`wang_magron_2021_example_6_1_mosek_n20_results.json`:

| Case | Lower bound | Paper | Abs. error | `mb` | Paper `mb` | Seconds |
|:-----------|-------------------:|------:|-----------:|-----:|-----------:|--------:|
| T3 dense   |  4.72e-4  | -1e-4 | 5.72e-4 | 59 | 59 | 11.14 |
| T3 sparse  | -2.45e-5  | -4e-4 | 3.75e-4 |  3 |  3 |  0.58 |
| T4 dense   |  1.000000 |  1.0  | 3.17e-12 | 40 | 40 |  0.06 |
| T4 sparse  |  1.000000 |  1.0  | 8.31e-8  |  3 |  3 |  0.01 |
| T5 dense   |  1.000000 |  1.0  | 2.66e-8  | 31 | 31 |  0.05 |
| T5 sparse  |  1.000000 |  1.0  | 1.11e-7  |  3 |  3 |  0.27 |
| T6 dense   |  1.04e-9  |  0.0  | 1.04e-9  | 41 | 41 |  0.05 |
| T6 sparse  | -5.96e-9  |  0.0  | 5.96e-9  |  5 |  5 |  0.01 |

Table 2 expectations (`paper = 0`, `mb = 61` dense / `mb = 15` sparse) are
tracked but no local Mosek solve has been recorded yet.

## Setup

When `MosekTools` is available, a silent Mosek optimizer is used by default.
Pass `optimizer=` to use a different backend or custom Mosek attributes.

````julia
using NCTSSoS, Printf, JSON3

const MOI = NCTSSoS.MOI
const EXAMPLE_6_1_DATA_DIR = joinpath(
    pkgdir(NCTSSoS), "docs", "src", "examples", "literate", "data"
)
const EXAMPLE_6_1_MOSEK_RESULTS_PATH = joinpath(
    EXAMPLE_6_1_DATA_DIR, "wang_magron_2021_example_6_1_mosek_n20_results.json"
)
const EXAMPLE_6_1_HAS_MOSEKTOOLS = Base.find_package("MosekTools") !== nothing

if EXAMPLE_6_1_HAS_MOSEKTOOLS
    import MosekTools
end
````

### Helpers

Internal utilities for optimizer resolution, expectation loading, result
verification, and summary formatting. These are used by the REPL entry
points defined further below.

````julia
function example_6_1_flatten_sizes(sizes)
    return reduce(vcat, sizes)
end
function example_6_1_max_block_size(result)
    return maximum(example_6_1_flatten_sizes(result.moment_matrix_sizes))
end
function example_6_1_default_optimizer()
    EXAMPLE_6_1_HAS_MOSEKTOOLS || error(
        "MosekTools is not available in the active environment. " *
        "Pass `optimizer=` explicitly to provide a compatible optimizer backend.",
    )
    return MOI.OptimizerWithAttributes(MosekTools.Optimizer, MOI.Silent() => true)
end
function example_6_1_resolve_optimizer(optimizer)
    return optimizer === nothing ? example_6_1_default_optimizer() : optimizer
end
function example_6_1_table_expectation_path(table::Int)
    2 <= table <= 6 || error("Only Tables 2--6 are fixture-backed in this example.")
    return joinpath(
        EXAMPLE_6_1_DATA_DIR, @sprintf("wang_magron_2021_example_6_1_table_%d.json", table)
    )
end

function example_6_1_load_table_expectation(table::Int)
    path = example_6_1_table_expectation_path(table)
    isfile(path) || error("Missing expectation JSON: $(repr(path))")
    return JSON3.read(read(path, String))
end

function example_6_1_load_mosek_n20_results()
    isfile(EXAMPLE_6_1_MOSEK_RESULTS_PATH) || error(
        "Missing tracked Mosek results JSON: $(repr(EXAMPLE_6_1_MOSEK_RESULTS_PATH))"
    )
    return JSON3.read(read(EXAMPLE_6_1_MOSEK_RESULTS_PATH, String))
end
function example_6_1_mode_key(mode::Symbol)
    mode in (:dense, :sparse) || error("Unsupported mode $(repr(mode)); use :dense or :sparse.")
    return String(mode)
end

function example_6_1_reference_row(table::Int, mode::Symbol, n::Int=20)
    data = example_6_1_load_table_expectation(table)["expected_values"]
    rows = data["$(example_6_1_mode_key(mode))_rows"]
    for row in rows
        Int(row["n"]) == n || continue
        return (
            table=table,
            mode=example_6_1_mode_key(mode),
            family=String(data["family"]),
            n=Int(row["n"]),
            paper=Float64(row["opt"]),
            paper_mb=Int(row["mb"]),
        )
    end
    error("No $(example_6_1_mode_key(mode)) expectation row for Table $(table) at n=$(n).")
end

function example_6_1_default_atol(table::Int, mode::Symbol=:sparse)
    example_6_1_mode_key(mode)
    return table == 6 ? 1e-6 : 1e-3
end
````

````
example_6_1_default_atol (generic function with 2 methods)
````

### Verification and Reporting

````julia
function example_6_1_verify_mosek_n20_results()
    return map(example_6_1_load_mosek_n20_results()["rows"]) do raw_row
        row = (
            table=Int(raw_row["table"]),
            family=String(raw_row["family"]),
            mode=String(raw_row["mode"]),
            n=Int(raw_row["n"]),
            objective=Float64(raw_row["objective"]),
            paper=Float64(raw_row["paper"]),
            abs_error=Float64(raw_row["abs_error"]),
            mb=Int(raw_row["mb"]),
            paper_mb=Int(raw_row["paper_mb"]),
            seconds=Float64(raw_row["seconds"]),
        )
        expected = example_6_1_reference_row(row.table, Symbol(row.mode), row.n)
        @assert row.family == expected.family (
            "Family mismatch for Table $(row.table) $(row.mode): got $(row.family), " *
            "expected $(expected.family)."
        )
        @assert row.paper == expected.paper (
            "Paper objective mismatch for Table $(row.table) $(row.mode): " *
            "got $(row.paper), expected $(expected.paper)."
        )
        @assert row.paper_mb == expected.paper_mb (
            "Paper block size mismatch for Table $(row.table) $(row.mode): " *
            "got $(row.paper_mb), expected $(expected.paper_mb)."
        )
        @assert row.mb == row.paper_mb (
            "Tracked Mosek block size mismatch for Table $(row.table) $(row.mode): " *
            "got $(row.mb), expected $(row.paper_mb)."
        )
        return row
    end
end

function example_6_1_summarize_rows(rows)
    io = IOBuffer()
    println(
        io,
        @sprintf(
            "%-7s %-8s %-7s %-4s %12s %10s %10s %5s %8s %8s",
            "table",
            "family",
            "mode",
            "n",
            "objective",
            "paper",
            "abs err",
            "mb",
            "nuniq",
            "sec",
        ),
    )
    println(io, repeat("-", 92))
    for row in rows
        println(
            io,
            @sprintf(
                "%-7s %-8s %-7s %-4d %12.6f %10.4f %10.3e %5d %8d %8.2f",
                "T$(row.table)",
                row.family,
                row.mode,
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
        mode=expected.mode,
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

function example_6_1_run_case(
    table::Int,
    mode::Symbol,
    pop,
    config;
    n::Int=20,
    atol::Real=example_6_1_default_atol(table, mode),
)
    expected = example_6_1_reference_row(table, mode, n)
    seconds = @elapsed result = cs_nctssos(pop, config)
    row = example_6_1_build_row(table, result, expected, seconds)
    example_6_1_verify_row(row, atol)
    return (row=row, result=result)
end
````

````
example_6_1_run_case (generic function with 1 method)
````

## Problem Builders

Each builder below encodes the polynomial optimisation problem for its
respective table, transcribed directly from the paper's formulations.

````julia
function build_table_2_problem(n::Int)
    n >= 2 || error("Table 2 requires n >= 2.")
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    P = typeof(x[1] + x[2])
    f = zero(P)

    for i in 1:n
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)
        coupling = isempty(jset) ? zero(P) : sum(x[j] + x[j]^2 for j in jset)
        g_i = 2.0 * x[i] + 5.0 * x[i]^3 + one(P) - coupling
        f += g_i * g_i
    end

    return polyopt(f, reg)
end

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

function example_6_1_problem_builder(table::Int)
    table == 2 && return build_table_2_problem
    table == 3 && return build_table_3_problem
    table == 4 && return build_table_4_problem
    table == 5 && return build_table_5_problem
    table == 6 && return build_table_6_problem
    error("Unsupported fixture-backed table $(table).")
end
````

````
example_6_1_problem_builder (generic function with 1 method)
````

## Solver Configuration

Each table uses a specific combination of relaxation order, basis, and
sparsity algorithms. The dense mode disables all elimination; the sparse
mode enables correlative sparsity (CS) and/or term sparsity (TS) via
`MF` and `MMD`.

````julia
function example_6_1_build_config(table::Int, mode::Symbol, pop; optimizer)
    if table == 2
        return SolverConfig(
            optimizer=optimizer,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=mode == :dense ? NoElimination() : MMD(),
        )
    elseif table == 3 || table == 4
        basis = newton_chip_basis(pop, 2)
        return SolverConfig(
            optimizer=optimizer,
            moment_basis=basis,
            cs_algo=mode == :dense ? NoElimination() : MF(),
            ts_algo=mode == :dense ? NoElimination() : MMD(),
        )
    elseif table == 5
        if mode == :dense
            return SolverConfig(
                optimizer=optimizer,
                moment_basis=newton_chip_basis(pop, 2),
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
            )
        end
        return SolverConfig(
            optimizer=optimizer,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=MMD(),
        )
    elseif table == 6
        return SolverConfig(
            optimizer=optimizer,
            moment_basis=newton_chip_basis(pop, 2),
            cs_algo=NoElimination(),
            ts_algo=mode == :dense ? NoElimination() : MMD(),
        )
    end
    error("Unsupported fixture-backed table $(table).")
end

function example_6_1_execute_case(
    table::Int,
    mode::Symbol;
    optimizer=nothing,
    n::Int=20,
    atol::Real=example_6_1_default_atol(table, mode),
    io::IO=stdout,
    print_summary::Bool=true,
)
    chosen_optimizer = example_6_1_resolve_optimizer(optimizer)
    pop = example_6_1_problem_builder(table)(n)
    report = example_6_1_run_case(
        table,
        mode,
        pop,
        example_6_1_build_config(table, mode, pop; optimizer=chosen_optimizer);
        n=n,
        atol=atol,
    )
    if print_summary
        print(io, example_6_1_summarize_rows([report.row]))
    end
    return report
end
````

````
example_6_1_execute_case (generic function with 1 method)
````

## REPL Entry Points

Each `test_table_<T>_<mode>_n20` function solves one case, verifies the
objective and block size against the JSON expectations, prints a summary
row, and returns `(row=..., result=...)`.

````julia
const EXAMPLE_6_1_N20_CASE_SPECS = [
    (table=2, mode=:dense),
    (table=2, mode=:sparse),
    (table=3, mode=:dense),
    (table=3, mode=:sparse),
    (table=4, mode=:dense),
    (table=4, mode=:sparse),
    (table=5, mode=:dense),
    (table=5, mode=:sparse),
    (table=6, mode=:dense),
    (table=6, mode=:sparse),
]

for spec in EXAMPLE_6_1_N20_CASE_SPECS
    fn = Symbol("test_table_", spec.table, "_", spec.mode, "_n20")
    table = spec.table
    mode = spec.mode
    @eval function $fn(;
        optimizer=nothing,
        atol::Real=example_6_1_default_atol($table, $(QuoteNode(mode))),
        io::IO=stdout,
    )
        return example_6_1_execute_case(
            $table,
            $(QuoteNode(mode));
            optimizer=optimizer,
            n=20,
            atol=atol,
            io=io,
        )
    end
end

function test_all_n20_reference_cases(; optimizer=nothing, io::IO=stdout)
    chosen_optimizer = example_6_1_resolve_optimizer(optimizer)
    reports = map(EXAMPLE_6_1_N20_CASE_SPECS) do spec
        example_6_1_execute_case(
            spec.table,
            spec.mode;
            optimizer=chosen_optimizer,
            n=20,
            atol=example_6_1_default_atol(spec.table, spec.mode),
            io=io,
            print_summary=false,
        )
    end
    print(io, example_6_1_summarize_rows([report.row for report in reports]))
    return reports
end

function test_all_fast_cases(; optimizer=nothing, io::IO=stdout)
    return test_all_n20_reference_cases(; optimizer=optimizer, io=io)
end
````

````
test_all_fast_cases (generic function with 1 method)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


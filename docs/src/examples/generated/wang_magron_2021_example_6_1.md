```@meta
EditURL = "../literate/wang_magron_2021_example_6_1.jl"
```

# [Wang-Magron 2021 Example 6.1: Size-Varying Eigenvalue Benchmarks](@id wang-magron-2021-example-6-1)

This page reproduces a small, docs-friendly slice of Section 6.1 from
Wang and Magron's NCTSSOS paper [wangExploitingTermSparsity2021](@cite).
It is intentionally separate from the tracial-polynomial "Example 6.1"
page: here we solve eigenvalue benchmarks for ordinary noncommutative
polynomials, not trace-polynomial problems.

Section 6.1 spans Tables 2--8:

- Table 2 (`f_Bb`, unconstrained Broyden banded): sparse optimum `0`.
- Table 3 (`f_cs`, unconstrained chained singular): sparse optimum near `0`.
- Table 4 (`f_gR`, unconstrained generalized Rosenbrock): sparse optimum `1`.
- Table 5 (`f_cW`, unconstrained chained Wood): sparse optimum near `1`.
- Table 6 (`f_Bt`, unconstrained Broyden tridiagonal): sparse optimum `0`.
- Table 7 (`f_Bb` over `D`): positive constrained optimum that grows with `n`.
- Table 8 (random quartics over multi-balls): large-scale random instances.

The docs page executes two representative families only:

1. Table 4 (`f_gR`) to show why some unconstrained benchmarks have known
   minimum `1`, and how correlative sparsity versus term sparsity changes cost.
2. Table 7 (`f_Bb` over `D`) to show how constraints move the optimum away
   from `0` and how runtime grows with problem size.

Tables 2, 3, 5, 6, and 8 are cited for context but omitted from execution to
keep the docs path short and deterministic.

## Setup

Mosek matches the workflow already used by the other executed example pages.

````julia
using NCTSSoS, MosekTools, Printf

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)

const TABLE_4_REFERENCE = Dict(20 => 1.0)
const TABLE_7_REFERENCE = Dict(5 => 3.113, 10 => 3.011, 20 => 9.658)

flatten_sizes(sizes) = reduce(vcat, sizes)
block_count(result) = length(flatten_sizes(result.moment_matrix_sizes))

function summarize_rows(rows)
    io = IOBuffer()
    println(
        io,
        @sprintf(
            "%-8s %-4s %-6s %12s %10s %10s %8s %8s %8s",
            "family", "n", "mode", "objective", "paper", "abs err", "nuniq", "blocks", "sec"
        ),
    )
    println(io, repeat("-", 84))
    for row in rows
        println(
            io,
            @sprintf(
                "%-8s %-4d %-6s %12.6f %10.3f %10.3e %8d %8d %8.2f",
                row.family,
                row.n,
                row.mode,
                row.objective,
                row.paper,
                row.abs_error,
                row.n_unique,
                row.blocks,
                row.seconds,
            ),
        )
    end
    return String(take!(io))
end

function run_case(family::AbstractString, n::Int, mode::AbstractString, pop, config, paper_ref::Real)
    seconds = @elapsed result = cs_nctssos(pop, config)
    return (
        family=family,
        n=n,
        mode=mode,
        objective=result.objective,
        paper=Float64(paper_ref),
        abs_error=abs(result.objective - paper_ref),
        n_unique=result.n_unique_moment_matrix_elements,
        blocks=block_count(result),
        seconds=round(seconds; digits=2),
    )
end

function generalized_rosenbrock_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    f = Float64(n) * one(typeof(x[1]))
    for i in 2:n
        f += 100.0 * x[i - 1]^4 - 200.0 * x[i - 1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
    end
    return polyopt(f, reg)
end

function constrained_broyden_banded_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    f = sum(1:n) do i
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)
        g = isempty(jset) ? 0.0 * x[1] : sum(x[j] + x[j]^2 for j in jset)
        (2.0 * x[i] + 5.0 * x[i]^3 + 1 - g)^2
    end
    ineq_constraints = [[1.0 - x[i]^2 for i in 1:n]; [x[i] - 1//3 for i in 1:n]]
    return polyopt(f, reg; ineq_constraints=ineq_constraints)
end
````

## Why the unconstrained minima are 0 or 1

The unconstrained Section 6.1 families are built from Hermitian-square
patterns [wangExploitingTermSparsity2021](@cite):

- Tables 2, 3, and 6 are minimized at the origin, so the benchmark minimum is `0`.
- Tables 4 and 5 add constant offsets that shift the minimum to `1`.

The constrained family in Table 7 breaks that trivial certificate by adding
`D = {1 - X_i^2, X_i - 1/3}`. The feasible set is pushed away from the origin,
so the optimum becomes positive and size-dependent.

## Unconstrained benchmark: Table 4 (`f_gR`)

Table 4 reports sparse optimum `1.0000` at `n = 20` for the generalized
Rosenbrock family. We solve the same `n = 20` instance twice:

- `CS`: correlative sparsity only
- `CS+TS`: correlative sparsity plus term sparsity

Dense relaxation is intentionally omitted here because the sparse comparison
already captures the practical point of Section 6.1 while keeping the docs
runtime modest.

````julia
gr20 = generalized_rosenbrock_problem(20)

table4_rows = [
    run_case(
        "f_gR",
        20,
        "CS",
        gr20,
        SolverConfig(
            optimizer=SILENT_MOSEK,
            order=2,
            cs_algo=MF(),
            ts_algo=NoElimination(),
        ),
        TABLE_4_REFERENCE[20],
    ),
    run_case(
        "f_gR",
        20,
        "CS+TS",
        gr20,
        SolverConfig(
            optimizer=SILENT_MOSEK,
            order=2,
            cs_algo=MF(),
            ts_algo=MMD(),
        ),
        TABLE_4_REFERENCE[20],
    ),
]

println(summarize_rows(table4_rows))
@assert all(row.abs_error < 1e-3 for row in table4_rows)
````

````
family   n    mode      objective      paper    abs err    nuniq   blocks      sec
------------------------------------------------------------------------------------
f_gR     20   CS         1.000000      1.000  7.553e-08      328       19     5.77
f_gR     20   CS+TS      1.000000      1.000  4.397e-07      117       94     0.14


````

Table 4's optimum is reproduced within `1e-3`, while the printed `nuniq` and
`blocks` columns show what term sparsity changes in practice: the answer stays
at `1`, but the SDP description becomes much smaller.

## Constrained benchmark: Table 7 (`f_Bb` over `D`)

Table 7 uses the noncommutative Broyden banded family with
`D = {1 - X_i^2, X_i - 1/3}`. The paper reports positive optima that grow with
`n`, which makes this a good "size-varying" benchmark even at modest sizes.

We keep the solver mode fixed at `CS+TS` and reproduce the table values for
`n = 5, 10, 20`.

````julia
table7_rows = [
    run_case(
        "f_Bb^D",
        n,
        "CS+TS",
        constrained_broyden_banded_problem(n),
        SolverConfig(
            optimizer=SILENT_MOSEK,
            order=3,
            cs_algo=MF(),
            ts_algo=MMD(),
        ),
        TABLE_7_REFERENCE[n],
    ) for n in (5, 10, 20)
]

println(summarize_rows(table7_rows))
@assert all(row.abs_error < 2e-3 for row in table7_rows)
````

````
family   n    mode      objective      paper    abs err    nuniq   blocks      sec
------------------------------------------------------------------------------------
f_Bb^D   5    CS+TS      3.113005      3.113  4.992e-06      233      145     0.18
f_Bb^D   10   CS+TS      3.011288      3.011  2.883e-04     1087     1544     1.58
f_Bb^D   20   CS+TS      9.658650      9.658  6.497e-04     2877     5404     6.66


````

The constrained objectives line up with Table 7, and the `sec` column makes
the scaling story concrete: even with `CS+TS`, the work rises quickly between
`n = 5` and `n = 20`. That is exactly why the paper pushes the same ideas to
much larger sparse instances instead of dense SDPs.

## Meaning

Three takeaways matter more than the raw numbers:

1. In the unconstrained families, the published optima are explained by the
   algebraic shape of the objective itself. Tables 2, 3, and 6 are anchored at
   `0`; Tables 4 and 5 are shifted to `1`.
2. Correlative sparsity and term sparsity mostly change cost, not the target
   optimum. On `f_gR`, both sparse runs stay at the Table 4 value, but
   `CS+TS` shrinks the SDP description substantially.
3. Constraints change the story qualitatively. Table 7 no longer has a trivial
   `0` or `1` minimum, and the positive optimum grows with `n`. The docs sweep
   already shows the trend; the full paper extends it to much larger sizes.

Table 8 is left out here because those random quartics come from external
instances rather than a compact closed-form family.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


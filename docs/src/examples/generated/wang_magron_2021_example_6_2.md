```@meta
EditURL = "../literate/wang_magron_2021_example_6_2.jl"
```

# Wang-Magron 2021 Example 6.2: Size-Varying Trace Benchmarks

This page recreates a small, docs-friendly slice of Section 6.2 from
Wang--Magron [wangExploitingTermSparsity2021](@cite): trace minimization
benchmarks whose size increases while the local interaction pattern stays
fixed.

The paper-scale takeaway is not that the optimum changes dramatically for the
unconstrained Broyden families. It is that sparse trace relaxations keep the
SDP blocks small enough to remain usable, while the dense trace relaxations
run out of memory almost immediately (Tables 9-10).

We execute three representative runs:

1. Unconstrained **Broyden banded** at `n = 10` (Table 9 reference).
2. Unconstrained **Broyden tridiagonal** at `n = 20` (Table 10 reference).
3. Constrained **Broyden banded over `D`** at `n = 5` (Table 11 reference).

Table 12 is included as a paper reference point only: those random multi-ball
instances are not shipped with this repository, so the deterministic Broyden
families are the stable in-tree reproductions.

## Setup

We use Mosek, solve the moment relaxation directly (`dualize=false`), and keep
the executed size sweep small enough for `make examples`.

````julia
using NCTSSoS, JuMP, MosekTools

const SILENT_MOSEK = optimizer_with_attributes(
    Mosek.Optimizer,
    "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
    "MSK_IPAR_LOG" => 0,
)

flatten_sizes(sizes) = reduce(vcat, sizes; init=Int[])

function summarize_result(result, elapsed_seconds)
    flat_sizes = flatten_sizes(result.moment_matrix_sizes)
    return (
        objective=Float64(result.objective),
        max_block=maximum(flat_sizes),
        nblocks=length(flat_sizes),
        n_unique_moments=result.n_unique_moment_matrix_elements,
        elapsed_seconds=round(elapsed_seconds; digits=2),
    )
end

function compare_to_paper(result_summary, paper_ref)
    return (
        reproduced_objective=result_summary.objective,
        paper_objective=paper_ref.objective,
        objective_error=abs(result_summary.objective - paper_ref.objective),
        reproduced_max_block=result_summary.max_block,
        paper_sparse_mb=paper_ref.paper_sparse_mb,
        nblocks=result_summary.nblocks,
        n_unique_moments=result_summary.n_unique_moments,
        elapsed_seconds=result_summary.elapsed_seconds,
    )
end

trace_lift(poly, exemplar) = tr(poly) * one(typeof(exemplar))

function broyden_banded_objective(x)
    poly_type = typeof(x[1] + x[min(2, length(x))])
    f = Float64(length(x)) * one(poly_type)

    for i in eachindex(x)
        jset = setdiff(max(1, i - 5):min(length(x), i + 1), i)
        f += 4.0 * x[i] + 4.0 * x[i]^2 + 10.0 * x[i]^3 + 20.0 * x[i]^4 + 25.0 * x[i]^6

        for j in jset
            f += -2.0 * x[j] - 2.0 * x[j]^2 - 4.0 * x[i] * x[j] -
                 4.0 * x[i] * x[j]^2 - 10.0 * x[i]^3 * x[j] - 10.0 * x[i]^3 * x[j]^2
        end

        for j in jset, k in jset
            f += x[j] * x[k] + 2.0 * x[j]^2 * x[k] + x[j]^2 * x[k]^2
        end
    end

    return f
end

function broyden_tridiagonal_objective(x)
    n = length(x)
    poly_type = typeof(x[1] + x[min(2, n)])
    f = Float64(n) * one(poly_type)

    f += 5.0 * x[1]^2 + 4.0 * x[1]^4 + 4.0 * x[2]^2 - 12.0 * x[1]^3 -
         12.0 * x[1] * x[2] + 6.0 * x[1] + 8.0 * x[1]^2 * x[2] - 4.0 * x[2]

    for i in 2:(n - 1)
        f += 5.0 * x[i]^2 + 4.0 * x[i]^4 + x[i - 1]^2 + 4.0 * x[i + 1]^2 -
             12.0 * x[i]^3 - 6.0 * x[i - 1] * x[i] - 12.0 * x[i] * x[i + 1] +
             6.0 * x[i] + 4.0 * x[i - 1] * x[i]^2 + 8.0 * x[i]^2 * x[i + 1] +
             4.0 * x[i - 1] * x[i + 1] - 2.0 * x[i - 1] - 4.0 * x[i + 1]
    end

    f += 5.0 * x[n]^2 + 4.0 * x[n]^4 + x[n - 1]^2 - 12.0 * x[n]^3 -
         6.0 * x[n - 1] * x[n] + 6.0 * x[n] + 4.0 * x[n - 1] * x[n]^2 - 2.0 * x[n - 1]

    return f
end

function run_trace_case(objective_builder, n, order; constrained=false, cs_algo=NoElimination(), ts_algo=MMD())
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    objective = trace_lift(objective_builder(x), x[1])

    ineq_constraints = if constrained
        vcat(
            [trace_lift(1.0 - x[i]^2, x[1]) for i in 1:n],
            [trace_lift(x[i] - 1.0 / 3, x[1]) for i in 1:n],
        )
    else
        typeof(objective)[]
    end

    pop = polyopt(objective, reg; ineq_constraints=ineq_constraints)
    solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order, cs_algo, ts_algo)

    elapsed_seconds = @elapsed result = cs_nctssos(pop, solver_config; dualize=false)
    return summarize_result(result, elapsed_seconds)
end

const TABLE_9_REF = (n=10, paper_sparse_mb=29, objective=0.0)
const TABLE_10_REF = (n=20, paper_sparse_mb=6, objective=0.0)
const TABLE_11_REF = (n=5, paper_sparse_mb=19, objective=3.113)
const TABLE_12_REF = (n=505, k=1, paper_sparse_mb=16, objective=-4.997)
````

## Paper Reference Points (Tables 9-12)

These are the paper values we use as compact reference markers:

| Paper table | Family | Selected row | Paper sparse `mb` | Paper objective |
|:------------|:-------|:-------------|------------------:|----------------:|
| Table 9 | Broyden banded | `n = 10` | 29 | 0 |
| Table 10 | Broyden tridiagonal | `n = 20` | 6 | 0 |
| Table 11 | Broyden banded over `D` | `n = 5` | 19 | 3.113 |
| Table 12 | Random multi-balls | `n = 505`, `k = 1` | 16 | -4.997 |

The exact block inventory produced by the current Julia implementation can
differ from the paper's `mb` column because elimination heuristics and clique
decompositions are implementation-dependent. The objective value is the main
reproducibility target here.

## Unconstrained Broyden Banded

Table 9 reports sparse trace minima of `0` for the banded family. We reproduce
the smallest paper row, `n = 10`, at relaxation order `3` and use term
sparsity only.

````julia
banded_n10 = compare_to_paper(
    run_trace_case(broyden_banded_objective, TABLE_9_REF.n, 3; ts_algo=MMD()),
    TABLE_9_REF,
)

@show banded_n10
@assert banded_n10.objective_error < 1e-4
````

````
banded_n10 = (reproduced_objective = 7.997448019381181e-7, paper_objective = 0.0, objective_error = 7.997448019381181e-7, reproduced_max_block = 15, paper_sparse_mb = 29, nblocks = 4404, n_unique_moments = 3708, elapsed_seconds = 22.82)

````

## Unconstrained Broyden Tridiagonal

Table 10 shows the same zero optimum for the tridiagonal family, but with much
smaller sparse blocks than the banded case. We again use term sparsity only.

````julia
tridiagonal_n20 = compare_to_paper(
    run_trace_case(broyden_tridiagonal_objective, TABLE_10_REF.n, 2; ts_algo=MMD()),
    TABLE_10_REF,
)

@show tridiagonal_n20
@assert tridiagonal_n20.objective_error < 1e-4
````

````
tridiagonal_n20 = (reproduced_objective = -1.3814031269987481e-6, paper_objective = 0.0, objective_error = 1.3814031269987481e-6, reproduced_max_block = 5, paper_sparse_mb = 6, nblocks = 1240, n_unique_moments = 1186, elapsed_seconds = 1.09)

````

## Constrained Broyden Banded over D

Section 6.2 also studies the same banded family over

```math
D = \{1 - X_i^2,\; X_i - \tfrac{1}{3} : i = 1,\ldots,n\}.
```

In the current tracial interface we lift these defining polynomials with
[`tr`](@ref), then solve the sparse relaxation with combined correlative and
term sparsity.

````julia
constrained_n5 = compare_to_paper(
    run_trace_case(
        broyden_banded_objective,
        TABLE_11_REF.n,
        3;
        constrained=true,
        cs_algo=MF(),
        ts_algo=MMD(),
    ),
    TABLE_11_REF,
)

@show constrained_n5
@assert constrained_n5.objective_error < 1e-2
````

````
constrained_n5 = (reproduced_objective = 3.1072928706826413, paper_objective = 3.113, objective_error = 0.005707129317358728, reproduced_max_block = 11, paper_sparse_mb = 19, nblocks = 626, n_unique_moments = 596, elapsed_seconds = 1.09)

````

## Meaning

The change from Section 6.1's eigenvalue benchmarks is mostly in the
**relaxation geometry**, not in the unconstrained optimum. For both Broyden
families the sparse trace relaxation still lands at `0`, so the practical
question becomes: how fast do the SDP blocks grow as `n` increases?

Tables 9-10 answer that directly. Dense trace runs are already omitted at the
smallest reported sizes, while the sparse trace relaxations keep small blocks
and continue scaling. The tridiagonal family is structurally easier than the
banded family, which is why its sparse `mb` column stays much smaller.

The constrained problem over `D` changes the story. Once the feasible set is
tightened, the optimum becomes strictly positive (Table 11), but the sparse
relaxation is still practical at sizes where dense relaxations are not. Table
12 shows the same scalability pattern on random multi-ball instances; we cite
it as context even though we only execute deterministic in-repo examples here.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


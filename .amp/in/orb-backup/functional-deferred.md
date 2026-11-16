# Deferred functional block work

Deferred by the owner until after the merge into main; do not merge.

Base: `579568d` on `redesign`. The source changes are an unfinished, unmerged draft. The three front-door guards remain unchanged on `redesign`. No finite-action routing or factor-merge rewrite was started.

## Completed observations

- State and trace CHSH: plain and dense-block routes both have one real PSD block of side 9, 20 free variables, and bound −2√2 within 1e-6.
- The old API reproduces the reviewed CHSH dense and sparse bounds. CS+MMD level 1 gives −4 on both old and new routes, not the dense −2√2.
- Legacy `src/optimization/interface.jl:500–506` replaces complete mixed rows with embedded operator-only rows for sparse pure-trace programs (five instead of nine CHSH rows).
- Legacy `src/optimization/sparsity.jl:1000–1013` omits mixed-row diagonal seeds for state-polynomial splitting. Next seeds those diagonals.
- Replay with the old row order does not remove the difference. Restoring legacy seeds reproduces the state sides `[1,1,1,1,1,3,3]`; keeping current seeds gives `[1,2,2,2,2,3,3]`. Selecting the old five trace rows reproduces `[1,3,3]`. All replay bounds agree with −2√2 within 1e-6.
- The prior sparse-route test revision passed 114/114 in 129 s. It checks both domains, all three rules, levels 1–3, complete starting-row equivalence and symbol localizers. The later ensemble/constraint-graph additions and the state-CHSH two-site adjustment have not been run. The ordinary functional checks passed 38/38 before those additions.

## Preserved scripts and data

All files below are in this directory. Run from the repository root with Julia 1.12.7 and the dependencies in the local solver environment (`Mosek`, `MosekTools`, `JuMP`, plus package dependencies). That environment and the license are not committed.

- `functional-sparse-pilot.jl` / `.log`: plain, block, CS, TS and CSTS CHSH measurements. The saved log ends at an earlier test-constructor error; the constructor was corrected later.
- `functional-sparse-legacy.jl` / `.log`: corresponding old-API solves, completed in 51 s.
- `functional-sparse-models.jl`: loads the unchanged bilocal and T2 reference functions from the source tree and constructs the same expressions through Next public variables. Retain access to the old source at the base commit when resuming after removal.
- `functional-sparse-replay.jl` / `.log`: old-order/old-seed controls, completed in 63 s. Also confirms the bilocal model has degree 6 and 44 terms.
- `functional-block-tests.log`: the earlier 114/114 sparse-route result.
- `functional-build-pilot.jl` / `.log`: cold compilation-inclusive sizing run, explicitly terminated when the owner deferred this work (exit 143, 61 s). It is incomplete and is not a benchmark result. No warmed distributions, product-share profile, bilocal solve comparison or T2 comparison was completed.

Remaining acceptance includes the final tests, independent ensembles, bilocal dense/MMD solves, covariance references, warmed build benchmarks, product-share profiles, the 242 nonfunctional fingerprints, quality checks, and review. Finite symmetry remains a separate deferred PR. Do not infer completion from the preserved pilot results.

# Current Task Context

## Mission

Design a **particle-number conservation public API** for `polyopt`,
built on top of the existing moment-equality-constraint mechanism, so users
can declare *"fix the total particle number to N"* for Fermionic or Bosonic
problems without hand-writing the number operator and wiring it into
`moment_eq_constraints` themselves.

## Why this worktree retracked off `origin/main`

A dedicated particle-conservation API does **not** exist on either:

- `origin/main` of `QuantumSOS/NCTSSoS.jl`, nor
- `origin/feat/h4-periodic-v2rdm-benchmark` (sibling worktree at
  `/Users/exaclior/QuantumSOS/NCTSSoS.jl-h4-periodic-v2rdm-benchmark`).

What **does** exist on the h4 branch — and is the natural starting point —
are worked examples that *use* `moment_eq_constraints` to impose a
number-operator equality by hand:

- `docs/src/examples/literate/bosonic_ground_state.jl` — fixes
  $\hat{N} = 2$ via a canonical particle-number constraint.
- `docs/src/examples/literate/hubbard_model.jl` — canonical sector with
  $N_\uparrow, N_\downarrow$ fixed.
- `docs/src/examples/literate/fermionic_ground_state.jl` — even-parity /
  particle-number sectors.
- `demos/h4_periodic_*.jl` and the V2RDM PQG demos — periodic H₄ / H₂
  at fixed filling.
- `src/optimization/moment.jl` (around the `moment_eq_constraints`
  docstring) — already names "particle-number fixing in fermionic systems"
  as the motivating use case.

For that reason, this worktree now tracks
`origin/feat/h4-periodic-v2rdm-benchmark` (HEAD `a07189a`) instead of
`origin/main`, so the API is designed *alongside* those examples rather
than ahead of them.

## End goal

A well-documented public API that lets a user say, in one line, *"constrain
this problem to total particle number $N$ (and optionally per-spin or
per-mode $N_i$)"*, for both Fermionic and Bosonic problems. Concretely:

- A user-facing constructor (working name TBD, e.g.
  `particle_number_constraint(...)`) that emits the right
  moment-equality constraint(s) and plugs into
  `polyopt(...; moment_eq_constraints = ...)`.
- Symmetric, principled support for:
  - `FermionicAlgebra` — total $N$, spin-resolved
    ($N_\uparrow$, $N_\downarrow$), mode-resolved.
  - `BosonicAlgebra` — total $N$, mode-resolved.
- Documented public exports, doctests, and at least one Literate example
  refactored to use the new API in place of the manual number-operator
  construction.
- Tests pinning behaviour against the existing canonical-sector benchmarks
  (Hubbard half-filling, bosonic ground state at $\hat{N} = 2$, fermionic
  even-parity sector) so the new API reproduces the values the manual
  `moment_eq_constraints` form already produces.

## Locked design decisions

| # | Decision | Choice |
|---|---|---|
| 1 | Public surface | Helper function only: `particle_number_constraint` |
| 2a | v1 scope | Total + per-group ship together |
| 2b | Group spelling | Pass monomial vectors (zero registry change; arbitrary subsets free) |
| 3 | Algebra dispatch | Infer from registry |
| 4 | CS interaction | Document only |
| 5 | Bosonic truncation | Out of scope — doc cross-reference to `ineq_constraints` |
| 6 | State/trace | Ordinary polyopt only; clean "not yet supported" error |

### Concrete public signature

```julia
# Total particle number across all modes in the registry
particle_number_constraint(registry::VariableRegistry{A,T}, N::Integer)
    where {A <: Union{FermionicAlgebra, BosonicAlgebra}, T}
    → Vector{Polynomial{A,T,Float64}}   # length 1: [N̂ - N·I]

# Per-group particle numbers (spin species, mode subsets, anything)
particle_number_constraint(registry, group => N, more_groups => more_Ns...)
    → Vector{Polynomial{A,T,Float64}}   # length = number of groups
```

User then splats into `polyopt`:

```julia
cons = particle_number_constraint(registry, c_up => 2, c_dn => 2)
pop  = polyopt(ham, registry; moment_eq_constraints = cons)
```

### Implementation plan, in order

1. **Survey** existing call sites — `docs/src/examples/literate/{bosonic_ground_state,hubbard_model,fermionic_ground_state}.jl`, `demos/h4_periodic_*.jl`, V2RDM PQG demos — to confirm every shape currently hand-rolled fits this signature. *Failure mode to look for: any example builds N̂ over a mode subset that isn't a slice of an existing monomial vector.*
2. **Extract registry-wide mode list** — read `idx_to_vars` for positive indices (annihilation operators). Verify against existing examples.
3. **Implement `particle_number_constraint`** in `src/optimization/`. Single file. ~50 lines.
4. **Public export** + **docstring** with cross-reference to `ineq_constraints` for per-mode bosonic caps.
5. **Tests** that pin against existing canonical-sector benchmarks: bosonic $\hat N = 2$ ground state, Hubbard $(N_\uparrow, N_\downarrow) = (2,2)$. The helper output must produce numerically identical SDP solutions to the manually-constructed `moment_eq_constraints`.
6. **Refactor one Literate example** (Hubbard is the loudest win) to use the new helper. Regenerate via `make examples`.
7. **Open issue** for Decision 6 follow-up: lift the state/trace `moment_eq_constraints` block at `src/optimization/problem.jl:227`.

## Current status

- Worktree retracked from `origin/main` →
  `origin/feat/h4-periodic-v2rdm-benchmark` (HEAD `a07189a`).
- Working branch: `feat/particle-number-conservation-api` (no upstream set).
- Design locked. Ready for implementation step 1 (call-site survey).

## Inherited context (still relevant from h4 branch)

- All heavy testing runs on the server via `easy-ssh`; do not run tests
  locally unless explicitly asked.
- Server / local path map and `uv` Python policy as in the previous
  TASK.md — preserved here only by reference. If you need them, recover
  the previous `TASK.md` from
  `git show origin/feat/h4-periodic-v2rdm-benchmark:TASK.md`.

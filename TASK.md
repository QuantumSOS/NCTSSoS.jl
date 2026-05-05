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

## Open design questions (resolve before adding any public symbols)

1. **Public surface** — standalone helper that returns a constraint, or a
   new keyword to `polyopt` (`particle_number = ...`)? Or both, with the
   keyword as sugar over the helper.
2. **Spin-resolved Fermionic spelling** — pass `(N_up, N_dn)`, a
   `Dict{Symbol,Int}` keyed by spin label, or a vector keyed by mode?
3. **Algebra dispatch** — infer from the `polyopt` call's variables, or
   require an explicit `algebra = ...` argument?
4. **Interaction with sparsity / CS** — number-operator equalities are
   dense across modes; document expected impact on clique structure and on
   the term-sparsity step.
5. **Bosonic truncation** — explicit cutoff vs. lazy vs. required user
   input.

Land these in a short design note (e.g. `docs/dev/particle_number_api.md`
or a section in this file) before introducing public symbols.

## Current status

- Worktree retracked from `origin/main` →
  `origin/feat/h4-periodic-v2rdm-benchmark` (HEAD `a07189a`).
- No new code yet. First step is the design note + a survey of existing
  manual `moment_eq_constraints` call sites in `docs/src/examples/literate/`
  and `demos/`, to make sure the API covers every shape currently
  hand-rolled.

## Inherited context (still relevant from h4 branch)

- All heavy testing runs on the server via `easy-ssh`; do not run tests
  locally unless explicitly asked.
- Server / local path map and `uv` Python policy as in the previous
  TASK.md — preserved here only by reference. If you need them, recover
  the previous `TASK.md` from
  `git show origin/feat/h4-periodic-v2rdm-benchmark:TASK.md`.

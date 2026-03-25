# nctssos-intermediate-sparsity-benchmarks

## Target Type
feature

## Target
Sparsity exploitation, mixed algebras, and benchmarking

## Use Case
Can a practitioner exploit correlative and term sparsity for performance, combine algebra types via `ComposedMonomial`, use `cs_nctssos_higher` for iterative refinement, and benchmark bounds/timing against TSSOS or ncpol2sdpa reference values? Covers the workflows in `sparsity_convergence.jl`, `bell.jl`, `mixed_algebras_tensor_products.jl`, and `certify_ground_state_property.jl`.

## Expected Outcome
- Toggles `cs_algo` between `MF()`, `MMD()`, `MaximalElimination()`, and `NoElimination()` and verifies that clique decomposition changes but final bounds stay consistent (within solver tolerance).
- Uses `ts_algo=MMD()` to exploit term sparsity; calls `cs_nctssos_higher(pop, prev_result, config)` and confirms monotonic bound improvement (as shown in `ground_state_energy.jl`).
- Combines Pauli + Fermionic variables via `ComposedMonomial((pauli_mono, fermionic_mono))` and runs `polyopt` on the composed problem without type errors.
- Supplies a custom `moment_basis` to `SolverConfig` and verifies the solver respects it (does not silently drop moments).
- Compares NCTSSoS CHSH bound (~2√2 ≈ 2.8284) against ncpol2sdpa output and finds agreement to solver tolerance.
- Recognizes the documented gap: "graph stabilization ≠ bound exactness" from `sparsity_convergence.jl` Examples 3.4 vs 3.8.

## Agent

### Background
Postdoc in quantum information theory (Rohan). Fluent in Julia 1.10+ and JuMP. Has used TSSOS (Python) and ncpol2sdpa for NPA-hierarchy computations over the past 3 years. Migrating existing Bell-inequality and ground-state workflows to NCTSSoS. Expects feature parity or explicit documentation of gaps. Has read the Burgdorf-Klep-Povh textbook.

### Experience Level
Intermediate

### Decision Tendencies
- Reads function signatures (`?cs_nctssos`, `?SolverConfig`) before looking at examples.
- First thing after install: checks `Project.toml` for solver deps and whether Mosek or COSMO is default.
- Will try all four elimination algorithms on the same problem and `@time` each one.
- When `cs_nctssos_higher` doesn't improve the bound, investigates whether term sparsity graph has stabilized rather than blindly increasing `order`.
- Expects `polyopt(StatePolynomial, reg)` to work for PBW algebras; will hit the `MonoidAlgebra`-only restriction and want a clear error message.
- Tries `NCTSSoS.get_term_sparsity_graph()` and `chordal_completion_edges()` from `sparsity_convergence.jl` — these are internal but documented in examples, creating API-boundary confusion.

### Quirks
- Wraps every solve in `@time` and `@allocated`; files performance regressions as issues.
- Expects TSSOS-style `nb=2, CS=true, TS="block"` keyword arguments; trips over NCTSSoS's `SolverConfig` struct pattern.
- Immediately substitutes `Mosek.Optimizer` for `COSMO.Optimizer` and adjusts `eps_abs` — frustrated if solver swap requires code changes beyond the optimizer field.
- When `certify_ground_state_property.jl` defines a local `cs_nctssos_with_entry()` helper (not public API), assumes it should be exported and opens an issue.
- Compares clique counts between NCTSSoS and TSSOS on the same problem; expects identical decompositions (they may differ due to elimination algorithm differences).

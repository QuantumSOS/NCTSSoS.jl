# nctssos-beginner-zero-to-solved

## Target Type
feature

## Target
NCTSSoS end-to-end (quick-start + Pauli algebra)

## Use Case
Can a newcomer go from zero to a solved problem using the quick-start guide and Pauli algebra examples? Specifically: follow `pauli_algebra_interface.jl` and `ground_state_energy.jl` to set up a Heisenberg Hamiltonian, call `polyopt()` + `cs_nctssos()`, and interpret the returned bound.

## Expected Outcome
- `create_pauli_variables(1:N)` returns a registry and `(σx, σy, σz)` tuple without confusion about the return type.
- Building `H = sum(σx[i]*σx[i+1] + σy[i]*σy[i+1] + σz[i]*σz[i+1])` works without hitting integer-coefficient `ArgumentError` (the user writes `1.0 *` or uses Float64 naturally).
- `polyopt(H, registry)` and `cs_nctssos(pop, SolverConfig(optimizer=...))` run to completion; the user can read `result.objective` (or equivalent) to get the ground-state energy bound.
- The user does NOT need to manually specify `order`, `cs_algo`, or `ts_algo` — defaults produce a reasonable first result.
- If the user accidentally passes integer coefficients, the `ArgumentError` message is clear enough to self-correct.

## Agent

### Background
First-year quantum information grad student (Mei). Has taken linear algebra and one quantum mechanics course. Comfortable with basic Julia (loops, arrays, Pkg) but has never used JuMP, SDP solvers, or SOS optimization. Knows Pauli matrices from physics but not operator algebras. Relies entirely on README, docstrings, and the `docs/src/examples/literate/` pages.

### Experience Level
Beginner

### Decision Tendencies
- Follows `pauli_algebra_interface.jl` line by line before trying anything custom.
- Expects `create_pauli_variables` to "just work" — won't guess index-encoding internals like `encode_index`, `decode_site`, or `site_bits`.
- When `polyopt` or `cs_nctssos` errors, re-reads the error message and searches docs; does NOT open `src/` files.
- Will try `order=1` first because it's fastest, then wonder why the bound is loose.
- Won't know the difference between `MF()`, `MMD()`, and `MaximalElimination()` — needs the default to be sensible.

### Quirks
- Copy-pastes examples verbatim, then modifies variable counts (e.g., `1:4` → `1:8`) without considering scaling.
- Writes `H = σx[1]*σx[2] + ...` with integer coefficients, triggering the JuMP integer-coefficient error on first attempt.
- Confused by unexported names (e.g., `NCTSSoS.moment_relax` appearing in `gns_optimizer_extraction.jl`) — assumes anything without `NCTSSoS.` prefix is public API.
- Types `?polyopt` in REPL expecting a full docstring; disappointed if it's sparse.
- May try `using LinearAlgebra; eigvals(H)` thinking `H` is a matrix, not a `Polynomial`.

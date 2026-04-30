# Repository Guidelines

## Current Task Context
- Read `TASK.md` before work for current task-specific constraints and environment notes.

## Project Layout
- `src/` — library code
  - `src/types/` — algebras, registries, monomials, polynomials
  - `src/simplification/` — algebra-specific rewrite rules (one file per algebra)
  - `src/optimization/` — sparsity + moment/SOS relaxations + JuMP model build
  - `src/states/` — state polynomials (quantum information)
- `test/` — curated suites (entry: `test/runtests.jl`)
- `test/data/` — reviewed expectation fixtures for solver-backed tests
- `docs/` — Documenter site (`docs/make.jl`, sources in `docs/src/`)
  - Literate sources: `docs/src/examples/literate/`
  - Generated (committed): `docs/src/examples/generated/`
- `test/oracles/` — optional reference outputs (requires external legacy `NCTSSOS` repo)

## Build, Test, and Docs Commands
CI baseline: Julia 1.11; solver: COSMO.
- `make init` — precompile root environment
- `make test` — full test suite (COSMO)
- `make coverage-ci` — CI-style coverage (`lcov.info`)
- `julia --project -e 'using Pkg; Pkg.test()'` — direct test run
- `make init-docs` — set up docs environment (`docs/Project.toml`)
- `make servedocs` — live docs server (LiveServer)
- `make examples` — execute Literate examples, regenerate `docs/src/examples/generated/`, and update `docs/examples_stamp.toml` (requires a Mosek license; CI only verifies the stamp via `julia docs/examples_stamp.jl`)
- `NCTSSOS_PATH=/path/to/NCTSSOS make oracle-chsh` — regenerate oracle values (see `Makefile` `oracle-%`)
- Documentation work: after any change under `docs/` or any docs-facing content, always run `make examples` to regenerate the example markdowns, then `make servedocs` to preview the docs locally before handoff.

## Architecture Overview
Type hierarchy:
```
AlgebraType
├── MonoidAlgebra        (NonCommutative/Projector/Unipotent)
├── TwistedGroupAlgebra  (Pauli)
└── PBWAlgebra           (Fermionic/Bosonic)
```
Core types: `NormalMonomial{A,T}` (immutable word), `Polynomial{A,T,C}` (mutable map), `VariableRegistry{A,T}` (symbol ↔ index).
Optimization flow: `polyopt()` → `cs_nctssos()` → `compute_sparsity()` → moment/SOS relaxation → JuMP model.

## Coding Style & Naming Conventions
- Indent: 4 spaces; no tabs.
- Naming: `CamelCase` types, `snake_case` functions, `UPPER_SNAKE_CASE` constants.
- Type params: keep algebra type `A` first for dispatch (e.g. `Polynomial{A,T,C}`).
- Prefer surgical changes; add regression tests for bug fixes.

## Testing Guidelines
- Canonical instructions: `TESTING.md`.
- Shared infra: `test/TestUtils.jl` defines `SOLVER` (COSMO) and helpers.
- Suites: `test/polynomials/` (no solver), `test/quality/` (Aqua/ExplicitImports/doctests), `test/relaxations/`, `test/state_poly/`, `test/correlated_sparsity/`, `test/trace_poly/`, `test/problems/`.
- Reviewed expectation fixtures live in `test/data/expectations/*.json`; keep stable case ids and follow `test/data/README.md`.
- Keep tests deterministic and solver-stable (COSMO in CI).

### Backlog-Driven Test Addition Workflow
- The canonical additional-test backlog currently lives outside the repo at:
  - `/Users/exaclior/notes/private-notes/test-expectations/test-examples.typ`
  - rendered companion: `/Users/exaclior/notes/private-notes/test-expectations/test-examples.pdf`
  - supporting paper summaries: `/Users/exaclior/notes/private-notes/test-expectations/references/md/`
- When working from that backlog, process one requested example id at a time unless the user explicitly asks for batching.
- Before adding a case, search existing coverage in `test/` and relevant docs examples under `docs/src/examples/literate/` to avoid duplicates. Prefer extending an existing nearby testset over creating a parallel near-copy.
- If a candidate overlaps an existing test, strengthen the existing test with the missing assertion/input instead of adding repetitive coverage.
- Route note examples to repo locations by feature, not by paper:
  - Algebra / rewrite / PBW / normal-ordering identities (including addendum items on `ComposedMonomial`, bosonic normal ordering, fermionic sign and zero rules) → `test/polynomials/composed_monomial.jl`, `test/polynomials/simplify.jl`, `test/polynomials/monomials.jl`, and related files under `test/polynomials/`.
  - Relaxation-interface or solver-construction guardrails → `test/relaxations/interface.jl`, `test/relaxations/sos.jl`, `test/relaxations/sparsity.jl`, `test/relaxations/gns.jl`, `test/relaxations/gns_pipeline.jl`, `test/relaxations/dualization.jl`.
  - Sparse-vs-dense NC polynomial cases and clique / RIP / sparsity-structure checks (notably `E6`, `E8`, `E9`, `E10`, `CS-failure`, and other CS/TS-focused examples) → `test/correlated_sparsity/` first; only use `test/problems/benchmarks/ncpop_benchmarks.jl` or `test/problems/nc_polynomial/*.jl` when the value being protected is the end-to-end objective rather than the sparsity structure.
  - Small NC eigenvalue benchmark families from the note (`E1`–`E5`, `E7`) → `test/problems/benchmarks/ncpop_benchmarks.jl` for benchmark-family regressions, or `test/problems/nc_polynomial/*.jl` for one-off end-to-end problem checks.
  - Bell / NPA / projector-algebra eigenvalue benchmarks (`E11`, `E12`, `B1`–`B12`) → `test/problems/bell_inequalities/`; if the same scientific example also has a trace or state formulation, add companion checks in `test/trace_poly/` or `test/state_poly/` rather than cloning the same assertion in multiple problem files.
  - Trace-polynomial examples from the note (`T1`–`T8`, Motzkin / trace-positive gap cases) → `test/trace_poly/` for hierarchy/API behavior and `test/problems/trace_polynomial/*.jl` for literature-style end-to-end regression cases.
  - State-polynomial examples from the note (`S1`–`S7`) → `test/state_poly/` for hierarchy behavior, `test/problems/state_polynomial/state_polynomial.jl` for consolidated literature examples, and `test/problems/quantum_networks/bilocal_networks.jl` for bilocal-network-specific cases.
  - Condensed-matter and many-body examples from the note (Heisenberg, Ising, PXP, Bose-Hubbard, XY, fermionic chain) → `test/problems/condensed_matter/*.jl` and `test/problems/fermionic/*.jl`.
- Use the note’s own priority buckets as the default placement heuristic:
  - Priority 1 (`E11`, `E7`, `T6`, `B1`, `T3`, `S6`, `E6`, `CS-failure`) are the best default candidates for CI coverage if they are COSMO-stable and reasonably fast.
  - Priority 2/3 large-scale, convergence-table, or hard Bell cases should usually stay out of always-on CI unless reduced to a cheap deterministic core.
- Verify each added case by running the narrowest relevant test file or suite first; run broader `make test` or `julia --project -e 'using Pkg; Pkg.test()'` when shared behavior or multiple suites are touched.
- For solver-backed numeric expectations meant to stay stable, prefer updating/adding reviewed fixtures in `test/data/expectations/*.json` instead of scattering magic numbers through tests.
- For each proposed case, explicitly decide and report whether it belongs in CI, docs, or both:
  - **Include in CI** when the case is deterministic, reasonably fast, COSMO-stable, and protects correctness or a regression.
  - **Include in docs/examples** when the case primarily teaches a public workflow, showcases an API, or is more useful as narrative runnable material than as repeated regression coverage.
  - **Include in both** only when the scenario is central and the CI version can stay minimal; avoid duplicating the same heavy example in two places.
- When promoting a backlog item to docs, route it to the nearest Literate page instead of inventing a new page unnecessarily:
  - Bell / CHSH / projector workflows → `docs/src/examples/literate/bell.jl`
  - trace/state workflows → `docs/src/examples/literate/trace_poly.jl`
  - sparsity / Newton-chip / GNS → `docs/src/examples/literate/sparsity_convergence.jl`, `newton_chip_method.jl`, `gns_construction_guide.jl`, `gns_optimizer_extraction.jl`
  - Pauli / ground-state / condensed-matter workflows → `docs/src/examples/literate/ground_state_energy.jl`, `certify_ground_state_property.jl`, `pauli_algebra_interface.jl`, `pauli_gns_construction.jl`
  - PBW / mixed-algebra showcase material → `docs/src/examples/literate/pbw_algebras_showcase.jl`, `mixed_algebras_tensor_products.jl`
- If a case is promoted to docs, add/update the Literate source, regenerate generated pages with `make examples`, and preview with `make servedocs` before handoff.

## Commit & Pull Request Guidelines
- Commit messages (recent): `feat: ...`, `fix(scope): ...`, `docs: ...`, `ci: ...` (often with `(#PR)` on merge/squash).
- PRs: small diffs, clear description, and a short “How to test” section (e.g., `make test`). Note any solver requirements (e.g., Mosek-only example regeneration).

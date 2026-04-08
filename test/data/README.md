# Test Data

This directory stores reviewed expectation values (fixtures) used by solver-dependent tests.

## Layout

- `expectations/*.toml`: numeric expectations keyed by stable `id` strings.

## TOML schema (v1)

Top-level keys:

- `schema_version` (Int): currently `1`
- `cases` (Array of tables): list of cases

Each case:

- `id` (String): stable identifier referenced by tests
- `expected` (Table):
  - `objective` (Float): expected objective value
  - `sides` (Array[Int], optional): expected flattened moment-matrix block sizes
  - `nuniq` (Int, optional): expected `n_unique_moment_matrix_elements`
  - suite-specific structural fields are allowed (for non-numeric oracle tests)
- `notes` (String, optional but strongly preferred): provenance / verification notes. State where the numbers came from: literature/table, exact analytic formula, exact diagonalization, or reviewed local solver run for implementation-specific structural values.

## Update workflow

- When a verified expectation changes, update the corresponding TOML file and keep `id`s stable.
- Every new fixture should explain provenance in `notes`.
- Be explicit about which values are mathematical/reference values (`objective` from literature, analytic formula, or ED) and which are implementation oracles (`sides`, `nuniq`, affine counts, etc.) captured from a reviewed local run.
- Prefer one TOML file per test suite (e.g. CHSH, NC examples, correlated sparsity).

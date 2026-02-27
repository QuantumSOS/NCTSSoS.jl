# Test Data

This directory stores reviewed expectation values (fixtures) used by solver-dependent tests.

## Layout

- `expectations/*.json`: numeric expectations keyed by stable `id` strings.

## JSON schema (v1)

Top-level object:

- `schema_version` (Int): currently `1`
- `cases` (Array): list of cases

Each case:

- `id` (String): stable identifier referenced by tests
- `expected` (Object):
  - `objective` (Number): expected objective value
  - `sides` (Array[Int], optional): expected flattened moment-matrix block sizes
  - `nuniq` (Int, optional): expected `n_unique_moment_matrix_elements`
- `notes` (String, optional): provenance / verification notes

## Update workflow

- When a verified expectation changes, update the corresponding JSON file and keep `id`s stable.
- Prefer one JSON file per test suite (e.g. CHSH, NC examples, correlated sparsity).


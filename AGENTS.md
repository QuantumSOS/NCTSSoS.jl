# Repository Guidelines

## Runs and validation
- Test commands, suite layout and CI solver policy: [TESTING.md](TESTING.md). Build and oracle commands: [Makefile](Makefile).
- Follow [test/data/README.md](test/data/README.md) for reviewed numeric fixtures.

### Run Policy
- Prefer Mosek for ad hoc solver-backed runs when available; keep COSMO as the CI baseline unless the task is explicitly Mosek-only. Before relying on Mosek, check for `~/mosek/mosek.lic` or a configured Mosek licensing environment variable and verify the license with a trivial solve; importing `MosekTools` alone is not sufficient. Orb setup writes the contents of the `MOSEK_LICENSE` workspace secret to the standard license path.
- Before launching any run, estimate the expected runtime and report it. After the run, report the actual elapsed time and whether it matched the estimate.

## Documentation
- Edit Literate sources in `docs/src/examples/literate/`; their generated pages in `docs/src/examples/generated/` are committed.
- After changes under `docs/` or to docs-facing content, run `make examples` (requires a Mosek license), then `make servedocs` before handoff. Generation updates `docs/examples_stamp.toml`; CI only verifies the stamp.

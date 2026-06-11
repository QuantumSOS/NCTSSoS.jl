# Wiring Changes and Verification

## Nav change: `docs/make.jl`

Add the new manual page in the Manual section, after "Extending Symmetry Support":

```julia
"Manual" => Any[
    ...
    "Extending Symmetry Support" => "manual/extending_symmetry.md",
    "Clifford Symmetry Detection" => "manual/clifford_symmetry_detection.md",  # NEW
    "SDP Relaxation" => "manual/sdp_relaxation.md",
    ...
]
```

## Bib entries: `docs/src/refs.bib`

See [source-material.md](source-material.md) for the list of entries to check/add.

## Execution order

1. Write `docs/src/manual/clifford_symmetry_detection.md` (Deliverable 1)
2. Update `docs/make.jl` nav + add bib entries
3. Update `docs/src/examples/literate/pauli_clifford_symmetry.jl` (Deliverable 2)
4. `make examples` + `make servedocs` + review
5. `make test` for safety

## Verification checklist

- [ ] `make examples` — all 22+ Literate pages regenerate without error
- [ ] `make servedocs` — new manual page renders: equations, step list, cross-links all work
- [ ] Timing comparison in the enhanced example shows meaningful speedup
- [ ] Summary table in the example includes the new 8-site row
- [ ] `make test` — no existing tests break
- [ ] New manual page is reachable from the sidebar nav
- [ ] Cross-links from the new page to SAB manual, extending-symmetry page, and API docs resolve

# Active Task: Match QMBCertify Performance Without Losing Generality

Close the 100–1000× symbolic construction gap against the hand-tuned QMBCertify
script for 1D Heisenberg chains, while keeping NCTSSoS.jl general across
algebras, lattices, and problem types. The bottleneck is symbolic matrix
materialization and generic Wedderburn decomposition — not the SDP solver.

## Plan

- [plan/qmbcertify_parity.md](plan/qmbcertify_parity.md) — full analysis,
  phased implementation plan, risk register, and benchmark targets.

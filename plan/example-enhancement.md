# Deliverable 2: Enhanced Example — Larger System with Timing

**File:** `docs/src/examples/literate/pauli_clifford_symmetry.jl` (modify in place)

## Goal

Keep the existing 2-site and 4-site sections. Append a new section that:
1. Builds a physically interesting Hamiltonian at $N = 8$ or $N = 10$ qubits — a periodic Heisenberg chain or XXZ ladder.
2. Times the dense baseline vs. symmetry-reduced solve using `@elapsed` (not `@benchmark` — no BenchmarkTools dependency in docs).
3. Prints a comparison table: system size, PSD block sizes, moment variable counts, wall-clock time, group order, number of generators.
4. Makes the payoff visceral — e.g. "dense: 3.2s with a single 25×25 block; symmetric: 0.4s with 8 blocks of size ≤5".

## Sizing decision

$N = 8$ periodic Heisenberg chain: 24 Pauli terms, order-1 basis has $3N + 1 = 25$ monomials. Dense moment matrix is $25 \times 25$. Symmetry group should be nontrivial (translation + reflection + spin-flip on a ring = dihedral × Z₂). Sweet spot — large enough for clear speedup, small enough to solve in <30s for docs regen.

If $N = 8$ is too fast to show timing differences, bump to $N = 10$.

## Structure of the new section

```
## Performance comparison: 8-site Heisenberg ring

### Build the Hamiltonian
### Dense baseline (timed)
### SympleQ detection + symmetry-reduced solve (timed)
### Comparison table
### Discussion: what makes the difference
```

## What NOT to change

- Existing 2-site manual SWAP and 4-site auto-detection sections stay intact (pedagogical purpose — incremental API introduction).
- The summary table at the bottom gets updated to include the new 8-site row.

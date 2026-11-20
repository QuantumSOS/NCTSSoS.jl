# Pauli relaxation prototype

`PauliWorkflowRedesign` is an experimental front end to NCTSSoS. It does not change the default solver or export a new API from `NCTSSoS`.

## Compile the relaxation before choosing the solver form

1. **Select the Gram rows.** The sparsity policy specifies which Pauli words can appear in a test operator. The graph-distance route builds them directly from the interaction graph. A sparse row pool changes the relaxation; it is not an exact reduction of the dense hierarchy.
2. **Verify symmetry on what will be solved.** Check candidate actions against the Hamiltonian, observable and selected rows before using them to identify moments or split cones. The bounded detector need not find every symmetry. In a benchmark that requires a periodic layout, a missing translation is a failed audit, not a successful smaller computation.
3. **Add strengthening in the same moment coordinates.** Record the actual stationarity equations, positive-semidefinite blocks and local density constraints. Requested flags are not a specification of the resulting model: adding Gram rows can also make previously unsupported stationarity equations available.
4. **Choose how to submit the plan.** `compile_relaxation` constructs the plan; `solve_ground_state` submits it. Explicit moment and sum-of-squares forms can have different numerical conditioning even for the same relaxation. Keep the model fingerprint, native status, primal/dual objectives and residuals in the original units. An accepted numerical estimate is not a rational certificate.

For a one-factor row-pool comparison, match the objective, symmetry actions, density blocks and stationarity coefficients by physical Pauli word, not just by count. On the 10×10 square, changing the pair stencil from a diamond to a box doubled the site pairs (2,000 to 4,000). Unrestricted assembly also added 115 stationarity equations. The isolated experiment retained the original 49 equations and five density blocks; it was **audited, not solved**.

See [Wang comparison and assumptions](WANG_MATCH.md), [benchmark protocol](benchmark/README.md), and [recorded results](benchmark/RESULTS.md).

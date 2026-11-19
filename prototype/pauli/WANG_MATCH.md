# Wang-matched Pauli workflow

Use `compile_relaxation(H; order=2, wang_matched=true, sparsity, strengthening)`.
The compiler receives no lattice, site order, shift, coordinate, or symmetry
generators. For direct rows use `PauliSparsitySpec(support_mode=:graph_distance,
power_graph_radius=r, ts_algo=:None, rewrite_complete_local_rdm=false)` or
`:local_products` with `local_product_depth`. Omit `sparsity` for the full-row
control. The route does not use MF, MMD, a chordal completion, or a dense Gram
followed by row pruning. Fourier assembly still visits triangles over reduced
row representatives; small row counts alone do not imply small assembly cost.

## Reference and provenance

The references are Jie Wang's [QMBCertify, CeCILL 2.1](https://github.com/wangjie212/QMBCertify/blob/9a3ae67ff2bb6da1c04b72318ecbfcec71db5495/LICENSE.en):

- Chain: [`0b152df` example](https://github.com/wangjie212/QMBCertify/blob/0b152dfcd9bcde3403aa4d967642ce5261c80712/examples/ground_state.jl),
  [`bound_gsp.jl`](https://github.com/wangjie212/QMBCertify/blob/0b152dfcd9bcde3403aa4d967642ce5261c80712/src/bound_gsp.jl),
  [`basic_function.jl`](https://github.com/wangjie212/QMBCertify/blob/0b152dfcd9bcde3403aa4d967642ce5261c80712/src/basic_function.jl).
  The first example uses a **negated** coefficient, order 4, `extra=4`,
  `lso=true`, `pso=3`, `rdm=10`; a frustrated example uses positive J1 and
  J2, order 4, `pso=3`, and `rdm=10`. The separate energy-window cases use
  `rdm=9`. Do not label those window bounds unrestricted ground-energy bounds.
- Square: [`9a3ae67` examples](https://github.com/wangjie212/QMBCertify/blob/9a3ae67ff2bb6da1c04b72318ecbfcec71db5495/examples/example.jl),
  [`sdp.jl`](https://github.com/wangjie212/QMBCertify/blob/9a3ae67ff2bb6da1c04b72318ecbfcec71db5495/src/sdp.jl),
  [`basicfunction.jl`](https://github.com/wangjie212/QMBCertify/blob/9a3ae67ff2bb6da1c04b72318ecbfcec71db5495/src/basicfunction.jl).
  The examples use order 4, `extra=0`: nearest-neighbor with no strengthening,
  and J2=0.3 with `rdm=8`, `pso=0` and default `lso=true`. That square objective
  encodes the diagonal neighbor as site 3 (`[1,7]`), but GSB's LSO/PSO loop
  interprets the same code as an **axial step of two**. This workflow uses the
  actual input `H` for commutators; it does not reproduce that mismatch.

`wang.jl` implements the mathematical conditions independently with NCTSSoS
word products, orbit coordinates and cone builders. No upstream implementation
is copied into NCTSSoS. The benchmark obtains clean detached upstream checkouts
and calls them as dependencies; it does not ship a fork, patch, or shim.

## What matches and what can differ

| Feature | This workflow |
| --- | --- |
| Moment and row symmetry | Verify coefficient/axis-colored graph actions against `H`, the observable, and the actual row family. On a recovered ring or square, verify the ring reflection or square mirrors and swaps from a regular bond-step orbit. Use that space group, not the extra automorphisms of the 4×4 interaction graph. Use verified uniform X/Z sign actions and signed global axis permutations in the little group, and complex conjugation as an antiunitary moment relation and, where possible, a real basis. |
| Momentum cones | Cyclic or product Fourier blocks plus the little group on each representative momentum. Equivalent sectors may merge or split differently from Wang's two explicit Pauli sectors and realifications. The report includes actual translation orders, point-group order and realification doublings. Continuous SU(2) is disabled in this lane, matching `SU2_symmetry=false` in GSB. |
| LSO | `commutator_linear_order=3`, `5`, or `7` selects the maximum perturbation degree. Use anchored odd-axis-count triples globally, quintics on the first `wang_lso_span` sites, and on a ring septics on the first `rdm_positivity_k-1` sites. Drop any commutator requiring an unsupported nonzero word after the moment, PSO and RDM pools are present. Deduplicate projected rows exactly up to scale rather than using Wang's seeded random hash. Report considered, unsupported, zero and emitted counts. |
| PSO | `kms_psd_order` selects nonconstant moment rows of degree at most that value. Apply `½(A†[H,B] + [A†,H]B)` on the same translation/point-group route and shared moment table, without the generic KMS null filter. Wang's source accumulates `Float16` coefficients and may omit diagonal-only moments from its support pool; this implementation does not deliberately round that way or omit such moments. Valid for global ground states, not arbitrary energy-window states. |
| RDM | `rdm_positivity_k=8`, `9`, or `10` uses one cluster: that many consecutive ring sites, or the expanding square-shell prefix (8 is 3×3 minus a corner). Project onto principal fixed-Hamming-weight blocks from weight 0 through `k÷2`. All terms have even counts of each axis after verified signs and conjugation. This matches the PSD domain in `posepsd8!/9!/10!`, not a full RDM PSD constraint without spin-charge projection. |

Spatial discovery is bounded. A nonperiodic, symmetry-broken, or observable-restricted
case may use a smaller verified subgroup without any lattice claim. The Wang
LSO/RDM path rejects cases without the needed site layout and global spin actions;
PSO requires a regular translation group on its rows. Layout recovery requires
an entire equal-coupling bond class that spans the ring/square, choosing the
largest magnitude class that fits. Equal-weight frustrated bonds that hide that
class are not promised a Wang layout. This route does not support Wang's optional
continuous-spin case, user-prescribed target sectors, or energy-window objective.

## Row pools and comparisons

- At order 2, after expanding the triplet sector by axis permutations and
  removing duplicates, chain graph radius `extra+1` matches Wang's words. On
  square-4 radius 4 matches all 1129 words. Square-6 Wang has 3997 words;
  radius 3 has 3673 (36 candidate-only / 72 Wang-only **pair supports**), and
  radius 4 has 4969 (108 candidate-only / 0 Wang-only). Neither is matched.
- `:local_products` retains 73/289/757 rows at depth 2 on chain-8/square-4/
  square-6. It can preserve the finite spin actions while dropping continuous
  rotations. In a separate isotropic three-spin check, graph rows preserved
  three continuous generators and local products none, though both gave −4.
- Graph-distance selects no ambient triples or plaquettes at higher order.
  Connected products can reach some of the same site supports, but not all
  axis assignments. This is **not** a reproduction of the order-4 GSB examples;
  it matches their applicable strengthening rules on the selected row pool.
- The benchmark's *matched-support* lane requires an actual word-set comparison
  and the expected finite symmetry checks, not equality of numerical objectives.
  `:local_products` and the square-6 radius choice are exploratory. Even a
  matched-support strengthened pair does not establish SDP equivalence: stock
  GSB does not return the emitted LSO rank, and its PSO arithmetic and pool
  bookkeeping differ. Record the active forms, not just requested flags.

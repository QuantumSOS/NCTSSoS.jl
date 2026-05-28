# [Extending Symmetry Support: The Fermionic Case](@id extending-symmetry)

The companion page [Symmetry-Adapted Basis](@ref symmetry-adapted-basis)
documents what the **dense, ordinary-polynomial symmetry path** does and
intentionally does not. This page is a contributor-facing audit of what remains
hard when extending that machinery beyond the narrow supported actions, using
fermionic problems (`FermionicAlgebra <: PBWAlgebra`) as the worked example.

The motivation for being concrete: every limitation listed below is enforced
by a fail-fast check in `src/optimization/symmetry.jl`. Lifting them is not
mysterious; it is a list of specific edits with specific failure modes. The
purpose of this page is to make that list legible so that the relative cost
of each step is clear before any code is written.

The discussion is anchored on fermions because that is the most-requested
extension, but most of the items below are not fermion-specific. Each item is
labeled as either **fermion-specific** or **algebra-general** so that the
ordering of work can be done with eyes open.

## [Audience and scope](@id ext-sym-audience)

This page is aimed at contributors planning to extend symmetry-adapted
basis support, and at users who want to understand why a fermionic symmetry
outside the supported mode-permutation / sector / spin path cannot use
`SolverConfig(; symmetry = ...)` today even though the underlying mathematics
is well-understood. It is not a tutorial; it assumes familiarity with the
[Symmetry-Adapted Basis](@ref symmetry-adapted-basis) page and with
`FermionicAlgebra`.

## [Two axes of limitation](@id ext-sym-axes)

Unsupported fermionic extensions fail for two different reasons:
*algebra-specific* normal-ordering issues and *algebra-general* limitations of
the dense symmetry path. Mixing the two makes "what's missing for fermions?"
sound bigger than it is.

- **Algebra-specific**: assumptions that hold for monoid algebras but not for
  PBW algebras such as `FermionicAlgebra`. Examples: "the action of a group
  element on a basis monomial returns a single basis monomial", "the
  symmetric canonicalization is the identity".
- **Algebra-general**: assumptions that the dense symmetry path makes for
  *any* algebra. Examples: "single clique", "no term sparsity", and "no GNS
  reconstruction from symmetry-reduced moments".

The supported fermionic MVP is deliberately the strict mode-permutation /
sector / spin subcase. Lifting the algebra-general items is the same problem
in any algebra and is strictly larger work.

## [The ten items](@id ext-sym-items)

### [1. Lift the algebra gate (algebra-specific, trivial)](@id ext-sym-1-gate)

Historically, every public entry point in `src/optimization/symmetry.jl` was
gated on `A <: MonoidAlgebra`. The current code has explicit gates for the
supported monoid, Pauli, and fermionic paths, but any new algebra/action pair
still has to audit the same call chain:

- `_check_symmetry_mvp_support` (in `src/optimization/interface.jl`)
- `moment_relax_symmetric`
- `_act_monomial`, `_act_polynomial`
- `_check_polynomial_invariance`, `_check_symmetry_invariance`
- `_check_basis_closure`
- `_build_orbit_reducer`, `_orbit_reduce_polynomial`, `_orbit_reduce_matrix`
- `_reduce_constraint_matrix_symmetric`
- `_collect_reduced_basis`, `_add_moment_eq_constraints_symmetric!`

Mechanically relaxing a gate is a few-line change. It is also **meaningless on
its own**: every item below describes something that would then be silently
wrong rather than loudly wrong.

### [2. The PBW expansion problem (algebra-specific, the real core change)](@id ext-sym-2-pbw)

Today:

```julia
function _act_monomial(g, mono::NormalMonomial{A,T}) where {A<:MonoidAlgebra,T}
    ...
    image = NormalMonomial{A,T}(simplify(A, word))   # ← single word
    return sign, image                               # ← (Int sign, single monomial)
end
```

For a monoid algebra, `simplify(A, word)` returns a single normal-form word.
For fermions, `simplify(FermionicAlgebra, word)` returns
`Vector{Tuple{Int, Vector{T}}}` — a *sum* of normal-ordered monomials with
integer coefficients (from anticommutation contractions). The "single
monomial out" return type is wired into:

- `_act_polynomial`, which pushes `(coef * sign, image)` and assumes one image,
- `_check_polynomial_invariance`, where `_act_polynomial(g, poly) == poly`
  presumes a polynomial-typed result,
- `_check_basis_closure`, where `image in lookup` only makes sense when
  `image` is a single basis element,
- `_build_orbit_reducer`, where orbits are equivalence classes of *monomials*
  with a single ±1 phase per element. If acting on a basis monomial produces
  a sum, "orbit" becomes "linear span" and the phase dictionary is no longer
  the right primitive.

There is one important escape hatch. If the user-supplied symmetry is restricted
to **pure mode permutations with optional ±1 phases per mode**, *and* the
action does not flip creation ↔ annihilation, then permuting a normal-ordered
fermionic word produces another single normal-ordered word: no two operators
on the same mode get colocated, so no contractions occur. In that strict
subcase `_act_monomial` *can* return a single monomial, and item 2 collapses
to "renormal-order after permuting and track the resulting sign carefully".

That subcase is the realistic fermionic MVP. Anything outside it (particle-hole
swaps, mixed mode/parity actions that cause re-contraction) requires a
polynomial-valued action. At that point the orbit reducer is no longer the
right primitive — you want an equivariant projector built from the group
representation directly.

### [3. The user-facing `SignedPermutation` cannot say what fermions need (algebra-specific)](@id ext-sym-3-signed-perm)

Fermionic operators have **signed integer indices**: `+i` is annihilation
`a_i`, `-i` is creation `a_i†`. The "sign" of the index is structural
(creation vs annihilation), not a free user-supplied phase.

A [`SignedPermutation`](@ref NCTSSoS.SignedPermutation) lets the user attach
an extra ±1 phase per index, so `1 => (-1, 1)` means "send `a_1` to `-a_1`" —
fine. But three things a fermionic user will reach for don't fit the current
type:

- **Mode permutation that auto-pairs `+i` and `-i`.** Today, swapping modes 1
  and 2 requires the user to write `1 => 2, 2 => 1, -1 => -2, -2 => -1`. Forget
  any of those and you get a non-bijection error or — worse — silent wrongness
  if `_validate_signed_permutation` only sees half the indices in the active
  domain. A `mode_permutation(1 => 2, 2 => 1)` constructor that materializes
  all four entries would carry its weight.
- **Particle-hole `a ↔ a†`.** Cannot be expressed at all in the current type.
  That image would be `+i ↔ -i` of the same magnitude. The current type
  permits the dictionary entry but the per-letter sign-tracking gets the wrong
  answer because re-normal-ordering after a creation/annihilation flip
  introduces residual signs that the current `_act_monomial` does not compute.
- **U(1) phase `a → e^{iθ} a`.** Continuous group; out of scope on purpose,
  but worth saying so explicitly so future contributors don't try to bolt it
  on top of `SignedPermutation`.

### [4. The decomposition is hardcoded to real characters (algebra-general, but bites fermions first)](@id ext-sym-4-real-chars)

In `_sw_decompose_half_basis`:

```julia
SymbolicWedderburn.symmetry_adapted_basis(Float64, group, action, basis)
```

`Float64` is baked in. For fermionic moment problems with complex coefficients
this is wrong. Two coordinated fixes are needed:

- Thread `MP_C` (or `Complex{Float64}` when `_is_complex_problem(A)`) through
  `_sw_decompose_half_basis` and the transformed-block reduction helpers.
- Switch the cone for the reduced blocks. The existing
  `_is_complex_problem(A) ? :HPSD : :PSD` already covers the constraint cone,
  but the local sign tracking in `_OrbitReducer.phase` is `Int`-valued. For
  genuinely complex-character actions (most fermionic point groups don't have
  these; spin SU(2) does), the orbit "phase" needs to become a complex
  number, not ±1.

### [5. Hermitian moments and conjugation (algebra-specific in degree, algebra-general in shape)](@id ext-sym-5-hermitian)

Fermionic moment matrices are complex Hermitian. The current orbit reducer
treats two monomials related by the group action as having the same moment
up to ±1. It does **not** know about `m ↔ m†`. For a real signed-permutation
action that *commutes* with the adjoint (the common case), the reduction is
still sound because `g \cdot m†` equals `(g \cdot m)†`, but you lose the
*further* reduction of identifying `m` with `m†*` (which the non-symmetric
path does separately via Hermitian conjugation of the moment matrix entries).

This is a *missed* reduction, not a *wrong* one. The fix lives in
`src/algorithms/canonicalization.jl`:

```julia
symmetric_canon(m::NormalMonomial{FermionicAlgebra}) = copy(m.word)   # identity
```

Replacing this with "pick the smaller of `m` and `m†` under some canonical
order, with the corresponding sign" lets the orbit reducer benefit from
conjugation identification on top of the group action. The same fix would
benefit `PauliAlgebra` and `BosonicAlgebra`; this is not symmetry-specific,
but it surfaces here.

### [6. `LinearAlgebra.adjoint(NormalMonomial{<:PBWAlgebra})` throws (algebra-specific, low-priority for symmetry)](@id ext-sym-6-adjoint)

In `src/types/monomial.jl`:

```julia
LinearAlgebra.adjoint(m::NormalMonomial{A,T}) where {A<:PBWAlgebra,T<:Signed} =
    throw(ArgumentError("adjoint not implemented for PBWAlgebra"))
```

The good news: the symmetry path does not currently call `adjoint` on a
`NormalMonomial`. The constraint-matrix construction (in both the regular and
the symmetric path) goes through `_neat_dot3(row_word, mono, col_word)` and
`_simplified_to_terms`, where the row word is already adjointed at the word
level. So this is **not a blocker for the symmetry path itself**, but it
will block any post-hoc verification or pretty-printing tooling that uses
`adjoint(::NormalMonomial)` for fermions. Worth fixing in passing because the
implementation is well-understood (reverse the word, flip index sign,
re-normal-order; for an already normal-ordered input the result is a single
monomial up to sign).

### [7. Domain enumeration assumes a flat index set (algebra-specific, ergonomic)](@id ext-sym-7-domain)

`_symmetry_domain` returns a `Vector{T}` of every index that appears anywhere.
For fermions, "every index" includes both `+i` and `-i` for every active mode,
and the bijection check on a `SignedPermutation` runs over that flat set.
That is sound but pushes a coordination burden onto the user (item 3 above).
The fix is the same `mode_permutation` helper proposed in item 3, plus a
slightly smarter `_symmetry_domain` that, for `FermionicAlgebra`, also
asserts that the user-supplied permutation respects the `+i ↔ -i` pairing.

### [8. Higher-dimensional symmetry blocks (algebra-general)](@id ext-sym-8-mult-free)

The dense symmetry path now keeps the block sizes returned by
`SymbolicWedderburn` and checks transformed off-block coupling. That lifted the
old scalar-only restriction for ordinary dense solves.

What remains hard is carrying those richer symmetry coordinates through the
parts of the package that still expect a dense unreduced moment table:

- GNS reconstruction and optimizer extraction from symmetry-reduced moments.
- Sparse/clique decompositions where each clique may only preserve a subgroup.
- Any future action whose representation has nontrivial complex phases rather
  than the current ±1 orbit bookkeeping.

This remaining work is not fermion-specific.

### [9. Composition with sparsity (algebra-general, practically critical for fermions)](@id ext-sym-9-sparsity)

Today symmetry requires `cs_algo = NoElimination()`,
`ts_algo = NoElimination()`, **single clique**, no term-sparsity block
splitting. Fermionic problems above ~6 sites are essentially unusable without
correlative sparsity, so even if everything else worked, the dense-only
restriction means the fermionic MVP would be a pedagogical demo rather than
a practical tool.

Lifting it requires per-clique symmetry groups (a subgroup that respects each
clique's variable subset) and per-clique decomposition, plus thinking
carefully about cross-clique constraints. Not algebra-specific, but a real
practical blocker for fermions.

### [10. Tests, fixtures, and an example (algebra-specific deliverables)](@id ext-sym-10-tests)

If items 1, 2 (escape-hatch form), 4, and 7 land, the contribution should
include:

- A `test/problems/fermionic/...` regression with a known ground-state energy
  under a small site-permutation symmetry. The 2-site Hubbard model at half
  filling is the obvious choice — translation/inversion symmetry is real,
  multiplicity-free in trivial cases, and the answer is known analytically.
- An expectation fixture in `test/data/expectations/`.
- A docs example, extending the existing
  [Fermionic Ground State](@ref fermionic-ground-state) Literate page rather
  than creating a new one (per the routing convention in `AGENTS.md`).

## [Honest verdict](@id ext-sym-verdict)

**For the supported fermionic MVP** — dense, single clique, pure mode
permutations plus optional sector/spin splitting — the critical restriction is
that acting on a basis monomial must still produce a single normal-ordered
monomial up to sign. That is the line that keeps the orbit reducer simple.

**For anything more ambitious** — particle-hole, continuous U(1) phases,
sparsity composition, or GNS reconstruction from reduced moments — you are no
longer in MVP territory. Those features require a richer intermediate
representation than a monomial orbit representative with a ±1 phase.

**The single highest-leverage piece of work that is not symmetry-specific**
is implementing `symmetric_canon` (item 5) and `LinearAlgebra.adjoint`
(item 6) properly for `NormalMonomial{FermionicAlgebra}`. Both fall out of
the same word-level adjoint helper, both have utility well beyond the
symmetry path, and both turn current "no-op or throw" stubs into something
the symmetry code can lean on.

## [Recommended ordering for one focused week](@id ext-sym-week)

If a contributor had one week and wanted to extend the current fermionic path
without making a mess:

1. Item 5 (real `symmetric_canon` for fermions) and item 6 (real
   `adjoint(::NormalMonomial{FermionicAlgebra})`) — both fall out of the same
   word-level helper, both are reusable beyond symmetry. **Day 1.**
2. Pick one new finite action and prove, with a property test, whether it maps
   every supported basis monomial to a single normal-ordered monomial up to
   sign. If not, stop; it needs a polynomial-valued action, not a patch.
   **Days 2–3.**
3. Thread complex character support only if the chosen action actually needs
   it. Do not make the whole path complex because it sounds general. **Day 4.**
4. Add one reviewed solver regression and one docs example in the nearest
   existing Literate page. **Day 5.**

What is **not** in that week: particle-hole/U(1) actions that expand basis
monomials, GNS/reconstruction from symmetry-reduced moments, and sparsity
composition. Each of those is its own project.

## [What this looks like for other algebras](@id ext-sym-other-algebras)

The same audit applies to the other non-monoid algebras with predictable
substitutions:

- **`PauliAlgebra` (`TwistedGroupAlgebra`)**: the supported path uses
  [`CliffordSymmetry`](@ref NCTSSoS.CliffordSymmetry), where Clifford
  conjugation maps each Hermitian Pauli word to another Pauli word with a real
  sign after simplification. General phase-valued or continuous actions remain
  outside the MVP for the same reason as fermionic U(1) phases.
- **`BosonicAlgebra` (`PBWAlgebra`)**: structurally analogous to fermions
  (PBW expansion is a sum of monomials), but bosonic normal ordering does
  not introduce signs from anticommutation, so item 2's sign bookkeeping is
  simpler. Modes are infinite-dimensional, which interacts badly with item 8
  (multi-dimensional irreps) for any non-trivial symmetry; the bosonic
  symmetry MVP is, in practice, even narrower than the fermionic one.

## See also

- The companion manual page describing what *is* implemented:
  [Symmetry-Adapted Basis](@ref symmetry-adapted-basis).
- The runnable reference example:
  [CHSH with Symmetry Reduction](@ref chsh-symmetry).
- Public types touched by this work: [`SignedPermutation`](@ref NCTSSoS.SignedPermutation),
  [`SymmetrySpec`](@ref NCTSSoS.SymmetrySpec),
  [`SymmetryReport`](@ref NCTSSoS.SymmetryReport).
- Fermionic algebra reference: the
  [Fermionic Ground State](@ref fermionic-ground-state) example, and the
  canonicalization stubs in `src/algorithms/canonicalization.jl`.

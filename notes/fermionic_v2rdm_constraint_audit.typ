#set page(
  paper: "a4",
  margin: (x: 2.2cm, y: 2.4cm),
  numbering: "1",
)
#set text(font: "New Computer Modern", size: 10.5pt)
#set par(justify: true, leading: 0.62em)
#set heading(numbering: "1.")
#show heading.where(level: 1): it => block(above: 1.4em, below: 0.8em)[
  #set text(size: 14pt, weight: "bold")
  #it
]
#show heading.where(level: 2): it => block(above: 1.0em, below: 0.6em)[
  #set text(size: 11.5pt, weight: "bold")
  #it
]
#show raw.where(block: false): box.with(
  fill: luma(240),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
)
#show raw.where(block: true): block.with(
  fill: luma(245),
  inset: 8pt,
  radius: 4pt,
  width: 100%,
)
#show link: it => underline(text(fill: rgb("#1f4ea1"), it))

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Fermionic Polyopt in NCTSSoS.jl: \
    What Variables and Constraints Actually Get Built?
  ]
  #v(0.3em)
  #text(size: 9.5pt, fill: gray)[
    Audit of the H#sub[4] periodic V2RDM constraint inflation. \
    Repository: `NCTSSoS.jl-h4-periodic-v2rdm-benchmark`
  ]
]

#v(0.5em)

#block(
  fill: luma(248),
  stroke: (left: 3pt + rgb("#1f4ea1")),
  inset: 10pt,
  radius: 2pt,
)[
  *Claim under audit.* The excessive number of constraints in NCTSSoS.jl's
  V2RDM-style formulation is caused by needing to establish equivalency of
  moment-matrix elements after normalization (Wick / normal ordering).

  *Verdict.* Partially correct. The 6.4#sym.times blow-up between distinct
  moments and JuMP scalar variables decomposes into two _structural_ factors
  (one unavoidable, one intrinsic to V2RDM's PQG conditions) plus one
  _hidden_ factor that really is the involution-equivalency issue the claim
  refers to. The honest decomposition is below.
]

= One picture of the pipeline

The build path for any fermionic polynomial optimization problem is:

#align(center)[
  `polyopt(...)` #sym.arrow.r `cs_nctssos(...)` #sym.arrow.r
  `compute_sparsity(...)` #sym.arrow.r \
  `moment_relax(...)` #sym.arrow.r `sos_dualize(...)` #sym.arrow.r solver
]

For `FermionicAlgebra`, `_is_complex_problem(FermionicAlgebra) == true`, so:

- moment matrices live in `:HPSD` (Hermitian PSD) cones,
- equality constraints live in `:Zero` cones,
- the symbolic moment problem is built over `Complex{Float64}`.

Two final passes inside `moment_relax` are fermionic-specific
(`src/optimization/moment.jl:421` and `:425`):

```julia
_add_parity_constraints!(mp)
_add_moment_eq_constraints!(mp, pop, moment_eq_row_bases, moment_eq_row_basis_degrees)
```

= The decision variables

There is _no explicit fermionic creation/annihilation operator_ in the SDP.
Variables are *moments*: complex numbers $y_alpha = lr(chevron.l alpha chevron.r)$,
one per canonical monomial #sym.alpha in the moment-problem basis.

== Primal (moment) form

`_solve_complex_moment_problem` (`src/optimization/moment.jl:614`):

```julia
basis = [symmetric_canon(expval(m)) for m in mp.total_basis]
sorted_unique!(basis)
@variable(model, y_re[1:n_basis], set_string_name=false)
@variable(model, y_im[1:n_basis], set_string_name=false)
@constraint(model, y_re[idx_one] == 1)        # normalization
@constraint(model, y_im[idx_one] == 0)
```

Each canonical basis monomial contributes *2 real JuMP variables*
(`y_re[k]`, `y_im[k]`).

== Dual (SOS) form

`_sos_dualize_hermitian` (`src/optimization/sos.jl:370`) creates one matrix
multiplier per constraint matrix:

- `:HPSD` cone of dim $n$ #sym.arrow.r real PSD lifted matrix of dim
  $2n times 2n$ (Goemans-style block embedding
  $[op("Re")(H), -op("Im")(H); op("Im")(H), op("Re")(H)]$).
- `:Zero` cone of dim $n$ #sym.arrow.r free symmetric multiplier of dim
  $2n times 2n$.

Plus one scalar $b$ (the SOS lower-bound variable). One
coefficient-matching constraint per canonical basis monomial — _real_ and
_imaginary_ parts each.

= The constraint zoo

== Moment-matrix blocks (HPSD)

For each clique $c$ and term-sparsity sub-basis $B_c$:

$ M_c [i, j] = op("simplify")(b_i^dagger dot 1 dot b_j) quad
  (op("HPSD"), thin op("dim") = abs(B_c)) $

After Wick normal ordering, *each entry is a polynomial* — a sum of
normal-ordered monomials with integer signs, possibly of multiple lower
degrees due to contractions. So one entry is _not_ one moment variable; it
is a coupling equation across many.

== Localizing-matrix blocks (HPSD or Zero)

For each constraint $g$ inside its clique, with sub-basis $B$:

$ L[i, j] = sum_t op("coef")_t dot
  op("simplify")(b_i^dagger dot g_t dot b_j) $

The cone is `:Zero` when $g in op("pop").op("eq_constraints")$, else
`:HPSD`.

== Complex Zero-cone splitter

`src/optimization/moment.jl:290`:

```julia
function _zero_constraint_components(mat)
    (!_is_complex_problem(A) || _is_hermitian_poly_matrix(mat)) && return [mat]
    hermitian_part     = (mat + adjoint(mat)) / 2
    skewhermitian_part = (mat - adjoint(mat)) / (2im)
    return [hermitian_part, skewhermitian_part]
end
```

Every non-Hermitian symbolic equality matrix is split into *two*
Hermitian matrices. So one equality on the $g$ side becomes up to 2
equality cones in the SDP.

== Parity superselection (`_add_parity_constraints!`)

`src/optimization/moment.jl:846` — the fermionic-only block:

```julia
A === FermionicAlgebra || return nothing
for (cone, mat) in mp.constraints
    for i in 1:dim, j in 1:dim
        poly = mat[i, j]
        if _has_odd_parity_only(poly)
            _append_constraint!(parity_constraints, :Zero, reshape([poly], 1, 1), P)
        end
    end
end
```

For *every* $(i, j)$ entry in *every* constraint matrix, if all surviving
Wick terms have odd parity, NCTSSoS adds an explicit *1#sym.times 1 zero
constraint*. For an order-2 V2RDM relaxation the moment matrix is roughly
half-filled with odd-parity entries, so this contributes $cal(O)(N^2)$
extra equality cones with $N$ = block size.

The comment in `moment_relax` (lines 376–381) acknowledges the trade-off:
keep odd-parity monomials in the basis _and_ add explicit zero constraints,
rather than filtering the basis. This is an _implementation choice that
grows the constraint count_ instead of shrinking the basis.

== Moment equality constraints (`_add_moment_eq_constraints!`)

`src/optimization/moment.jl:449` — for each state-sector constraint $g$
(typical example: fix particle number $sum_i a_i^dagger a_i = n$), this
adds *one 1#sym.times 1 zero constraint per row basis element $b_i$*. With
$K$ such constraints and a row basis of size $R$ within budget, you get on
the order of $K dot R$ extra equality entries.

= Why the headline numbers look the way they do

This is the part that addresses your original question directly. The
H#sub[4]-Nk2 demo log
(`demos/results/h4_periodic_term_sparsity_cosmo_benchmark.md`) reports:

#table(
  columns: (auto, auto),
  align: (left, right),
  inset: 6pt,
  stroke: 0.5pt + gray,
  [*Quantity*], [*Value*],
  [Distinct moments in moment-matrix entries], [26,817],
  [Basis monomials kept (`mp.total_basis`)], [85,941],
  [JuMP scalar variables (real)], [171,882],
  [JuMP scalar constraints (all types)], [16,831],
  [PSD-block scalar entries], [498,845],
  [Symbolic `moment_relax` peak memory], [10.12 GiB],
  [Primal moment build peak memory], [106.4 GiB],
)

The 6.4#sym.times gap between distinct moments and JuMP scalar variables is
*not* one phenomenon. It factors as:

#align(center)[
  $ underbrace(171{,}882, "JuMP scalar vars") quad
    arrow.l.long.bar^(div 2) quad
    underbrace(85{,}941, "complex moments")
    quad arrow.l.long.bar^(div 3.2) quad
    underbrace(26{,}817, "moments in moment matrix only") $
]

== Step 1 — the $times 2$: HPSD-to-real lift

`FermionicAlgebra` uses Hermitian PSD constraints. JuMP/MOI on the standard
real backend doesn't carry complex variables natively, so
`_solve_complex_moment_problem` creates *two real variables per complex
moment*:

```julia
@variable(model, y_re[1:n_basis], set_string_name=false)
@variable(model, y_im[1:n_basis], set_string_name=false)
```

This factor is structural to the embedding $H in op("HPSD") arrow.l.r.double
[op("Re")(H), -op("Im")(H); op("Im")(H), op("Re")(H)] in op("PSD")$. It is
*not* waste — each complex moment genuinely carries two real degrees of
freedom — but it is also where the Hermitian symmetry "force every moment to
obey real / imaginary structure" is paid for at the variable level.

_Avoidable?_ Only by switching to a backend with native
`HermitianPSDCone` (e.g. Mosek) and short-circuiting the lift in
`_solve_complex_moment_problem` and `_sos_dualize_hermitian`.

== Step 2 — the $times 3.2$: basis is the union over _all_ constraint matrices

This is the part most people get wrong on first reading. Compare the two
basis-construction functions:

`_moment_matrix_basis` (`src/optimization/moment.jl:144`) — used to compute
the headline 26,817:

```julia
for basis in term_sparsities[1].block_bases    # only term_sparsities[1] = moment matrix
    ...
end
```

`_polynomial_total_basis` (`src/optimization/moment.jl:166`) — used to size
the JuMP variables:

```julia
for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities)
    # poly = 1   → moment matrix
    # poly = g_k → localizing matrix for each constraint g_k
    ...
end
# plus, if pop.moment_eq_constraints is nonempty:
#   add monomials from b_i† · g for one-sided moment equalities
```

So `mp.total_basis` is the union of monomials reachable from:

+ $op("simplify")(b_i^dagger dot 1 dot b_j)$ — the moment matrix (this
  contributes the 26,817).
+ $op("simplify")(b_i^dagger dot g_k dot b_j)$ for *every* localizing /
  equality polynomial $g_k$. PQG adds three such matrix conditions
  ($D, Q, G$); each one drags in fresh higher-order monomials when Wick
  expansion meets the bilinear form.
+ $op("simplify")(b_i^dagger dot g_(op("meq")))$ for state-sector
  moment-equality constraints (particle number, $k$-momentum). See
  `_collect_moment_eq_row_bases` and `_truncate_moment_eq_row_bases`.

The 59,124 "extra" monomials in
`total_basis ∖ moment_matrix_basis` are *not duplication*. They are
genuinely new degrees of freedom introduced by Wick-expanding the
localizing constraints. Each is a free moment variable until the localizing
PSD blocks and the moment-eq scalar zero-cones tie it down.

This 3.2#sym.times is *not* about missing involution canonicalization. It
is intrinsic to V2RDM's PQG having many localizing matrices, each one
dragging in fresh products.

_Avoidable?_ No — not without rewriting the relaxation. PQG's three
N-representability matrices are the entire physical content of the V2RDM
bound.

== Step 3 — the hidden $times 2$: the involution-canonicalization gap

This is where the original claim _is_ correct, but it lives _inside_ both
of the numbers above, not as a separate factor between them.

`src/algorithms/canonicalization.jl:19–21`:

```julia
symmetric_canon(m::NormalMonomial{PauliAlgebra})    = copy(m.word)
symmetric_canon(m::NormalMonomial{FermionicAlgebra}) = copy(m.word)   # ← identity!
symmetric_canon(m::NormalMonomial{BosonicAlgebra})  = copy(m.word)
```

vs. the real `MonoidAlgebra` version a few lines below:

```julia
function symmetric_canon(m::NormalMonomial{A}) where {A<:MonoidAlgebra}
    word = copy(m.word)
    word_rev = simplify!(A, reverse(word))
    word_rev < word && return word_rev
    return word
end
```

For non-commuting projector strings the basis is collapsed under
$m arrow.l.r m^dagger$ (i.e. word $arrow.l.r$ reversed word) — so
$lr(chevron.l A B C chevron.r)$ and $lr(chevron.l C B A chevron.r)$ become *the same*
moment variable.

*For `FermionicAlgebra` this collapse is disabled.* `symmetric_canon` is
the identity, so the basis used by both `_solve_complex_moment_problem`
and `_sos_dualize_hermitian` keeps $alpha$ and $alpha^dagger$ (after Wick
normal ordering) as *two distinct complex moments* $y_alpha$ and
$y_(alpha^dagger)$.

*Counting impact.* Inside _each_ of the 85,941 and 26,817 numbers,
fermionic monomials are not collapsed under $m arrow.l.r m^dagger$. Almost
every degree-$gt.eq 2$ fermionic word has a distinct adjoint, so a proper
non-trivial `symmetric_canon` would shrink each number by roughly 2#sym.times.

#table(
  columns: (auto, auto, auto),
  align: (left, right, right),
  inset: 6pt,
  stroke: 0.5pt + gray,
  [*Quantity*], [*Current*], [*With involution dedup*],
  [Distinct moments in moment matrix], [26,817], [#sym.tilde 13,500],
  [`mp.total_basis` length], [85,941], [#sym.tilde 43,000],
  [JuMP real scalar variables], [171,882], [#sym.tilde 86,000],
)

*Mechanism, precisely.* NCTSSoS does *not* add explicit
$y_alpha - overline(y_(alpha^dagger)) = 0$ rows. The equivalence between
$y_alpha$ and $overline(y_(alpha^dagger))$ is enforced *implicitly* by
two structural pieces:

+ The HPSD cone forces $M[i,j] = overline(M[j,i])$ within each block.
  Symbolically $M[i, j]$ simplifies to a polynomial in $y_beta$'s and
  $M[j, i]$ to a polynomial in $y_(beta^dagger)$'s; the Hermitian PSD
  constraint thus implicitly couples them.
+ The SOS dual emits one real and one imaginary $f_alpha = 0$ coefficient
  equation per canonical basis monomial. Because $alpha$ and $alpha^dagger$
  are not merged, both equations are emitted independently; their
  consistency is established only via the dual variable structure of the
  point above.

So the user-facing symptom is exactly what you were told:

- The dual SOS has roughly *2#sym.times more $f_alpha = 0$ equations* than
  the `MonoidAlgebra` formulation would for the same underlying physics.
- Each equation also tends to be denser, because Wick contractions inside
  one $M[i, j]$ polynomial spread one coefficient match across many
  low-degree moments.

This is the precise sense in which "we need to establish equivalency of
moment-matrix elements after normalization" is the dominant _avoidable_
cost source.

= Two cheap experiments to verify this

== Patch test

Locally redefine the fermionic `symmetric_canon`:

```julia
function symmetric_canon(m::NormalMonomial{FermionicAlgebra})
    w = copy(m.word)
    # Use the same involution as MonoidAlgebra, going through the fermionic
    # adjoint convention (negate signs and reverse) — adapt the existing
    # adjoint(::NormalMonomial{FermionicAlgebra}) implementation here.
    w_adj = ...
    return min(w, w_adj)
end
```

then re-run any small fermionic benchmark and compare:

- `length(_sorted_symmetric_basis(mp.total_basis))`
- `JuMP.num_variables(sos.model)`
- `JuMP.num_constraints(sos.model, ...)`

Expectation: roughly a 2#sym.times reduction in JuMP variables and in
`fα_constraints_re/_im` counts.

== Audit script

For an existing `MomentProblem` `mp`, count distinct entries that already
obey $alpha = alpha^dagger$ after `simplify(adjoint(...))`:

```julia
mb = NCTSSoS._sorted_symmetric_basis(mp.total_basis)
adjoint_basis = unique(simplify(NCTSSoS.adjoint(m)) for m in mp.total_basis)
length(mb), length(adjoint_basis), length(intersect(Set(mb), Set(adjoint_basis)))
```

The intersection size relative to `length(mb)` tells you exactly how many
moment variables are duplicated by the missing involution canonicalization
— and hence how much the implicit Hermiticity machinery is doing for
nothing.

= Verdict matrix

#table(
  columns: (1fr, auto),
  align: (left, center),
  inset: 7pt,
  stroke: 0.5pt + gray,
  [*Sub-claim*], [*Status*],
  [Constraint count is dominated by something tied to normalization],
    [*Yes*],
  [The cause is per-element equality $y_alpha = overline(y_(alpha^dagger))$
   style rows],
    [*No* — those rows are not emitted explicitly],
  [The cause is failure to canonicalize moments under the
   $m arrow.l.r m^dagger$ involution after Wick, so the SDP keeps both
   copies and the SOS dual emits #sym.tilde 2#sym.times more
   $f_alpha = 0$ equations],
    [*Yes — main avoidable cause*],
  [Per-entry parity superselection adds $cal(O)(N^2)$ extra zero
   constraints],
    [*Yes — second contributor*],
  [Per-row state-sector moment equalities
   (`_add_moment_eq_constraints!`) add $cal(O)(K dot R)$ extra zero
   constraints],
    [*Yes — third contributor*],
  [Hermitian #sym.arrow.r real PSD lift doubles every PSD block
   dimension],
    [*Yes — affects size, not count*],
)

= Decomposition summary

#table(
  columns: (auto, auto, 1fr, 1fr),
  align: (center, center, left, left),
  inset: 7pt,
  stroke: 0.5pt + gray,
  [*Factor*], [*Ratio*], [*Source*], [*Avoidable?*],
  [#sym.times 2], [171{,}882 / 85{,}941],
   [Real / imag split for HPSD on real-only JuMP backends],
   [Only with native complex SDP cones (e.g. Mosek's `HermitianPSDCone`)],
  [#sym.times 3.2], [85{,}941 / 26{,}817],
   [Localizing matrices (PQG: $D, Q, G$) and moment-eq constraints add new
    monomials beyond what the moment matrix indexes],
   [No — fundamental to V2RDM's N-representability conditions],
  [(hidden #sym.times 2 inside both)], [#sym.minus],
   [Fermionic `symmetric_canon` is the identity, so $alpha$ and
    $alpha^dagger$ are kept as two separate complex moments],
   [*Yes* — local fix in `src/algorithms/canonicalization.jl:20`],
)

= Leverage points for a real fix

If you want to attack the inflation rather than just account for it, the
three local interventions are:

+ Implement a non-trivial
  `symmetric_canon(::NormalMonomial{FermionicAlgebra})`
  (`src/algorithms/canonicalization.jl:20`). This kills the hidden
  #sym.times 2 across the entire build.
+ Replace `_add_parity_constraints!` with basis filtering of odd-parity
  monomials. The comment at `src/optimization/moment.jl:376` explicitly
  describes the trade-off you'd be reversing.
+ Vectorize / batch the moment-equality cones in
  `_add_moment_eq_constraints!` so $K dot R$ scalar 1#sym.times 1 cones
  become a single length-$K dot R$ zero cone.

Each of those is local — no algorithmic change to the relaxation itself.

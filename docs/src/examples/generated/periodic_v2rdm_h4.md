```@meta
EditURL = "../literate/periodic_v2rdm_h4.jl"
```

# [Periodic V2RDM Benchmark (H₄ Chain)](@id periodic-v2rdm-h4)

This page is a **specification-first walkthrough** of the periodic H₄
chain benchmark from Schouten, Ewing & Mazziotti (2025, arXiv:2504.02861) —
the first ab initio V2RDM test case NCTSSoS.jl targets.

It has two jobs:

1. Pin down every coefficient, index convention, and spin structure
   needed to build the problem correctly — so that when you rebuild it,
   you have no silent convention drift.
2. Report honestly what the package can and cannot do today: the
   small Hubbard-dimer sanity check runs end-to-end with COSMO, but
   the faithful 32-mode periodic solve is not yet CI-ready.

**Status in one line.** This page is a *spec + smoke test*, not yet
a reproducible periodic run.  Section 9 names the two blockers;
Section 11 is the piece that does run today.

!!! note "Prerequisites"
    This page assumes familiarity with fermionic creation / annihilation
    operators and the canonical anticommutation relations.  If those are
    new, start with
    [Fermionic Ground State (XY Model)](@ref fermionic-ground-state).
    The
    [Hubbard Model](@ref hubbard-model)
    example introduces spin degrees of freedom and the
    `moment_eq_constraints` mechanism used below.

---

## Opening — why bother with periodic V2RDM?

Most quantum-materials calculations run on density functional theory
(DFT), which works well for simple metals and semiconductors but
fails for strongly correlated systems (Mott insulators, high-Tc
superconductors, many 2D materials).

An alternative is to drop the wave function entirely and solve
directly for the **2-electron reduced density matrix**
${}^{2}\!D^{p_1 p_2}_{p_3 p_4} = \langle\Psi|\,a^\dagger_{p_1} a^\dagger_{p_2} a_{p_4} a_{p_3}\,|\Psi\rangle$
(the reversed $p_3 \leftrightarrow p_4$ in the annihilation pair is a
convention — Section 2 explains why).
Because the electronic Hamiltonian has only one- and two-body terms,
the energy is linear in ${}^{2}\!D$:
$E = \text{Tr}({}^{2}\!K\, {}^{2}\!D)$.
The hard part is forcing ${}^{2}\!D$ to come from a real
$N$-electron state — the **$N$-representability** problem.
Approximating it by the necessary conditions
${}^{2}\!D, {}^{2}\!Q, {}^{2}\!G \succeq 0$ (the **PQG conditions** of
Garrod–Percus, 1964) turns the problem into a semidefinite program
(SDP).  This is the **variational 2-RDM** (V2RDM) method.

In the language of polynomial optimisation over a fermionic algebra,
V2RDM is an **order-2 moment relaxation**: the moment matrix at order
2 contains every degree-$\le 4$ expectation, which includes the full
${}^{2}\!D$, ${}^{2}\!Q$, and ${}^{2}\!G$ blocks.

Schouten et al. (2025) extend V2RDM to periodic crystals by
exploiting crystal-momentum conservation to block-diagonalise the
2-RDM.  The H₄ chain is their simplest benchmark: real Coulomb
integrals from a Gaussian basis set, not a lattice model.

---

## Challenge — three things that must line up

Getting a V2RDM benchmark right is mostly about **conventions**.
Three independent choices can all silently flip signs or factors:

1. **Operator ordering** in the two-body term.  The paper writes
   $\hat{a}^{\dagger}_{p_1}\hat{a}^{\dagger}_{p_2}\hat{a}_{p_4}\hat{a}_{p_3}$
   with the two annihilation indices reversed vs. the superscripts
   of ${}^{2}\!V$.  Get it backwards and you invert the sign
   convention of the whole 2-RDM tensor (via
   $a_{p_4}a_{p_3} = -a_{p_3}a_{p_4}$), which silently breaks energy
   evaluation.
2. **Chemist vs. physicist notation** for the Coulomb integrals.
   Extracted data files use chemist notation `(pr|qs)`; the standard
   second-quantised Hamiltonian wants physicist notation
   $\langle pq|rs\rangle$.  The two are related by a single index
   swap, and mixing them gives wrong correlation energies.
3. **Per-unit-cell normalisation.**  A crystal repeats forever, so
   one- and two-electron integrals get divided by $N_k$ and $N_k^2$
   respectively.  Skipping this step produces energies off by the
   number of $k$-points — a 100% error at $N_k = 2$.

Where each is pinned down: § 2 fixes the operator ordering;
§ 8 fixes the chemist↔physicist translation; § 5 fixes the
per-unit-cell normalisation.  § 3 (spin), § 4 (momentum), § 6–7
(reduced Hamiltonian + counting) fill in the surrounding structure.

---

## 1. The physical system

| Parameter | Value |
|:----------|:------|
| Unit cell | 4 hydrogen atoms, equally spaced |
| Bond distance | 1.0 Å |
| Lattice constant | $a = 4.0$ Å (along $x$) |
| Vacuum padding | 10 Å in $y$ and $z$ (models a 1D chain) |
| Basis set | gth-tzvp (Gaussian crystalline basis) |
| Pseudopotential | gth-pade |
| Density fitting | Gaussian density fitting (GDF) |
| Active space | $[4, 8]$ = 4 electrons in 8 spatial MOs per unit cell |

Hydrogen has **no core electrons** — its single $1s$ electron is a
valence electron.  Every electron in the system enters the active
space, so $E_\text{core} = 0$ and the total energy is
$E_\text{tot} = E_\text{active} + E_\text{nuc} + E_\text{Madelung}$,
where $E_\text{Madelung}$ is the classical Madelung correction that
handles the divergent periodic Coulomb sum.

!!! warning "Energy bookkeeping"
    Any SDP run on the active-space Hamiltonian only returns
    $E_\text{active}$.  The paper's figures report the **total**
    energy $E_\text{tot}$.  To compare, you must add the constant
    shift $E_\text{nuc} + E_\text{Madelung}$ (stored as
    `hf_constant_shift` in the vendored expectation file).  The
    regression test `test/problems/fermionic/h4_periodic_v2rdm.jl`
    verifies
    $E^\text{HF}_\text{active} + \text{shift} = E^\text{HF}_\text{total}$
    for the Nk = 2 asset.

---

## 2. The Hamiltonian

The electronic Hamiltonian in $k$-space (paper, Eq. 2) is

```math
\hat{H} = \sum_{\substack{p_1 p_2 \\ k_1 k_2}}
{}^{1}\!H^{p_1 k_1}_{p_2 k_2}\;
\hat{a}^{\dagger}_{p_1 k_1}\hat{a}_{p_2 k_2}
\;+\;
\sum_{\substack{p_1 p_2 p_3 p_4 \\ k_1 k_2 k_3 k_4}}
{}^{2}\!V^{p_1 k_1,\, p_2 k_2}_{p_3 k_3,\, p_4 k_4}\;
\hat{a}^{\dagger}_{p_1 k_1}
\hat{a}^{\dagger}_{p_2 k_2}
\hat{a}_{p_4 k_4}
\hat{a}_{p_3 k_3}
```

Each index $p_i$ is a **spin-orbital**: a spatial orbital $p$, a spin
$\sigma\in\{\uparrow,\downarrow\}$, and a crystal momentum $k$.

!!! warning "Reversed annihilation indices"
    The annihilation operators appear as
    $\hat{a}_{p_4 k_4}\, \hat{a}_{p_3 k_3}$ — with indices
    **swapped** relative to the superscripts of ${}^{2}\!V$.  This
    ordering is chosen so that the 2-RDM matrix-element reads
    naturally:
    ```math
    {}^{2}\!D^{p_1 k_1,\, p_2 k_2}_{p_3 k_3,\, p_4 k_4}
    = \langle \Psi |\,
      \hat{a}^{\dagger}_{p_1 k_1}
      \hat{a}^{\dagger}_{p_2 k_2}
      \hat{a}_{p_4 k_4}
      \hat{a}_{p_3 k_3}
    \,| \Psi \rangle.
    ```
    Swapping $p_3 \leftrightarrow p_4$ is the same as swapping the two
    annihilation operators; by the CAR
    $a_{p_4}\,a_{p_3} = -a_{p_3}\,a_{p_4}$, that flips the sign of the
    whole ${}^{2}\!D$ tensor and produces wrong energies.

**Sanity check — the trace identity.** A fast catch for an index
mistake is the particle-number trace

```math
\sum_{p<q} {}^{2}\!D^{pq}_{pq}
\;=\; \tfrac{1}{2}\, N(N-1).
```

You do not need a wave function to check this.  It is an *operator*
identity from the CAR (worked out in detail below):

```math
\sum_{p<q} a^{\dagger}_p a^{\dagger}_q a_q a_p
  \;=\; \tfrac{1}{2}(\hat N^2 - \hat N)
  \;=\; \tfrac{1}{2}\hat N(\hat N - 1).
```

So on any state with $\hat N\,|\Psi\rangle = N\,|\Psi\rangle$ the
expectation is $N(N-1)/2$ automatically.

A **cheap check in code** uses a Hartree–Fock Slater determinant.
One step of the CAR rewrites the diagonal summand as

```math
a^\dagger_p a^\dagger_q a_q a_p \;=\; n_p n_q \;-\; \delta_{pq}\, n_p,
\qquad n_p \equiv a^\dagger_p a_p,
```

so ${}^{2}\!D^{pq}_{pq} = \langle n_p n_q\rangle$ whenever $p\ne q$.
The remaining claim is that on a Slater determinant the right-hand
side is just a product of two 0/1 occupation bits.  That rests on
one fact worth seeing spelled out.

#### A Slater det is a simultaneous eigenstate of every $n_p$

Write a Slater determinant with occupied set
$\mathcal O = \{i_1, \ldots, i_N\}$ (distinct indices) in second
quantisation:

```math
|\Phi\rangle \;=\; a^\dagger_{i_1}\, a^\dagger_{i_2}\, \cdots\,
  a^\dagger_{i_N}\, |0\rangle.
```

One commutator from the CAR:

```math
[n_p,\, a^\dagger_q]
\;=\; a^\dagger_p\{a_p, a^\dagger_q\} - \{a^\dagger_p, a^\dagger_q\}a_p
\;=\; \delta_{pq}\, a^\dagger_q,
```

which reads: "creating a particle in mode $q$ raises the $n_p$
count by $\delta_{pq}$." Propagate $n_p$ rightwards through all $N$
creation operators, then use $n_p|0\rangle = 0$:

```math
n_p\,|\Phi\rangle
\;=\; \Bigl(\sum_{k=1}^{N} \delta_{p,\,i_k}\Bigr)\,|\Phi\rangle.
```

The indicator sum is 1 if $p \in \mathcal O$ and 0 otherwise —
never 2, because the $i_k$ are distinct ($(a^\dagger_p)^2 = 0$ is
Pauli exclusion).  So every $n_p$ acts as the scalar
$n_p \in \{0,1\}$ on $|\Phi\rangle$.

The number operators also commute pairwise ($[n_p, n_q] = 0$), so
$|\Phi\rangle$ simultaneously diagonalises all of them.  That is
what lets expectations of *products* collapse into products of
scalars:

```math
\langle\Phi|\, n_p n_q\,|\Phi\rangle \;=\; n_p\, n_q
\quad (\text{plain bit product}).
```

This step — "expectation of a product = product of expectations" —
is specific to Slater determinants.  A correlated
$|\Psi\rangle = \sum_I c_I\,|\Phi_I\rangle$ is a superposition of
Slater determinants with *different* occupation sets $\mathcal O_I$,
so it is no longer an $n_p$ eigenstate and the factorisation
breaks; the size of that gap is literally what electron correlation
measures.

#### Back to the trace identity

With $\langle n_p n_q\rangle = n_p n_q$ in hand,

```math
\sum_{p<q}\,{}^{2}\!D^{pq}_{pq}
\;=\; \sum_{p<q} n_p n_q
\;=\; \tfrac{1}{2}\Bigl[\bigl(\textstyle\sum_p n_p\bigr)^2
                         - \sum_p n_p^2\Bigr]
\;=\; \tfrac{1}{2}\bigl(N^2 - N\bigr),
```

using $n_p^2 = n_p$ (bits square to themselves) in the last step.
Trivial to evaluate from an occupation vector — no SDP required.

If your ${}^{2}\!D$ has the wrong sign, $\sum_{p<q} {}^{2}\!D^{pq}_{pq}$
will come out as $-N(N-1)/2$, not $+N(N-1)/2$.  No large Hilbert space
or MPS needed.

---

## 3. Spin structure

### Where the selection rule comes from

The non-relativistic Coulomb operator $\hat V = \sum_{i<j} 1/r_{ij}$
is a purely *spatial* operator — it commutes with $\hat S^2$ and
$\hat S_z$ and acts as the identity on spin.

The key structural fact is that the one-electron Hilbert space is a
**tensor product** of a spatial factor and a spin factor,
$\mathcal{H}_1 = \mathcal{H}_\text{space} \otimes \mathcal{H}_\text{spin}$,
so a spin-orbital is literally a tensor product of a spatial wave
function and a spin ket:

```math
|\phi_{p\sigma}\rangle
\;=\; |\phi_p\rangle \otimes |\sigma\rangle,
\qquad
\phi_{p\sigma}(x) \;=\; \phi_p(\mathbf r)\,|\sigma\rangle,
\qquad x = (\mathbf r, \text{spin}).
```

Every matrix element of an operator of the form
$\hat O = \hat O_\text{space} \otimes \hat{\mathbb 1}_\text{spin}$
therefore factorises into a spatial integral times a spin inner
product — the two factors never talk to each other.  For the
two-electron Coulomb operator:

```math
V^{(p,\sigma_1),(q,\sigma_2)}_{(r,\sigma_3),(s,\sigma_4)}
\;=\;
\underbrace{\bigl\langle p,q \,\bigl|\, 1/r_{12} \,\bigr|\, r,s \bigr\rangle_\text{spatial}}_{V^\text{spatial}_{pq,rs}}
\;\cdot\;
\underbrace{\langle \sigma_1 | \sigma_3 \rangle}_{\text{vertex 1}}
\;\cdot\;
\underbrace{\langle \sigma_2 | \sigma_4 \rangle}_{\text{vertex 2}}
\;=\;
V^\text{spatial}_{pq,rs}\;\delta_{\sigma_1\sigma_3}\;\delta_{\sigma_2\sigma_4}.
```

The Kronecker deltas are just the orthonormality of the spin states
($\langle\uparrow|\downarrow\rangle = 0$).  The same argument for the
one-body Hamiltonian $\hat h$ — which is kinetic energy plus a
spin-independent external potential — gives:

```math
{}^{1}\!H^{(p,\sigma,k)}_{(q,\sigma',k')}
= \delta_{\sigma\sigma'}\;\delta_{kk'}\; h^\text{spatial}_{pq}(k),
\qquad
{}^{2}\!V^{(p,\sigma_1,k_1),\,(q,\sigma_2,k_2)}_{(r,\sigma_3,k_3),\,(s,\sigma_4,k_4)}
= \delta_{\sigma_1\sigma_3}\;\delta_{\sigma_2\sigma_4}\;
V^\text{spatial}_{pq,rs}(k_1 k_2 k_3 k_4).
```

(The $k$-index structure is handled separately in Section 4.)

### Allowed vs. forbidden spin patterns

Read the rule vertex by vertex: **creator and annihilator at the
same vertex must have the same spin**.  In the
$(\sigma_1\sigma_2\sigma_3\sigma_4)$ ordering of ${}^{2}\!V$ — creators
$\sigma_1, \sigma_2$, annihilators $\sigma_3, \sigma_4$ — vertex 1 is
$(\sigma_1, \sigma_3)$ and vertex 2 is $(\sigma_2, \sigma_4)$:

| $\sigma_1\sigma_2\sigma_3\sigma_4$ | Allowed? | What it is |
|:----------|:--------:|:-----------|
| $\uparrow\uparrow\uparrow\uparrow$, $\downarrow\downarrow\downarrow\downarrow$ | ✓ | same-spin pair Coulomb + exchange |
| $\uparrow\downarrow\uparrow\downarrow$, $\downarrow\uparrow\downarrow\uparrow$ | ✓ | opposite-spin pair Coulomb |
| $\uparrow\downarrow\downarrow\uparrow$, $\downarrow\uparrow\uparrow\downarrow$ | ✗ | opposite-spin exchange ($\sigma_1\ne\sigma_3$) |
| $\uparrow\uparrow\downarrow\downarrow$, $\downarrow\downarrow\uparrow\uparrow$ | ✗ | both vertices spin-flip |
| any with $\sigma_1 \ne \sigma_3$ or $\sigma_2 \ne \sigma_4$ | ✗ | single-vertex spin-flip |

Two things that *look* like counterexamples but are not:

- **Opposite-spin repulsion is included.** The
  $\uparrow\downarrow\uparrow\downarrow$ integral is non-zero and
  encodes Coulomb repulsion between a spin-up and a spin-down
  electron.  “Spin-diagonal per vertex” does not mean electrons of
  opposite spin ignore each other.
- **Same-spin exchange is included.** The $\uparrow\uparrow$ channel
  contains both the Coulomb integral $V_{pq,rs}$ and the exchange
  integral $V_{pq,sr}$; antisymmetry of the fermionic creation
  operators builds the exchange term automatically when you sum
  over ordered $(p_1,p_2,p_3,p_4)$.  What’s forbidden is
  *opposite-spin* exchange — that’s the spin-flip at each vertex.

### Counting

Each spatial one-electron integral spawns 2 spin-orbital copies
($\uparrow\uparrow$, $\downarrow\downarrow$).  Each spatial
two-electron integral spawns 4 spin-orbital copies
($\uparrow\uparrow\uparrow\uparrow$,
$\uparrow\downarrow\uparrow\downarrow$,
$\downarrow\uparrow\downarrow\uparrow$,
$\downarrow\downarrow\downarrow\downarrow$,
with the spin ordering $\sigma_1\sigma_2\sigma_3\sigma_4$).

### References

Standard textbook derivations of spin integration and the selection
rule:

- Szabo & Ostlund, *Modern Quantum Chemistry* (Dover, 1996),
  §§2.3–2.4 — spin orbitals and reduction of spin-orbital integrals
  to spatial integrals.
- Helgaker, Jørgensen, Olsen, *Molecular Electronic-Structure
  Theory* (Wiley, 2000), §1.8 — spin-independent Hamiltonians and
  spin-adapted second quantisation.

---

## 4. Momentum conservation

!!! warning "Notation — three uses of subscripted $k$"
    Three things share the letter $k$ on this page.  Quick glossary:

    | Term | Meaning |
    |:-----|:--------|
    | **BZ mesh** | The $N_k$ sample points of the Brillouin zone ($-\pi/a \le k < \pi/a$ in 1D) used by the periodic SCF. |
    | **Mesh values** $\kappa_0, \kappa_1, \ldots$ | Actual points on the mesh.  For $N_k = 2$ on a 1D chain: $\kappa_0 \equiv 0$ (Γ) and $\kappa_1 \equiv \pi/a$ (X). |
    | **Slot labels** $k_1, k_2, k_3, k_4$ | “The momentum of the $i$-th spin-orbital in the four-index integral.”  Each $k_i$ ranges over the mesh values. |
    | **Tuple entries** $(0,1,\ldots)$ | Mesh-point *indices* used in tables below (0 for $\kappa_0$, 1 for $\kappa_1$); matches the `k0`, `k1` labels in code. |

The two-electron integrals vanish unless crystal momentum is
conserved up to a reciprocal lattice vector $\mathbf{G}$:

```math
k_1 + k_2 - k_3 - k_4 \;\equiv\; \mathbf{G} \pmod{\mathbf{G}_\text{lattice}}.
```

This is not imposed by hand — it is what happens when you evaluate a
Coulomb integral between four Bloch orbitals on an infinite lattice.
Each Bloch orbital carries a phase $e^{ik\cdot R}$ under lattice
translation $R$; the integrand picks up
$e^{i(-k_1-k_2+k_3+k_4)\cdot R}$; summing over all $R$ forces this
phase to 1.

For $N_k = 2$ on a 1D chain, the two mesh values give Bloch phases
$e^{i\kappa_0 a} = +1$ and $e^{i\kappa_1 a} = -1$.  So a tuple
$(k_1, k_2, k_3, k_4)$ satisfies conservation iff the product of the
four corresponding signs is $+1$, i.e. iff an **even** number of the
four slots sit at $\kappa_1$.

**All 8 allowed $(k_1, k_2, k_3, k_4)$ tuples**, with each entry a
mesh-point index (0 for $\kappa_0$, 1 for $\kappa_1$).  The
regression test `test/problems/fermionic/h4_periodic_v2rdm.jl` locks
this list in.

| Sector $K = (k_1 + k_2) \bmod N_k$ | $(k_1, k_2, k_3, k_4)$ |
|:----------|:------------------------|
| $K = 0$: even–even | $(0,0,0,0),\ (0,0,1,1),\ (1,1,0,0),\ (1,1,1,1)$ |
| $K = 1$: odd–odd | $(0,1,0,1),\ (0,1,1,0),\ (1,0,0,1),\ (1,0,1,0)$ |

Counting: there are 2 ways to pick $(k_1, k_2)$ with even sum and 2
with odd sum; given that, exactly 2 of the 4 $(k_3, k_4)$ choices
conserve momentum; so the allowed count is
$2 \times 2 + 2 \times 2 = 8$ out of $2^4 = 16$.  Every other
combination — e.g. $(0,0,0,1)$ with a lone $\kappa_1$ — has an odd
number of $-1$ phases and evaluates to zero.

This is a **block-diagonal** structure in the 2-RDM: entries with
different total momentum $K = (k_1 + k_2) \bmod N_k$ live in
separate blocks, and the SDP could in principle optimise each
independently.
The paper's implementation exploits this explicitly to reduce the
SDP size; how much of it NCTSSoS.jl's generic correlative-sparsity
detector recovers from the Bloch Hamiltonian is a current-limitation
issue discussed in Section 9.

---

## 5. Normalisation (per unit cell)

The paper normalises integrals so the SDP objective is energy **per
unit cell**:

```math
h_\text{norm} = \frac{h^\text{spatial}}{N_k}, \qquad
V_\text{norm} = \frac{V^\text{spatial}}{N_k^2}.
```

!!! danger "The repo asset stores un-normalised integrals"
    `test/data/assets/h4_chain_nk2_integrals.txt` stores the raw
    PySCF output (un-normalised).  Forgetting to divide by $N_k$ and
    $N_k^2$ rescales the energy by $N_k$ — at $N_k = 2$ that is a
    factor-of-two error.

---

## 6. The reduced Hamiltonian ${}^{2}\!K$

V2RDM does not optimise using $\hat{H}$ directly.  Instead, the
one-body term is folded into a reduced Hamiltonian so the whole
objective is a single trace:

**Notation bridge.** From here on we switch from the paper's
labels $(p_1, p_2, p_3, p_4)$ to the shorter $(p, q, r, s)$ under
the identification $p_1 \to p,\ p_2 \to q,\ p_3 \to r,\ p_4 \to s$.
Under this relabelling the paper's operator ordering
$a^\dagger_{p_1} a^\dagger_{p_2} a_{p_4} a_{p_3}$ becomes
$a^\dagger_p a^\dagger_q a_s a_r$, and the tensor indexing
${}^{2}\!D^{p_1 p_2}_{p_3 p_4}$ becomes ${}^{2}\!D_{pq,\,rs}$
(same object, same Section-2 index reversal, shorter letters).

```math
{}^{2}\!K_{pq,\,rs}
\;=\; \frac{h_{pr}\,\delta_{qs}}{N-1}
\;+\; \tfrac{1}{2}\, V_{pq,\,rs},
\qquad
E \;=\; \sum_{pqrs}{}^{2}\!K_{pq,\,rs}\; {}^{2}\!D_{pq,\,rs}.
```

The $1/(N-1)$ comes from the contraction
$^{1}\!D_{pr} = \tfrac{1}{N-1}\sum_q {}^{2}\!D_{pq,\,rq}$, which
lets the one-body term be written as a linear functional of
${}^{2}\!D$.  The factor $1/2$ avoids double-counting of pairs under
the unrestricted sum over $p,q,r,s$ (Mazziotti, PRA 57, 4219,
1998).  For the H₄ chain with $N_k = 2$ and 4 electrons per cell,
$N = 8$ so $1/(N-1) = 1/7$.

Equivalent antisymmetrised-pair forms using the Grassmann wedge
product appear in the paper; they are equivalent on antisymmetric
pair indices.  The regression helper `build_k2` in
`test/H4PeriodicAssets.jl` uses the full-tensor form above.

!!! tip "You usually do not build ²K directly"
    NCTSSoS.jl accepts the full operator form

    ```math
    \hat H \;=\; \sum_{ps} h_{ps}\, a^\dagger_p a_s
      \;+\; \tfrac{1}{2}\sum_{pqrs} V_{pq,\,rs}\,
      a^\dagger_p a^\dagger_q a_s a_r,
    ```

    and computes its moments internally.  The ${}^{2}\!K$ form is
    only needed when you compare against a hand-rolled V2RDM code or
    against the paper's tables.
    (Annihilator order $a_s a_r$ matches the index reversal fixed in
    Section 2.)

---

## 7. Counting for $N_k = 2$

| Quantity | Formula | Value |
|:---------|:--------|------:|
| Spatial orbitals / cell | $r$ | 8 |
| $k$-points | $N_k$ | 2 |
| Spatial orbitals total | $r \times N_k$ | 16 |
| **Spin-orbitals total** ($M$) | $2r \times N_k$ | **32** |
| **Electrons total** ($N$) | $4 \times N_k$ | **8** |
| Antisymmetric pairs | $\binom{M}{2}$ | 496 |
| Non-zero ERI blocks | — | 8 of 16 |

### Momentum blocking of ${}^{2}\!D$

With $R = 2r = 16$ spin-orbitals per mesh point ($\kappa_0$ and
$\kappa_1$ from Section 4):

| Sector | Pair content | Block size |
|:-------|:-------------|:-----------|
| $K = 0$ | $(\kappa_0, \kappa_0)$ antisym. + $(\kappa_1, \kappa_1)$ antisym. | $\binom{16}{2} + \binom{16}{2} = 240$ |
| $K = 1$ | $(\kappa_0, \kappa_1)$ all combinations | $16 \times 16 = 256$ |

The two blocks sum to $240 + 256 = 496 = \binom{32}{2}$, matching
the total pair count. ${}^{2}\!G$, indexed by ordered pairs and
blocking by $\Delta\kappa = \kappa_\text{create} - \kappa_\text{annihilate}$,
splits into two $512 \times 512$ blocks (down from
$1024\times 1024$).

---

## 8. Integral format and chemist–physicist mapping

The pre-extracted data in `h4_chain_nk2_integrals.txt` uses
**chemist notation** for the two-electron integrals:

```math
(p_1 k_1,\, p_3 k_3 \mid p_2 k_2,\, p_4 k_4)
= \int\!\!\int
\phi^*_{p_1 k_1}(1)\, \phi_{p_3 k_3}(1)\;
\frac{1}{r_{12}}\;
\phi^*_{p_2 k_2}(2)\, \phi_{p_4 k_4}(2)\;
d1\,d2.
```

The standard physicist-convention integral is related by a single
index swap:

```math
\langle p\, q \mid r\, s \rangle_\text{physicist}
= (p\, r \mid q\, s)_\text{chemist}.
```

To assemble the physicist-convention Hamiltonian

```math
\hat H \;=\; \sum_{ps} h_{ps}\, a^\dagger_p a_s
\;+\; \tfrac{1}{2}\sum_{pqrs}\langle pq|rs\rangle\,
a^\dagger_p a^\dagger_q a_s a_r,
```

read the ERI array as
$\langle pq \mid rs \rangle = \texttt{eri}[(k_p, k_q, k_r, k_s)][p, r, q, s]$.
(Annihilator order $a_s a_r$ enforces Section 2's index reversal,
and the storage indexing matches the loader in `test/H4PeriodicAssets.jl`,
which reads the chemist-style record $(p\,k_1, r\,k_3\,|\,q\,k_2, s\,k_4)$
into `eri[(k1,k2,k3,k4)][p, r, q, s]`.)

### Observed degeneracies at $\kappa_1 = \pi/a$

From the vendored asset, the diagonal one-electron integrals (in
Hartree, un-normalised) are:

| Mesh point | $h_{11}$ | $h_{22}$ | $h_{33}$ | $h_{44}$ | $h_{55}$ | $h_{66}$ | $h_{77}$ | $h_{88}$ |
|:----------|:---------|:---------|:---------|:---------|:---------|:---------|:---------|:---------|
| $\kappa_0 = 0$ | $-1.563$ | $-1.216$ | $-1.201$ | $-0.398$ | $-0.053$ | $-0.069$ | $-0.349$ | $-0.248$ |
| $\kappa_1 = \pi/a$ | $-1.474$ | $-1.474$ | $-0.625$ | $-0.625$ | $-0.459$ | $-0.459$ | $-0.155$ | $-0.155$ |

The pair-by-pair match at $\kappa_1$ holds only to $\sim 10^{-6}$ Ha
(see `h1e_diag_k1` in the pinned expectation file).  It reflects an
additional symmetry of the Bloch Hamiltonian at the zone boundary;
see Schouten et al. (2025) for details.

---

## 9. Why this page stops at a sketch

Before walking through the sketch in Section 10 or the smoke test
in Section 11, the honest thing to report is **why the faithful
periodic solve is not CI-ready**.  Plugging the Section-10 sketch
into the current `cs_nctssos` path does **not** reproduce the
paper's $E \approx -2.188$ Ha/cell.  Two concrete blockers are
pinned in `test/data/expectations/h4_periodic_v2rdm.toml`.

If you want the explicit 32-mode Hamiltonian and a runnable proof that
its support graph is already complete, see
[Periodic V2RDM: Why Correlative Sparsity Misses the $k$-Blocks](@ref periodic-v2rdm-cs-k-blocks).

The short version is that the paper blocks a **pair-indexed 2-RDM**,
while generic correlative sparsity here sees only **single-mode
co-occurrence in Hamiltonian monomials**.

### Blocker 1 — Correlative sparsity collapses to a single clique

On the Bloch Hamiltonian with 32 spin-orbital modes, the spatial
support graph comes out as the **complete graph** on 16 spatial
orbitals (120 edges).  Every pair of compound orbitals appears in
*some* monomial of the ab-initio Hamiltonian, so `cs_algo = MF()`
finds exactly one clique of size 32.  The generic
correlative-sparsity detector does not, today, distinguish
momentum-conserving monomials from merely-allowed ones, so the
paper's momentum blocking of Section 7 is not recovered
automatically.

### Blocker 2 — The order-2 moment problem does not fit in memory

Absent the block decomposition, the full 32-mode order-2 lift
produces:

| Quantity | Value |
|:---------|------:|
| Order-2 basis size | 2081 |
| Unique order-2 moments | 679121 |
| Moment-problem constraints | 508611 |

A direct Mosek solve on this formulation aborts with
`MSK_RES_ERR_SPACE` (code 1051) after symbolic model construction.
COSMO fares no better on a 2081-dimensional PSD block at this
constraint count.

### What is missing

A faithful $N_k = 2$ periodic V2RDM solve on this asset requires
**either** (a) a spin-adapted / $k$-blocked variable creator that
exposes momentum sectors as distinct operator species, **or** (b) a
sparsity-construction hook that accepts user-provided block
structure instead of rediscovering it from the Hamiltonian's
support graph.  Neither is currently in the public API.  Until
one of them lands, the full periodic solve is a spec, not a run.

!!! note "Paper targets (for eventual comparison)"
    Energy per unit cell at $N_k = 2$ from Schouten et al. (2025),
    Fig. 1: V2RDM $[4,8]$ $\approx -2.188$ Ha/cell
    (HF $-2.107$, MP2 $-2.182$, CCSD $-2.203$, CCSD(T) $-2.206$;
    paper figure uncertainty $\pm 0.002$ Ha).  The full $N_k$
    convergence table and $k$-convergence notes live in the
    vendored reference `test/data/assets/h4_chain_reference.toml`.
    V2RDM sits above CCSD(T) at $[4,8]$ — not a contradiction:
    CCSD(T) is not variational (non-Hermitian similarity transform
    plus perturbative triples), and active-space V2RDM is not a
    pure lower bound because the excluded virtual orbitals carry
    correlation the SDP cannot see.  Full-orbital V2RDM with PQG
    would be a rigorous lower bound; $[4,8]$ is not.

---

## 10. What would the full H₄ chain ($N_k = 2$) look like?

Sections 2–8 give you everything needed to write the Hamiltonian
for the full benchmark — it just won't run end-to-end until the
blockers in Section 9 lift.

!!! tip "Want the executable version of Sections 10–11's recipe?"
    The companion page
    [Periodic V2RDM: Building the Paper's PQG Relaxation](@ref periodic-v2rdm-pqg-construction)
    carries out the *symbolic* construction end-to-end — paper's
    PQG basis (2017 monomials), custom $(\Delta N, K)$ block
    partition matching the paper's $\{240, 240, 256, 256, 512, 513\}$
    sizes, and a fully-built `MomentProblem` whose sizes are
    printed. It stops short of `optimize!`.

- 32 spin-orbital modes — 8 spatial orbitals × 2 spins × 2 mesh
  points ($\kappa_0, \kappa_1$).
- One-body term from `h1e[ik]` (per Section 3, diagonal in spin and
  $k$), divided by $N_k$ (Section 5).
- Two-body term from `eri[(k1,k2,k3,k4)]` for the 8 allowed
  $k$-tuples from Section 4, divided by $N_k^2$, translated from
  chemist to physicist notation (Section 8), and summed over the
  four spin-conserving combinations (Section 3).
- A particle-number constraint, either on the total
  $N_\text{total} = 8$ or on each spin sector
  ($N_\uparrow = N_\downarrow = 4$).

Concretely (the code labels `k0`, `k1` correspond to the mesh
values $\kappa_0 = 0$ and $\kappa_1 = \pi/a$ from Section 4):

```julia
nk = 2
norb = 8  # spatial orbitals per mesh point

registry, ((c_up_k0, c_up_k0_dag),
           (c_dn_k0, c_dn_k0_dag),
           (c_up_k1, c_up_k1_dag),
           (c_dn_k1, c_dn_k1_dag)) = create_fermionic_variables([
    ("c_up_k0", 1:norb),
    ("c_dn_k0", 1:norb),
    ("c_up_k1", 1:norb),
    ("c_dn_k1", 1:norb),
])
# ... build one- and two-body terms as above ...
```

`test/H4PeriodicAssets.jl` ships a loader for the vendored integral
text file, so you can bolt the asset onto your own Hamiltonian
builder without re-running PySCF.

---

## 11. Smoke test: the fermionic SDP path

What *does* run end-to-end today is the smallest fermionic system
that exercises the same pipeline: a **Hubbard dimer** — 2 sites,
4 spin-orbitals, 2 electrons (equivalent to H₂ in a minimal
basis).  Within the $N_\uparrow = N_\downarrow = 1$ singlet sector
the Hamiltonian reduces to a $2\times 2$ block on
$\{|D_+\rangle, |S\rangle\}$ (symmetric doubly-occupied combination
$|D_+\rangle$ and the spin-singlet $|S\rangle$):
$H^{(+)} = \begin{pmatrix} U & -2t \\ -2t & 0 \end{pmatrix}$,
whose smaller eigenvalue is

```math
E_0 \;=\; \frac{U - \sqrt{U^2 + 16t^2}}{2}
      \;\stackrel{t=1,\,U=2}{=}\; 1 - \sqrt{5} \;\approx\; -1.2361.
```

At $N = 2$ the PQG conditions are *exact*, so the order-2 moment
relaxation must reproduce $E_0$ to solver tolerance; any larger
gap is a construction error in the fermionic SDP path.

This is deliberately minimal — a smoke test, not a tutorial.  For a
full Hubbard walkthrough, see [Hubbard Model](@ref hubbard-model).

````julia
using NCTSSoS, COSMO, JuMP

t_hop, U_int = 1.0, 2.0
E_exact = (U_int - sqrt(U_int^2 + 16 * t_hop^2)) / 2

nsites = 2
registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
    ("c_up", 1:nsites),
    ("c_dn", 1:nsites),
])

hopping = -t_hop * (
    c_up_dag[1]*c_up[2] + c_up_dag[2]*c_up[1] +
    c_dn_dag[1]*c_dn[2] + c_dn_dag[2]*c_dn[1]
)
interaction = U_int * sum(
    (c_up_dag[i]*c_up[i]) * (c_dn_dag[i]*c_dn[i]) for i in 1:nsites
)
ham = hopping + interaction

# Half-filled: fix $N_\uparrow = N_\downarrow = 1$ on the state.
n_up = 1.0 * sum(c_up_dag[i]*c_up[i] for i in 1:nsites)
n_dn = 1.0 * sum(c_dn_dag[i]*c_dn[i] for i in 1:nsites)

pop = polyopt(ham, registry;
    moment_eq_constraints = [n_up - 1.0*one(ham),
                             n_dn - 1.0*one(ham)])

SOLVER = optimizer_with_attributes(COSMO.Optimizer,
    "max_iter" => 50_000, "eps_abs" => 1e-7, "eps_rel" => 1e-7,
    "verbose" => false)
result = cs_nctssos(pop, SolverConfig(optimizer = SOLVER, order = 2))

println("SDP objective: $(round(result.objective; digits=6))")
println("Exact:         $(round(E_exact; digits=6))")
````

````
SDP objective: -1.236068
Exact:         -1.236068

````

COSMO typically reaches the closed-form $E_0$ to $\sim 10^{-5}$ on
this problem; the `1e-3` assertion is a deliberately loose CI
guard.  What this confirms is that the fermionic SDP *path*
— variable creation → Hamiltonian arithmetic → `polyopt` with
moment constraints → `cs_nctssos` order-2 relaxation → SDP solve
— is wired correctly.  The full periodic problem uses the same
code path, but scaling up is not the only missing piece: as
Section 9 spells out, it also needs explicit momentum-block
structure that the generic sparsity detector does not recover.

---

## Summary

### What runs today

The Section-11 smoke test: build a fermionic Hamiltonian with
[`create_fermionic_variables`](@ref), impose particle number with
`moment_eq_constraints` in [`polyopt`](@ref), solve order-2 with
[`cs_nctssos`](@ref) + COSMO.  Matches $1 - \sqrt 5$ to solver
tolerance on the dimer.

### Locked-down specification for the full benchmark (Sections 2–8)

| Fact | Where |
|:-----|:-----|
| Operator ordering $a^\dagger a^\dagger a_{p_4} a_{p_3}$ | Section 2 |
| Spin is diagonal per vertex, opposite spins still repel | Section 3 |
| 8 allowed $k$-tuples out of 16 | Section 4 |
| $h/N_k$, $V/N_k^2$ normalisation | Section 5 |
| ${}^{2}\!K = h\delta/(N-1) + V/2$ | Section 6 |
| Chemist `(pr|qs)` vs. physicist $\langle pq|rs\rangle$ | Section 8 |

### What does not run yet (Section 9)

| Blocker | Source |
|:--------|:-------|
| `cs_algo = MF()` collapses 32-mode Bloch Hamiltonian into one clique of 32 | `h4_periodic_v2rdm.toml` → `spin_orbital_order2_blocker` |
| Order-2 moment problem (508611 constraints) runs Mosek out of memory | `h4_periodic_v2rdm.toml` → `correlative_sparsity_only_attempt` |
| Spin-adapted / $k$-blocked formulation | not yet in public API |

### See also

- [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
  creation/annihilation operators, CAR, and parity constraints.
- [Hubbard Model](@ref hubbard-model) — spin degrees of freedom,
  particle-number constraints, `moment_eq_constraints`.
- [Kitaev Chain](@ref kitaev-chain) — pairing terms and topological
  phases in fermionic systems.

### References

- A. O. Schouten, S. Ewing, D. A. Mazziotti, "Bootstrapping the
  Electronic Structure of Quantum Materials," arXiv:2504.02861
  (2025). — The paper this benchmark comes from; extends V2RDM to
  $k$-space and applies it to H₄ chains, MoS₂, and NiO.
- D. A. Mazziotti, "Contracted Schrödinger equation: Determining
  quantum energies and two-particle density matrices without wave
  functions," *Phys. Rev. A* **57**, 4219 (1998). — Defines the
  reduced Hamiltonian ${}^{2}\!K$.
- C. Garrod, J. K. Percus, "Reduction of the $n$-particle variational
  problem," *J. Math. Phys.* **5**, 1756 (1964). — The original PQG
  2-positivity conditions.
- S. Ewing, D. A. Mazziotti, "Correlation-driven phenomena in
  periodic molecular systems from variational two-electron reduced
  density matrix theory," *J. Chem. Phys.* **154**, 214106 (2021).
  — First periodic V2RDM implementation; precursor to the 2025
  $k$-space generalisation.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


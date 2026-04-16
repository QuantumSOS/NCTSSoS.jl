```@meta
EditURL = "../literate/periodic_v2rdm_h4.jl"
```

# [Periodic V2RDM Benchmark (H₄ Chain)](@id periodic-v2rdm-h4)

This page documents the **Hamiltonian specification** for the periodic
hydrogen chain benchmark from Schouten, Ewing & Mazziotti (2025) —
the first ab initio test case for periodic V2RDM in NCTSSoS.jl.

It pins down every coefficient, index convention, and spin structure
so there is no ambiguity when building the SDP.  We then verify the
full V2RDM construction on a **two-site toy model** (Hubbard dimer)
where the exact answer is known.

!!! note "Prerequisites"
    This page assumes familiarity with fermionic creation/annihilation
    operators and the canonical anticommutation relations.  If those are
    new, start with the
    [Fermionic Ground State (XY Model)](@ref fermionic-ground-state)
    example.  The
    [Hubbard Model](@ref hubbard-model)
    example introduces spin degrees of freedom and particle-number
    constraints, which we also use here.

## Why this benchmark?

The **variational 2-electron reduced density matrix** (V2RDM) method
solves for the ground-state energy of a many-electron system by
semidefinite programming — minimise $E = \text{Tr}({}^2\!K\; {}^2\!D)$
subject to N-representability constraints (${}^2\!D \succeq 0$,
${}^2\!Q \succeq 0$, ${}^2\!G \succeq 0$).  This is exactly an
order-2 moment relaxation in the language of NCTSSoS.jl.

Schouten et al. (2025) extend V2RDM to **periodic crystals** by
exploiting crystal momentum conservation to block-diagonalise the
2-RDM — precisely the same trick as correlative sparsity in
polynomial optimisation.  The H₄ chain is their simplest benchmark:
real Coulomb integrals from a Gaussian basis, not a lattice model.

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
space, so $E_\text{core} = 0$.  The total energy is simply
$E_\text{tot} = E_\text{active} + E_\text{nuc} + E_\text{Madelung}$.

---

## 2. The Hamiltonian

The electronic Hamiltonian in $k$-space (paper's Eq. 2) is:

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

Each index $p_i$ is a **spin-orbital**: it encodes a spatial orbital
$p$, a spin $\sigma \in \{\uparrow, \downarrow\}$, and a crystal
momentum $k$.

!!! warning "Operator ordering in the two-body term"
    The annihilation operators appear as
    $\hat{a}_{p_4 k_4}\, \hat{a}_{p_3 k_3}$ — indices **reversed**
    from the superscripts of ${}^{2}\!V$.  This is the standard
    quantum chemistry convention, chosen so that the 2-RDM
    definition reads naturally:
    ```math
    {}^{2}\!D^{p_1 k_1,\, p_2 k_2}_{p_3 k_3,\, p_4 k_4}
    = \langle \Psi |\,
      \hat{a}^{\dagger}_{p_1 k_1}
      \hat{a}^{\dagger}_{p_2 k_2}
      \hat{a}_{p_4 k_4}
      \hat{a}_{p_3 k_3}
    \,| \Psi \rangle
    ```
    Mixing up $p_3 \leftrightarrow p_4$ flips the sign of every
    off-diagonal block and produces wrong energies.

---

## 3. Spin structure

The Coulomb interaction does not depend on spin.  Consequently, both
the one-electron and two-electron integrals are diagonal in spin:

```math
{}^{1}\!H^{(p,\sigma,k)}_{(q,\sigma',k')}
= \delta_{\sigma\sigma'}\;\delta_{kk'}\; h^\text{spatial}_{pq}(k)
```

```math
{}^{2}\!V^{(p,\sigma_1,k_1),\,(q,\sigma_2,k_2)}_{(r,\sigma_3,k_3),\,(s,\sigma_4,k_4)}
= \delta_{\sigma_1\sigma_3}\;\delta_{\sigma_2\sigma_4}\;
V^\text{spatial}_{pq,rs}(k_1 k_2 k_3 k_4)
```

In words: the spin of a created electron must match the spin of the
corresponding annihilated electron.  Spin-up and spin-down electrons
are coupled only through the **Pauli exclusion principle**
(anticommutation relations), not through the integrals themselves.

This means each spatial one-electron integral spawns **2** spin-orbital
copies ($\uparrow\uparrow$ and $\downarrow\downarrow$), and each
spatial two-electron integral spawns **4** spin-orbital copies
($\uparrow\uparrow\uparrow\uparrow$,
$\uparrow\downarrow\uparrow\downarrow$,
$\downarrow\uparrow\downarrow\uparrow$,
$\downarrow\downarrow\downarrow\downarrow$, where the spin
ordering follows $\sigma_1 \sigma_2 \sigma_3 \sigma_4$).

---

## 4. Momentum conservation

The two-electron integrals are **zero** unless crystal momentum is
conserved:

```math
k_1 + k_2 - k_3 - k_4 = \mathbf{G}
\quad\text{(reciprocal lattice vector)}
```

For $N_k = 2$ with $k_0 = 0$ and $k_1 = \pi/a$, each Bloch orbital
picks up a phase $e^{ik \cdot a}$ under a lattice translation:

```math
e^{i k_0 \cdot a} = +1, \qquad e^{i k_1 \cdot a} = -1
```

A two-electron integral picks up a net phase
$e^{i(-k_1 - k_2 + k_3 + k_4) \cdot a}$.  This must equal $+1$
for the integral to survive.  Since each factor is $\pm 1$, the
product of the four signs must be $+1$.

**Allowed $k$-combinations** (8 out of 16 survive):

| $(k_1, k_2, k_3, k_4)$ | Total $K = k_1{+}k_2$ | Status |
|:------------------------|:-----------------------|:-------|
| $(0, 0, 0, 0)$ | $0$ | ✓ |
| $(0, 1, 0, 1)$ | $1$ | ✓ |
| $(0, 1, 1, 0)$ | $1$ | ✓ |
| $(1, 0, 0, 1)$ | $1$ | ✓ |
| $(1, 0, 1, 0)$ | $1$ | ✓ |
| $(1, 1, 0, 0)$ | $0\;\text{mod}\;G$ | ✓ |
| $(1, 1, 1, 1)$ | $0\;\text{mod}\;G$ | ✓ |
| $(0, 0, 0, 1)$ | — | ✗ momentum violated |

This creates **block-diagonal structure** in the 2-RDM: entries
with different total momentum $K = k_1 + k_2 \;\text{mod}\; G$ live
in separate blocks.  The SDP optimises each block independently.

This is **exactly correlative sparsity** in polynomial optimisation.
When the Hamiltonian is written in the Bloch (k-space) basis,
NCTSSoS.jl's `cs_algo` detects the sparsity pattern automatically.

---

## 5. Normalization

The paper normalises integrals so that the energy is **per unit cell**:

```math
h_\text{normalised} = \frac{h_\text{spatial}}{N_k}, \qquad
V_\text{normalised} = \frac{V_\text{spatial}}{N_k^2}
```

!!! danger "Don't forget this step"
    The PySCF-extracted `.npz` files store **un-normalised** integrals.
    Omitting the $1/N_k$ and $1/N_k^2$ factors gives energies that are
    off by factors of $N_k$ — wrong by 100% or more.

---

## 6. The reduced Hamiltonian ${}^{2}\!K$

V2RDM does not optimise with $\hat{H}$ directly.  Instead, it folds
the one-body part into the two-body part using the **reduced
Hamiltonian** (paper's Eq. 5–6):

```math
{}^{2}\!K^{p_1 k_1,\, p_2 k_2}_{p_3 k_3,\, p_4 k_4}
= \frac{4}{N-1}\;
{}^{1}\!H^{p_1 k_1}_{p_3 k_3} \wedge \delta^{p_2 k_2}_{p_4 k_4}
\;+\;
{}^{2}\!V^{p_1 k_1,\, p_2 k_2}_{p_3 k_3,\, p_4 k_4}
```

where $\wedge$ is the **Grassmann wedge product** (antisymmetric
tensor product):

```math
(A \wedge B)^{ij}_{kl}
= \tfrac{1}{2}\bigl(A^i_k\, B^j_l - A^i_l\, B^j_k\bigr)
```

The energy then becomes a single trace:
$E = \text{Tr}({}^{2}\!K\; {}^{2}\!D)$.

For $N_k = 2$: $N_\text{total} = 4 \times 2 = 8$ electrons,
so the prefactor is $4/(N{-}1) = 4/7 \approx 0.5714$.

!!! tip "NCTSSoS.jl handles this automatically"
    When you pass a Hamiltonian with one-body and two-body terms to
    `polyopt`, the library builds the moment matrix from the full
    operator — the ${}^{2}\!K$ folding is implicit.  You only need
    ${}^{2}\!K$ if you're verifying against a hand-built V2RDM code.

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

| Sector | Pair types | Block size |
|:-------|:-----------|:-----------|
| $K = 0$ | $(k_0, k_0)$ antisym + $(k_1, k_1)$ antisym | $\binom{16}{2} + \binom{16}{2} = 240$ |
| $K = 1$ | $(k_0, k_1)$ all combinations | $16 \times 16 = 256$ |

Without blocking: one $496 \times 496$ matrix for ${}^{2}\!D$.
With blocking: two blocks ($240 \times 240$ and $256 \times 256$).

The ${}^{2}\!G$ matrix blocks by momentum *difference*
$\Delta k = k_\text{creation} - k_\text{annihilation}$ and splits
into two $512 \times 512$ blocks (down from $1024 \times 1024$).

---

## 8. Integral format and notation mapping

The pre-extracted data (in `h4_chain_nk2.npz`) uses **chemist
notation** for the two-electron integrals:

```math
(p_1 k_1,\, p_3 k_3 \mid p_2 k_2,\, p_4 k_4)
= \int\!\!\int
\phi^*_{p_1 k_1}(1)\, \phi_{p_3 k_3}(1)\;
\frac{1}{r_{12}}\;
\phi^*_{p_2 k_2}(2)\, \phi_{p_4 k_4}(2)\;
d1\,d2
```

The physicist-notation conversion is:

```math
\langle p\, q \mid r\, s \rangle_\text{physicist}
= (p\, r \mid q\, s)_\text{chemist}
```

The paper's ${}^{2}\!V$ uses:

```math
{}^{2}\!V^{p_1 k_1,\, p_2 k_2}_{p_3 k_3,\, p_4 k_4}
= (p_1 k_1,\, p_3 k_3 \mid p_2 k_2,\, p_4 k_4)_\text{chemist}
= \langle p_1 k_1,\, p_2 k_2 \mid p_3 k_3,\, p_4 k_4 \rangle_\text{physicist}
```

To build the standard physicist-convention Hamiltonian
``H = \sum h_{ps}\, a^\dagger_p a_s
+ \tfrac{1}{2}\sum \langle pq|st\rangle\, a^\dagger_p a^\dagger_q a_t a_s``,
read the ERI array as:

```math
\langle pq \mid rs \rangle = \texttt{eri}[(k_p, k_q, k_r, k_s)][p, r, q, s]
```

One-electron integral diag values (un-normalised, in Hartree):

| $k$-point | $h_{11}$ | $h_{22}$ | $h_{33}$ | $h_{44}$ | $h_{55}$ | $h_{66}$ | $h_{77}$ | $h_{88}$ |
|:----------|:---------|:---------|:---------|:---------|:---------|:---------|:---------|:---------|
| $k_0 = 0$ | $-1.563$ | $-1.216$ | $-1.201$ | $-0.398$ | $-0.053$ | $-0.069$ | $-0.349$ | $-0.248$ |
| $k_1 = \pi/a$ | $-1.474$ | $-1.474$ | $-0.625$ | $-0.625$ | $-0.459$ | $-0.459$ | $-0.155$ | $-0.155$ |

Note the pairwise degeneracies at $k_1$ — a consequence of sublattice
symmetry at the zone boundary.

---

## 9. Verification: Hubbard dimer (2 sites, 4 spin-orbitals)

Before tackling the full H₄ chain, we verify the entire V2RDM
construction on the simplest possible system: **2 electrons in
2 spatial orbitals** — a Hubbard dimer equivalent to H₂ in a
minimal basis.

| Parameter | Value |
|:----------|:------|
| Spatial orbitals | 2 |
| Spin-orbitals ($M$) | 4 |
| Electrons ($N$) | 2 |
| Hopping ($t$) | 1.0 |
| On-site repulsion ($U$) | 2.0 |

The exact ground-state energy is $E_0 = (U - \sqrt{U^2 + 16t^2})/2
= 1 - \sqrt{5} \approx -1.2361$.  For $N = 2$ electrons, PQG
conditions are **exact** — the SDP must reproduce this to solver
tolerance.

### Exact diagonalisation

````julia
using LinearAlgebra

t_hop = 1.0
U_int = 2.0

# Singlet Hamiltonian in the basis {|1↑1↓⟩, |2↑2↓⟩}
H_singlet = [U_int -2t_hop; -2t_hop 0.0]
E_exact = eigmin(H_singlet)
println("Exact ground-state energy: $(round(E_exact; digits=6))")
````

````
Exact ground-state energy: -1.236068

````

### Building the Hamiltonian with NCTSSoS.jl

````julia
using NCTSSoS

nsites = 2
registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
    ("c_up", 1:nsites),
    ("c_dn", 1:nsites),
]);
println("Created $(length(registry)) variables: $(2nsites) spin-orbitals")
````

````
Created 8 variables: 4 spin-orbitals

````

Hopping — kinetic energy

````julia
hopping = -t_hop * sum(
    c_up_dag[i] * c_up[j] + c_up_dag[j] * c_up[i] +
    c_dn_dag[i] * c_dn[j] + c_dn_dag[j] * c_dn[i]
    for (i, j) in [(1, 2)]
)
````

````
-c_dn⁺₂c_dn₁ + -c_dn⁺₁c_dn₂ + -c_up⁺₂c_up₁ + -c_up⁺₁c_up₂
````

On-site repulsion — the quartic Coulomb term

````julia
interaction = U_int * sum(
    (c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
    for i in 1:nsites
)

ham = hopping + interaction
````

````
-c_dn⁺₂c_dn₁ + -c_dn⁺₁c_dn₂ + -c_up⁺₂c_up₁ + -c_up⁺₁c_up₂ + 2.0 * c_dn⁺₂c_up⁺₂c_up₂c_dn₂ + 2.0 * c_dn⁺₁c_up⁺₁c_up₁c_dn₁
````

### Fixing particle number

The half-filled sector has $N_\uparrow = 1$, $N_\downarrow = 1$.
We impose this as a **moment equality constraint**
$(N_\sigma - 1)|\psi\rangle = 0$:

````julia
n_up = 1.0 * sum(c_up_dag[i] * c_up[i] for i in 1:nsites)
n_dn = 1.0 * sum(c_dn_dag[i] * c_dn[i] for i in 1:nsites)

pop = polyopt(ham, registry;
    moment_eq_constraints = [n_up - 1.0 * one(ham),
                             n_dn - 1.0 * one(ham)])
````

````
Optimization Problem (FermionicAlgebra)
────────────────────────────────────
Objective:
    -c_dn⁺₂c_dn₁ + -c_dn⁺₁c_dn₂ + -c_up⁺₂c_up₁ + -c_up⁺₁c_up₂ + 2.0 * c_dn⁺₂c_up⁺₂c_up₂c_dn₂ + 2.0 * c_dn⁺₁c_up⁺₁c_up₁c_dn₁

Equality constraints (0):
    (none)

Inequality constraints (0):
    (none)

Moment equality constraints (2):
    -1.0 + c_up⁺₂c_up₂ + c_up⁺₁c_up₁ = 0
            -1.0 + c_dn⁺₂c_dn₂ + c_dn⁺₁c_dn₁ = 0

Variables (8):
    c_dn⁺₂, c_dn⁺₁, c_up⁺₂, c_up⁺₁, c_up₁, c_up₂, c_dn₁, c_dn₂

````

### Solving the SDP

At order 2, the moment matrix contains all monomials up to degree 4.
For 4 spin-orbitals this includes the full 2-RDM
$\langle a^\dagger_p a^\dagger_q a_t a_s \rangle$ as well as the
related Q and G subblocks.  If the CAR relations fully propagate
through the moment matrix, PQG positivity is automatic.

````julia
using COSMO, JuMP

SOLVER = optimizer_with_attributes(COSMO.Optimizer,
    "max_iter" => 50_000,
    "eps_abs" => 1e-7,
    "eps_rel" => 1e-7,
    "verbose" => false)

config = SolverConfig(optimizer = SOLVER, order = 2)
result = cs_nctssos(pop, config)

println("SDP objective:    $(round(result.objective; digits=6))")
println("Exact energy:     $(round(E_exact; digits=6))")
println("Gap:              $(round(abs(result.objective - E_exact); digits=6))")
````

````
SDP objective:    -1.236068
Exact energy:     -1.236068
Gap:              0.0

````

If the gap is $\lesssim 10^{-5}$, the V2RDM construction is correct.
PQG at $N = 2$ is exact, so any residual gap is solver tolerance.

---

## 10. Template: H₄ chain at $N_k = 2$

The full periodic benchmark follows the same pattern, with two
additions: (1) k-point labels on every mode, and (2) ab initio
integrals from PySCF instead of Hubbard parameters.

### Step 1 — Create 32 fermionic modes

Each mode is labeled by (spatial orbital $p$, spin $\sigma$,
$k$-point $k$).  With 8 spatial orbitals, 2 spins, and 2 $k$-points:

```julia
nk = 2
norb = 8  # spatial orbitals per k-point

# One species per (spin, k-point) combination → 4 groups
registry, ((c_up_k0, c_up_k0_dag),
           (c_dn_k0, c_dn_k0_dag),
           (c_up_k1, c_up_k1_dag),
           (c_dn_k1, c_dn_k1_dag)) = create_fermionic_variables([
    ("c_up_k0", 1:norb),
    ("c_dn_k0", 1:norb),
    ("c_up_k1", 1:norb),
    ("c_dn_k1", 1:norb),
])
# Total: 32 spin-orbitals
```

### Step 2 — Load and normalise integrals

The integrals come from PySCF (see `extract_h4_integrals.py`).
Normalise before building the Hamiltonian:

```julia
using NPZ
data = npzread("h4_chain_nk2.npz")
h1e_k0 = ComplexF64.(data["h1e_k0"]) / nk
h1e_k1 = ComplexF64.(data["h1e_k1"]) / nk
# ... similarly for eri blocks, divided by nk²
```

### Step 3 — Build the Hamiltonian

One-body terms (diagonal in spin and $k$):

```julia
ham = zero(typeof(c_up_k0[1]))

# k₀ block, spin-up
for p in 1:norb, q in 1:norb
    h1e_k0[p, q] == 0 && continue
    ham += h1e_k0[p, q] * c_up_k0_dag[p] * c_up_k0[q]
end
# ... repeat for (k₀, ↓), (k₁, ↑), (k₁, ↓)
```

Two-body terms (with momentum conservation and spin selection):

```julia
for (ik1, ik2, ik3, ik4) in allowed_k_combinations
    eri = data["eri_$(ik1)_$(ik2)_$(ik3)_$(ik4)"] / nk^2
    for p in 1:norb, q in 1:norb, r in 1:norb, s in 1:norb
        # physicist convention: ⟨pq|rs⟩ = eri[p, r, q, s]
        v = eri[p, r, q, s]
        abs(v) < 1e-12 && continue
        # Sum over spin-conserving combinations:
        # (σ₁=σ₃, σ₂=σ₄) → 4 terms: ↑↑↑↑, ↑↓↑↓, ↓↑↓↑, ↓↓↓↓
        ham += 0.5 * v * c_dag[p,↑,k1] * c_dag[q,↑,k2] * c[s,↑,k4] * c[r,↑,k3]
        ham += 0.5 * v * c_dag[p,↑,k1] * c_dag[q,↓,k2] * c[s,↓,k4] * c[r,↑,k3]
        # ... etc.
    end
end
```

The factor $1/2$ avoids double-counting when summing over all
$p, q, r, s$ independently.

### Step 4 — Fix particle number and solve

```julia
N_total = 8  # 4 electrons/cell × 2 k-points
N_op = sum of all number operators across all modes

pop = polyopt(ham, registry;
    moment_eq_constraints = [N_op - N_total * one(ham)])

config = SolverConfig(optimizer = SOLVER, order = 2, cs_algo = MF())
result = cs_nctssos(pop, config)
# Expected: E ≈ -2.188 Ha/cell (paper, Nk=2)
```

Setting `cs_algo = MF()` enables correlative sparsity detection,
which automatically discovers the momentum-conservation block
structure from Section 4.

---

## 11. Expected results

Energy per unit cell in Hartree (digitised from paper Fig. 1,
$\pm 0.002$ Ha):

| $N_k$ | HF | MP2 | CCSD | CCSD(T) | V2RDM $[4,8]$ |
|:------|:------|:------|:------|:--------|:--------------|
| 2 | $-2.108$ | $-2.182$ | $-2.203$ | $-2.206$ | $-2.188$ |
| 4 | $-2.137$ | $-2.215$ | $-2.233$ | $-2.236$ | $-2.221$ |
| 6 | $-2.141$ | $-2.221$ | $-2.237$ | $-2.241$ | $-2.225$ |
| 8 | $-2.142$ | $-2.223$ | $-2.238$ | $-2.243$ | $-2.225$ |

V2RDM converges by $N_k = 6$ (energy changes $< 10^{-4}$ Ha beyond
that).  At $N_k = 8$, V2RDM is within $\sim 0.75\%$ of CCSD(T) —
recovering $\sim 82\%$ of the correlation energy.

!!! note "Why is V2RDM above CCSD(T)?"
    Two facts that are easy to confuse:
    (1) CCSD(T) is **not variational** — it can go below the true
    ground-state energy.
    (2) Active-space V2RDM is **not a pure lower bound** — the
    excluded virtual orbitals contribute correlation that the SDP
    cannot capture.
    So $E_\text{CCSD(T)} < E_\text{V2RDM}$ does not contradict
    theory.

---

## Summary

### The code recipe

| Step | What | API / tool |
|:-----|:-----|:-----------|
| 1 | Create fermionic modes for each $(p, \sigma, k)$ | [`create_fermionic_variables`](@ref) with one group per (spin, k-point) |
| 2 | Extract integrals from PySCF | `extract_h4_integrals.py` → `.npz` file |
| 3 | Normalise: $h/N_k$, $V/N_k^2$ | Manual division before Hamiltonian build |
| 4 | Build $\hat{H}$ with spin selection and momentum conservation | Standard Julia arithmetic on fermionic polynomials |
| 5 | Fix particle number $\hat{N} = N_\text{total}$ | `moment_eq_constraints` in [`polyopt`](@ref) |
| 6 | Solve with correlative sparsity for k-blocking | [`cs_nctssos`](@ref) with `cs_algo = MF()` |
| 7 | Compare energy to paper's Table / Fig. 1 | $E \approx -2.188$ Ha/cell at $N_k = 2$ |

### Key concepts

| Concept | One-liner |
|:--------|:----------|
| **V2RDM** | Minimise $E = \text{Tr}({}^{2}\!K\; {}^{2}\!D)$ subject to PQG semidefinite constraints — equivalent to order-2 moment relaxation |
| **2-RDM** (${}^{2}\!D$) | $\binom{M}{2} \times \binom{M}{2}$ matrix encoding all pairwise electron correlations |
| **PQG conditions** | ${}^{2}\!D \succeq 0$, ${}^{2}\!Q \succeq 0$, ${}^{2}\!G \succeq 0$ — necessary N-representability constraints |
| **Momentum conservation** | $k_1 + k_2 = k_3 + k_4 + \mathbf{G}$ — kills half the ERI blocks and block-diagonalises the 2-RDM |
| **Correlative sparsity** | NCTSSoS.jl's `cs_algo` detects k-blocking automatically from the Hamiltonian's support graph |
| **Normalization** | $h \to h/N_k$, $V \to V/N_k^2$ for per-unit-cell energies |
| **Chemist vs. physicist** | $(pr\|qs)_\text{chem} = \langle pq\|rs \rangle_\text{phys}$ — always check which convention your integrals use |

### See also

- [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
  creation/annihilation operators, CAR, and parity constraints
- [Hubbard Model](@ref hubbard-model) — spin degrees of freedom,
  particle-number constraints, `moment_eq_constraints`
- [Kitaev Chain](@ref kitaev-chain) — pairing terms and topological
  phases in fermionic systems

### References

- A. O. Schouten, S. Ewing, and D. A. Mazziotti, "Bootstrapping the
  Electronic Structure of Quantum Materials," arXiv:2504.02861 (2025).
- D. A. Mazziotti, "Contracted Schrödinger equation: Determining
  quantum energies and two-particle density matrices without wave
  functions," *Phys. Rev. A* **57**, 4219 (1998).
- C. Garrod and J. K. Percus, "Reduction of the n-particle variational
  problem," *J. Math. Phys.* **5**, 1756 (1964).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


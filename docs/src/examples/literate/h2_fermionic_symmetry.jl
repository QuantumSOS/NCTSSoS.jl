# # [H₂/STO-3G with Fermionic Symmetry Adaptation](@id h2-fermionic-symmetry)
#
# This page walks through the **full hydrogen pipeline** — from
# molecular-orbital integrals to a certified electronic ground-state
# energy — using the **fermionic + symmetry** path of `NCTSSoS.jl`.
#
# H₂ is the smallest non-trivial molecule.  In the minimal STO-3G basis
# it has just **two spatial orbitals** ($\sigma_g$, $\sigma_u$) and
# **four spin orbitals** ($\sigma_g\!\uparrow$, $\sigma_g\!\downarrow$,
# $\sigma_u\!\uparrow$, $\sigma_u\!\downarrow$).  Two electrons live on
# this stage.  That is small enough to:
#
# 1. write down the *exact* electronic Hamiltonian as a polynomial in
#    fermionic creation/annihilation operators;
# 2. solve it exactly by full diagonalisation, restricted to the physical
#    $N = 2$ sector, as a verification oracle;
# 3. solve the same problem as an SDP relaxation, **once dense** and
#    **once symmetry-adapted**, and recover the same energy with much
#    smaller PSD blocks.
#
# The point of using H₂ here is **not** chemistry.  The point is that
# this is the smallest setting where every piece of the fermionic
# symmetry machinery — sector splitting (parity / particle number /
# $S_z$ / Abelian orbital irrep) and SU(2) spin adaptation — does
# something visible and easy to interpret.
#
# Prerequisites:
#
# - Comfort with creation/annihilation operators and the fermionic CAR;
#   see [Fermionic Ground State (XY Model)](@ref fermionic-ground-state)
#   for a from-scratch build.
# - The [Hubbard Model](@ref hubbard-model) example, which introduces
#   `moment_eq_constraints` for fixing particle number.
# - The [CHSH with Symmetry Reduction](@ref chsh-symmetry) example and
#   the manual page [Symmetry-Adapted Basis](@ref symmetry-adapted-basis),
#   which describe the design of NCTSSoS's symmetry path.
#
# Throughout this page, the **reference energy is electronic only** —
# we never add the nuclear-repulsion constant.  Tacking that constant on
# at the end gives the usual textbook total energy, but it is a
# scalar that does not affect any of the operator-algebra discussion.

# ---
#
# ## Why H₂ deserves a fermionic + symmetry treatment
#
# Even at this size, H₂ already has *all* of the structure that makes
# molecular electronic structure hard:
#
# | Structure | What it gives us |
# |:----------|:------------------|
# | **Fermionic CAR** | electrons obey Pauli exclusion; the algebra is non-commutative |
# | **Particle-number conservation** | only the $N = 2$ sector is physical |
# | **Total $S_z$ conservation** | the Hamiltonian is spin-rotation invariant |
# | **Total $S^2$ conservation (SU(2))** | physical eigenstates are singlets ($S = 0$) or triplets ($S = 1$) |
# | **Inversion symmetry $\hat i$** | reflecting through the bond midpoint leaves $H$ unchanged, so every operator carries a parity $\pm 1$ under $\hat i$: *gerade* ($\sigma_g$, label `:Ag`) or *ungerade* ($\sigma_u$, label `:B1u`) |
#
# That last row is the one that usually reads as jargon, so here is the
# one-paragraph unpacking. H₂ has a geometric symmetry: put the origin
# at the midpoint of the H–H bond and send every point
# $\mathbf r \mapsto -\mathbf r$. The two protons swap, but because the
# nuclei are identical the electronic Hamiltonian is unchanged — i.e.
# $[H, \hat i] = 0$. Since $\hat i^2 = 1$, its eigenvalues are $\pm 1$,
# so every operator has a definite *parity* under $\hat i$: even
# (*gerade*, `g`, label `:Ag`) or odd (*ungerade*, `u`, label `:B1u`).
# The two spatial orbitals $\sigma_g$ and $\sigma_u$ are already
# eigenorbitals of $\hat i$, which is why the integral matrix $h$
# above is diagonal.
#
# **Why that matters for the SDP.** A symmetry of $H$ is a free
# block-diagonalisation. Moments of the form
# $\langle \text{even}\cdot\text{odd} \rangle$ must vanish, so the one
# big moment matrix splits as
# $M = M_\text{even} \oplus M_\text{odd}$, and the single big PSD
# constraint becomes two smaller ones. In NCTSSoS this is the
# `split_irrep = true` piece of [`FermionicSectorSpec`](@ref) further
# down this page.
#
# **Where the names come from.** H₂'s full point group is the
# continuous $D_{\infty h}$; quantum-chemistry codes can't handle the
# continuous part, so they keep the largest *finite Abelian* subgroup
# $D_{2h}$ (8 one-dimensional irreps) and label orbitals with its
# Mulliken names. For H₂/STO-3G only two of those names appear in
# practice — `:Ag` for the gerade orbital and `:B1u` for the
# ungerade one — and only the $\hat i$-parity bit of $D_{2h}$ is
# actually doing work here. The rest is inherited naming.
#
# A *dense* moment relaxation ignores all of that: it builds one large PSD
# block over every monomial up to some order.  The fermionic **sector
# split** discards monomials whose quantum numbers cannot survive in the
# physical sector, and the **spin adaptation** further block-diagonalises
# the surviving piece into total-$S^2$ blocks.
#
# Concretely, the order-2 dense moment matrix on these four spin orbitals
# is a single PSD constraint of side close to $50$.  The symmetry-adapted
# relaxation in this page returns the *same* electronic ground-state
# energy, but the largest reduced PSD block has side $3$.

# ---
#
# ## Exact diagonalisation oracle (reviewed integrals)
#
# We need an exact answer to compare against.  Build the four-mode
# Hamiltonian via the **Jordan–Wigner** mapping (fermion operators →
# tensor products of Pauli matrices) and take the smallest eigenvalue
# *restricted to* the physically relevant $N = 2$ sector.

using LinearAlgebra

## Single-mode Jordan–Wigner pieces (occupation-basis convention).
const PAULI_Z   = ComplexF64[1 0; 0 -1]
const PAULI_I   = Matrix{ComplexF64}(I, 2, 2)
const JW_RAISE  = ComplexF64[0 1; 0 0]   ## σ⁺ = |0⟩⟨1|, the annihilator on |1⟩
const JW_LOWER  = ComplexF64[0 0; 1 0]   ## σ⁻ = |1⟩⟨0|, the creator on |0⟩

"""
    jw_op(mode, kind, nmodes)

Jordan–Wigner image of a fermion operator on `nmodes` sites.  `kind` is
`:annihilate` or `:create`.
"""
function jw_op(mode::Int, kind::Symbol, nmodes::Int)
    mats = [site < mode ? PAULI_Z :
            site == mode ? (kind === :annihilate ? JW_RAISE : JW_LOWER) :
            PAULI_I for site in 1:nmodes]
    return reduce(kron, mats)
end

# ### Where these integrals come from
#
# The numbers below are the **standard minimal-basis H₂/STO-3G**
# molecular-orbital integrals at the equilibrium bond length
# $R = 1.4\,\text{a.u.}\;(\approx 0.7414\,\text{Å})$, in chemists'
# notation:
#
# - $h_{pq} = \int \phi_p(\mathbf r)\,\bigl(-\tfrac12\nabla^2 + V_\text{ne}(\mathbf r)\bigr)\,\phi_q(\mathbf r)\,d\mathbf r$ — one-body MO matrix.
# - $v_{pqrs} = (pq|rs) = \iint \phi_p(\mathbf r_1)\phi_q(\mathbf r_1)\,\frac{1}{r_{12}}\,\phi_r(\mathbf r_2)\phi_s(\mathbf r_2)\,d\mathbf r_1 d\mathbf r_2$ — two-body MO tensor.
#
# **How they were generated.** They were produced once with
# [PySCF](https://pyscf.org/) (v2.x), printed at machine precision,
# reviewed, and then frozen so this page stays Python-free. The
# reproducer is short:
#
# ```python
# from pyscf import gto, scf, ao2mo
# import numpy as np
#
# mol = gto.M(atom="H 0 0 0; H 0 0 1.4", unit="Bohr",
#             basis="sto-3g", verbose=0)
# mf  = scf.RHF(mol).run()
# C    = mf.mo_coeff                                      # AO -> MO transform
# h_mo = C.T @ mf.get_hcore() @ C                         # 1-body, MO basis
# v_mo = ao2mo.restore(1, ao2mo.kernel(mol, C), C.shape[1])  # (pq|rs), chemists'
# np.set_printoptions(precision=12, suppress=True)
# print("h =", h_mo); print("v =", v_mo)
# ```
#
# **How to verify them without running PySCF.** The same minimal-basis
# H₂ example is worked out by hand in Szabo & Ostlund, *Modern Quantum
# Chemistry* (Dover, 1996), §3.5 — every constant below matches that
# table to its 4-figure printed precision; PySCF just gives more
# digits. The committed test fixture
# `test/data/expectations/fermionic_symmetry.toml` (case
# `h2_sto3g_n2_spin_adapted`) and the regression test
# `test/problems/fermionic/fermionic_symmetry.jl` use *these exact
# numbers*.

function h2_sto3g_integrals()
    h = [
        -1.252477303982   0.0;
         0.0            -0.475934275367;
    ]
    v = zeros(2, 2, 2, 2)
    v[1, 1, 1, 1] = 0.674493166046
    v[2, 2, 2, 2] = 0.69739794957
    v[1, 2, 1, 2] = 0.663472101055
    v[2, 1, 2, 1] = 0.663472101055
    v[1, 2, 2, 1] = 0.181287518306
    v[2, 1, 1, 2] = 0.181287518306
    v[1, 1, 2, 2] = 0.181287518306
    v[2, 2, 1, 1] = 0.181287518306
    return h, v
end

# Mode layout for the 4-mode Jordan–Wigner construction:
# modes `1, 2` are $\sigma_g\!\uparrow, \sigma_u\!\uparrow$,
# modes `3, 4` are $\sigma_g\!\downarrow, \sigma_u\!\downarrow$.

function h2_exact_n2_energy()
    h, v = h2_sto3g_integrals()
    nmodes = 4
    a   = [jw_op(i, :annihilate, nmodes) for i in 1:nmodes]
    a_d = [jw_op(i, :create,     nmodes) for i in 1:nmodes]

    up = (1, 2)   ## (σ_g↑, σ_u↑)
    dn = (3, 4)   ## (σ_g↓, σ_u↓)
    dim = 2^nmodes
    H = zeros(ComplexF64, dim, dim)

    for p in 1:2, q in 1:2
        H .+= h[p, q] * (a_d[up[p]] * a[up[q]] + a_d[dn[p]] * a[dn[q]])
    end
    for p in 1:2, q in 1:2, r in 1:2, s in 1:2
        coef = 0.5 * v[p, q, r, s]
        coef == 0 && continue
        for (cd_σ, c_σ) in ((a_d[up[p]], a[up[r]]), (a_d[dn[p]], a[dn[r]])),
            (cd_τ, c_τ) in ((a_d[up[q]], a[up[s]]), (a_d[dn[q]], a[dn[s]]))
            ## chemists' notation: ½ Σ v_{pqrs} c†_{pσ} c†_{qτ} c_{sτ} c_{rσ}
            ## with the up/down assignment matching (cd_σ, c_σ) / (cd_τ, c_τ).
            H .+= coef * cd_σ * cd_τ * c_τ * c_σ
        end
    end

    ## Restrict to the N = 2 occupation sector.
    keep = [idx for idx in 1:dim if count_ones(UInt(idx - 1)) == 2]
    return real(eigmin(Hermitian(H[keep, keep])))
end

E_exact = h2_exact_n2_energy()
@assert abs(E_exact - (-1.8510462258674936)) < 1e-10   #src
println("Exact electronic N = 2 ground state:  E₀ = $E_exact")

# This is the **electronic-only** ground-state energy of H₂/STO-3G with
# two electrons.  Add the nuclear-repulsion constant separately if you
# want the textbook total energy; the polynomial we feed the SDP does
# not contain it, and neither does the value above.
#
# ### Sanity-check against the literature
#
# Two independent ways to convince yourself the number above is right:
#
# 1. **Comparison to a tabulated total energy.** In a minimal basis,
#    full configuration interaction (FCI) is exact, and the FCI/STO-3G
#    total energy of H₂ at $R = 1.4\,\text{a.u.}$ is the textbook value
#    $E_\text{total} \approx -1.1368\,\text{hartree}$ (e.g. Szabo &
#    Ostlund, *Modern Quantum Chemistry*, §3.5 / §4.3, and reproduced
#    by every standard QC code). Adding the nuclear-repulsion constant
#    $1/R = 1/1.4 = 0.71428\!\dots$ to our electronic value gives
#
#    ```
#    E_total = E_exact + 1/R
#            = -1.85104622586… + 0.71428571428…
#            = -1.13676051158… hartree
#    ```
#
#    matching the tabulated FCI/STO-3G value to all printed digits.
#
# 2. **Internal Jordan–Wigner cross-check.** The `E_exact` above *is*
#    the Jordan–Wigner exact diagonalisation of the same fermionic
#    polynomial we hand to the SDP (restricted to the $N = 2$ sector),
#    so the FCI comparison is what `h2_exact_n2_energy()` already does.
#    The reviewed fixture
#    `test/data/expectations/fermionic_symmetry.toml`
#    (case `h2_sto3g_n2_spin_adapted`) records the same value
#    $-1.8510462258674936$ to machine precision, and CI checks both
#    the exact diagonalisation and the SDP against it.

# ---
#
# ## Building H₂ as a fermionic polynomial
#
# We now build the *same* Hamiltonian via the polynomial interface.
# [`create_fermionic_variables`](@ref) gives us creation and annihilation
# operators that already obey every CAR — including cross-species
# (different spin) anticommutation.  Two species, two spatial orbitals:

using NCTSSoS

registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
    ("c_up", 1:2),   ## (σ_g↑, σ_u↑)
    ("c_dn", 1:2),   ## (σ_g↓, σ_u↓)
]);
length(registry)

# The four spin-orbital modes, in the order NCTSSoS assigned them:

up_modes = [Int(op.word[1]) for op in c_up]
dn_modes = [Int(op.word[1]) for op in c_dn]
(up_modes, dn_modes)

# !!! note "What are `op.word[1]` indices?"
#     Every fermionic operator returned by
#     [`create_fermionic_variables`](@ref) is a one-letter
#     `NormalMonomial`; `op.word[1]` is the underlying registry index of
#     that letter.  We use these indices below to tell the symmetry layer
#     which physical mode is spin-up vs. spin-down, and which orbital is
#     gerade vs. ungerade.

# The electronic Hamiltonian is built term by term from the integrals.

function h2_hamiltonian(c_up, c_up_dag, c_dn, c_dn_dag)
    h, v = h2_sto3g_integrals()
    H = zero(1.0 * one(c_up[1]))

    for p in 1:2, q in 1:2
        H += h[p, q] * (c_up_dag[p] * c_up[q] + c_dn_dag[p] * c_dn[q])
    end

    for p in 1:2, q in 1:2, r in 1:2, s in 1:2
        coef = 0.5 * v[p, q, r, s]
        coef == 0 && continue
        for (cσ, cσ_dag) in ((c_up, c_up_dag), (c_dn, c_dn_dag)),
            (cτ, cτ_dag) in ((c_up, c_up_dag), (c_dn, c_dn_dag))
            H += coef * (cσ_dag[p] * cτ_dag[q] * cτ[s] * cσ[r])
        end
    end
    return H
end

ham = h2_hamiltonian(c_up, c_up_dag, c_dn, c_dn_dag);

# ---
#
# ## The fermionic mode layout
#
# Sector splitting and spin adaptation both need to know **which
# physical mode is what** — its orbital, its spin label, and which
# Abelian irrep its orbital carries. Before we fill in that table,
# here is what each of those words means in this concrete H₂ setting.
#
# ### Background: orbitals, spin labels, and orbital irreps
#
# - **Spatial orbital.** A one-electron wavefunction in 3D space —
#   the *spatial* part of where an electron lives. In the minimal
#   STO-3G basis for H₂ there are exactly two spatial orbitals:
#   $\sigma_g$ (the *bonding* combination of the two atomic 1s
#   orbitals, with extra amplitude between the nuclei) and $\sigma_u$
#   (the *antibonding* combination, with a node between the nuclei).
#
# - **Spin orbital / fermionic mode.** A spatial orbital times a spin
#   component $\uparrow$ (spin up) or $\downarrow$ (spin down).
#   H₂/STO-3G therefore has $2 \times 2 = 4$ spin orbitals, and each
#   spin orbital is one *fermionic mode* that an electron can either
#   occupy or not. Those four modes are exactly the four operators
#   `c_up[1]`, `c_up[2]`, `c_dn[1]`, `c_dn[2]` we built above.
#
# - **The number `2`$S_z$.** $S_z$ is the $z$-component of an
#   electron's intrinsic spin and only takes two values: $+\tfrac12$
#   for spin up, $-\tfrac12$ for spin down. We store **twice** that
#   value, $2 S_z \in \{+1, -1\}$, so dictionary keys stay integer
#   instead of half-integer. There is no physics in the factor of two
#   — it is purely a bookkeeping convention.
#
# - **Inversion symmetry and the labels `:Ag` / `:B1u`.** H₂ has an
#   *inversion symmetry* $i$: reflect every coordinate through the
#   midpoint of the bond. Acting on a spatial orbital, $i$ either
#   leaves it alone or flips its sign:
#
#     - $i\,\sigma_g = +\sigma_g$ — *gerade* (German for "even");
#     - $i\,\sigma_u = -\sigma_u$ — *ungerade* ("odd").
#
#   That single $\pm 1$ bit is the only orbital-symmetry information
#   we are using here. It is what physicists call an *irreducible
#   representation* (irrep) of the inversion group $\{1, i\} \cong
#   \mathbb{Z}_2$: a label that records how the orbital transforms
#   under the symmetry. The names `:Ag` and `:B1u` are inherited from
#   the larger $D_{2h}$ point group that quantum-chemistry codes use
#   for linear molecules; for our purposes they are just shorthand for
#
#     - `:Ag`  ↔ even under inversion (no sign flip), and
#     - `:B1u` ↔ odd under inversion (sign flip).
#
#   The label belongs to the **spatial** orbital, not to the spin —
#   $\sigma_g\!\uparrow$ and $\sigma_g\!\downarrow$ both inherit
#   `:Ag`, because spin-flip and spatial inversion are unrelated.
#
# - **Multiplying irreps.** When you multiply two operators, their
#   orbital irreps multiply. With only two labels and the inversion
#   rule "two evens or two odds give even, one of each gives odd",
#   the full multiplication table is:
#
#   | $\times$ | `:Ag`  | `:B1u` |
#   |:--------:|:------:|:------:|
#   | `:Ag`    | `:Ag`  | `:B1u` |
#   | `:B1u`   | `:B1u` | `:Ag`  |
#
#   This is just the $\mathbb{Z}_2$ group table with `:Ag` as the
#   identity. NCTSSoS reads it from the `multiply` callback we hand
#   to [`AbelianIrrepTable`](@ref) below.
#
# ### The H₂/STO-3G layout
#
# With those terms defined, the four-mode layout reads:
#
# | Mode          | Spatial orbital | Spin           | $2 S_z$ | Orbital irrep |
# |:--------------|:----------------|:---------------|:-------:|:--------------|
# | `up_modes[1]` | $\sigma_g$      | $\uparrow$     | `+1`    | `:Ag`         |
# | `up_modes[2]` | $\sigma_u$      | $\uparrow$     | `+1`    | `:B1u`        |
# | `dn_modes[1]` | $\sigma_g$      | $\downarrow$   | `-1`    | `:Ag`         |
# | `dn_modes[2]` | $\sigma_u$      | $\downarrow$   | `-1`    | `:B1u`        |
#
# We pass that information to [`FermionicModeLayout`](@ref): the
# orbital each mode belongs to, the $2 S_z$ of each mode, the irrep
# of each *orbital*, and the irrep multiplication table.

h2_irrep_multiply(a, b) = a === b ? :Ag : :B1u
nothing #hide

layout = FermionicModeLayout(
    Dict(vcat(
        [up_modes[i] => i for i in eachindex(up_modes)],
        [dn_modes[i] => i for i in eachindex(dn_modes)],
    ));
    spin2_of = Dict(vcat(
        [up_modes[i] => 1  for i in eachindex(up_modes)],
        [dn_modes[i] => -1 for i in eachindex(dn_modes)],
    )),
    irrep_of    = Dict(1 => :Ag, 2 => :B1u),
    irrep_table = AbelianIrrepTable(:Ag; multiply = h2_irrep_multiply, dual = identity),
)

# One API note worth flagging: the **first positional argument** to
# [`FermionicModeLayout`](@ref) is the `orbital_of` map
# `mode => orbital_index`. That is what tells the symmetry layer which
# *spin* orbital (mode) belongs to which *spatial* orbital, so that
# both spin copies of $\sigma_g$ correctly inherit the `:Ag` label and
# both copies of $\sigma_u$ inherit `:B1u`. The `spin2_of` and
# `irrep_of` keyword arguments then attach the spin and the irrep on
# top.

# ---
#
# ## Sector split + SU(2) spin adaptation
#
# A [`FermionicSectorSpec`](@ref) declares which conserved quantum
# numbers should be used to **partition** monomials into independent
# sectors.  We turn on all four for H₂:

sector = FermionicSectorSpec(
    mode_layout    = layout,
    split_parity   = true,    ## fermionic (-1)^N parity, always exact
    split_number   = true,    ## total electron count N
    split_spin     = true,    ## total 2 S_z
    split_irrep    = true,    ## Abelian orbital-irrep label (here Ag/B1u)
)

# A [`FermionicSpinAdaptationSpec`](@ref) layers on top of that: within
# each sector that has well-defined $S_z$, it diagonalises the SU(2)
# Casimir $\hat S^2$ on the basis to identify singlets ($S = 0$,
# `total_spin2 = 0`), triplets ($S = 1$, `total_spin2 = 2`), etc.

spin = FermionicSpinAdaptationSpec(mode_layout = layout)

# The two are wrapped in a single [`SymmetrySpec`](@ref).  No raw
# [`SignedPermutation`](@ref) generators are needed for this example —
# all the structure comes from the conserved quantum numbers and from
# the Casimir, not from a finite group action on operators.

symspec = SymmetrySpec(sector = sector, spin_adaptation = spin)

# ---
#
# ## Fixing $N = 2$ as a state constraint
#
# H₂ is a two-electron problem.  Any solver should look only at the
# $N = 2$ sector of Fock space.  We do that with `moment_eq_constraints`,
# which imposes $g\,|\psi\rangle = 0$ on the **target state** rather than
# as an operator identity (cf. the warning in [Hubbard
# Model](@ref hubbard-model)).

n_total = sum(c_up_dag[i] * c_up[i] + c_dn_dag[i] * c_dn[i] for i in 1:2)

pop = polyopt(ham, registry; moment_eq_constraints = [n_total - 2.0 * one(ham)]);

# ---
#
# ## Solver setup

# !!! note "SDP Solver"
#     These examples use [Mosek](https://www.mosek.com/) via `MosekTools`.
#     Any SDP-capable solver works: replace `Mosek.Optimizer` with
#     `COSMO.Optimizer` or `Clarabel.Optimizer` for open-source
#     alternatives.

using MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG"         => 0,
    "MSK_IPAR_NUM_THREADS" => 0)
nothing #hide

# ---
#
# ## Run 1 — Dense order-2 baseline
#
# No symmetry, no sparsity: just the full order-2 dense moment relaxation.

dense_config = SolverConfig(
    optimizer = SOLVER,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
)

dense_result = cs_nctssos(pop, dense_config);

println("Dense  : objective = $(dense_result.objective)")
println("Exact  : objective = $E_exact")
println("Gap    : ", abs(dense_result.objective - E_exact))
@assert abs(dense_result.objective - E_exact) < 5e-6   #src

# A single PSD block; its size is the dense order-2 baseline we will
# shrink:

dense_block_sizes = only(dense_result.moment_matrix_sizes)

# (Here `only(...)` extracts the per-clique vector for the single
# correlative-sparsity clique that the relaxation has by default.)

dense_result.symmetry === nothing

# ---
#
# ## Run 2 — Symmetry-adapted relaxation
#
# The only change is the `symmetry` keyword.  Same Hamiltonian, same
# constraint, same order, same solver.

sym_config = SolverConfig(
    optimizer = SOLVER,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = symspec,
)

sym_result = cs_nctssos(pop, sym_config);

println("Sym    : objective = $(sym_result.objective)")
println("Dense  : objective = $(dense_result.objective)")
println("|Δ|    : ", abs(sym_result.objective - dense_result.objective))
@assert abs(sym_result.objective - E_exact) < 5e-6           #src
@assert abs(sym_result.objective - dense_result.objective) < 5e-6  #src

# Same answer — to solver tolerance.  The reduction is **exact**, not an
# additional relaxation: every isotypic component is a separate PSD
# block, and the original PSD constraint is the direct sum of those
# blocks.

# ---
#
# ## Inspect the [`SymmetryReport`](@ref)
#
# `result.symmetry` is no longer `nothing`; it is a small struct
# summarising what the symmetry path did.

report = sym_result.symmetry

# The pre-symmetry half-basis (i.e. the basis the dense path would have
# used for each side of its single moment block):

(report.basis_half_size, report.basis_full_size)

# The single dense PSD block of side $\approx 50$ has been replaced by
# many small blocks, the largest of side $3$:

report.psd_block_sizes

# Comparison with the pre-symmetry single-block layout:

(maximum(report.psd_block_sizes), maximum(dense_block_sizes))

# Each block carries a **label** describing where it came from.  Most
# come from the spin-adaptation step and are
# [`FermionicSpinBlockLabel`](@ref); the rest are
# [`FermionicSectorLabel`](@ref) for sectors with no spin structure
# (e.g. the trivial / parity-only blocks).

unique(typeof.(report.block_labels))

# A quick distribution of block provenance — which step in the pipeline
# produced each block:

unique(report.block_provenance)

# ---
#
# ## What the labels mean for H₂
#
# The two structural facts you should be able to read off the report:
#
# **(a) Orbital irreps `:Ag` vs. `:B1u`.**
# Densities like $a^\dagger_{\sigma_g \uparrow}\,a_{\sigma_g \uparrow}$
# or $a^\dagger_{\sigma_u \downarrow}\,a_{\sigma_u \downarrow}$ are
# *irrep-neutral* — they multiply to `:Ag`.  Transitions like
# $a^\dagger_{\sigma_g \uparrow}\,a_{\sigma_u \uparrow}$ flip gerade $↔$
# ungerade and carry the lone odd irrep `:B1u`.  The set of irrep labels
# that show up among the spin-adapted blocks is exactly $\{ \mathrm{Ag},
# \mathrm{B1u} \}$:

irreps_present = Set(label.sector.irrep
    for label in report.block_labels
    if label isa FermionicSpinBlockLabel)

@assert irreps_present == Set([:Ag, :B1u])   #src
irreps_present

# **(b) Total spin $S^2$ values.**
# Within each sector that has a well-defined $S_z$, the Casimir-based
# spin adaptation splits the basis into singlet ($S = 0$,
# `total_spin2 = 0`) and triplet ($S = 1$, `total_spin2 = 2`) blocks.
# Both show up:

spins_present = Set(label.total_spin2
    for label in report.block_labels
    if label isa FermionicSpinBlockLabel)

@assert 0 in spins_present && 2 in spins_present   #src
spins_present

# !!! note "`total_spin2` is `2S`, twice the spin quantum number"
#     The convention is the same as for `spin2_of` in the layout: we
#     store **twice** the conserved spin quantum number, so everything
#     stays integer.  So `total_spin2 = 0` is $S = 0$ (singlet),
#     `total_spin2 = 2` is $S = 1$ (triplet), `total_spin2 = 1` would
#     be a doublet, and so on.  The associated Casimir eigenvalue is
#     $S(S+1) =$ `total_spin2 * (total_spin2 + 2) / 4`.
#
# That the H₂ ground state is a singlet is well known.  The relaxation
# does not need that knowledge — both singlet and triplet blocks live
# in the SDP, and the $S^2$ adaptation is what keeps them from being
# tangled together inside one block.

# ---
#
# ## Summary
#
# ### Pipeline at a glance
#
# | Step | What we did | API |
# |:-----|:------------|:----|
# | 1 | Build the exact electronic Hamiltonian as a polynomial in $a$, $a^\dagger$ from reviewed STO-3G integrals | [`create_fermionic_variables`](@ref) + standard Julia arithmetic |
# | 2 | Tag each spin-orbital mode with its spin and orbital irrep (Ag/B1u) | [`FermionicModeLayout`](@ref), [`AbelianIrrepTable`](@ref) |
# | 3 | Declare which conserved quantum numbers should split the moment basis into sectors | [`FermionicSectorSpec`](@ref) |
# | 4 | Layer Casimir-based SU(2) spin adaptation on top | [`FermionicSpinAdaptationSpec`](@ref) |
# | 5 | Wrap (3) + (4) in a single symmetry spec | [`SymmetrySpec`](@ref) |
# | 6 | Fix the physical sector $N = 2$ as a state constraint | `moment_eq_constraints` of [`polyopt`](@ref) |
# | 7 | Solve dense and symmetry-adapted relaxations | [`cs_nctssos`](@ref) with and without `symmetry = symspec` |
# | 8 | Read the symmetry report; check Ag/B1u and singlet/triplet labels | `result.symmetry` ([`SymmetryReport`](@ref)) |
#
# ### Numerical headline
#
# | Quantity | Dense order 2 | Symmetry-adapted order 2 |
# |:---------|:--------------|:--------------------------|
# | objective (electronic only)   | matches $-1.85104622\!\dots$ | matches the same value |
# | largest PSD block side        | dense (single big block)     | $3$ |
# | block labels available?       | no                            | yes — Ag/B1u + singlet/triplet |
# | output type                   | [`PolyOptResult`](@ref NCTSSoS.PolyOptResult) | same, plus `result.symmetry` populated |
#
# ### Caveats
#
# - The energy reported here is **electronic only**.  If you want the
#   physical total energy, add the constant nuclear-repulsion term
#   yourself; it does not interact with the operator algebra and so
#   it would only shift the SDP objective by a constant.
# - This page exercises the *currently supported* fermionic-symmetry
#   path: sector splitting and Casimir-based spin adaptation, with no
#   raw [`SignedPermutation`](@ref) generators in the same spec.  See
#   [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for the
#   precise scope and fail-fast guards.
# - The integrals used are reviewed offline molecular-orbital
#   constants for STO-3G H₂; they live alongside the regression test
#   `test/problems/fermionic/fermionic_symmetry.jl` and the reviewed
#   expectation fixture `test/data/expectations/fermionic_symmetry.toml`
#   (case `h2_sto3g_n2_spin_adapted`).  Swapping in your own integrals
#   for a different molecule means changing only `h2_sto3g_integrals`
#   and the orbital → irrep assignment.
#
# ### See also
#
# - [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
#   from-scratch build of the fermionic interface.
# - [Hubbard Model](@ref hubbard-model) — particle-number constraints
#   via `moment_eq_constraints`.
# - [CHSH with Symmetry Reduction](@ref chsh-symmetry) — the original
#   walk-through of the symmetry path, on a non-fermionic problem.
# - [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) — design and
#   supported scope of the symmetry path.
# - The companion regression suite
#   `test/problems/fermionic/fermionic_symmetry.jl` (case
#   `h2_sto3g_n2_spin_adapted`) keeps the numbers in this page honest.

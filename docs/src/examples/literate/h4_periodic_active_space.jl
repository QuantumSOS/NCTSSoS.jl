# # [H₄ Periodic Active-Space Workflow](@id h4-periodic-active-space)
#
# This page **inspects** the vendored data for a periodic hydrogen chain and
# explains what it contains — without running an SDP solve.  It reconstructs
# the Hartree–Fock bookkeeping, builds a derived tensor, and shows why the
# current solver API is not yet the right interface for the full periodic
# V2RDM benchmark.
#
# If you want to understand the background concepts (what H₄ is, what
# "periodic" means, how SDP helps, what a moment relaxation is), start with
# the companion page:
# [H₄ Chain — k=0 Moment Relaxation](@ref h4-chain-energy-benchmark).
#
# ## What this page covers
#
# 1. Load the reviewed Nk=2 assets shipped with the repository.
# 2. Inspect the **k-blocked one- and two-body integrals** (the numbers
#    that encode how electrons interact).
# 3. Reconstruct the **Hartree–Fock energy** and the **constant shift**
#    (see below).
# 4. Rebuild the derived $K_2$ tensor.
# 5. Quantify why a naive 32-mode order-2 lift becomes a dense relaxation
#    in disguise.
#
# What we **do not** do here is claim that `NCTSSoS.jl` already reproduces the
# full periodic V2RDM solve from the paper.  The k-blocked / spin-adapted
# formulation behind that result is still future work for the package.
#
# !!! warning "Current scope"
#     Exact on this page: Nk=2 asset parsing, HF reconstruction, additive
#     constant recovery, and ``K_2`` construction.  Context only:
#     figure-digitized HF / MP2 / CCSD / CCSD(T) / V2RDM values.  Not yet in
#     the API: a faithful periodic V2RDM SDP solve on the original k-blocked
#     formulation.
#
# If you want a **runnable reduced solve** on the same assets, see
# [H₄ Chain — k=0 Moment Relaxation](@ref h4-chain-energy-benchmark).  That
# companion page keeps only the k=0 intra-sector Hamiltonian and is explicit
# about what it does — and does not — certify.
#
# If you want the **first honest full-Nk=2 attempt** with the current solver
# interface — keep the full Bloch-basis Hamiltonian, turn on only correlative
# sparsity, and see where Mosek runs out of road — see
# [H₄ Periodic Chain — Correlative-Sparsity Attempt](@ref h4-periodic-correlative-sparsity).
#
# If you then want a **runnable full-Nk=2 compromise solve** that keeps the
# whole periodic Hamiltonian but uses a hand-built incomplete order-2 basis,
# see [H₄ Periodic Chain — Incomplete Momentum-Basis Relaxation](@ref h4-periodic-incomplete-momentum-basis).
#
# ## Key concepts for this page
#
# ### What is Hartree–Fock (HF)?
#
# **Hartree–Fock** is the simplest approximation to the electronic structure
# problem.  It assumes each electron moves independently in an average field
# created by all the others — like assuming each party guest reacts only to
# the average mood, not to specific conversations.  HF gives a baseline
# energy, but misses **correlation energy** — the energy lowering that comes
# from electrons actually avoiding each other.
#
# ### What is an active space?
#
# Real materials have many electrons, most of which sit in low-energy "core"
# orbitals that barely participate in interesting physics.  The **active-space
# approximation** freezes core electrons at the HF level and runs the
# expensive calculation only on the "active" electrons and orbitals near the
# energy frontier.  The notation **[4, 8]** means 4 active electrons in 8
# active orbitals per unit cell.
#
# The energy of the frozen core is accounted for by a **constant shift**
# (also called the "additive constant") that gets added back after
# optimization.
#
# ### What are "integrals" in this context?
#
# The Hamiltonian (the operator encoding all physics) is built from numbers
# called **integrals** — computed from the overlap of electron orbitals:
#
# - **One-body integrals** ($h_{pq}$): encode kinetic energy and
#   electron–nucleus attraction.  One matrix per k-point.
# - **Two-body integrals** (ERI = electron repulsion integrals,
#   $V_{pqrs}$): encode electron–electron Coulomb repulsion.  Stored as
#   blocks labeled by four k-point indices $(k_1, k_2, k_3, k_4)$.
#
# Crystal momentum conservation means most k-combinations are zero —
# only blocks where $k_1 + k_2 \equiv k_3 + k_4 \pmod{N_k}$ survive.
#
# !!! info "What is crystal momentum, physically?"
#     The momentum label $k$ describes the **phase pattern between unit
#     cells** — not a velocity.  Think of a standing wave on a guitar
#     string: $k=0$ means every cell has the same sign (all in sync);
#     $k=\pi/a$ means the sign alternates (plus, minus, plus, minus).
#     Momentum conservation then says: the phase patterns of a pair of
#     interacting electrons must be consistent before and after the
#     interaction.  Combinations that violate this give zero integrals,
#     which is why most ERI blocks vanish.

# ## Setup
#
# The raw data lives under `test/data/assets/`.  To keep the docs example and
# the regression test in lockstep, we reuse the same small helper module in both
# places.

using NCTSSoS
using LinearAlgebra, Printf

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: H4_INTEGRALS_PATH,
                         build_k2,
                         complete_graph_edge_count,
                         compound_spatial_edge_count,
                         fermionic_order2_basis_size,
                         fermionic_order2_nuniq,
                         hf_energy,
                         load_nk2_asset

count_nonzero_entries(block; atol = 1e-10) = count(x -> abs(x) > atol, block)
nothing #hide

# ## Physical setup and vendored assets
#
# The reviewed H₄ asset corresponds to a periodic hydrogen chain in **k-space**
# (see the [companion page](@ref h4-chain-energy-benchmark) for what "k-space"
# means):
#
# - **4 H atoms per unit cell** — the repeating unit of the crystal,
# - bond distance **1.0 Å** — the spacing between adjacent atoms,
# - lattice constant **4.0 Å** — the length of one full repeating unit,
# - active space **[4,8] per unit cell** — 4 electrons, 8 orbitals,
# - **Nk = 2** k-points — two discrete "frequency channels" ($k=0$ and $k=1$).
#
# For this Nk=2 asset that means:
#
# - **4 active electrons per cell**,
# - **8 active spatial orbitals per k-point** — each orbital is a function
#   $\phi(x,y,z) \to \mathbb{C}$ whose $|\phi|^2$ gives the probability
#   of finding the electron at that point.  The 8 per k-point come from
#   the 4 hydrogen atoms each contributing basis functions that combine
#   into crystal orbitals of different shapes and energies.  The k-point
#   labels the inter-cell phase pattern; the orbital index labels which
#   shape/energy level within that phase pattern —
# - **16 spatial orbitals total** (8 per k-point × 2 k-points),
# - **32 spin-orbital modes** when we account for spin (each spatial orbital
#   splits into spin-up and spin-down).

asset = load_nk2_asset()
(; reference, preflight, blocker, figure, nk, n_active_orb, n_active_elec,
   total_active_electrons, total_spatial_orbitals, total_spin_orbital_modes,
   h1e, eri) = asset

integral_path = relpath(H4_INTEGRALS_PATH, pkgdir(NCTSSoS))

println("Asset file: ", integral_path)
println("atoms/cell = ", reference["system"]["atoms_per_cell"])
println("lattice constant = ", reference["system"]["lattice_constant_angstrom"], " Å")
println("active space per cell = [", n_active_elec, ", ", n_active_orb, "]")
println("Nk = ", nk)
println("total spatial orbitals = ", total_spatial_orbitals)
println("naive spin-orbital modes = ", total_spin_orbital_modes)

@assert reference["system"]["atoms_per_cell"] == 4 #src
@assert reference["system"]["lattice_constant_angstrom"] == 4.0 #src
@assert reference["system"]["active_space_electrons"] == n_active_elec #src
@assert reference["system"]["active_space_orbitals"] == n_active_orb #src

# ## One-body block structure: ``h_{pq}^k``
#
# The one-body integrals (kinetic energy + nuclear attraction) are stored as
# separate **Hermitian** blocks for each k-point.  "Hermitian" means the
# matrix equals its conjugate transpose — a physical requirement ensuring
# energies are real numbers.
#
# For Nk=2 we expect exactly two blocks, `k = 0` and `k = 1`.
# We summarize each by shape, number of nonzeros, and diagonal values.
# The diagonal values represent the orbital energies at each k-point.

h1e_keys = sort(collect(keys(h1e)))
@assert h1e_keys == [0, 1] #src
@assert all(h1e[k] ≈ h1e[k]' for k in h1e_keys) #src
@assert real.(diag(h1e[0])) ≈ preflight["h1e_diag_k0"] atol = 1e-10 #src
@assert real.(diag(h1e[1])) ≈ preflight["h1e_diag_k1"] atol = 1e-10 #src

for k in h1e_keys
    diagonal = round.(real.(diag(h1e[k])); digits = 6)
    println("k = ", k,
        ": size = ", size(h1e[k]),
        ", nnz = ", count_nonzero_entries(h1e[k]),
        ", diag = ", diagonal)
end

# The diagonal values line up with the reviewed fixture.  The off-diagonal part
# is sparse, but it is not negligible — this is already an ab initio active-space
# Hamiltonian, not a toy nearest-neighbour model.

# ## Two-body momentum-conserving ERI blocks
#
# The two-electron repulsion integrals are stored as blocks keyed by four
# k-point indices `(k₁, k₂, k₃, k₄)`.  **Momentum conservation** —
# $k_1 + k_2 \equiv k_3 + k_4 \pmod{N_k}$ — means most combinations
# vanish.  For Nk=2 that leaves eight nonzero blocks.
#
# This is the same idea as the [companion page](@ref h4-chain-energy-benchmark)
# restricting to the $(0,0,0,0)$ block — but here we see *all* surviving
# blocks, including the cross-k ones that the reduced model ignores.

eri_keys = sort(collect(keys(eri)))
@assert length(eri_keys) == preflight["eri_blocks"] #src
@assert all(mod(key[1] + key[2] - key[3] - key[4], nk) == 0 for key in eri_keys) #src
@assert maximum(maximum(abs, block) for block in values(eri)) ≈ preflight["max_abs_eri"] atol = 1e-12 #src

println("Momentum-conserving ERI block keys:")
for key in eri_keys
    println("  ", key, "  nnz = ", count_nonzero_entries(eri[key]))
end
println(@sprintf("max |ERI| = %.12f", maximum(maximum(abs, block) for block in values(eri))))

# Translational symmetry removes many formally possible blocks, but the blocks
# that remain are still dense enough to define the problem's real structure.

# ## Reconstructing the Hartree–Fock energy bookkeeping
#
# The integral file represents the **active-space** electronic Hamiltonian —
# it describes only the active electrons, not the frozen core.  To recover
# the total HF energy that PySCF reports, we must add back the **constant
# shift**: the combined energy of (a) frozen-core electrons, (b) nuclear
# repulsion, and (c) core-active interactions.
#
# This bookkeeping matters: comparing the active-space energy directly
# against the total PySCF number would be mixing two different conventions.

hf_electronic = hf_energy(h1e, eri; nk, n_active_elec)
hf_total = preflight["hf_total_energy"]
hf_shift = hf_total - hf_electronic

@assert hf_electronic ≈ preflight["hf_active_space_electronic_energy"] atol = 1e-10 #src
@assert hf_shift ≈ preflight["hf_constant_shift"] atol = 1e-12 #src
@assert hf_total ≈ reference["pyscf_nk2"]["HF_energy"] atol = 1e-12 #src

println(@sprintf("HF active-space electronic energy = %.13f Ha", hf_electronic))
println(@sprintf("HF additive constant shift        = %.13f Ha", hf_shift))
println(@sprintf("HF total energy per cell          = %.13f Ha", hf_total))

# That distinction is not cosmetic.  Comparing the raw active-space Hamiltonian
# directly against the total PySCF number would be mixing two different energy
# conventions.

# ## Rebuilding the derived ``K_2`` tensor
#
# The **reduced Hamiltonian** ${}^2K$ is a 4-index tensor that combines the
# one-body and two-body integrals into a single object.  The ground-state
# energy is then simply $E = \operatorname{Tr}({}^2K \cdot {}^2D)$ — a
# linear function of the 2-RDM.  (See the
# [companion page](@ref h4-chain-energy-benchmark) for what a 2-RDM is.)
#
# Flattening the compound index `(k, p)` gives a `(16, 16, 16, 16)` tensor
# for this Nk=2 asset.  Reconstructing it validates the asset-processing
# pipeline.

k2 = build_k2(h1e, eri; nk, n_active_orb, n_total_electrons = total_active_electrons)
k2_nnz = count_nonzero_entries(k2)

@assert size(k2) == (nk * n_active_orb, nk * n_active_orb, nk * n_active_orb, nk * n_active_orb) #src
@assert maximum(abs, k2) ≈ preflight["max_abs_k2"] atol = 1e-12 #src

println("K₂ tensor size = ", size(k2))
println("K₂ nonzeros    = ", k2_nnz, " / ", length(k2))
println(@sprintf("max |K₂|       = %.12f", maximum(abs, k2)))
println(@sprintf("sample K₂[1,1,1,1] = %.12f%+.12fim", real(k2[1, 1, 1, 1]), imag(k2[1, 1, 1, 1])))

# Reconstructing ``K_2`` validates the asset-processing pipeline.  It does **not**
# mean the current package is already solving the periodic V2RDM SDP.

# ## Why a naive order-2 spin-orbital lift is the wrong interface
#
# A tempting shortcut: flatten the 16 spatial orbitals into 32 spin orbitals,
# feed that into the existing fermionic API, and call order-2.  Why doesn't
# this work?
#
# **Correlative sparsity** is the technique `NCTSSoS.jl` uses to break large
# SDPs into smaller, cheaper blocks.  It works by finding groups of variables
# that don't interact — if orbital A never appears in the same term as
# orbital B, their moment-matrix blocks can be solved independently.
#
# But for real molecules, the electron repulsion integrals couple essentially
# **every** orbital to every other.  The coupling graph is *complete* — there
# are no independent groups to split off.  Correlative sparsity collapses to
# a single dense clique, and you're back to solving one giant SDP.

spatial_edges = compound_spatial_edge_count(h1e, eri; n_active_orb)
complete_edges = complete_graph_edge_count(total_spatial_orbitals)
order2_basis_size = Int(fermionic_order2_basis_size(total_spin_orbital_modes))
order2_nuniq = Int(fermionic_order2_nuniq(total_spin_orbital_modes))

@assert spatial_edges == blocker["spatial_graph_edges"] #src
@assert spatial_edges == complete_edges #src
@assert blocker["spatial_graph_complete"] == true #src
@assert order2_basis_size == blocker["order2_basis_size"] #src
@assert order2_nuniq == blocker["order2_nuniq"] #src

println("compound spatial orbitals = ", total_spatial_orbitals)
println("spatial graph edges       = ", spatial_edges, " / ", complete_edges, " (complete graph)")
println("naive spin-orbital modes  = ", total_spin_orbital_modes)
println("order-2 basis size        = ", order2_basis_size)
println("unique order-2 moments    = ", order2_nuniq)

# 2081 basis elements and 679121 unique moments — this is not a lightweight
# docs or CI example.  The right approach is a formulation that keeps the
# periodic structure explicit — **k-blocked**, **spin-adapted**, or both —
# rather than collapsing everything into one dense matrix.  That is future
# work for `NCTSSoS.jl`.

# ## Figure values are context, not oracles
#
# The reference TOML also stores energies visually read off ("digitized")
# from the paper's figure.  These are useful for orientation — they show
# where different methods land — but they are **not exact numbers**.
# Treat them like plot labels, not regression targets.
#
# The methods shown below form a hierarchy of increasing accuracy (and cost):
#
# - **HF** (Hartree–Fock): mean-field, no correlation energy.
# - **MP2** (2nd-order perturbation theory): cheapest correlation correction.
# - **CCSD** (coupled-cluster singles & doubles): systematic wave-function method.
# - **CCSD(T)**: CCSD with perturbative triples — the "gold standard" of
#   quantum chemistry.
# - **V2RDM [4,8]**: the SDP-based method from the paper, with active space
#   [4, 8].

figure_uncertainty = reference["figure_reference"]["energy_uncertainty_Ha"]
@assert abs(hf_total - figure["HF"]) > figure_uncertainty #src

println("Approximate Nk=2 figure values (context only, ±", figure_uncertainty, " Ha):")
println("  HF      = ", figure["HF"])
println("  MP2     = ", figure["MP2"])
println("  CCSD    = ", figure["CCSD"])
println("  CCSD(T) = ", figure["CCSD_T"])
println("  V2RDM   = ", figure["V2RDM_4_8"])
println(@sprintf("Exact PySCF HF differs from the figure HF by %.6f Ha", abs(hf_total - figure["HF"])))

# So the clean reading of the asset is:
#
# - **exact / reviewed today**: Nk=2 metadata, integral parsing, HF reconstruction,
#   additive shift, and ``K_2`` construction,
# - **context only**: figure-digitized HF / MP2 / CCSD / CCSD(T) / V2RDM values,
# - **future implementation target**: a periodic V2RDM solve on a formulation
#   that preserves the k-blocked structure instead of collapsing immediately to
#   one dense 32-mode order-2 relaxation.

# ## Summary
#
# This page gives an honest entry point for the periodic H₄ benchmark: inspect
# the vendored assets, verify the physics that is actually encoded there, and be
# explicit about where the current solver interface stops.
#
# See also:
#
# - [H₄ Chain — k=0 Moment Relaxation](@ref h4-chain-energy-benchmark) for a
#   reduced but genuinely runnable SDP on the same Nk=2 asset,
# - [H₄ Periodic Chain — Correlative-Sparsity Attempt](@ref h4-periodic-correlative-sparsity)
#   for the first full-Nk=2, no-custom-basis Mosek attempt and the point where
#   current correlative sparsity stops helping,
# - [H₄ Periodic Chain — Incomplete Momentum-Basis Relaxation](@ref h4-periodic-incomplete-momentum-basis)
#   for a full-Nk=2 compromise solve with a hand-built incomplete basis,
# - [Hubbard Model](@ref hubbard-model) for a fermionic example that the current
#   API *does* solve directly,
# - [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) for the
#   creation/annihilation-operator basics,
# - [PBW Algebra Showcase](@ref pbw-algebras-showcase) for CAR simplification rules.

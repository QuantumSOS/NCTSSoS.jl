# # [H₄ Periodic Active-Space Workflow](@id h4-periodic-active-space)
#
# This page is an **asset-backed workflow example**.  It inspects the vendored
# Nk=2 periodic H₄ active-space data, reconstructs the Hartree–Fock bookkeeping
# already encoded there, and shows why the current solver API is not yet the
# right interface for the full periodic V2RDM benchmark.
#
# Concretely, we will:
#
# 1. load the reviewed Nk=2 assets shipped with the repository,
# 2. inspect the k-blocked one- and two-body integrals,
# 3. reconstruct the active-space HF energy, total HF energy, and constant shift,
# 4. rebuild the derived ``K_2`` tensor,
# 5. quantify why a naive 32-mode order-2 lift becomes a dense relaxation in disguise.
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
# The reviewed H₄ asset corresponds to a periodic hydrogen chain in **k-space**:
#
# - **4 H atoms per unit cell**,
# - bond distance **1.0 Å**,
# - lattice constant **4.0 Å**,
# - active space **[4,8] per unit cell**,
# - **Nk = 2** k-points.
#
# For this Nk=2 asset that means:
#
# - **4 active electrons per cell**,
# - **8 active spatial orbitals per k-point**,
# - **16 spatial orbitals total**,
# - **32 spin-orbital modes** under a naive spin lift.
#
# We load the reviewed metadata plus the parsed integral blocks.  The printed
# summary stays compact on purpose; the full tensors remain in the vendored text
# file rather than being dumped into the docs page.

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
# The file stores one-electron integrals as separate Hermitian blocks for each
# k-point.  For Nk=2 we expect exactly two blocks, `k = 0` and `k = 1`.
# We summarize each block by its shape, number of nonzeros, and diagonal values.

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
# The two-electron integrals are stored as blocks keyed by
# `(k₁, k₂, k₃, k₄)`.  Only momentum-conserving combinations are present.
# For Nk=2 that leaves eight blocks.  Again we print only compact summaries.

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
# The integral dump represents the **active-space electronic Hamiltonian**.  To
# recover the total PySCF Hartree–Fock energy per cell, we must add back the
# constant shift that is not explicit in the second-quantized active-space data.

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
# The research workflow also constructs a derived ``K_2`` tensor from the
# normalized one- and two-body blocks.  Flattening the compound index `(k, p)`
# gives a `(16, 16, 16, 16)` tensor for this Nk=2 asset.

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
# A tempting idea is to flatten the 16 spatial orbitals into 32 spin orbitals,
# feed that Hamiltonian into the current fermionic API, and call order-2.  For
# this asset, that would mostly hide dense work behind a familiar interface.
#
# The extracted compound-orbital coupling graph is effectively **complete**:
# once the ab initio ERIs are included, every spatial orbital couples to every
# other one.  Correlative sparsity therefore collapses to a single dense clique.

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

# `2081` basis elements and `679121` unique moments are not a lightweight docs
# or CI example.  Without a formulation that keeps the periodic structure explicit
# — k-blocked, spin-adapted, or both — the current API is doing dense work in
# disguise.

# ## Figure values are context, not oracles
#
# The reference TOML also stores energies visually digitized from the paper's
# figure.  Those numbers are useful for orientation, but they are **not exact**.
# Treat them like plot labels, not regression targets.

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
# - [Hubbard Model](@ref hubbard-model) for a fermionic example that the current
#   API *does* solve directly,
# - [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) for the
#   creation/annihilation-operator basics,
# - [PBW Algebra Showcase](@ref pbw-algebras-showcase) for CAR simplification rules.

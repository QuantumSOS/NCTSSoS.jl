```@meta
EditURL = "../literate/h4_periodic_correlative_sparsity.jl"
```

# [H₄ Periodic Chain — Correlative-Sparsity Attempt](@id h4-periodic-correlative-sparsity)

This page tests the simplest serious periodic-V2RDM experiment on the full
**Nk = 2** H₄ active-space asset:

- keep the **full Bloch-basis Hamiltonian**,
- keep the **full order-2 fermionic basis**,
- turn on **correlative sparsity only** (`cs_algo = MF()`), and
- use **Mosek** as the intended SDP backend.

The hope — spelled out in the periodic-V2RDM note — is that the **zero
Hamiltonian coefficients** forced by crystal-momentum conservation might be
enough for `NCTSSoS.jl` to discover the periodic block structure
automatically.

Short answer: **not with the current variable-level API.**

We can get through the expensive symbolic stages:

1. build the full 32-mode Hamiltonian,
2. compute correlative sparsity,
3. build the full order-2 moment relaxation.

But the correlative graph still collapses to **one 32-mode clique**, so the
resulting order-2 relaxation is essentially dense.  In a manual direct Mosek
run on this formulation, the solve eventually failed with
`MSK_RES_ERR_SPACE` (`1051`, out of space).

That is not failure of the paper's physics.  It is a statement about the
current formulation.  The periodic V2RDM blocking lives in **pair-momentum
sectors**; the present correlative-sparsity graph is built on **single
spin-orbital modes**.  Those are not the same thing.

If you want the asset inspection and bookkeeping first, see
[H₄ Periodic Active-Space Workflow](@ref h4-periodic-active-space).
If you want the later custom-basis compromise after this page, see
[H₄ Periodic Chain — Incomplete Momentum-Basis Relaxation](@ref h4-periodic-incomplete-momentum-basis).

## Setup

We reuse the vendored Nk=2 asset and the reviewed expectation fixtures from
`test/` so the docs stay tied to the same data as the regression suite.

````julia
using NCTSSoS
using MosekTools
using LinearAlgebra

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: load_expectation,
                         load_nk2_asset

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(
    Mosek.Optimizer,
    MOI.Silent() => true,
    "MSK_IPAR_NUM_THREADS" => 0,
)
const RUN_FULL_MOSEK = get(ENV, "NCTSSOS_RUN_HEAVY_H4_PERIODIC_CS", "0") == "1"

asset = load_nk2_asset();
cs_only = load_expectation("correlative_sparsity_only_attempt");

(; nk, n_active_orb, h1e, eri) = asset;

n_spatial = nk * n_active_orb         # 16 spatial orbitals total
n_modes = 2 * n_spatial               # 32 spin-orbital modes total

println("Nk = ", nk)
println("spatial orbitals = ", n_spatial)
println("spin-orbital modes = ", n_modes)
println("momentum-conserving ERI blocks = ", length(eri))

````

````
Nk = 2
spatial orbitals = 16
spin-orbital modes = 32
momentum-conserving ERI blocks = 8

````

## Full Bloch-basis Hamiltonian

We flatten each spin-orbital mode as `(k, p, σ)`:

- `k ∈ {0, 1}` is crystal momentum,
- `p ∈ {0, …, 7}` is the active orbital inside that k-sector,
- `σ ∈ {0, 1}` is spin.

The integrals are normalized exactly as in the other full-Nk=2 page:
one-body terms by `1 / Nk`, two-body terms by `1 / Nk^2`.

The only zeros we exploit are the zeros already present in the Bloch-basis
integral file.  Nothing is deleted by hand.

````julia
mode(k::Int, orb::Int, spin::Int) = spin * n_spatial + k * n_active_orb + orb + 1

function build_full_periodic_hamiltonian(
    a,
    a_dag,
    h1e,
    eri;
    nk::Int,
    n_orb::Int,
)
    h = zero(a_dag[1] * a[1])
    h1e_scale = 1 / nk
    eri_scale = 1 / nk^2

    for k in 0:nk-1, p in 0:n_orb-1, q in 0:n_orb-1
        hpq = h1e_scale * h1e[k][p + 1, q + 1]
        abs(hpq) < 1e-14 && continue
        for σ in 0:1
            h += hpq * a_dag[mode(k, p, σ)] * a[mode(k, q, σ)]
        end
    end

    for ((k1, k2, k3, k4), block) in eri
        for σ in 0:1, τ in 0:1
            for p in 0:n_orb-1, q in 0:n_orb-1, r in 0:n_orb-1, s in 0:n_orb-1
                v = eri_scale * block[p + 1, r + 1, q + 1, s + 1]
                abs(v) < 1e-14 && continue
                h += 0.5v *
                     a_dag[mode(k1, p, σ)] *
                     a_dag[mode(k2, q, τ)] *
                     a[mode(k4, s, τ)] *
                     a[mode(k3, r, σ)]
            end
        end
    end

    return simplify((h + adjoint(h)) / 2)
end

registry, (a, a_dag) = create_fermionic_variables(1:n_modes)
h = build_full_periodic_hamiltonian(a, a_dag, h1e, eri; nk, n_orb=n_active_orb)


println("Hamiltonian terms = ", length(h.terms))
````

````
Hamiltonian terms = 23752

````

So the data is doing what it should: only the momentum-conserving ERI blocks
survive.  Now we ask the actual question.

## Correlative sparsity only

We keep the physically relevant canonical sector — 4 spin-up and 4 spin-down
electrons in the full active space — but we do **not** turn on term sparsity
and we do **not** supply a custom basis.

This is the pure “let the zero Bloch-basis coefficients do the work” test.

````julia
I_poly = one(h)
n_up = 1.0 * sum(a_dag[i] * a[i] for i in 1:n_spatial)
n_dn = 1.0 * sum(a_dag[i] * a[i] for i in (n_spatial + 1):n_modes)

pop = polyopt(h, registry;
    moment_eq_constraints=[n_up - 4.0 * I_poly,
                           n_dn - 4.0 * I_poly])

config = SolverConfig(
    optimizer=SILENT_MOSEK,
    order=2,
    cs_algo=MF(),
    ts_algo=NoElimination(),
)

sparsity = compute_sparsity(pop, config)
clique_sizes = sort(length.(sparsity.corr_sparsity.cliques); rev=true)
moment_basis_sizes = length.(sparsity.corr_sparsity.clq_mom_mtx_bases)


println("correlative cliques = ", length(clique_sizes))
println("clique sizes        = ", clique_sizes)
println("order-2 basis sizes = ", moment_basis_sizes)
println("unique order-2 moments = ", cs_only["order2_nuniq"])
````

````
correlative cliques = 1
clique sizes        = [32]
order-2 basis sizes = [2081]
unique order-2 moments = 679121

````

The key result is brutally simple: **one clique**.  So pure correlative
sparsity does **not** recover a useful periodic decomposition here.

Why not?  Because the present graph is built on **single modes**.  The full
periodic quartic Hamiltonian contains enough allowed scattering terms that,
at the single-mode level, every mode becomes connected to every other mode
somewhere.  The paper's periodic blocking is finer: it lives on **degree-2
pair sectors** labeled by total momentum `K`, not on individual mode
vertices.

## Building the full order-2 moment relaxation

Even though the clique structure is useless, we can still ask how far the
current stack gets before numerical solving.  The symbolic relaxation does
build.

This stage is already heavy.  On the manual run used to write this page, the
combination of `compute_sparsity` and `moment_relax` peaked at roughly
**4–5 GB** of RSS before Mosek even started solving.

````julia
moment_problem = NCTSSoS.moment_relax(
    pop,
    sparsity.corr_sparsity,
    sparsity.cliques_term_sparsities,
)


println("moment constraints      = ", length(moment_problem.constraints))
println("unique moment variables = ", moment_problem.n_unique_moment_matrix_elements)
````

````
moment constraints      = 508611
unique moment variables = 679121

````

That is the farthest stage the current formulation reaches comfortably:
full Hamiltonian, full order-2 symbolic relaxation, one dense clique.

## Optional: direct Mosek solve

A **direct** solve (`dualize = false`) is the closest match to the raw
order-2 moment SDP we are trying to study here.

!!! warning "Manual run result"
    On the reference machine used for this rewrite, the direct Mosek solve of
    this formulation eventually failed with `Mosek.MosekError(1051)`, i.e.
    `MSK_RES_ERR_SPACE` (“out of space”), after memory use grew to roughly
    **46 GB**.  So this block is opt-in, not part of the normal docs build.

````julia
if RUN_FULL_MOSEK
    println("Attempting the full direct Mosek solve …")
    try
        result = NCTSSoS.solve_sdp(moment_problem, SILENT_MOSEK; dualize=false)
        println("termination status      = ", result.status)
        println("electronic lower bound  = ", result.objective)
    catch err
        println("Mosek solve failed with:")
        println(err)
    end
else
    println("Skipping the full Mosek solve.")
    println("Set ENV[\"NCTSSOS_RUN_HEAVY_H4_PERIODIC_CS\"] = \"1\" to attempt it.")
end
````

````
Skipping the full Mosek solve.
Set ENV["NCTSSOS_RUN_HEAVY_H4_PERIODIC_CS"] = "1" to attempt it.

````

## Reading this result correctly

This page answers a narrow but important question:

- **Do the Bloch-basis zero coefficients exist in the shipped H₄ asset?** Yes.
- **Does `cs_algo = MF()` alone recover the periodic V2RDM blocking?** No.
- **Can the current API still build the full order-2 symbolic relaxation?** Yes.
- **Does a full direct Mosek solve look practical in this formulation?** Not on
  the tested machine; it runs out of space.

So the next step should not be “turn the same knob harder.”

The next step should be a formulation that exposes the right structure:

1. **k-blocked / pair-momentum** variables or basis elements,
2. possibly **spin adaptation**,
3. or, as a pragmatic stopgap, a **smaller custom basis**.

The last option is exactly what the companion page below does.

See also:

- [H₄ Periodic Active-Space Workflow](@ref h4-periodic-active-space) — asset
  inspection, HF bookkeeping, and the dense-order-2 blocker,
- [H₄ Periodic Chain — Incomplete Momentum-Basis Relaxation](@ref h4-periodic-incomplete-momentum-basis) — a solvable full-Nk=2 compromise using a custom basis,
- [H₄ Chain — k=0 Moment Relaxation](@ref h4-chain-energy-benchmark) — the
  smaller reduced model that already gives a meaningful certified bound.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


```@meta
EditURL = "../literate/hubbard_filling_scan.jl"
```

# [Hubbard Model (Filling-Sector Scan)](@id hubbard-filling-scan)

This companion example answers one very specific question:
**what happens if you remove `moment_eq_constraints` from the canonical
Hubbard example?**

You no longer solve the half-filled problem.  You solve the
**grand-canonical** one, so the SDP is free to drop into whichever
particle-number sector has the lowest energy.

For the same $N = 4$ periodic Hubbard ring used in
[Hubbard Model](@ref hubbard-model), we will:

1. **exactly diagonalise every** $(N_\uparrow, N_\downarrow)$ sector,
2. identify the true minimum across all fillings, and
3. show that the **unconstrained SDP lower bound** lands on that sector,
   not on the half-filled $(2,2)$ sector.

The punchline is simple: if you omit `moment_eq_constraints`, the lower
bound drops from the half-filled energy $\approx -2.103$ to the global
minimum $\approx -3.419$.

---

## Exact diagonalisation across all fillings

Because this example is explicitly about **filling numbers**, we use the
occupation-basis Jordan–Wigner convention $|0\rangle =$ empty and
$|1\rangle =$ occupied.  Then bit counts in the computational basis are
literal particle counts.

````julia
using LinearAlgebra, Printf

const PAULI_Z = ComplexF64[1 0; 0 -1]
const PAULI_I = Matrix{ComplexF64}(I, 2, 2)
const JW_CREATE = ComplexF64[0 0; 1 0]   # |1⟩⟨0|
const JW_DESTROY = ComplexF64[0 1; 0 0]  # |0⟩⟨1|

function jw_fermion_op(mode::Int, kind::Symbol, nmodes::Int)
    mats = [site < mode ? PAULI_Z :
            site == mode ? (kind === :annihilate ? JW_DESTROY : JW_CREATE) :
            PAULI_I for site in 1:nmodes]
    return reduce(kron, mats)
end

function hubbard_exact_matrix(nsites::Int; t::Real = 1.0, U::Real = 4.0)
    nmodes = 2 * nsites
    a  = [jw_fermion_op(i, :annihilate, nmodes) for i in 1:nmodes]
    ad = [jw_fermion_op(i, :create,     nmodes) for i in 1:nmodes]
    dim = 2^nmodes
    H = zeros(ComplexF64, dim, dim)
    for i in 1:nsites
        j = mod1(i + 1, nsites)
        H .-= t .* (ad[i] * a[j] + ad[j] * a[i])
        H .-= t .* (ad[i + nsites] * a[j + nsites] + ad[j + nsites] * a[i + nsites])
    end
    for i in 1:nsites
        H .+= U .* (ad[i] * a[i]) * (ad[i + nsites] * a[i + nsites])
    end
    return Hermitian(H)
end

function sector_ground_state_energy(H::AbstractMatrix, nsites::Int; nup::Int, ndn::Int)
    nmodes = 2 * nsites
    upmask = (UInt(1) << nsites) - 1
    keep = [s + 1 for s in 0:(2^nmodes - 1)
            if count_ones(UInt(s) & upmask) == nup &&
               count_ones(UInt(s) >> nsites) == ndn]
    return eigmin(Hermitian(H[keep, keep]))
end

function scan_sector_ground_states(H::AbstractMatrix, nsites::Int)
    return [(; nup, ndn, total_filling = nup + ndn,
             energy = sector_ground_state_energy(H, nsites; nup, ndn))
            for nup in 0:nsites for ndn in 0:nsites]
end

function print_sector_matrix(sector_energy::Dict{Tuple{Int, Int}, Float64}, nsites::Int)
    println("Exact sector energies E₀(N↑, N↓):")
    println(rpad("N↑ \\ N↓", 8) * join([lpad(string(ndn), 12) for ndn in 0:nsites], ""))
    for nup in 0:nsites
        formatted = join([lpad(@sprintf("%.6f", sector_energy[(nup, ndn)]), 12) for ndn in 0:nsites], "")
        println(rpad(string(nup), 8) * formatted)
    end
end

function print_total_filling_summary(rows)
    println("\nBest exact sector at fixed total filling Nₑ = N↑ + N↓:")
    for ntotal in minimum(row.total_filling for row in rows):maximum(row.total_filling for row in rows)
        candidates = [row for row in rows if row.total_filling == ntotal]
        _, best_idx = findmin([row.energy for row in candidates])
        best = candidates[best_idx]
        println("Nₑ = $(best.total_filling): sector ($(best.nup), $(best.ndn)) -> E₀ = $( @sprintf("%.6f", best.energy) )")
    end
end

N = 4
t = 1.0
U = 4.0

H_mat = hubbard_exact_matrix(N; t, U)
sector_rows = scan_sector_ground_states(H_mat, N)
sector_energy = Dict((row.nup, row.ndn) => row.energy for row in sector_rows)
E_full = eigmin(H_mat)
E_half = sector_energy[(2, 2)]
E_best, best_idx = findmin([row.energy for row in sector_rows])
best_sector = sector_rows[best_idx]


print_sector_matrix(sector_energy, N)
print_total_filling_summary(sector_rows)
println("\nGlobal exact minimum:      sector ($(best_sector.nup), $(best_sector.ndn)) with E₀ = $( @sprintf("%.6f", E_best) )")
println("Half-filled exact sector:  sector (2, 2) with E₀ = $( @sprintf("%.6f", E_half) )")
println("Missing-constraint shift:  $(round(E_half - E_best; digits=6))")
````

````
Exact sector energies E₀(N↑, N↓):
N↑ \ N↓            0           1           2           3           4
0           0.000000   -2.000000   -2.000000   -2.000000    0.000000
1          -2.000000   -3.418551   -2.752158   -1.806424    2.000000
2          -2.000000   -2.752158   -2.102748    1.247842    6.000000
3          -2.000000   -1.806424    1.247842    4.581449   10.000000
4           0.000000    2.000000    6.000000   10.000000   16.000000

Best exact sector at fixed total filling Nₑ = N↑ + N↓:
Nₑ = 0: sector (0, 0) -> E₀ = 0.000000
Nₑ = 1: sector (0, 1) -> E₀ = -2.000000
Nₑ = 2: sector (1, 1) -> E₀ = -3.418551
Nₑ = 3: sector (1, 2) -> E₀ = -2.752158
Nₑ = 4: sector (2, 2) -> E₀ = -2.102748
Nₑ = 5: sector (2, 3) -> E₀ = 1.247842
Nₑ = 6: sector (3, 3) -> E₀ = 4.581449
Nₑ = 7: sector (3, 4) -> E₀ = 10.000000
Nₑ = 8: sector (4, 4) -> E₀ = 16.000000

Global exact minimum:      sector (1, 1) with E₀ = -3.418551
Half-filled exact sector:  sector (2, 2) with E₀ = -2.102748
Missing-constraint shift:  1.315802

````

The exact scan says it plainly: the **lowest** sector is
$(N_\uparrow, N_\downarrow) = (1,1)$ with energy
$E_0 \approx -3.418551$.  The half-filled sector $(2,2)$ sits higher at
$E_0 \approx -2.102748$.

So if you drop the canonical constraints, the SDP should target the
$(1,1)$ sector.  Anything else would be nonsense.

---

## Build the same Hamiltonian with NCTSSoS.jl

This is exactly the same polynomial model as in the canonical Hubbard
example.  The only difference is that we will **not** pass
`moment_eq_constraints` to [`polyopt`](@ref).

````julia
using NCTSSoS

registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
    ("c_up", 1:N),
    ("c_dn", 1:N),
]);

bonds = [(i, mod1(i + 1, N)) for i in 1:N]

hopping = -t * sum(
    c_up_dag[i] * c_up[j] + c_up_dag[j] * c_up[i] +
    c_dn_dag[i] * c_dn[j] + c_dn_dag[j] * c_dn[i]
    for (i, j) in bonds
)

interaction = U * sum(
    (c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
    for i in 1:N
)

ham = hopping + interaction;
````

---

## Solve the unconstrained SDP

No `moment_eq_constraints`.  No fixed $(N_\uparrow, N_\downarrow)$.
Just the Hamiltonian over the full Fock space.

````julia
using MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG" => 0,
    "MSK_IPAR_NUM_THREADS" => 0);

config = SolverConfig(optimizer = SOLVER, order = 3, ts_algo = MMD())

pop_unconstrained = polyopt(ham, registry)
result_first = cs_nctssos(pop_unconstrained, config)
result_unconstrained = cs_nctssos_higher(pop_unconstrained, result_first, config)


println("Unconstrained SDP (first sparse pass):  $( @sprintf("%.6f", result_first.objective) )")
println("Unconstrained SDP (refined):            $( @sprintf("%.6f", result_unconstrained.objective) )")
println("Exact minimum over all sectors:         $( @sprintf("%.6f", E_best) )")
println("Exact half-filled sector (2,2):         $( @sprintf("%.6f", E_half) )")
````

````
Unconstrained SDP (first sparse pass):  -4.962560
Unconstrained SDP (refined):            -3.418551
Exact minimum over all sectors:         -3.418551
Exact half-filled sector (2,2):         -2.102748

````

After one refinement step, the unconstrained SDP lands on the exact
minimum over **all** filling sectors.  That is exactly what it should do.
Without `moment_eq_constraints`, there is nothing in the problem telling
the solver to stay in the half-filled sector.

---

## Results

| Quantity | Value |
|:--------|------:|
| Exact minimum over all $(N_\uparrow, N_\downarrow)$ sectors | $-3.419$ at $(1,1)$ |
| Exact half-filled sector $(2,2)$ | $-2.103$ |
| Unconstrained SDP (order 3 TS + higher) | $\approx -3.419$ |
| Shift from half-filling when constraints are omitted | $1.316$ |

The lower bound you get **without** `moment_eq_constraints` is therefore
the **grand-canonical** value, not the half-filled canonical one.

---

## Summary

If you remove `moment_eq_constraints`, you are solving the wrong problem.
The solver is then free to pick the lowest filling sector in the entire
Fock space, and for this $N = 4$, $U/t = 4$ ring that sector is
$(N_\uparrow, N_\downarrow) = (1,1)$ with energy
$E_0 \approx -3.418551$.

The half-filled sector $(2,2)$ instead has
$E_0 \approx -2.102748$, so omitting the canonical constraints lowers the
answer by about $1.3158$.  That is not a numerical quirk.  It is the SDP
correctly solving a different optimisation problem.

### See also

- [Hubbard Model](@ref hubbard-model) — canonical half-filling with
  `moment_eq_constraints`
- [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
  fermionic creation/annihilation operators and CAR

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


```@meta
EditURL = "../literate/strengthening.jl"
```

# [Strengthening and the States a Bound Applies To](@id strengthening)

A stronger relaxation adds conditions. Before using its bound, check which
states satisfy those conditions. `result.valid_for` prints that domain;
`result.bound_kind` distinguishes a numerical optimum from a certificate.

| Kind | What it adds |
|:-----|:-------------|
| `DensityMatrix(sites)` | Positivity of the normalized density matrix on Pauli sites |
| `InvariantCone(sites, r)` | A named Kramers-invariant cone on three or four Pauli sites |
| `Stationarity()` | Commutator equalities for moment words in the selected blocks |
| `StateOptimality()` | Positivity of the objective's second variation on each selected row block |
| `FermionicConditions(modes)` | Fermionic positivity conditions, without a fixed particle number |
| `Sector(Q, q)` | A restriction to states supported on the eigenspace `Q = q` |

Scalar expectation constraints are different: `ExpectationPositive(p)` and
`ExpectationZero(p)` impose one scalar inequality or equality. They do not
create an operator localizing matrix or declare a sector.

The examples use the public [polynomial optimization API](@ref polynomial-optimization)
and absolute tolerance `1e-6`. The comparison values are numerical results,
not exact solver certificates.

````julia
using NCTSSoS, Mosek, LinearAlgebra, Markdown
import NCTSSoS
solver = MosekBackend();
````

Print the complete validity report, including all combined conditions:

````julia
function report_bound(label, result)
    println(label, ": ", result.bound)
    println("bound_kind: ", result.bound_kind)
    println("valid_for: ", replace(sprint(show, result.valid_for),
        "; " => ";\n  ", "eigenspace " => "eigenspace\n    "))
end
````

## Compare the four Pauli kinds on one problem

Use three sites and the Hamiltonian

```math
H = -(X_1Y_2 + X_1Y_3 + X_2X_3 + Y_2Y_3 + Z_2Z_3).
```

Its exact minimum is `-3`. The ordinary matrix below provides a separate
reference. Its spectrum is `[-3,-3,-1,-1,1,1,3,3]`.

````julia
_, (x, y, z) = pauli_variables(3);
H = -(x[1]*y[2] + x[1]*y[3] + x[2]*x[3] + y[2]*y[3] + z[2]*z[3]);
problem = Problem(H);
X = [0 1; 1 0];
Y = Complex{Int}[0 -im; im 0];
Z = [1 0; 0 -1];
I2 = Matrix{Int}(I, 2, 2);
Hmat = -(kron(X,Y,I2) + kron(X,I2,Y) + kron(I2,X,X) +
         kron(I2,Y,Y) + kron(I2,Z,Z));
eigvals(Hermitian(Hmat))
````

````
8-element Vector{Float64}:
 -3.0
 -2.999999999999999
 -1.0000000000000002
 -1.0
  0.9999999999999993
  0.9999999999999994
  2.999999999999999
  2.9999999999999996
````

Keep order 1 fixed. Only the added condition changes. The table is built
from these executed results, rather than copied from a previous run.

````julia
kinds = ((), Stationarity(), StateOptimality(),
         InvariantCone(1:3, "YII"), DensityMatrix(1:3));
labels = ("Plain", "Stationarity", "StateOptimality", "InvariantCone", "DensityMatrix");
results = [solve(problem, Relaxation(order=1, strengthen=kind); solver) for kind in kinds];
expected = (-5, -5, -5, -sqrt(15), -3);
table_rows = ["| $(label) | $(round(r.bound; digits=9)) | -3 |"
              for (label, r) in zip(labels, results)];
Markdown.parse(join(vcat(["| Condition | Numerical bound | Exact ground energy |",
                          "|:----------|----------------:|--------------------:|"],
                         table_rows), "\n"))
````

```@raw html
<div class="markdown"><table><tr><th align="left">Condition</th><th align="right">Numerical bound</th><th align="right">Exact ground energy</th></tr><tr><td align="left">Plain</td><td align="right">-5.000000004</td><td align="right">-3</td></tr><tr><td align="left">Stationarity</td><td align="right">-5.000000004</td><td align="right">-3</td></tr><tr><td align="left">StateOptimality</td><td align="right">-4.999999999</td><td align="right">-3</td></tr><tr><td align="left">InvariantCone</td><td align="right">-3.872983359</td><td align="right">-3</td></tr><tr><td align="left">DensityMatrix</td><td align="right">-3.000000003</td><td align="right">-3</td></tr></table>
</div>
```

```@raw html
<br>
```

`Stationarity` and `StateOptimality` give **no gain on this problem**.
The invariant cone improves the bound but does not reach `-3`; its value
near `-sqrt(15)` is a checked numerical value, not an exact certificate of
the relaxation optimum. Density-matrix positivity reaches the reference
within the stated tolerance.

The domains below are part of the answer. Density-matrix positivity and
the invariant cone are valid for each state. Stationarity is valid for
eigenstates of the objective; optimality is valid for its ground space.
A ground state is an eigenstate, so both kinds still give lower bounds on
the ground energy of an **unconstrained** problem.

````julia
for (label, result) in zip(labels, results)
    report_bound(label, result)
end
````

````
Plain: -5.000000004400016
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions
Stationarity: -5.000000004400017
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions;
  Stationarity: valid for eigenstates of the objective
StateOptimality: -4.999999998864657
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions;
  StateOptimality: valid for the ground space of the objective
InvariantCone: -3.872983358787062
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions;
  InvariantCone: valid for each state
DensityMatrix: -3.000000002704607
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions;
  DensityMatrix: valid for each state

````

### Why the invariant cone is valid for each state

Let ``\Theta=P_r K``, where ``K`` is complex conjugation. The Pauli label
``r`` has an odd number of Y factors, so ``\Theta^2=-1``. The cone keeps
non-identity Pauli words with ``\Theta w\Theta^{-1}=w``.

Let ``E`` be the sum of squared expectations over the kept words, and
``O`` the sum over the other non-identity words. For a pure state,
``\langle\psi|\Theta\psi\rangle=0``, so
``\operatorname{tr}(\rho\Theta\rho\Theta^{-1})=0``. Together with purity,
this gives

```math
1+E-O=0,\qquad 1+E+O=2^k,
\qquad E=2^{k-1}-1.
```

The sum of squares ``E`` is convex, so mixed states satisfy
``E\leq 2^{k-1}-1``. Thus the cone radius is
``\sqrt{2^{k-1}-1}``. This argument applies to both supported sizes,
three and four sites. It does not require the input state to be invariant
under ``\Theta``.

## A strict stationarity gain at order 2

On four sites, use four Pauli words. Stationarity improves this bound
from about `-4.477216` to `-3.889855`. It does not reach the ground energy
`1-2sqrt(5)`, about `-3.472136`. The independent matrix check below gives
that reference; all solver comparisons use absolute tolerance `1e-6`.

````julia
_, (x4, y4, z4) = pauli_variables(4);
H4 = -2*x4[1]*z4[2]*x4[3]*x4[4] + y4[2]*x4[3]*z4[4] +
     2*y4[1]*x4[2]*y4[4] - 2*y4[1]*z4[2]*x4[3]*x4[4];
H4mat = -2*kron(X,Z,X,X) + kron(I2,Y,X,Z) +
        2*kron(Y,X,I2,Y) - 2*kron(Y,Z,X,X);
reference4 = 1-2sqrt(5);
eigvals(Hermitian(H4mat))
````

````
16-element Vector{Float64}:
 -3.472135954999581
 -3.47213595499958
 -3.472135954999579
 -3.472135954999577
 -3.0000000000000018
 -3.0000000000000004
 -3.0
 -2.999999999999999
  0.9999999999999994
  1.0000000000000013
  1.0000000000000016
  1.0000000000000044
  5.4721359549995725
  5.472135954999576
  5.47213595499958
  5.47213595499958
````

Compare the two order-2 relaxations with that matrix reference.

````julia
plain4 = solve(Problem(H4), Relaxation(order=2); solver);
stationary4 = solve(Problem(H4), Relaxation(order=2,
    strengthen=Stationarity()); solver);
report_bound("Four-site plain", plain4)
report_bound("Four-site stationary", stationary4)
````

````
Four-site plain: -4.477216145767562
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions
Four-site stationary: -3.8898545235604844
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions;
  Stationarity: valid for eigenstates of the objective

````

## Fermionic conditions do not fix particle number

The default conditions are P, Q and G. The next example requests all six
supported kinds. P/Q use annihilator/creator pairs; G uses particle-hole
rows and identity; T1 uses annihilator triples; T2 uses creator-annihilator-
annihilator triples. T2prime augments T2 with annihilator rows and is not
the same condition as T2.

For three modes, the exact minimum of `-N` is `-3`, attained by the fully
occupied state. This example checks the conditions and their state domain;
it does not claim a gain over the plain relaxation.

````julia
_, (a, ad) = fermionic_variables(3);
N = sum(ad .* a);
fermi = solve(Problem(-N), Relaxation(order=1,
    strengthen=FermionicConditions(1:3;
        conditions=(:P, :Q, :G, :T1, :T2, :T2prime))); solver);
report_bound("Fermionic conditions", fermi)
````

````
Fermionic conditions: -3.000000000950703
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions;
  FermionicConditions: valid for each state

````

A particle-number sector is a separate declaration. On the `N=1`
eigenspace, the exact minimum of `-N` is `-1`, not the unrestricted `-3`.
`Sector(N,1)` restricts states; it is not the operator identity `Zero(N-1)`.

````julia
number_sector = solve(Problem(-N), Relaxation(order=1,
    strengthen=Sector(N, 1)); solver);
report_bound("One-particle sector", number_sector)
````

````
One-particle sector: -1.0000000032551957
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions;
  Sector: valid for states supported on the eigenspace
    (ad1 a1 + ad2 a2 + ad3 a3) = 1

````

## A singlet sector changes the reference energy

On two Pauli sites, use `H = -(X₁X₂ + Y₁Y₂ + Z₁Z₂)` and total spin squared
`S²`. The exact unrestricted minimum is `-1`. On the singlet the energy
is instead `+3`. Check both with ordinary matrices before solving.

````julia
_, (xs, ys, zs) = pauli_variables(2);
Hs = -(xs[1]*xs[2] + ys[1]*ys[2] + zs[1]*zs[2]);
S2 = 0.25*((xs[1]+xs[2])^2 + (ys[1]+ys[2])^2 + (zs[1]+zs[2])^2);
Hs_matrix = -(kron(X,X) + kron(Y,Y) + kron(Z,Z));
singlet = [0, 1, -1, 0] / sqrt(2);
(eigmin(Hermitian(Hs_matrix)), real(dot(singlet, Hs_matrix*singlet)))
````

````
(-1.0, 2.9999999999999996)
````

The plain order-1 relaxation gives `-3`, a non-tight lower bound on `-1`.
The selected singlet gives `+3`, the minimum only in that sector.

````julia
unrestricted = solve(Problem(Hs), Relaxation(order=1); solver);
in_singlet = solve(Problem(Hs), Relaxation(order=1,
    strengthen=Sector(S2, 0)); solver);
report_bound("Unrestricted order 1", unrestricted)
report_bound("Singlet sector", in_singlet)
````

````
Unrestricted order 1: -3.0000000000000466
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions
Singlet sector: 2.9999999912045983
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states;
  no additional state conditions;
  Sector: valid for states supported on the eigenspace
    (1.5 + (0.5) Z1 Z2 + (0.5) X1 X2 + (0.5) Y1 Y2) = 0

````

**Do not use the sector value `+3` as a lower bound on the unrestricted
minimum `-1`.** Read the complete `valid_for` report each time.
Bosonic sectors are outside the scope of this page.

## Scalar constraints are not operator constraints

`ExpectationPositive(q)` imposes only `⟨q⟩ ≥ 0`. In contrast, `Positive(q)`
uses a localizing matrix for operator positivity. For involutions `q` and
`v`, minimizing `v*q*v` distinguishes them: the exact minima are `-1`
and `+1`, respectively. The scalar reference is `q=Z`, `v=X`, state `|0⟩`.

````julia
_, (u,) = site_variables(Unipotent(), 1, 2);
q, v = u;
scalar = solve(Problem(v*q*v; constraints=(ExpectationPositive(q),)),
    Relaxation(order=3); solver);
operator = solve(Problem(v*q*v; constraints=(Positive(q),)),
    Relaxation(order=3); solver);
report_bound("Scalar positivity", scalar)
report_bound("Operator positivity", operator)
````

````
Scalar positivity: -1.0000000046945017
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states satisfying: ⟨s1_1⟩ ≥ 0 (scalar)
Operator positivity: 0.99999998385029
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states satisfying: s1_1 ⪰ 0 (operator)

````

`ExpectationZero(q-c)` fixes the single expectation `⟨q⟩ = c`. It does
not set all localized products to zero and does not declare an eigenspace.
Here the objective itself is fixed, so the exact answer is `1/4`.

````julia
equal = solve(Problem(q; constraints=(ExpectationZero(q-0.25),)),
    Relaxation(order=1); solver);
report_bound("Scalar equality", equal)
````

````
Scalar equality: 0.24999999513611726
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states satisfying: ⟨-0.25 + s1_1⟩ = 0 (scalar)

````

**Scalar constraint with strengthening.** A scalar expectation constraint
adds one scalar row to the program and combines with symmetry, grouping,
split and strengthening.

````julia
scalar_density = NCTSSoS.build(Problem(z[1]; constraints=(ExpectationZero(z[1]),)),
    Relaxation(order=1, strengthen=DensityMatrix(1:1)));
scalar_density_result = solve(scalar_density; solver);
report_bound("Scalar equality with density matrix", scalar_density_result)
````

````
Scalar equality with density matrix: -1.1292437644944657e-17
bound_kind: numerical_optimum_not_certificate
valid_for: Normalized states satisfying: ⟨Z1⟩ = 0 (scalar);
  DensityMatrix: valid for each state

````

## Five rejected combinations and how to correct them

These messages are produced by executed calls. Each block expects an
`ArgumentError`; any other failure is rethrown. They describe the current
API, not a complete compatibility table.

````julia
function show_rejection(f)
    try
        f()
        error("Expected the call to be rejected")
    catch err
        err isa ArgumentError || rethrow()
        println(replace(sprint(showerror, err),
            ": " => ":\n", "; " => ";\n", ", " => ",\n", ". " => ".\n"))
    end
end
````

**Wrong algebra.** Density-matrix positivity requires Pauli variables.

````julia
show_rejection() do
    NCTSSoS.build(Problem(-N), Relaxation(order=1, strengthen=DensityMatrix(1:3)))
end
````

````
ArgumentError:
DensityMatrix and InvariantCone require Pauli variables;
for bosonic variables use Sector,
and for fermionic variables use FermionicConditions or Sector.

````

**Stationarity with a constraint.** A constrained minimizer need not be an
eigenstate of the objective. Remove stationarity for this problem.

````julia
show_rejection() do
    NCTSSoS.build(Problem(z[1]; constraints=(Positive(x[1]),)),
        Relaxation(order=1, strengthen=Stationarity()))
end
````

````
ArgumentError:
Stationarity and StateOptimality require an unconstrained problem:
a constraint restricts the state set,
so its minimizer need not be an eigenstate of the objective.
Remove these strengthening kinds.

````

**Sector with StateOptimality.** This combination is not supported.

````julia
show_rejection() do
    NCTSSoS.build(Problem(z[1]), Relaxation(order=1,
        strengthen=(Sector(z[1], 1), StateOptimality())))
end
````

````
ArgumentError:
StateOptimality with a sector is not supported in this delivery;
remove StateOptimality or the sector restriction.

````

**Symmetry moves the sector.** The swap preserves the objective but does
not preserve the selected `Z₁=1` sector. Remove it or select an invariant sector.

````julia
show_rejection() do
    NCTSSoS.build(Problem(z[1]+z[2]; symmetry=(SitePermutation([2,1,3]),)),
        Relaxation(order=1, strengthen=Sector(z[1], 1)))
end
````

````
ArgumentError:
A supplied symmetry moves the selected sector;
remove the generator or choose an invariant sector.

````

**Non-Hermitian sector operator.** Supply the intended self-adjoint
observable; the API does not silently replace it with its Hermitian part.

````julia
show_rejection() do
    Sector(im*x[1], 0)
end
````

````
ArgumentError:
Sector operator ((0.0 + 1.0im) X1) is not Hermitian;
supply a self-adjoint polynomial.
If the Hermitian part is intended,
explicitly use (p + p')/2 via 0.5*(p + p').

````

## Read the domain before comparing bounds

Keep `bound_kind` and all of `valid_for` with each reported number.
Compare bounds only for the state domain they describe. In particular,
a sector restriction changes that domain, while scalar and operator
constraints impose different conditions even when they use the same polynomial.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


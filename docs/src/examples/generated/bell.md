```@meta
EditURL = "../literate/bell.jl"
```

# [Bell inequalities](@id bell-inequalities)

Bell inequalities test whether quantum mechanics can be explained by local hidden variable
theories. They are linear combinations of expectation values with bounds that differ between
classical and quantum theories.

The general form of a Bell inequality is:

```math
\sum_{i,j} c_{ij} \langle A_i B_j \rangle \leq C
```

where $A_i$ and $B_j$ are observables measured by Alice and Bob, $c_{ij}$ are coefficients,
and $C$ is the classical bound. Quantum mechanics can exceed this bound.

## Setup

We use `NCTSSoS.jl` for polynomial optimization and `Mosek` as the SDP solver backend.
Any SDP solver works — replace `Mosek.Optimizer` with `COSMO.Optimizer` or
`Clarabel.Optimizer` for open-source alternatives.

````julia
using NCTSSoS, MosekTools
````

````
Precompiling packages...
   1535.0 ms  ✓ ChordalGraph
   1997.4 ms  ✓ CliqueTrees → AMDExt
  40158.6 ms  ✓ MathOptInterface
  23276.0 ms  ✓ JuMP
  44074.6 ms  ✓ MathOptInterface → MathOptInterfaceCliqueTreesExt
  96259.6 ms  ✓ Clarabel
  31265.3 ms  ✓ NCTSSoS
  7 dependencies successfully precompiled in 168 seconds. 92 already precompiled.
Precompiling packages...
   3410.0 ms  ✓ MosekTools
  1 dependency successfully precompiled in 4 seconds. 59 already precompiled.

````

## Key Concepts: Unipotent and Projector Variables

Bell inequalities use two types of measurement operators:

1. **Unipotent operators** ($U^2 = I$): Model ±1-valued observables like Pauli matrices
2. **Projector operators** ($P^2 = P$): Model projection measurements like $|0\rangle\langle 0|$

Let's demonstrate both:

### Unipotent Variables (U² = I)

Create operators that square to identity:

````julia
reg_unip, (A, B) = create_unipotent_variables([("A", 1:2), ("B", 1:2)])
````

````
(VariableRegistry with 4 variables: A₁, A₂, B₁, B₂, (NCTSSoS.NormalMonomial{NCTSSoS.UnipotentAlgebra, UInt8}[A₁, A₂], NCTSSoS.NormalMonomial{NCTSSoS.UnipotentAlgebra, UInt8}[B₁, B₂]))
````

reg_unip: registry storing variable names and algebra type
A: Alice's measurement operators [A₁, A₂] on site 1
B: Bob's measurement operators [B₁, B₂] on site 2

Verify the unipotent property (U² = I):

````julia
A[1] * A[1]  # should simplify to identity
````

````
1
````

Check that operators on different sites commute:

````julia
A[1] * B[1] == B[1] * A[1]  # true: different sites commute
````

````
true
````

### Projector Variables (P² = P)

Create operators that are idempotent:

````julia
reg_proj, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 1:2)])
````

````
(VariableRegistry with 4 variables: P₁, P₂, Q₁, Q₂, (NCTSSoS.NormalMonomial{NCTSSoS.ProjectorAlgebra, UInt8}[P₁, P₂], NCTSSoS.NormalMonomial{NCTSSoS.ProjectorAlgebra, UInt8}[Q₁, Q₂]))
````

reg_proj: registry for projector algebra
P: Alice's projectors [P₁, P₂] on site 1
Q: Bob's projectors [Q₁, Q₂] on site 2

Verify the idempotent property (P² = P):

````julia
monomials(P[1] * P[1])  # should be [P[1]]
````

````
1-element Vector{NCTSSoS.NormalMonomial{NCTSSoS.ProjectorAlgebra, UInt8}}:
 P₁
````

---
## Linear Bell Inequalities
---

### CHSH Inequality

The CHSH inequality involves two parties with two ±1-valued observables each.
The objective function is:

```math
f(A_1, A_2, B_1, B_2) = \langle A_1 B_1 \rangle + \langle A_1 B_2 \rangle + \langle A_2 B_1 \rangle - \langle A_2 B_2 \rangle
```

Classical bound: $f \leq 2$. Quantum bound (Tsirelson): $f \leq 2\sqrt{2} \approx 2.828$.

#### Step 1: Create unipotent variables for CHSH

````julia
registry, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
````

````
(VariableRegistry with 4 variables: x₁, x₂, y₁, y₂, (NCTSSoS.NormalMonomial{NCTSSoS.UnipotentAlgebra, UInt8}[x₁, x₂], NCTSSoS.NormalMonomial{NCTSSoS.UnipotentAlgebra, UInt8}[y₁, y₂]))
````

registry: variable registry encoding U² = I constraint
x: Alice's observables [x₁, x₂] = [A₁, A₂]
y: Bob's observables [y₁, y₂] = [B₁, B₂]

#### Step 2: Define the CHSH objective function

````julia
f = 1.0 * x[1] * y[1] +  # ⟨A₁B₁⟩ term
    1.0 * x[1] * y[2] +  # ⟨A₁B₂⟩ term
    1.0 * x[2] * y[1] -  # ⟨A₂B₁⟩ term
    1.0 * x[2] * y[2]    # -⟨A₂B₂⟩ term
````

````
x₁y₁ + x₁y₂ + x₂y₁ + -x₂y₂
````

f: polynomial representing the CHSH Bell operator

Inspect the polynomial structure:

````julia
(monomials(f),      # list of monomials in f
 coefficients(f))   # corresponding coefficients
````

````
(NCTSSoS.NormalMonomial{NCTSSoS.UnipotentAlgebra, UInt8}[x₁y₁, x₁y₂, x₂y₁, x₂y₂], [1.0, 1.0, 1.0, -1.0])
````

#### Step 3: Create the optimization problem

`polyopt` creates a **minimization** problem. To find the maximum quantum
violation we minimize `-f` and negate the result (same pattern as I₃₃₂₂ below).

````julia
pop = polyopt(-f, registry)
````

````
Optimization Problem (UnipotentAlgebra)
────────────────────────────────────
Objective:
    -x₁y₁ + -x₁y₂ + -x₂y₁ + x₂y₂

Equality constraints (0):
    (none)

Inequality constraints (0):
    (none)

Moment equality constraints (0):
    (none)

Variables (4):
    x₁, x₂, y₁, y₂

````

pop: minimize -f ≡ maximize f, subject to algebraic constraints (U² = I)

#### Step 4: Configure and run the SDP solver

````julia
solver_config = SolverConfig(
    optimizer = Mosek.Optimizer,  # SDP solver backend
    order = 1                      # relaxation order (hierarchy level)
)
````

````
NCTSSoS.SolverConfig(Mosek.Optimizer, 1, nothing, NCTSSoS.NoElimination(), NCTSSoS.NoElimination(), nothing)
````

solver_config: specifies solver and relaxation parameters

````julia
result = cs_nctssos(pop, solver_config)
````

````
Objective: -2.828427124561678
Correlative Sparsity (UnipotentAlgebra): 

   maximum clique size: 4
   number of cliques: 1
   Clique 1: 
       Variables: [:x₁, :x₂, :y₁, :y₂]
       Bases length: 5
       Constraints: 
   Global Constraints: 
Term Sparsity:
Clique 1:
   Moment Matrix Block Sizes: [5]
   Moment Matrix:
Number of Activated supp:   11
Number of Bases Activated in each sub-block[5]

   Localizing Matrix:
Unique Moment Matrix Elements: 11

````

result: optimization result containing objective value and solver info

#### Step 5: Extract the upper bound

````julia
chsh_bound = -result.objective
````

````
2.828427124561678
````

chsh_bound: upper bound on maximal quantum violation (negate since we minimized -f)

Compare with Tsirelson's bound:

````julia
tsirelson_bound = 2 * sqrt(2)
````

````
2.8284271247461903
````

tsirelson_bound: theoretical maximum = 2√2 ≈ 2.828

````julia
abs(chsh_bound - tsirelson_bound)  # difference (should be ~1e-7)
````

````
1.8451240535455327e-10
````

!!! tip "Going further: shrink this SDP with symmetry"
    The CHSH operator is invariant under a 16-element symmetry group. The
    [CHSH with Symmetry Reduction](@ref chsh-symmetry) example shows how to
    use that group to replace the dense `5\times 5` PSD block with three
    independent `1\times 1` PSD blocks while still recovering `2\sqrt{2}`.
    The architectural picture is on the [Symmetry-Adapted Basis](@ref symmetry-adapted-basis)
    manual page.

---
### $I_{3322}$ Inequality

The $I_{3322}$ inequality uses **projector** observables (P² = P) with three measurements
per party [pal2010maximal](@cite).

```math
f = \langle A_1(B_1+B_2+B_3) \rangle + \langle A_2(B_1+B_2-B_3) \rangle + \langle A_3(B_1-B_2) \rangle - \langle A_1 \rangle - 2\langle B_1 \rangle - \langle B_2 \rangle
```

Classical bound: $f \leq 0$. Quantum bound: $f \leq 0.25$.

#### Step 1: Create projector variables

````julia
registry, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
````

````
(VariableRegistry with 6 variables: x₁, x₂, x₃, y₁, y₂, y₃, (NCTSSoS.NormalMonomial{NCTSSoS.ProjectorAlgebra, UInt8}[x₁, x₂, x₃], NCTSSoS.NormalMonomial{NCTSSoS.ProjectorAlgebra, UInt8}[y₁, y₂, y₃]))
````

registry: variable registry encoding P² = P constraint
x: Alice's projectors [x₁, x₂, x₃] = [A₁, A₂, A₃]
y: Bob's projectors [y₁, y₂, y₃] = [B₁, B₂, B₃]

#### Step 2: Define the I₃₃₂₂ objective function

````julia
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) +  # A₁(B₁+B₂+B₃)
    1.0 * x[2] * (y[1] + y[2] - y[3]) +  # A₂(B₁+B₂-B₃)
    1.0 * x[3] * (y[1] - y[2]) -         # A₃(B₁-B₂)
    1.0 * x[1] -                          # -A₁
    2.0 * y[1] -                          # -2B₁
    1.0 * y[2]                            # -B₂
````

````
-x₁ + -2.0 * y₁ + -y₂ + x₁y₁ + x₁y₂ + x₁y₃ + x₂y₁ + x₂y₂ + -x₂y₃ + x₃y₁ + -x₃y₂
````

f: I₃₃₂₂ Bell polynomial

Check the number of terms:

````julia
length(monomials(f))  # number of monomials
````

````
11
````

#### Step 3: Solve (minimizing -f to find maximum of f)

````julia
pop = polyopt(-f, registry)
````

````
Optimization Problem (ProjectorAlgebra)
────────────────────────────────────
Objective:
    x₁ + 2.0 * y₁ + y₂ + -x₁y₁ + -x₁y₂ + -x₁y₃ + -x₂y₁ + -x₂y₂ + x₂y₃ + -x₃y₁ + x₃y₂

Equality constraints (0):
    (none)

Inequality constraints (0):
    (none)

Moment equality constraints (0):
    (none)

Variables (6):
    x₁, x₂, x₃, y₁, y₂, y₃

````

pop: minimize -f (equivalent to maximize f)

````julia
solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=2)
````

````
NCTSSoS.SolverConfig(Mosek.Optimizer, 2, nothing, NCTSSoS.NoElimination(), NCTSSoS.NoElimination(), nothing)
````

order=2: second level of the moment hierarchy

````julia
result = cs_nctssos(pop, solver_config)
i3322_bound = -result.objective
````

````
0.25093972390866476
````

i3322_bound: upper bound on I₃₃₂₂ violation (negate since we minimized -f)

````julia
i3322_bound  # should be close to 0.25
````

````
0.25093972390866476
````

---
### Exploiting Sparsity for Larger Problems

Higher relaxation orders improve bounds but increase SDP size.
**Sparsity exploitation** reduces computational cost:

1. **Correlative Sparsity (CS)**: Decomposes problem by variable interactions
2. **Term Sparsity (TS)**: Removes unnecessary monomials from moment matrices

Let's solve I₃₃₂₂ at order=6 using correlative sparsity:

#### Without sparsity (for comparison, order=3)

````julia
registry, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]
pop = polyopt(-f, registry)

solver_config_dense = SolverConfig(optimizer=Mosek.Optimizer, order=3)
````

````
NCTSSoS.SolverConfig(Mosek.Optimizer, 3, nothing, NCTSSoS.NoElimination(), NCTSSoS.NoElimination(), nothing)
````

solver_config_dense: no sparsity exploitation

````julia
@time result_dense = cs_nctssos(pop, solver_config_dense)
bound_dense = -result_dense.objective
````

````
0.25087567011299783
````

bound_dense: bound without sparsity

````julia
bound_dense
````

````
0.25087567011299783
````

#### With correlative sparsity (order=6)

````julia
solver_config_sparse = SolverConfig(
    optimizer = Mosek.Optimizer,
    order = 6,             # higher order for better bound
    cs_algo = MF()         # use MaxFlow algorithm for correlative sparsity
)
````

````
NCTSSoS.SolverConfig(Mosek.Optimizer, 6, nothing, CliqueTrees.MF(), NCTSSoS.NoElimination(), nothing)
````

cs_algo=MF(): enables correlative sparsity via chordal graph decomposition

````julia
@time result_sparse = cs_nctssos(pop, solver_config_sparse)
bound_sparse = -result_sparse.objective
````

````
0.25087541281436
````

bound_sparse: improved bound using sparsity

````julia
bound_sparse  # closer to theoretical 0.25
````

````
0.25087541281436
````

Improvement in bound:

````julia
bound_dense - bound_sparse  # positive = improvement
````

````
2.5729863784018647e-7
````

---
## Nonlinear Bell Inequalities

Nonlinear Bell inequalities involve polynomial functions of expectation values,
not just linear combinations. They can detect non-locality where linear inequalities fail.

### Covariance Bell Inequality

The covariance between observables A and B is:

```math
\text{Cov}(A, B) = \langle AB \rangle - \langle A \rangle \langle B \rangle
```

This is **nonlinear** because it involves products of expectation values.

The covariance Bell inequality [pozsgay2017Covariance](@cite):

```math
f = \sum_{i,j} s_{ij} \text{Cov}(A_i, B_j)
```

with signs $s_{ij} \in \{+1, -1\}$. Classical bound: $f \leq 4.5$. Quantum bound: $f = 5$.

#### Step 1: Create unipotent variables

````julia
registry, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
````

````
(VariableRegistry with 6 variables: x₁, x₂, x₃, y₁, y₂, y₃, (NCTSSoS.NormalMonomial{NCTSSoS.UnipotentAlgebra, UInt8}[x₁, x₂, x₃], NCTSSoS.NormalMonomial{NCTSSoS.UnipotentAlgebra, UInt8}[y₁, y₂, y₃]))
````

x: Alice's observables [A₁, A₂, A₃]
y: Bob's observables [B₁, B₂, B₃]

#### Step 2: Define the identity monomial

The identity monomial `ID` serves a dual purpose: it anchors the algebra type
and promotes `StatePolynomial` → `NCStatePolynomial` (what `polyopt` expects).
Multiplying by `ID` is mathematically a no-op (×𝟙 = identity) but tells the
type system which algebra the non-commutative layer belongs to.

````julia
ID = one(typeof(x[1]))
````

````
𝟙
````

ID: identity monomial (𝟙) inferred from the variable type

````julia
ID  # display the identity
````

````
𝟙
````

#### Step 3: Define the covariance function using state polynomials

State polynomials use `ς(·)` (type `\varsigma` + Tab, or use the ASCII alias
`varsigma`) to denote expectation values ⟨·⟩.

`ς` applied to a `Polynomial` returns a `StatePolynomial`; applied to a
single `NormalMonomial` it returns a `StateWord` (a single expectation factor).

````julia
cov(a, b) = 1.0 * ς(x[a] * y[b]) * ID -  # ⟨AᵢBⱼ⟩
            1.0 * ς(x[a]) * ς(y[b]) * ID  # -⟨Aᵢ⟩⟨Bⱼ⟩
````

````
cov (generic function with 1 method)
````

cov(a,b): covariance Cov(Aₐ, Bᵦ) as a state polynomial
ς (varsigma): expectation value operator, type \varsigma + Tab

Example: Cov(A₁, B₁)

````julia
cov(1, 1)
````

````
-⟨x₁⟩⟨y₁⟩ + ⟨x₁y₁⟩
````

#### Step 4: Build the objective function

````julia
sp = cov(1,1) + cov(1,2) + cov(1,3) +  # Cov(A₁, B₁) + Cov(A₁, B₂) + Cov(A₁, B₃)
     cov(2,1) + cov(2,2) - cov(2,3) +  # Cov(A₂, B₁) + Cov(A₂, B₂) - Cov(A₂, B₃)
     cov(3,1) - cov(3,2)               # Cov(A₃, B₁) - Cov(A₃, B₂)
````

````
-⟨x₁⟩⟨y₁⟩ - ⟨x₁⟩⟨y₂⟩ - ⟨x₁⟩⟨y₃⟩ - ⟨x₂⟩⟨y₁⟩ - ⟨x₂⟩⟨y₂⟩ + ⟨x₂⟩⟨y₃⟩ - ⟨x₃⟩⟨y₁⟩ + ⟨x₃⟩⟨y₂⟩ + ⟨x₁y₁⟩ + ⟨x₁y₂⟩ + ⟨x₁y₃⟩ + ⟨x₂y₁⟩ + ⟨x₂y₂⟩ - ⟨x₂y₃⟩ + ⟨x₃y₁⟩ - ⟨x₃y₂⟩
````

sp: state polynomial for covariance Bell inequality

#### Step 5: Create optimization problem and solve

````julia
spop = polyopt(sp, registry)
````

````
Optimization Problem (UnipotentAlgebra)
────────────────────────────────────
Objective:
    -⟨x₁⟩⟨y₁⟩ - ⟨x₁⟩⟨y₂⟩ - ⟨x₁⟩⟨y₃⟩ - ⟨x₂⟩⟨y₁⟩ - ⟨x₂⟩⟨y₂⟩ + ⟨x₂⟩⟨y₃⟩ - ⟨x₃⟩⟨y₁⟩ + ⟨x₃⟩⟨y₂⟩ + ⟨x₁y₁⟩ + ⟨x₁y₂⟩ + ⟨x₁y₃⟩ + ⟨x₂y₁⟩ + ⟨x₂y₂⟩ - ⟨x₂y₃⟩ + ⟨x₃y₁⟩ - ⟨x₃y₂⟩

Equality constraints (0):
    (none)

Inequality constraints (0):
    (none)

Moment equality constraints (0):
    (none)

Variables (6):
    x₁, x₂, x₃, y₁, y₂, y₃

````

spop: state polynomial optimization problem

````julia
solver_config = SolverConfig(
    optimizer = Mosek.Optimizer,
    order = 2
)

result = cs_nctssos(spop, solver_config)
cov_bound = -result.objective
````

````
5.000000000001876
````

cov_bound: upper bound on covariance Bell violation

````julia
cov_bound  # should be close to 5.0
````

````
5.000000000001876
````

Compare with known quantum value:

````julia
abs(cov_bound - 5.0)  # difference from theoretical value
````

````
1.8758328224066645e-12
````

#### Step 6: Improve bound using term sparsity and higher-order iteration

````julia
solver_config_ts = SolverConfig(
    optimizer = Mosek.Optimizer,
    order = 3,
    ts_algo = MF()  # term sparsity
)
````

````
NCTSSoS.SolverConfig(Mosek.Optimizer, 3, nothing, NCTSSoS.NoElimination(), CliqueTrees.MF(), nothing)
````

ts_algo=MF(): enables term sparsity exploitation

````julia
result_ts = cs_nctssos(spop, solver_config_ts)
````

````
State Optimization Result
Objective: -5.000000021561767
Correlative Sparsity (UnipotentAlgebra, Arbitrary): 

   maximum clique size: 6
   number of cliques: 1
   Clique 1: 
       Variables: [:x₁, :x₂, :x₃, :y₁, :y₂, :y₃]
       Bases length: 690
       Constraints: 
   Global Constraints: 
Term Sparsity:
Clique 1:
   Moment Matrix Block Sizes: [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 9, 6, 5, 5, 5, 8, 7, 6, 5, 5, 8, 5, 7, 7, 5, 5, 8, 5, 7, 7, 5, 5, 5, 8, 7, 6, 5, 7, 7, 7, 7, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 9, 6, 7, 5, 9, 9, 5, 5, 9, 8, 5, 5, 5, 8, 7, 8, 11, 5, 5, 5, 9, 8, 7, 5, 9, 9, 5, 5, 9, 8, 5, 5, 5, 8, 7, 8, 11, 5, 5, 5, 8, 7, 5, 5, 5, 8, 7, 5, 5, 9, 8, 8, 9, 5, 5, 9, 6, 5, 5, 5, 8, 7, 5, 5, 9, 8, 8, 9, 13, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 9, 8, 5, 5, 5, 8, 7, 5, 5, 5, 8, 7, 5, 5, 9, 8, 8, 9, 5, 5, 5, 8, 7, 5, 5, 9, 8, 8, 5, 5, 9, 8, 5, 5, 9, 8, 8, 11, 14, 14, 17, 15, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   Moment Matrix:
Number of Activated supp:   447
Number of Bases Activated in each sub-block[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 9, 6, 5, 5, 5, 8, 7, 6, 5, 5, 8, 5, 7, 7, 5, 5, 8, 5, 7, 7, 5, 5, 5, 8, 7, 6, 5, 7, 7, 7, 7, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 9, 6, 7, 5, 9, 9, 5, 5, 9, 8, 5, 5, 5, 8, 7, 8, 11, 5, 5, 5, 9, 8, 7, 5, 9, 9, 5, 5, 9, 8, 5, 5, 5, 8, 7, 8, 11, 5, 5, 5, 8, 7, 5, 5, 5, 8, 7, 5, 5, 9, 8, 8, 9, 5, 5, 9, 6, 5, 5, 5, 8, 7, 5, 5, 9, 8, 8, 9, 13, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 9, 8, 5, 5, 5, 8, 7, 5, 5, 5, 8, 7, 5, 5, 9, 8, 8, 9, 5, 5, 5, 8, 7, 5, 5, 9, 8, 8, 5, 5, 9, 8, 5, 5, 9, 8, 8, 11, 14, 14, 17, 15, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

   Localizing Matrix:
Unique Moment Matrix Elements: 375

````

result_ts: first iteration with term sparsity

````julia
result_higher = cs_nctssos_higher(spop, result_ts, solver_config_ts)
````

````
State Optimization Result
Objective: -5.0000000016276065
Correlative Sparsity (UnipotentAlgebra, Arbitrary): 

   maximum clique size: 6
   number of cliques: 1
   Clique 1: 
       Variables: [:x₁, :x₂, :x₃, :y₁, :y₂, :y₃]
       Bases length: 690
       Constraints: 
   Global Constraints: 
Term Sparsity:
Clique 1:
   Moment Matrix Block Sizes: [3, 3, 3, 3, 4, 7, 3, 3, 4, 7, 11, 11, 13, 11, 16, 3, 3, 4, 4, 3, 3, 7, 7, 11, 11, 4, 3, 3, 7, 7, 11, 17, 3, 3, 7, 7, 3, 4, 7, 11, 2, 3, 7, 10, 3, 3, 4, 3, 3, 7, 7, 11, 3, 4, 7, 11, 16, 16, 20, 21, 19, 2, 2, 2, 2, 3, 3, 3, 4, 5, 4, 5, 2, 3, 3, 4, 3, 4, 2, 4, 5, 4, 5, 6, 7, 11, 13, 3, 10, 3, 10, 3, 10, 4, 11, 9, 9, 10, 10, 10, 22, 4, 2, 7, 14, 21, 22, 19, 28, 31, 33, 3, 10, 3, 10, 3, 10, 4, 11, 9, 9, 10, 10, 10, 22, 4, 2, 7, 14, 21, 22, 19, 3, 9, 3, 10, 3, 10, 4, 11, 4, 11, 9, 9, 9, 10, 10, 21, 17, 20, 18, 33, 36, 54, 3, 10, 3, 10, 3, 10, 4, 11, 9, 9, 10, 10, 10, 22, 4, 2, 7, 14, 21, 22, 19, 3, 4, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 11, 15, 3, 4, 4, 5, 5, 7, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 11, 15, 22, 28, 35, 35, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 11, 15, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 11, 15, 3, 3, 3, 3, 4, 4, 5, 4, 4, 5, 6, 10, 13, 16, 28, 37, 3, 4, 5, 10, 4, 5, 4, 5, 6, 7, 3, 4, 5, 10, 13, 17, 3, 3, 3, 3, 4, 4, 5, 4, 5, 5, 6, 10, 13, 16, 30, 3, 9, 3, 10, 3, 10, 4, 11, 4, 11, 9, 9, 9, 10, 10, 21, 17, 20, 18, 3, 3, 3, 3, 4, 4, 5, 4, 5, 5, 6, 10, 13, 16, 29, 32, 48, 66, 2, 3, 3, 3, 4, 3, 4, 3, 4, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 3, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 11, 15, 22, 4, 5, 4, 5, 5, 7, 3, 4, 5, 10, 4, 5, 4, 4, 4, 7, 3, 4, 5, 10, 13, 17, 23, 25, 39, 3, 3, 4, 3, 4, 4, 5, 4, 5, 5, 7, 11, 13, 3, 4, 4, 5, 6, 7, 3, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 11, 15, 22, 3, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 11, 15, 3, 10, 3, 10, 3, 10, 4, 11, 9, 9, 10, 10, 10, 22, 4, 2, 7, 14, 21, 22, 19, 29, 35, 3, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 11, 15, 3, 3, 3, 3, 4, 4, 5, 4, 5, 5, 6, 10, 13, 16, 28, 49, 62, 70, 2]
   Moment Matrix:
Number of Activated supp:   3747
Number of Bases Activated in each sub-block[3, 3, 3, 3, 4, 7, 3, 3, 4, 7, 11, 11, 13, 11, 16, 3, 3, 4, 4, 3, 3, 7, 7, 11, 11, 4, 3, 3, 7, 7, 11, 17, 3, 3, 7, 7, 3, 4, 7, 11, 2, 3, 7, 10, 3, 3, 4, 3, 3, 7, 7, 11, 3, 4, 7, 11, 16, 16, 20, 21, 19, 2, 2, 2, 2, 3, 3, 3, 4, 5, 4, 5, 2, 3, 3, 4, 3, 4, 2, 4, 5, 4, 5, 6, 7, 11, 13, 3, 10, 3, 10, 3, 10, 4, 11, 9, 9, 10, 10, 10, 22, 4, 2, 7, 14, 21, 22, 19, 28, 31, 33, 3, 10, 3, 10, 3, 10, 4, 11, 9, 9, 10, 10, 10, 22, 4, 2, 7, 14, 21, 22, 19, 3, 9, 3, 10, 3, 10, 4, 11, 4, 11, 9, 9, 9, 10, 10, 21, 17, 20, 18, 33, 36, 54, 3, 10, 3, 10, 3, 10, 4, 11, 9, 9, 10, 10, 10, 22, 4, 2, 7, 14, 21, 22, 19, 3, 4, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 11, 15, 3, 4, 4, 5, 5, 7, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 11, 15, 22, 28, 35, 35, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 11, 15, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 11, 15, 3, 3, 3, 3, 4, 4, 5, 4, 4, 5, 6, 10, 13, 16, 28, 37, 3, 4, 5, 10, 4, 5, 4, 5, 6, 7, 3, 4, 5, 10, 13, 17, 3, 3, 3, 3, 4, 4, 5, 4, 5, 5, 6, 10, 13, 16, 30, 3, 9, 3, 10, 3, 10, 4, 11, 4, 11, 9, 9, 9, 10, 10, 21, 17, 20, 18, 3, 3, 3, 3, 4, 4, 5, 4, 5, 5, 6, 10, 13, 16, 29, 32, 48, 66, 2, 3, 3, 3, 4, 3, 4, 3, 4, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 7, 3, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 11, 15, 22, 4, 5, 4, 5, 5, 7, 3, 4, 5, 10, 4, 5, 4, 4, 4, 7, 3, 4, 5, 10, 13, 17, 23, 25, 39, 3, 3, 4, 3, 4, 4, 5, 4, 5, 5, 7, 11, 13, 3, 4, 4, 5, 6, 7, 3, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 11, 15, 22, 3, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 11, 15, 3, 10, 3, 10, 3, 10, 4, 11, 9, 9, 10, 10, 10, 22, 4, 2, 7, 14, 21, 22, 19, 29, 35, 3, 3, 4, 4, 5, 3, 4, 4, 5, 6, 7, 11, 15, 3, 3, 3, 3, 4, 4, 5, 4, 5, 5, 6, 10, 13, 16, 28, 49, 62, 70, 2]

   Localizing Matrix:
Unique Moment Matrix Elements: 2580

````

result_higher: higher-order iteration refining the bound

````julia
improved_bound = -result_higher.objective
````

````
5.0000000016276065
````

improved_bound: refined upper bound

````julia
(improved_bound,               # closer to 5.0
 abs(improved_bound - 5.0))    # very small difference from theoretical value
````

````
(5.0000000016276065, 1.6276064940257129e-9)
````

---
## Summary

| Inequality | Operator Type | Classical Bound | Quantum Bound | API |
|:-----------|:--------------|:----------------|:--------------|:----|
| CHSH | Unipotent (U²=I) | 2 | 2√2 ≈ 2.828 | `create_unipotent_variables` |
| I₃₃₂₂ | Projector (P²=P) | 0 | 0.25 | `create_projector_variables` |
| Covariance | Unipotent + State | 4.5 | 5 | `ς(·)` state polynomials |

Key functions:
- `create_unipotent_variables`: Create U² = I operators
- `create_projector_variables`: Create P² = P operators
- `polyopt`: Create optimization problem
- `SolverConfig`: Configure solver and sparsity options
- `cs_nctssos`: Solve using moment-SOS hierarchy
- `ς(·)`: Expectation value for state polynomials

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


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

We use `NCTSSoS.jl` for polynomial optimization and Mosek as the SDP solver backend.

````julia
using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true);
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
(VariableRegistry with 4 variables: A₁, A₂, B₁, B₂, (NormalMonomial{UnipotentAlgebra, UInt8}[UInt8[0x05], UInt8[0x09]], NormalMonomial{UnipotentAlgebra, UInt8}[UInt8[0x0e], UInt8[0x12]]))
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
(VariableRegistry with 4 variables: P₁, P₂, Q₁, Q₂, (NormalMonomial{ProjectorAlgebra, UInt8}[UInt8[0x05], UInt8[0x09]], NormalMonomial{ProjectorAlgebra, UInt8}[UInt8[0x0e], UInt8[0x12]]))
````

reg_proj: registry for projector algebra
P: Alice's projectors [P₁, P₂] on site 1
Q: Bob's projectors [Q₁, Q₂] on site 2

Verify the idempotent property (P² = P):

````julia
monomials(P[1] * P[1])  # should be [P[1]]
````

````
1-element Vector{NormalMonomial{ProjectorAlgebra, UInt8}}:
 UInt8[0x05]
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
(VariableRegistry with 4 variables: x₁, x₂, y₁, y₂, (NormalMonomial{UnipotentAlgebra, UInt8}[UInt8[0x05], UInt8[0x09]], NormalMonomial{UnipotentAlgebra, UInt8}[UInt8[0x0e], UInt8[0x12]]))
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
UInt8[0x05, 0x0e] + UInt8[0x05, 0x12] + UInt8[0x09, 0x0e] + -UInt8[0x09, 0x12]
````

f: polynomial representing the CHSH Bell operator

Inspect the polynomial structure:

````julia
(monomials(f),      # list of monomials in f
 coefficients(f))   # corresponding coefficients
````

````
(NormalMonomial{UnipotentAlgebra, UInt8}[UInt8[0x05, 0x0e], UInt8[0x05, 0x12], UInt8[0x09, 0x0e], UInt8[0x09, 0x12]], [1.0, 1.0, 1.0, -1.0])
````

#### Step 3: Create the optimization problem

````julia
pop = polyopt(f, registry)
````

````
Optimization Problem (UnipotentAlgebra)
────────────────────────────────────
Objective:
    x₁y₁ + x₁y₂ + x₂y₁ + -x₂y₂

Equality constraints (0):
    (none)

Inequality constraints (0):
    (none)

Variables (4):
    x₁, x₂, y₁, y₂

````

pop: polynomial optimization problem maximizing f
     subject to the algebraic constraints in registry (U² = I)

#### Step 4: Configure and run the SDP solver

````julia
solver_config = SolverConfig(
    optimizer = SILENT_MOSEK,  # SDP solver backend (silent mode)
    order = 1                    # relaxation order (hierarchy level)
)
````

````
SolverConfig(MathOptInterface.OptimizerWithAttributes(Mosek.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[MathOptInterface.Silent() => true]), 1, NoElimination(), NoElimination())
````

solver_config: specifies solver and relaxation parameters

````julia
result = cs_nctssos(pop, solver_config)
````

````
Objective: -2.82842713216232
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
chsh_bound = result.objective
````

````
-2.82842713216232
````

chsh_bound: upper bound on maximal quantum violation

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
5.65685425690851
````

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
(VariableRegistry with 6 variables: x₁, x₂, x₃, y₁, y₂, y₃, (NormalMonomial{ProjectorAlgebra, UInt8}[UInt8[0x05], UInt8[0x09], UInt8[0x0d]], NormalMonomial{ProjectorAlgebra, UInt8}[UInt8[0x12], UInt8[0x16], UInt8[0x1a]]))
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
-UInt8[0x05] + -2.0 * UInt8[0x12] + -UInt8[0x16] + UInt8[0x05, 0x12] + UInt8[0x05, 0x16] + UInt8[0x05, 0x1a] + UInt8[0x09, 0x12] + UInt8[0x09, 0x16] + -UInt8[0x09, 0x1a] + UInt8[0x0d, 0x12] + -UInt8[0x0d, 0x16]
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

Variables (6):
    x₁, x₂, x₃, y₁, y₂, y₃

````

pop: minimize -f (equivalent to maximize f)

````julia
solver_config = SolverConfig(optimizer=SILENT_MOSEK, order=2)
````

````
SolverConfig(MathOptInterface.OptimizerWithAttributes(Mosek.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[MathOptInterface.Silent() => true]), 2, NoElimination(), NoElimination())
````

order=2: second level of the moment hierarchy

````julia
result = cs_nctssos(pop, solver_config)
i3322_bound = -result.objective
````

````
0.2509397976352535
````

i3322_bound: upper bound on I₃₃₂₂ violation (negate since we minimized -f)

````julia
i3322_bound  # should be close to 0.25
````

````
0.2509397976352535
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

solver_config_dense = SolverConfig(optimizer=SILENT_MOSEK, order=3)
````

````
SolverConfig(MathOptInterface.OptimizerWithAttributes(Mosek.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[MathOptInterface.Silent() => true]), 3, NoElimination(), NoElimination())
````

solver_config_dense: no sparsity exploitation

````julia
@time result_dense = cs_nctssos(pop, solver_config_dense)
bound_dense = -result_dense.objective
````

````
0.2508757549955246
````

bound_dense: bound without sparsity

````julia
bound_dense
````

````
0.2508757549955246
````

#### With correlative sparsity (order=6)

````julia
solver_config_sparse = SolverConfig(
    optimizer = SILENT_MOSEK,
    order = 6,             # higher order for better bound
    cs_algo = MF()         # use MaxFlow algorithm for correlative sparsity
)
````

````
SolverConfig(MathOptInterface.OptimizerWithAttributes(Mosek.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[MathOptInterface.Silent() => true]), 6, MF(), NoElimination())
````

cs_algo=MF(): enables correlative sparsity via chordal graph decomposition

````julia
@time result_sparse = cs_nctssos(pop, solver_config_sparse)
bound_sparse = -result_sparse.objective
````

````
0.2508754080902566
````

bound_sparse: improved bound using sparsity

````julia
bound_sparse  # closer to theoretical 0.25
````

````
0.2508754080902566
````

Improvement in bound:

````julia
bound_dense - bound_sparse  # positive = improvement
````

````
3.4690526801162136e-7
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
(VariableRegistry with 6 variables: x₁, x₂, x₃, y₁, y₂, y₃, (NormalMonomial{UnipotentAlgebra, UInt8}[UInt8[0x05], UInt8[0x09], UInt8[0x0d]], NormalMonomial{UnipotentAlgebra, UInt8}[UInt8[0x12], UInt8[0x16], UInt8[0x1a]]))
````

x: Alice's observables [A₁, A₂, A₃]
y: Bob's observables [B₁, B₂, B₃]

#### Step 2: Define the identity monomial

````julia
ID = one(NormalMonomial{UnipotentAlgebra, UInt8})
````

````
𝟙
````

ID: identity element (𝟙) needed for state polynomial arithmetic

````julia
ID  # display the identity
````

````
𝟙
````

#### Step 3: Define the covariance function using state polynomials

State polynomials use `ς(·)` to denote expectation values ⟨·⟩.

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
-⟨UInt8[0x05]⟩⟨UInt8[0x12]⟩ + ⟨UInt8[0x05, 0x12]⟩
````

#### Step 4: Build the objective function

````julia
sp = cov(1,1) + cov(1,2) + cov(1,3) +  # Cov(A₁, B₁) + Cov(A₁, B₂) + Cov(A₁, B₃)
     cov(2,1) + cov(2,2) - cov(2,3) +  # Cov(A₂, B₁) + Cov(A₂, B₂) - Cov(A₂, B₃)
     cov(3,1) - cov(3,2)               # Cov(A₃, B₁) - Cov(A₃, B₂)
````

````
-⟨UInt8[0x05]⟩⟨UInt8[0x12]⟩ - ⟨UInt8[0x05]⟩⟨UInt8[0x16]⟩ - ⟨UInt8[0x05]⟩⟨UInt8[0x1a]⟩ - ⟨UInt8[0x09]⟩⟨UInt8[0x12]⟩ - ⟨UInt8[0x09]⟩⟨UInt8[0x16]⟩ + ⟨UInt8[0x09]⟩⟨UInt8[0x1a]⟩ - ⟨UInt8[0x0d]⟩⟨UInt8[0x12]⟩ + ⟨UInt8[0x0d]⟩⟨UInt8[0x16]⟩ + ⟨UInt8[0x05, 0x12]⟩ + ⟨UInt8[0x05, 0x16]⟩ + ⟨UInt8[0x05, 0x1a]⟩ + ⟨UInt8[0x09, 0x12]⟩ + ⟨UInt8[0x09, 0x16]⟩ - ⟨UInt8[0x09, 0x1a]⟩ + ⟨UInt8[0x0d, 0x12]⟩ - ⟨UInt8[0x0d, 0x16]⟩
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
    -⟨UInt8[0x05]⟩⟨UInt8[0x12]⟩ - ⟨UInt8[0x05]⟩⟨UInt8[0x16]⟩ - ⟨UInt8[0x05]⟩⟨UInt8[0x1a]⟩ - ⟨UInt8[0x09]⟩⟨UInt8[0x12]⟩ - ⟨UInt8[0x09]⟩⟨UInt8[0x16]⟩ + ⟨UInt8[0x09]⟩⟨UInt8[0x1a]⟩ - ⟨UInt8[0x0d]⟩⟨UInt8[0x12]⟩ + ⟨UInt8[0x0d]⟩⟨UInt8[0x16]⟩ + ⟨UInt8[0x05, 0x12]⟩ + ⟨UInt8[0x05, 0x16]⟩ + ⟨UInt8[0x05, 0x1a]⟩ + ⟨UInt8[0x09, 0x12]⟩ + ⟨UInt8[0x09, 0x16]⟩ - ⟨UInt8[0x09, 0x1a]⟩ + ⟨UInt8[0x0d, 0x12]⟩ - ⟨UInt8[0x0d, 0x16]⟩

Equality constraints (0):
    (none)

Inequality constraints (0):
    (none)

Variables (6):
    x₁, x₂, x₃, y₁, y₂, y₃

````

spop: state polynomial optimization problem

````julia
solver_config = SolverConfig(
    optimizer = SILENT_MOSEK,
    order = 2
)

result = cs_nctssos(spop, solver_config)
cov_bound = -result.objective
````

````
4.999999999618061
````

cov_bound: upper bound on covariance Bell violation

````julia
cov_bound  # should be close to 5.0
````

````
4.999999999618061
````

Compare with known quantum value:

````julia
abs(cov_bound - 5.0)  # difference from theoretical value
````

````
3.8193892493154635e-10
````

#### Step 6: Improve bound using term sparsity and higher-order iteration

````julia
solver_config_ts = SolverConfig(
    optimizer = SILENT_MOSEK,
    order = 3,
    ts_algo = MF()  # term sparsity
)
````

````
SolverConfig(MathOptInterface.OptimizerWithAttributes(Mosek.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[MathOptInterface.Silent() => true]), 3, NoElimination(), MF())
````

ts_algo=MF(): enables term sparsity exploitation

````julia
result_ts = cs_nctssos(spop, solver_config_ts)
````

````
State Optimization Result
Objective: -4.999999997729253
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
Objective: -4.999999998288536
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
   Moment Matrix Block Sizes: [3, 3, 3, 3, 3, 3, 3, 4, 3, 11, 18, 4, 4, 4, 11, 4, 11, 13, 14, 13, 15, 22, 28, 4, 4, 3, 3, 4, 4, 4, 11, 4, 11, 13, 13, 4, 4, 11, 4, 4, 11, 22, 24, 4, 4, 4, 4, 11, 4, 11, 13, 15, 26, 30, 32, 29, 28, 3, 2, 3, 3, 3, 3, 2, 2, 3, 3, 7, 2, 2, 3, 4, 3, 4, 3, 3, 3, 3, 3, 3, 3, 7, 3, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 7, 4, 6, 3, 3, 4, 2, 3, 4, 2, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 10, 10, 10, 10, 10, 24, 22, 23, 30, 4, 6, 3, 4, 4, 5, 6, 7, 15, 22, 4, 6, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 6, 7, 8, 15, 23, 32, 45, 49, 49, 55, 3, 4, 6, 3, 4, 6, 4, 6, 4, 6, 8, 10, 4, 6, 4, 6, 7, 8, 16, 4, 6, 4, 6, 7, 8, 4, 6, 4, 6, 7, 8, 18, 22, 28, 29, 34, 41, 51, 76, 73, 73, 73, 76, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 9, 10, 10, 10, 10, 23, 22, 20, 27, 2, 7, 3, 4, 4, 6, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 6, 7, 8, 15, 23, 32, 3, 4, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 5, 5, 7, 15, 22, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 10, 10, 10, 10, 10, 24, 22, 23, 30, 41, 52, 50, 57, 64, 72, 2, 2, 3, 4, 2, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 9, 10, 10, 10, 10, 23, 22, 20, 27, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 10, 10, 10, 10, 10, 24, 22, 23, 30, 4, 6, 3, 4, 4, 5, 6, 7, 15, 22, 4, 6, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 6, 7, 8, 15, 23, 32, 45, 49, 49, 61, 64, 2, 3, 3, 4, 4, 6, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 5, 6, 7, 15, 22, 28, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 10, 10, 10, 10, 10, 24, 22, 23, 29, 36, 44, 45, 3, 4, 6, 3, 4, 6, 4, 6, 4, 6, 8, 10, 19, 4, 6, 4, 6, 7, 8, 16, 4, 6, 4, 6, 7, 8, 4, 6, 4, 6, 7, 8, 18, 22, 28, 35, 44, 59, 75, 72, 72, 72, 93, 93, 4]
   Moment Matrix:
Number of Activated supp:   4977
Number of Bases Activated in each sub-block[3, 3, 3, 3, 3, 3, 3, 4, 3, 11, 18, 4, 4, 4, 11, 4, 11, 13, 14, 13, 15, 22, 28, 4, 4, 3, 3, 4, 4, 4, 11, 4, 11, 13, 13, 4, 4, 11, 4, 4, 11, 22, 24, 4, 4, 4, 4, 11, 4, 11, 13, 15, 26, 30, 32, 29, 28, 3, 2, 3, 3, 3, 3, 2, 2, 3, 3, 7, 2, 2, 3, 4, 3, 4, 3, 3, 3, 3, 3, 3, 3, 7, 3, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 7, 4, 6, 3, 3, 4, 2, 3, 4, 2, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 10, 10, 10, 10, 10, 24, 22, 23, 30, 4, 6, 3, 4, 4, 5, 6, 7, 15, 22, 4, 6, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 6, 7, 8, 15, 23, 32, 45, 49, 49, 55, 3, 4, 6, 3, 4, 6, 4, 6, 4, 6, 8, 10, 4, 6, 4, 6, 7, 8, 16, 4, 6, 4, 6, 7, 8, 4, 6, 4, 6, 7, 8, 18, 22, 28, 29, 34, 41, 51, 76, 73, 73, 73, 76, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 9, 10, 10, 10, 10, 23, 22, 20, 27, 2, 7, 3, 4, 4, 6, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 6, 7, 8, 15, 23, 32, 3, 4, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 5, 5, 7, 15, 22, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 10, 10, 10, 10, 10, 24, 22, 23, 30, 41, 52, 50, 57, 64, 72, 2, 2, 3, 4, 2, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 9, 10, 10, 10, 10, 23, 22, 20, 27, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 10, 10, 10, 10, 10, 24, 22, 23, 30, 4, 6, 3, 4, 4, 5, 6, 7, 15, 22, 4, 6, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 6, 7, 8, 15, 23, 32, 45, 49, 49, 61, 64, 2, 3, 3, 4, 4, 6, 3, 4, 4, 6, 7, 8, 15, 4, 6, 3, 4, 4, 5, 6, 7, 15, 22, 28, 3, 10, 3, 10, 3, 10, 4, 12, 4, 12, 10, 10, 10, 10, 10, 24, 22, 23, 29, 36, 44, 45, 3, 4, 6, 3, 4, 6, 4, 6, 4, 6, 8, 10, 19, 4, 6, 4, 6, 7, 8, 16, 4, 6, 4, 6, 7, 8, 4, 6, 4, 6, 7, 8, 18, 22, 28, 35, 44, 59, 75, 72, 72, 72, 93, 93, 4]

   Localizing Matrix:
Unique Moment Matrix Elements: 3184

````

result_higher: higher-order iteration refining the bound

````julia
improved_bound = -result_higher.objective
````

````
4.999999998288536
````

improved_bound: refined upper bound

````julia
(improved_bound,               # closer to 5.0
 abs(improved_bound - 5.0))    # very small difference from theoretical value
````

````
(4.999999998288536, 1.7114638595217002e-9)
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


# Bell inequalities

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

We use `NCTSSoS.jl` for polynomial optimization and `COSMO` as the SDP solver backend.

````julia
using NCTSSoS, COSMO
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

reg_unip: registry storing variable names and algebra type
A: Alice's measurement operators [A₁, A₂] on site 1
B: Bob's measurement operators [B₁, B₂] on site 2

Verify the unipotent property (U² = I):

````julia
A[1] * A[1]  # should simplify to identity
````

Check that operators on different sites commute:

````julia
A[1] * B[1] == B[1] * A[1]  # true: different sites commute
````

### Projector Variables (P² = P)

Create operators that are idempotent:

````julia
reg_proj, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 1:2)])
````

reg_proj: registry for projector algebra
P: Alice's projectors [P₁, P₂] on site 1
Q: Bob's projectors [Q₁, Q₂] on site 2

Verify the idempotent property (P² = P):

````julia
monomials(P[1] * P[1])  # should be [P[1]]
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

f: polynomial representing the CHSH Bell operator

Inspect the polynomial structure:

````julia
(monomials(f),      # list of monomials in f
 coefficients(f))   # corresponding coefficients
````

#### Step 3: Create the optimization problem

````julia
pop = polyopt(f, registry)
````

pop: polynomial optimization problem maximizing f
     subject to the algebraic constraints in registry (U² = I)

#### Step 4: Configure and run the SDP solver

````julia
solver_config = SolverConfig(
    optimizer = Mosek.Optimizer,  # SDP solver backend
    order = 1                      # relaxation order (hierarchy level)
)
````

solver_config: specifies solver and relaxation parameters

````julia
result = cs_nctssos(pop, solver_config)
````

result: optimization result containing objective value and solver info

#### Step 5: Extract the upper bound

````julia
chsh_bound = result.objective
````

chsh_bound: upper bound on maximal quantum violation

Compare with Tsirelson's bound:

````julia
tsirelson_bound = 2 * sqrt(2)
````

tsirelson_bound: theoretical maximum = 2√2 ≈ 2.828

````julia
abs(chsh_bound - tsirelson_bound)  # difference (should be ~1e-7)
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

f: I₃₃₂₂ Bell polynomial

Check the number of terms:

````julia
length(monomials(f))  # number of monomials
````

#### Step 3: Solve (minimizing -f to find maximum of f)

````julia
pop = polyopt(-f, registry)
````

pop: minimize -f (equivalent to maximize f)

````julia
solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=2)
````

order=2: second level of the moment hierarchy

````julia
result = cs_nctssos(pop, solver_config)
i3322_bound = -result.objective
````

i3322_bound: upper bound on I₃₃₂₂ violation (negate since we minimized -f)

````julia
i3322_bound  # should be close to 0.25
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

solver_config_dense: no sparsity exploitation

````julia
@time result_dense = cs_nctssos(pop, solver_config_dense)
bound_dense = -result_dense.objective
````

bound_dense: bound without sparsity

````julia
bound_dense
````

#### With correlative sparsity (order=6)

````julia
solver_config_sparse = SolverConfig(
    optimizer = Mosek.Optimizer,
    order = 6,             # higher order for better bound
    cs_algo = MF()         # use MaxFlow algorithm for correlative sparsity
)
````

cs_algo=MF(): enables correlative sparsity via chordal graph decomposition

````julia
@time result_sparse = cs_nctssos(pop, solver_config_sparse)
bound_sparse = -result_sparse.objective
````

bound_sparse: improved bound using sparsity

````julia
bound_sparse  # closer to theoretical 0.25
````

Improvement in bound:

````julia
bound_dense - bound_sparse  # positive = improvement
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

x: Alice's observables [A₁, A₂, A₃]
y: Bob's observables [B₁, B₂, B₃]

#### Step 2: Define the identity monomial

````julia
ID = one(NormalMonomial{UnipotentAlgebra, UInt8})
````

ID: identity element (𝟙) needed for state polynomial arithmetic

````julia
ID  # display the identity
````

#### Step 3: Define the covariance function using state polynomials

State polynomials use `ς(·)` to denote expectation values ⟨·⟩.

````julia
cov(a, b) = 1.0 * ς(x[a] * y[b]) * ID -  # ⟨AᵢBⱼ⟩
            1.0 * ς(x[a]) * ς(y[b]) * ID  # -⟨Aᵢ⟩⟨Bⱼ⟩
````

cov(a,b): covariance Cov(Aₐ, Bᵦ) as a state polynomial
ς (varsigma): expectation value operator, type \varsigma + Tab

Example: Cov(A₁, B₁)

````julia
cov(1, 1)
````

#### Step 4: Build the objective function

````julia
sp = cov(1,1) + cov(1,2) + cov(1,3) +  # Cov(A₁, B₁) + Cov(A₁, B₂) + Cov(A₁, B₃)
     cov(2,1) + cov(2,2) - cov(2,3) +  # Cov(A₂, B₁) + Cov(A₂, B₂) - Cov(A₂, B₃)
     cov(3,1) - cov(3,2)               # Cov(A₃, B₁) - Cov(A₃, B₂)
````

sp: state polynomial for covariance Bell inequality

#### Step 5: Create optimization problem and solve

````julia
spop = polyopt(sp, registry)
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

cov_bound: upper bound on covariance Bell violation

````julia
cov_bound  # should be close to 5.0
````

Compare with known quantum value:

````julia
abs(cov_bound - 5.0)  # difference from theoretical value
````

#### Step 6: Improve bound using term sparsity and higher-order iteration

````julia
solver_config_ts = SolverConfig(
    optimizer = Mosek.Optimizer,
    order = 3,
    ts_algo = MF()  # term sparsity
)
````

ts_algo=MF(): enables term sparsity exploitation

````julia
result_ts = cs_nctssos(spop, solver_config_ts)
````

result_ts: first iteration with term sparsity

````julia
result_higher = cs_nctssos_higher(spop, result_ts, solver_config_ts)
````

result_higher: higher-order iteration refining the bound

````julia
improved_bound = -result_higher.objective
````

improved_bound: refined upper bound

````julia
(improved_bound,               # closer to 5.0
 abs(improved_bound - 5.0))    # very small difference from theoretical value
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


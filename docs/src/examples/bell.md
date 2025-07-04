# [Bell inequalities](@id bell-inequalities)

Bell inequalities are mathematical expressions that test whether the predictions of quantum mechanics can be explained by local hidden variable theories. They were first introduced by John Stewart Bell in 1964 and have since become fundamental tools in quantum information theory and quantum foundations.
A Bell inequality is typically expressed as a linear combination of expectation values of observables, with bounds that differ between classical and quantum theories. In the classical case, these inequalities must be satisfied if the system can be described by local hidden variables. However, quantum mechanics can violate these inequalities, demonstrating the non-local nature of quantum correlations.
The general form of a Bell inequality can be written as:

$$\sum_{i,j} c_{ij} \langle A_i B_j \rangle \leq C$$

where $A_i$ and $B_j$ are observables measured by two parties (traditionally called Alice and Bob), $c_{ij}$ are real coefficients, and $C$ is the classical bound. Quantum mechanics can violate this inequality, with the maximum violation known as the quantum bound.

## Linear Bell inequalities

### CHSH inequality
The most famous Bell inequality is the CHSH (Clauser-Horne-Shimony-Holt) inequality, which involves two parties, each measuring two observables. For unipotent (square to 1) observables $A_1, A_2$ measured by Alice and $B_1, B_2$ measured by Bob. We define the objective function as:

$$f(A_1, A_2, B_1, B_2) = \langle A_1B_1 \rangle + \langle A_1B_2 \rangle + \langle A_2B_1 \rangle - \langle A_2B_2 \rangle$$

The CHSH inequality is then given by $$f(A_1, A_2, B_1, B_2) \leq 2$$, which must be satisfied by any local hidden variable theory. However, quantum mechanics can violate this inequality up to the value $2\sqrt{2}$, known as the Tsirelson bound. This violation demonstrates that quantum mechanics cannot be described by any local hidden variable theory.
The CHSH inequality is particularly important because it is the simplest non-trivial Bell inequality and has been experimentally verified numerous times, providing strong evidence for the non-local nature of quantum mechanics.

The upper bound of the CHSH inequality can be computed using the following code:

```julia chsh
using NCTSSoS, MosekTools

@ncpolyvar x[1:2]  # x = (A_1, A_2)
@ncpolyvar y[1:2]  # y = (B_1, B_2)
f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]  # objective function

pop = polyopt(             # the optimization problem
        f;
        comm_gps= [x, y],   # commutative group
        is_unipotent=true  # the variables are unipotent
    )

solver_config = SolverConfig(;
    optimizer=Mosek.Optimizer,  # the solver backend
    order=1                    # the order of the moment matrix
)
result = cs_nctssos(pop, solver_config)
result.objective  # the upper bound of the CHSH inequality
```

```julia
-2.8284271321623202
```

Here, we first declare some operators as non-commutative variables, and then construct the optimization problem. In `PolyOpt` constructor,
- `comm_gps` argument specifies the commutative group of the variables, which means variables in different commutative groups commute with each other.
- `is_unipotent` argument specifies that the variables are unipotent, which means they square to 1 (e.g. Pauli operators).

Here, since the variables on different qubits commute with each other, we can group them into different commutative groups.

In the solver configuration, we use [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl) as the semidefinite programming solver backend. It is an open-source solver for conic programs with quadratic objectives, and it uses the interior point method to solve the problem [Clarabel_2024](@cite).

The resulting upper bound is very close to the theoretical exact value $2\sqrt{2} \approx 2.8284271247461903$ (accurate up to 7 decimals!!).

### $I_{3322}$ inequality

The $I_{3322}$ inequality is a more complex inequality that involves three parties, each measuring three observables. Let $A_1, A_2, A_3$ be the projective (square to itself) observables measured by Alice and $B_1, B_2, B_3$ be the projective observables measured by Bob. We define the objective function as [pal2010maximal](@cite):

```math
f(A_1, A_2, A_3, B_1, B_2, B_3) = \langle A_1(B_1+B_2+B_3) \rangle + \langle A_2(B_1+B_2-B_3) \rangle\\
+ \langle A_3(B_1-B_2) \rangle
- \langle A_1 \rangle - 2\langle B_1 \rangle - \langle B_2 \rangle
```

In classical mechanics, the inequality $f(A_1, A_2, A_3, B_1, B_2, B_3) \leq 0$ must be satisfied. However, quantum mechanics can violate this inequality up to the value $0.25$, known as the Tsirelson bound. This violation demonstrates that quantum mechanics cannot be described by any local hidden variable theory.

The upper bound of the $I_{3322}$ inequality can be computed using the following code:

```julia i3322
using NCTSSoS, MosekTools

@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]  # objective function

pop = polyopt(-f; comm_gps= [x, y], is_projective=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer; order=2)

result = cs_nctssos(pop, solver_config)
result.objective
```

```julia
-0.25093972222278366
```

Here, the `is_projective` argument specifies that the variables are projective, which means they square to themselves (e.g. $|0\rangle\langle 0|$ and $|1\rangle\langle 1|$).

The resulting upper bound is close to the theoretical exact value $0.25$. By increasing the order of the moment matrix, this upper bound can be improved.

#### Reducing SDP Problem Size with Sparsity

To reach the theoretical exact value of $0.25$, we can increase the order of the moment matrix [magronSparsePolynomialOptimization2023](@cite).

```julia
using NCTSSoS, MosekTools

@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]  # objective function

pop = polyopt(-f; comm_gps= [x,y], is_projective=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer; order=3)

@time result = cs_nctssos(pop, solver_config)
@show result.objective
```

```julia
1.923037 seconds (43.04 M allocations: 2.210 GiB, 20.41% gc time)
Objective: -0.2508755502587585
```

Indeed, by increasing the order of the moment matrix to 3, have improved the lower bound from $-0.25093972222278366$ to $-0.2508755502587585$.

However, keep increase the order can lead to a large semidefinite programming (SDP) problem size, which can be computationally expensive. To reduce the problem size, we may exploit the sparsity of the problem [magronSparsePolynomialOptimization2023](@cite). There are two sparsity patterns that can be used to reduce the problem size:

1. **Correlation Sparsity**: exploits the fact that few variable products exists in the objective function. Therefore, we could break down the objective function into smaller parts, each involving fewer variables. This reduces the moment matrix size and the number of constraints in the SDP problem, making it more tractable.

2. **Term Sparsity**: exploits the fact that not all monomials in the moment matrix are needed to represent the objective function. By identifying and removing unnecessary monomials, we can further reduce the size of the moment matrix and the SDP problem.

To take advantage of these sparsity patterns,

```julia i3322_sparsity
using NCTSSoS, MosekTools

@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]  # objective function

pop = polyopt(-f; comm_gps= [x,y], is_projective=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer; order=6, cs_algo=MF())

@time result = cs_nctssos(pop, solver_config)
@show result.objective
```

```julia
46.996790 seconds (14.14 M allocations: 1.579 GiB, 0.90% gc time, 0.16% compilation time)
-0.2508753195677618
```

Using almost half of the time, we are able to improve the $7$-th digit of the upper bound!

## Nonlinear Bell Inequalities

Non-linear Bell inequalities are extensions of the standard linear Bell inequalities. Instead of being linear combinations of expectation values, they involve polynomial functions of these expectation values. These inequalities arise naturally when considering more complex scenarios, such as multi-party settings or when the parties can perform sequences of measurements.

The significance of non-linear Bell inequalities in quantum information lies in their ability to detect non-locality in situations where linear inequalities might fail. They can provide tighter bounds on classical correlations and reveal quantum non-locality in a broader range of experimental setups. Furthermore, studying non-linear Bell inequalities helps in understanding the structure of quantum correlations and the boundary between classical and quantum physics more deeply. They are also relevant in the context of quantum cryptography and communication complexity, where understanding the limits of classical and quantum correlations is crucial.

### Covariance Bell Inequality

The covariance Bell inequality is a non-linear Bell inequality that involves the covariance of measurements. It can be expressed as:

$$\text{Cov}(A, B) = \langle A B \rangle - \langle A \rangle \langle B \rangle$$

where $A$ and $B$ are observables measured by two parties. Comparing with the linear Bell inequality, the covariance Bell inequality is non-linear because it involves the product of two observables.

Let us define the objective function as:
```math
f(A_1,A_2,A_3, B_1,B_2,B_3) = \text{Cov}(A_1, B_1) + \text{Cov}(A_1, B_2) + \text{Cov}(A_1,B_3)  + \\ \text{Cov}(A_2, B_1) + \text{Cov}(A_2, B_2) - \text{Cov}(A_2, B_3) + \text{Cov}(A_3, B_1) - \text{Cov}(A_3,B_2)
```

it was shown that $f(A_1,A_2,A_3,B_1,B_2,B_3) \leq \frac{9}{2}$ in classical models, while it can attain a maximum value of $5$ in spatial quantum model of qubits and a maximally entangled state [pozsgay2017Covariance](@cite).

An *open question* was whether a higher bound can be attained in a spatial quantum model of qudits, i.e., systems with more than two levels. Using State Polynomial Optimization [klep2024State](@cite) , we can certify the upper bound of this inequality:

```julia covariance
using NCTSSoS, MosekTools, NCTSSoS.FastPolynomials

@ncpolyvar x[1:3] y[1:3]  # x = (A_1, A_2, A_3), y = (B_1, B_2, B_3)

# covariance function
cov(a, b) = 1.0 * ς(x[a] * y[b]) * one(Monomial) -
            1.0 * ς(x[a]) * ς(y[b]) * one(Monomial)

# objective function
sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)


spop = polyopt(
        sp;                                 # the optimization problem
        is_unipotent=true,                  # the variables are unipotent
        comm_gps=[x[1:3], y[1:3]]           # the commutative groups of the variables
        )

solver_config = SolverConfig(
    optimizer=Mosek.Optimizer;           # the solver backend
    order=2                             # the order of the moment matrix
)

result = cs_nctssos(spop, solver_config)
result
```

```julia
Objective: -5.000271541108556
```

!!! note "Typing Unicodes"
    You can type the unicode characters in the code by using `\varsigma` and pressing `Tab` to get the unicode character `ς`.

The resulting upper bound is very close to the previously known best value of $5$ (accurate up to 3 decimals!!). It accertains the value of $5$ for any system size.

We can use sparsity to improve the performance of the algorithm.

```julia
using NCTSSoS, MosekTools, NCTSSoS.FastPolynomials

@ncpolyvar x[1:3] y[1:3]  # x = (A_1, A_2, A_3), y = (B_1, B_2, B_3)

# covariance function
cov(a, b) = 1.0 * ς(x[a] * y[b]) * one(Monomial) -
            1.0 * ς(x[a]) * ς(y[b]) * one(Monomial)

# objective function
sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)


spop = polyopt(
        sp;                                 # the optimization problem
        is_unipotent=true,                  # the variables are unipotent
        comm_gps=[x[1:3], y[1:3]]           # the commutative groups of the variables
        )

solver_config = SolverConfig(
    optimizer=Mosek.Optimizer;              # the solver backend
    order=3,                            # the order of the moment matrix
    cs_algo = MF()                         # term sparse algorithm
)

result = cs_nctssos(spop, solver_config)

result_higher = cs_nctssos_higher(spop, result,solver_config)
result_higher
```

```julia
Objective: -4.999999981821947
```

This is accurate up to $10$ decimals.

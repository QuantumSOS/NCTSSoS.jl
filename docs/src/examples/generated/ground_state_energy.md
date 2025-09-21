```@meta
EditURL = "../literate/ground_state_energy.jl"
```

# Obtaining Ground State Energy Lower Bound

Finding the ground state of a quantum system is a fundamental problem in quantum
mechanics [wang2024Certifying](@cite). Variational methods are commonly used to
approximate the ground state. Due to the variational nature of these methods,
only an upper bound can be obtained [kull2024Lower](@cite). Polynomial
optimization techniques provides a way to find the lower bound of the ground
state energy.

In general, we consider the following Hamiltonian:
```math
H = \frac{1}{4} \sum_{i \u003c j} J_{ij} \sum_{ a \in \{x,y,z\}} \sigma_i^a \sigma_j^a
```

## 1D Heisenberg Model with Nearest Neighbor Interaction

Firstly, let's consider the simplest case of 1D Heisenberg chain with nearest
neighbor interaction and periodic boundary condition.

```julia

````@example ground_state_energy
using NCTSSoS, MosekTools
N = 6
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = cpolyopt(
            ham;                                        # the Hamiltonian
            eq_constraints=eq_cons,                     # anti-commutation relation between Pauli Operators
            comm_gps=[[x[i], y[i], z[i]] for i in 1:N], # commutation relation between Pauli Operators
            is_unipotent=true                           # Pauli operators square to identity
            )

solver_config = SolverConfig(
                    optimizer=Mosek.Optimizer,          # the solver backend
                    order=2,                        # moment matrix order
                    ts_algo = MMD(),                    # term sparsity algorithm
                    )

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(
            pop,                                        # Polynomial Optimization Problem
            res,                                        # Solution of First Order Term Sparsity Iteration
            solver_config                               # Solver Configuration
        )
res.objective / N
````

```

The returned result matches the actual ground state energy $-0.467129$ to $6$
digits. [wang2024Certifying](@cite)

## 1D Heisenberg Model with next nearest neighbor interaction

Polynomial Optimization framework is quite general. Almost no modification is
required to handle more complex Hamiltonian. 1D Heisenberg Model with geometric
frustration induced by next nearest neighbor interaction can be solved as:

```julia

````@example ground_state_energy
using NCTSSoS, MosekTools
N = 6
J1 = 1.0                            # Nearest Neighbor Interaction
J2 = 0.2                            # Next Nearest Neighbor Interaction
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] + ComplexF64(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=2, ts_algo = MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)
res.objective/N
````

```

We are able to obtain the ground state energy of $-0.4270083225302217$, accurate
to $6$ digits!

## 2D Square Lattice

Extending Heisenberg model to $2$-D case is also straightforward. However `NCTSSoS.jl` is not efficient enough to handle system at this size.

```julia

````@example ground_state_energy
using NCTSSoS, MosekTools
Nx = 3
Ny = 3
N = Nx * Ny
J1 = 1.0
J2 = 0.0
@ncpolyvar x[1:N] y[1:N] z[1:N]

LI = LinearIndices((1:Nx, 1:Ny))

ham = sum(ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(i, mod1(j + 1, Ny))]] + ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), j)]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j + 1, Ny))]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j - 1, Ny))]] for op in [x, y, z] for i in 1:Nx for j in 1:Ny)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=3, cs_algo=MF(), ts_algo=MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)

res.objective / N
````

```

## Next step

With such lower bounds, estimates of properties of the ground
state, correlations functions, structure factors and order parameters, can also
be obtained. We provide examples in another [section](@ref certify-property).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


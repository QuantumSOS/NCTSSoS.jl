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

````julia
using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)
N = 6
````

````
6
````

Create Pauli variables using the typed algebra system
This automatically encodes all Pauli commutation relations

````julia
registry, (σx, σy, σz) = create_pauli_variables(1:N)

ham = sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [σx, σy, σz] for i in 1:N)
````

````
0.25 + 0.0im * UInt8[0x01, 0x04] + 0.25 + 0.0im * UInt8[0x01, 0x10] + 0.25 + 0.0im * UInt8[0x02, 0x05] + 0.25 + 0.0im * UInt8[0x02, 0x11] + 0.25 + 0.0im * UInt8[0x03, 0x06] + 0.25 + 0.0im * UInt8[0x03, 0x12] + 0.25 + 0.0im * UInt8[0x04, 0x07] + 0.25 + 0.0im * UInt8[0x05, 0x08] + 0.25 + 0.0im * UInt8[0x06, 0x09] + 0.25 + 0.0im * UInt8[0x07, 0x0a] + 0.25 + 0.0im * UInt8[0x08, 0x0b] + 0.25 + 0.0im * UInt8[0x09, 0x0c] + 0.25 + 0.0im * UInt8[0x0a, 0x0d] + 0.25 + 0.0im * UInt8[0x0b, 0x0e] + 0.25 + 0.0im * UInt8[0x0c, 0x0f] + 0.25 + 0.0im * UInt8[0x0d, 0x10] + 0.25 + 0.0im * UInt8[0x0e, 0x11] + 0.25 + 0.0im * UInt8[0x0f, 0x12]
````

No need to manually specify constraints - they're encoded in the algebra type!

````julia
pop = polyopt(ham, registry)

solver_config = SolverConfig(
                    optimizer=SILENT_MOSEK,          # the solver backend
                    order=2,                            # moment matrix order
                    ts_algo = MMD(),                    # term sparsity algorithm
                    )

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(
            pop,                                        # Polynomial Optimization Problem
            res,                                        # Solution of First Order Term Sparsity Iteration
            solver_config                               # Solver Configuration
        )
energy_per_site = res.objective / N
@show energy_per_site
````

````
-0.4671292719994368
````

The returned result matches the actual ground state energy $-0.467129$ to $6$
digits. [wang2024Certifying](@cite)

## 1D Heisenberg Model with next nearest neighbor interaction

Polynomial Optimization framework is quite general. Almost no modification is
required to handle more complex Hamiltonian. 1D Heisenberg Model with geometric
frustration induced by next nearest neighbor interaction can be solved as:

````julia
using NCTSSoS, MosekTools
N = 6
J1 = 1.0                            # Nearest Neighbor Interaction
J2 = 0.2                            # Next Nearest Neighbor Interaction

registry, (σx, σy, σz) = create_pauli_variables(1:N)

ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] + ComplexF64(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [σx, σy, σz] for i in 1:N)

pop = polyopt(ham, registry)

solver_config = SolverConfig(optimizer=SILENT_MOSEK, order=2, ts_algo = MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop, res, solver_config)
energy_per_site = res.objective / N
@show energy_per_site
````

````
-0.4270083217760736
````

The literature value for this $J_1=1$, $J_2=0.2$ chain is
$-0.4270083225302217$, and the output above matches it to $6$ digits
[wang2024Certifying](@cite).

## 2D Square Lattice

Extending Heisenberg model to $2$-D case is also straightforward. However `NCTSSoS.jl` is not efficient enough to handle system at this size.

````julia
using NCTSSoS, MosekTools
Nx = 3
Ny = 3
N = Nx * Ny
J1 = 1.0
J2 = 0.0

registry, (σx, σy, σz) = create_pauli_variables(1:N)

LI = LinearIndices((1:Nx, 1:Ny))

ham = sum(ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(i, mod1(j + 1, Ny))]] + ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), j)]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j + 1, Ny))]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j - 1, Ny))]] for op in [σx, σy, σz] for i in 1:Nx for j in 1:Ny)

pop = polyopt(ham, registry)

solver_config = SolverConfig(optimizer=SILENT_MOSEK, order=3, cs_algo=MF(), ts_algo=MMD())
````

````
NCTSSoS.SolverConfig(MathOptInterface.OptimizerWithAttributes(Mosek.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[MathOptInterface.Silent() => true]), 3, CliqueTrees.MF(), CliqueTrees.MMD(0))
````

## Next step

With such lower bounds, estimates of properties of the ground
state, correlations functions, structure factors and order parameters, can also
be obtained. We provide examples in another section.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


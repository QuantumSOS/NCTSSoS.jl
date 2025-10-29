```@meta
EditURL = "../literate/pauli_algebra_interface.jl"
```

# Simplified Quantum Spin Models with Pauli Algebra Interface

NCTSSoS provides a convenient interface for working with common quantum algebras,
eliminating the need to manually specify commutation relations and constraints.
This tutorial demonstrates the `pauli_algebra` constructor for quantum spin systems.

## Basic Usage: Heisenberg XXX Model

The traditional approach requires manually defining Pauli commutation relations:

````julia
using NCTSSoS, MosekTools
N = 6
@ncpolyvar x[1:N] y[1:N] z[1:N]
````

````
(NCTSSoS.FastPolynomials.Variable[x₁, x₂, x₃, x₄, x₅, x₆], NCTSSoS.FastPolynomials.Variable[y₁, y₂, y₃, y₄, y₅, y₆], NCTSSoS.FastPolynomials.Variable[z₁, z₂, z₃, z₄, z₅, z₆])
````

Construct the XXX Heisenberg Hamiltonian

````julia
ham = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)] for op in [x, y, z] for i in 1:N)
````

````
0.25 + 0.0im * x₁¹x₂¹ + 0.25 + 0.0im * x₂¹x₃¹ + 0.25 + 0.0im * x₃¹x₄¹ + 0.25 + 0.0im * x₄¹x₅¹ + 0.25 + 0.0im * x₅¹x₆¹ + 0.25 + 0.0im * x₆¹x₁¹ + 0.25 + 0.0im * y₁¹y₂¹ + 0.25 + 0.0im * y₂¹y₃¹ + 0.25 + 0.0im * y₃¹y₄¹ + 0.25 + 0.0im * y₄¹y₅¹ + 0.25 + 0.0im * y₅¹y₆¹ + 0.25 + 0.0im * y₆¹y₁¹ + 0.25 + 0.0im * z₁¹z₂¹ + 0.25 + 0.0im * z₂¹z₃¹ + 0.25 + 0.0im * z₃¹z₄¹ + 0.25 + 0.0im * z₄¹z₅¹ + 0.25 + 0.0im * z₅¹z₆¹ + 0.25 + 0.0im * z₆¹z₁¹
````

**Old way**: Manual constraint specification

````julia
eq_cons = reduce(vcat, [
    [x[i] * y[i] - im * z[i],
     y[i] * x[i] + im * z[i],
     y[i] * z[i] - im * x[i],
     z[i] * y[i] + im * x[i],
     z[i] * x[i] - im * y[i],
     x[i] * z[i] + im * y[i]]
    for i in 1:N
])

pop_old = cpolyopt(ham;
    eq_constraints=eq_cons,
    comm_gps=[[x[i], y[i], z[i]] for i in 1:N],
    is_unipotent=true
)
````

````
    obj: 

        0.25 + 0.0im * x₁¹x₂¹ + 0.25 + 0.0im * x₂¹x₃¹ + 0.25 + 0.0im * x₃¹x₄¹ + 0.25 + 0.0im * x₄¹x₅¹ + 0.25 + 0.0im * x₅¹x₆¹ + 0.25 + 0.0im * x₆¹x₁¹ + 0.25 + 0.0im * y₁¹y₂¹ + 0.25 + 0.0im * y₂¹y₃¹ + 0.25 + 0.0im * y₃¹y₄¹ + 0.25 + 0.0im * y₄¹y₅¹ + 0.25 + 0.0im * y₅¹y₆¹ + 0.25 + 0.0im * y₆¹y₁¹ + 0.25 + 0.0im * z₁¹z₂¹ + 0.25 + 0.0im * z₂¹z₃¹ + 0.25 + 0.0im * z₃¹z₄¹ + 0.25 + 0.0im * z₄¹z₅¹ + 0.25 + 0.0im * z₅¹z₆¹ + 0.25 + 0.0im * z₆¹z₁¹ 

    constraints: 

        0.0 - 1.0im * z₁¹ + 1.0 + 0.0im * x₁¹y₁¹ = 0 	0.0 + 1.0im * z₁¹ + 1.0 + 0.0im * y₁¹x₁¹ = 0 	0.0 - 1.0im * x₁¹ + 1.0 + 0.0im * y₁¹z₁¹ = 0 	0.0 + 1.0im * x₁¹ + 1.0 + 0.0im * z₁¹y₁¹ = 0 	0.0 - 1.0im * y₁¹ + 1.0 + 0.0im * z₁¹x₁¹ = 0 	0.0 + 1.0im * y₁¹ + 1.0 + 0.0im * x₁¹z₁¹ = 0 	0.0 - 1.0im * z₂¹ + 1.0 + 0.0im * x₂¹y₂¹ = 0 	0.0 + 1.0im * z₂¹ + 1.0 + 0.0im * y₂¹x₂¹ = 0 	0.0 - 1.0im * x₂¹ + 1.0 + 0.0im * y₂¹z₂¹ = 0 	0.0 + 1.0im * x₂¹ + 1.0 + 0.0im * z₂¹y₂¹ = 0 	0.0 - 1.0im * y₂¹ + 1.0 + 0.0im * z₂¹x₂¹ = 0 	0.0 + 1.0im * y₂¹ + 1.0 + 0.0im * x₂¹z₂¹ = 0 	0.0 - 1.0im * z₃¹ + 1.0 + 0.0im * x₃¹y₃¹ = 0 	0.0 + 1.0im * z₃¹ + 1.0 + 0.0im * y₃¹x₃¹ = 0 	0.0 - 1.0im * x₃¹ + 1.0 + 0.0im * y₃¹z₃¹ = 0 	0.0 + 1.0im * x₃¹ + 1.0 + 0.0im * z₃¹y₃¹ = 0 	0.0 - 1.0im * y₃¹ + 1.0 + 0.0im * z₃¹x₃¹ = 0 	0.0 + 1.0im * y₃¹ + 1.0 + 0.0im * x₃¹z₃¹ = 0 	0.0 - 1.0im * z₄¹ + 1.0 + 0.0im * x₄¹y₄¹ = 0 	0.0 + 1.0im * z₄¹ + 1.0 + 0.0im * y₄¹x₄¹ = 0 	0.0 - 1.0im * x₄¹ + 1.0 + 0.0im * y₄¹z₄¹ = 0 	0.0 + 1.0im * x₄¹ + 1.0 + 0.0im * z₄¹y₄¹ = 0 	0.0 - 1.0im * y₄¹ + 1.0 + 0.0im * z₄¹x₄¹ = 0 	0.0 + 1.0im * y₄¹ + 1.0 + 0.0im * x₄¹z₄¹ = 0 	0.0 - 1.0im * z₅¹ + 1.0 + 0.0im * x₅¹y₅¹ = 0 	0.0 + 1.0im * z₅¹ + 1.0 + 0.0im * y₅¹x₅¹ = 0 	0.0 - 1.0im * x₅¹ + 1.0 + 0.0im * y₅¹z₅¹ = 0 	0.0 + 1.0im * x₅¹ + 1.0 + 0.0im * z₅¹y₅¹ = 0 	0.0 - 1.0im * y₅¹ + 1.0 + 0.0im * z₅¹x₅¹ = 0 	0.0 + 1.0im * y₅¹ + 1.0 + 0.0im * x₅¹z₅¹ = 0 	0.0 - 1.0im * z₆¹ + 1.0 + 0.0im * x₆¹y₆¹ = 0 	0.0 + 1.0im * z₆¹ + 1.0 + 0.0im * y₆¹x₆¹ = 0 	0.0 - 1.0im * x₆¹ + 1.0 + 0.0im * y₆¹z₆¹ = 0 	0.0 + 1.0im * x₆¹ + 1.0 + 0.0im * z₆¹y₆¹ = 0 	0.0 - 1.0im * y₆¹ + 1.0 + 0.0im * z₆¹x₆¹ = 0 	0.0 + 1.0im * y₆¹ + 1.0 + 0.0im * x₆¹z₆¹ = 0
        
    variables:
        x₁ x₂ x₃ x₄ x₅ x₆ y₁ y₂ y₃ y₄ y₅ y₆ z₁ z₂ z₃ z₄ z₅ z₆ 

    is_unipotent:
        true 

    is_projective:
        false 


````

**New way**: Using the algebra constructor

````julia
sys = pauli_algebra(N)
x_new, y_new, z_new = sys.variables
````

````
(NCTSSoS.FastPolynomials.Variable[x₁, x₂, x₃, x₄, x₅, x₆], NCTSSoS.FastPolynomials.Variable[y₁, y₂, y₃, y₄, y₅, y₆], NCTSSoS.FastPolynomials.Variable[z₁, z₂, z₃, z₄, z₅, z₆])
````

Construct the same Hamiltonian with new variables

````julia
ham_new = sum(ComplexF64(1/4) * op[i] * op[mod1(i+1, N)]
              for op in [x_new, y_new, z_new] for i in 1:N)

pop_new = cpolyopt(ham_new, sys)
````

````
    obj: 

        0.25 + 0.0im * x₁¹x₂¹ + 0.25 + 0.0im * x₂¹x₃¹ + 0.25 + 0.0im * x₃¹x₄¹ + 0.25 + 0.0im * x₄¹x₅¹ + 0.25 + 0.0im * x₅¹x₆¹ + 0.25 + 0.0im * x₆¹x₁¹ + 0.25 + 0.0im * y₁¹y₂¹ + 0.25 + 0.0im * y₂¹y₃¹ + 0.25 + 0.0im * y₃¹y₄¹ + 0.25 + 0.0im * y₄¹y₅¹ + 0.25 + 0.0im * y₅¹y₆¹ + 0.25 + 0.0im * y₆¹y₁¹ + 0.25 + 0.0im * z₁¹z₂¹ + 0.25 + 0.0im * z₂¹z₃¹ + 0.25 + 0.0im * z₃¹z₄¹ + 0.25 + 0.0im * z₄¹z₅¹ + 0.25 + 0.0im * z₅¹z₆¹ + 0.25 + 0.0im * z₆¹z₁¹ 

    constraints: 

        0.0 - 1.0im * z₁¹ + 1.0 + 0.0im * x₁¹y₁¹ = 0 	0.0 + 1.0im * z₁¹ + 1.0 + 0.0im * y₁¹x₁¹ = 0 	0.0 - 1.0im * x₁¹ + 1.0 + 0.0im * y₁¹z₁¹ = 0 	0.0 + 1.0im * x₁¹ + 1.0 + 0.0im * z₁¹y₁¹ = 0 	0.0 - 1.0im * y₁¹ + 1.0 + 0.0im * z₁¹x₁¹ = 0 	0.0 + 1.0im * y₁¹ + 1.0 + 0.0im * x₁¹z₁¹ = 0 	0.0 - 1.0im * z₂¹ + 1.0 + 0.0im * x₂¹y₂¹ = 0 	0.0 + 1.0im * z₂¹ + 1.0 + 0.0im * y₂¹x₂¹ = 0 	0.0 - 1.0im * x₂¹ + 1.0 + 0.0im * y₂¹z₂¹ = 0 	0.0 + 1.0im * x₂¹ + 1.0 + 0.0im * z₂¹y₂¹ = 0 	0.0 - 1.0im * y₂¹ + 1.0 + 0.0im * z₂¹x₂¹ = 0 	0.0 + 1.0im * y₂¹ + 1.0 + 0.0im * x₂¹z₂¹ = 0 	0.0 - 1.0im * z₃¹ + 1.0 + 0.0im * x₃¹y₃¹ = 0 	0.0 + 1.0im * z₃¹ + 1.0 + 0.0im * y₃¹x₃¹ = 0 	0.0 - 1.0im * x₃¹ + 1.0 + 0.0im * y₃¹z₃¹ = 0 	0.0 + 1.0im * x₃¹ + 1.0 + 0.0im * z₃¹y₃¹ = 0 	0.0 - 1.0im * y₃¹ + 1.0 + 0.0im * z₃¹x₃¹ = 0 	0.0 + 1.0im * y₃¹ + 1.0 + 0.0im * x₃¹z₃¹ = 0 	0.0 - 1.0im * z₄¹ + 1.0 + 0.0im * x₄¹y₄¹ = 0 	0.0 + 1.0im * z₄¹ + 1.0 + 0.0im * y₄¹x₄¹ = 0 	0.0 - 1.0im * x₄¹ + 1.0 + 0.0im * y₄¹z₄¹ = 0 	0.0 + 1.0im * x₄¹ + 1.0 + 0.0im * z₄¹y₄¹ = 0 	0.0 - 1.0im * y₄¹ + 1.0 + 0.0im * z₄¹x₄¹ = 0 	0.0 + 1.0im * y₄¹ + 1.0 + 0.0im * x₄¹z₄¹ = 0 	0.0 - 1.0im * z₅¹ + 1.0 + 0.0im * x₅¹y₅¹ = 0 	0.0 + 1.0im * z₅¹ + 1.0 + 0.0im * y₅¹x₅¹ = 0 	0.0 - 1.0im * x₅¹ + 1.0 + 0.0im * y₅¹z₅¹ = 0 	0.0 + 1.0im * x₅¹ + 1.0 + 0.0im * z₅¹y₅¹ = 0 	0.0 - 1.0im * y₅¹ + 1.0 + 0.0im * z₅¹x₅¹ = 0 	0.0 + 1.0im * y₅¹ + 1.0 + 0.0im * x₅¹z₅¹ = 0 	0.0 - 1.0im * z₆¹ + 1.0 + 0.0im * x₆¹y₆¹ = 0 	0.0 + 1.0im * z₆¹ + 1.0 + 0.0im * y₆¹x₆¹ = 0 	0.0 - 1.0im * x₆¹ + 1.0 + 0.0im * y₆¹z₆¹ = 0 	0.0 + 1.0im * x₆¹ + 1.0 + 0.0im * z₆¹y₆¹ = 0 	0.0 - 1.0im * y₆¹ + 1.0 + 0.0im * z₆¹x₆¹ = 0 	0.0 + 1.0im * y₆¹ + 1.0 + 0.0im * x₆¹z₆¹ = 0
        
    variables:
        x₁ x₂ x₃ x₄ x₅ x₆ y₁ y₂ y₃ y₄ y₅ y₆ z₁ z₂ z₃ z₄ z₅ z₆ 

    is_unipotent:
        true 

    is_projective:
        false 


````

Both approaches produce identical optimization problems, but the new interface is
more concise and less error-prone.

## Solving the Optimization Problem

````julia
solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=2)
res = cs_nctssos(pop_new, solver_config)
energy_per_site = res.objective / N
````

````
-0.46712927253210834
````

The result matches the known ground state energy per site of approximately -0.467129
for the 6-site XXX Heisenberg model.

## Combining Algebra with Custom Constraints

## Summary

The `pauli_algebra` interface provides:
1. **Automatic constraint generation**: No need to manually write commutation relations
2. **Error prevention**: Correct algebra structure guaranteed
3. **Code clarity**: Physics intent is clear from `pauli_algebra(N)`
4. **Flexibility**: Can still add custom constraints when needed
5. **Consistency**: Same interface works for different system sizes

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


using NCTSSoS
using MosekTools
using JuMP

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
                    )

res = cs_nctssos(pop, solver_config)

model = res.model

cons = all_constraints(model, include_variable_in_set_constraints=true)

value.(cons[1])

primal_var = Matrix(dual(cons[1]))

"""
    convert_X_to_H(X)

Convert a Hermitian matrix X in the block form [[X₁, X₃†], [X₃, X₂]] 
to H = H_R + im * H_I.

The conversion follows:
- H_R = X₁ + X₂
- H_I = (X₃ - X₃†) / (1im)
- H = H_R + im * H_I

Returns the complex Hermitian matrix H.
"""
function convert_X_to_H(X)
    n = size(X, 1)
    @assert iseven(n) "Matrix size must be even"
    
    # Split X into blocks
    half = n ÷ 2
    X₁ = X[1:half, 1:half]
    X₂ = X[half+1:end, half+1:end]
    X₃ = X[half+1:end, 1:half]
    X₃_dagger = X[1:half, half+1:end]
    
    # Compute H_R and H_I
    H_R = X₁ + X₂
    H_I = (X₃ - X₃_dagger) 
    
    # Return H = H_R + im * H_I
    return H_R + im * H_I
end

# Convert primal_var from X form to H
H = convert_X_to_H(primal_var)

using NCTSSoS.FastPolynomials: get_basis, neat_dot

get_basis([x; y; z], 2)

# size mismatch because the original problem is complex valued and realification comes with a factor of 2 
paulis = reconstruct(H, [x; y; z], 2, 2)

# how do I limit re-constructed paulis to proper size?
2^6

paulis[1]

res.objective / N
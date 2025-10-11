using NCTSSoS
using COSMO
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
                    optimizer=COSMO.Optimizer,          # the solver backend
                    order=2,                        # moment matrix order
                    )

res = cs_nctssos(pop, solver_config)

model = res.model

cons = all_constraints(model, include_variable_in_set_constraints=true)

value.(cons[1])

hankel = Matrix(dual(cons[1]))


using NCTSSoS.FastPolynomials: get_basis, neat_dot

get_basis([x; y; z], 2)

# size mismatch because the original problem is complex valued and realification comes with a factor of 2 
reconstruct(hankel, [x; y; z], 2, 2)

res.objective / N
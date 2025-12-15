using BenchmarkTools
using MosekTools, NCTSSoS
using Profile
using JET

using NCTSSoS.FastPolynomials: neat_dot, star, simplify

order = 3
n = 20
reg, (x,) = create_noncommutative_variables([("x", 1:n)])
f = 0.0
for i = 1:n
    jset = max(1, i - 5):min(n, i + 1)
    jset = setdiff(jset, i)
    f += (2x[i] + 5 * x[i]^3 + 1)^2
    f -= sum([
        4x[i] * x[j] +
        10x[i]^3 * x[j] +
        2x[j] +
        4x[i] * x[j]^2 +
        10x[i]^3 * x[j]^2 +
        2x[j]^2 for j in jset
    ])
    f += sum([
        x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
    ])
end

cons = vcat([(1 - x[i]^2) for i = 1:n], [(x[i] - 1 / 3) for i = 1:n])

pop = polyopt(f, reg; ineq_constraints=cons)
solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=order,
    cs_algo=MF(), ts_algo=MMD())


br = @benchmark result = cs_nctssos($pop, $solver_config)

result = cs_nctssos(pop, solver_config)
Profile.clear()
@profile result = cs_nctssos(pop, solver_config)
Profile.print(mincount=300)


# Test simplification with typed algebras
reg_uni, (xu,) = create_unipotent_variables([("x", 1:10)])
reg_proj, (xp,) = create_projector_variables([("x", 1:10)])

# Simplification benchmarks
@benchmark simplify(a) setup = (
    a = xu[6] * xu[3] * xu[6] * xu[1] * xu[2] * xu[5] * xu[4]
)

a_uni = xu[6] * xu[3] * xu[6] * xu[1] * xu[2] * xu[5] * xu[4]
Profile.clear()
@profile for _ in 1:5000000 simplify(a_uni) end
Profile.print(mincount=50)


@benchmark simplify(a) setup = (
    a = xp[6] * xp[3] * xp[6] * xp[1] * xp[2] * xp[5] * xp[4]
)
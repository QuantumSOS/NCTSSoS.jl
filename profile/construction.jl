using BenchmarkTools
using MosekTools, NCTSSoS
using Profile
using Test

n = 10
reg, (x,) = create_noncommutative_variables([("x", 1:n)])
f = 0.0
for i = 1:n
    jset = max(1, i - 5):min(n, i + 1)
    jset = setdiff(jset, i)
    f += (2x[i] + 5 * x[i]^3 + 1)^2
    f -= sum([4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] + 4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset])
    f += sum([x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset])
end

cons = vcat([(1 - x[i]^2) for i = 1:n], [(x[i] - 1 / 3) for i = 1:n])
pop = polyopt(f, reg; ineq_constraints=cons)
solver_config = SolverConfig(optimizer=Mosek.Optimizer; order=3, cs_algo=MF(), ts_algo=MMD())

@btime cs_nctssos($pop, $solver_config);
cs_nctssos(pop, solver_config)

Profile.clear()
@profile result_cs_ts = cs_nctssos(pop, solver_config)
Profile.print(mincount=100, format=:tree)


using NCTSSoS.FastPolynomials: ς, Monomial

reg2, (x2, y2) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
cov(a, b) = 1.0 * ς(x2[a] * y2[b]) - 1.0 * ς(x2[a]) * ς(y2[b])
sp = cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2)

spop = polyopt(sp * one(Monomial), reg2)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)

@btime cs_nctssos($spop, $solver_config)

Profile.clear()
@profile for _ in 1:20 result_cs_ts = cs_nctssos(spop, solver_config) end
Profile.print(mincount=100, format=:tree)

using BenchmarkTools
using MosekTools, NCTSSoS
using Profile
using NCTSSoS.FastPolynomials: ς, Monomial

reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
cov(a, b) = 1.0 * ς(x[a] * y[b]) - 1.0 * ς(x[a]) * ς(y[b])
sp = cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2)

spop = polyopt(sp * one(Monomial), reg)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)

@btime cs_nctssos($spop, $solver_config)

Profile.clear()
@profile for _ in 1:5 result_cs_ts = cs_nctssos(spop, solver_config) end
Profile.print(mincount=200, format=:flat)
Profile.print(mincount=100, format=:tree, recur=:flat)
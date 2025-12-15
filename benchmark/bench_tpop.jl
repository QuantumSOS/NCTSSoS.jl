module BenchTracePolyOpt
using BenchmarkTools
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials: tr, Monomial

const SUITE = BenchmarkGroup()

reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])

cov(i, j) = tr(x[i] * y[j]) - tr(x[i]) * tr(y[j])
p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2))

SUITE["Covariance Example"] = @benchmarkable cs_nctssos(tpop, solver_config) setup = (tpop = polyopt(p * one(Monomial), reg);
solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)
)


end

BenchTracePolyOpt.SUITE
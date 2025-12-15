module BenchStatePolyOpt
using BenchmarkTools
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials: ς, Monomial

const SUITE = BenchmarkGroup()

reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

d = 3

SUITE["State Poly SOS Problem"] = @benchmarkable cs_nctssos(spop, solver_config) setup = (
    spop = polyopt(sp * one(Monomial); is_unipotent=true, comm_gps=[x, y]);
    solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=d)
)

end

BenchStatePolyOpt.SUITE
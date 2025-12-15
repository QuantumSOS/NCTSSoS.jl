module BenchComplexPolyOpt
using BenchmarkTools
using NCTSSoS, MosekTools

const SUITE = BenchmarkGroup()
N = 3
reg, (x, y, z) = create_pauli_variables(1:N)

J = 1.0
h = 2.0
ham = sum(-complex(J / 4) * z[i] * z[mod1(i + 1, N)] for i in 1:N) + sum(-h / 2 * x[i] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

SUITE["Complex Poly Opt Problem"] = @benchmarkable cs_nctssos(cpop, solver_config) setup = (
    cpop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true); solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)
)

end

BenchComplexPolyOpt.SUITE
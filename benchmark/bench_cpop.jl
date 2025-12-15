module BenchComplexPolyOpt
using BenchmarkTools
using NCTSSoS, MosekTools

const SUITE = BenchmarkGroup()
const N = 3
const reg, (x, y, z) = create_pauli_variables(1:N)

const J = 1.0
const h = 2.0
const ham = sum(-complex(J / 4) * z[i] * z[mod1(i + 1, N)] for i in 1:N) + sum(-h / 2 * x[i] for i in 1:N)

# Using Pauli algebra, constraints are handled automatically by the registry
SUITE["Pauli Poly Opt Problem"] = @benchmarkable cs_nctssos(pop, solver_config) setup = (
    pop = polyopt(ham, reg);
    solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)
)

end

BenchComplexPolyOpt.SUITE
using NCTSSoS, MosekTools
using NCTSSoS: tr, Monomial

registry, (x, y) = create_unipotent_variables([(x, 1:3), (y, 1:3)])

cov(i, j) = tr(x[i] * y[j]) - tr(x[i]) * tr(y[j])
p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2))
tpop = polyopt(p * one(Monomial), registry)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)

result = cs_nctssos(tpop, solver_config)

@info "Result objective: $(result.objective)"
@assert isapprox(result.objective, -5.0, atol=1e-5) "Expected -5.0, got $(result.objective)"

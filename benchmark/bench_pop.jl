module BenchPolyOpt
using BenchmarkTools
using NCTSSoS, MosekTools

const SUITE = BenchmarkGroup()

function make_poly(x, n)
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
    return f
end

SUITE["Example 1"] = @benchmarkable result = cs_nctssos(pop, solver_config) setup = (
    begin
        n = 20
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])
        f = make_poly(x, n)
        order = 3
        cons = vcat([(1 - x[i]^2) for i = 1:n], [(x[i] - 1 / 3) for i = 1:n])
        pop = polyopt(f; ineq_constraints=cons)
        solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=order,
            cs_algo=MF(), ts_algo=MMD())
    end
)


SUITE["Example 2"] = @benchmarkable result = cs_nctssos(pop, solver_config) setup = (
    begin
        reg, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])

        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]  # objective function
        pop = polyopt(-f, reg)
        solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=3,
            cs_algo=MF(), ts_algo=MMD())
    end
)

SUITE["Example 3"] = @benchmarkable result = cs_nctssos(pop, solver_config) setup = (
    begin
        num_sites = 6
        J1_interactions =
            unique!([tuple(sort([i, mod1(i + 1, num_sites)])...) for i = 1:num_sites])
        J2_interactions =
            unique!([tuple(sort([i, mod1(i + 2, num_sites)])...) for i = 1:num_sites])

        J1 = 2.0
        J2 = 1.0

        true_ans = -num_sites / 4 * 6

        ij2idx_dict = Dict(
            zip(
                [(i, j) for i in 1:num_sites, j in 1:num_sites if j > i],
                1:(num_sites*(num_sites-1)÷2),
            ),
        )
        reg, (pij,) = create_projector_variables([("pij", 1:(num_sites*(num_sites-1)÷2))])

        objective = (
            sum([J1 * pij[ij2idx_dict[(i, j)]] for (i, j) in J1_interactions]) + sum([J2 * pij[ij2idx_dict[(i, j)]] for (i, j) in J2_interactions])
        )

        gs = unique!([
            (
                pij[ij2idx_dict[tuple(sort([i, j])...)]] *
                pij[ij2idx_dict[tuple(sort([j, k])...)]] +
                pij[ij2idx_dict[tuple(sort([j, k])...)]] *
                pij[ij2idx_dict[tuple(sort([i, j])...)]] -
                0.5 * (
                    pij[ij2idx_dict[tuple(sort([i, j])...)]] +
                    pij[ij2idx_dict[tuple(sort([j, k])...)]] -
                    pij[ij2idx_dict[tuple(sort([i, k])...)]]
                )
            ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
            (i != j && j != k && i != k)
        ])

        pop = polyopt(-objective, reg; eq_constraints=gs)

        solver_config = SolverConfig(optimizer=Mosek.Optimizer; order=1)
    end
)

end

BenchPolyOpt.SUITE

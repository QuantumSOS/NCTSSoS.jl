using Test, NCTSSoS, NCTSSoS.FastPolynomials
if haskey(ENV, "LOCAL_TESTING") 
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end
using COSMO
const QUICK_SOLVER = COSMO.Optimizer
using JuMP
using NCTSSoS:
    get_state_basis,
    neat_dot,
    NCStateWord,
    NCStatePolynomial,
    constrain_moment_matrix!,
    substitute_variables,
    moment_relax
using NCTSSoS.FastPolynomials: expval, terms, symmetric_canonicalize, monomials

using NCTSSoS:
    correlative_sparsity,
    iterate_term_sparse_supp,
    sorted_union,
    MinimalChordal,
    NoElimination


@testset "State Polynomial Opt 7.2.0" begin
    @ncpolyvar x[1:2] y[1:2]
    sp =
        -1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) +
        1.0 * ς(x[2] * y[2])
    spop = polyopt(sp * one(Monomial); is_unipotent = true, comm_gps = [x, y])

    d = 1

    solver_config = SolverConfig(; optimizer = SOLVER, order = d)

    result_mom = cs_nctssos(spop, solver_config; dualize=false)
    result_sos = cs_nctssos(spop, solver_config)


    @test isapprox(result_mom.objective, -2.8284271321623202, atol=1e-5)

    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-5)
end

@testset "State Polynomial Opt 7.2.1" begin
    @ncpolyvar x[1:2] y[1:2]
    sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
    sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
    sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

    spop = polyopt(sp * one(Monomial); is_unipotent=true, comm_gps=[x, y])

    d = 3
    cr = correlative_sparsity(spop, d, NoElimination())

    cliques_objective = [
        reduce(
            +,
            [
                issubset(variables(t[2]), clique) ? t[1] * t[2] : zero(t[2]) for
                t in terms(spop.objective)
            ],
        ) for clique in cr.cliques
    ]

    initial_activated_supp = [
        sorted_union(
            symmetric_canonicalize.(monomials(obj_part)),
            mapreduce(
                a -> monomials(a),
                vcat,
                spop.constraints[cons_idx];
                init = typeof(monomials(spop.objective)[1])[],
            ),
        ) for (obj_part, cons_idx, idcs_bases) in
        zip(cliques_objective, cr.cliques_cons, cr.cliques_idcs_bases)
    ]

    cliques_term_sparsities = map(
        zip(initial_activated_supp, cr.cliques_cons, cr.cliques_idcs_bases),
    ) do (activated_supp, cons_idx, idcs_bases)
        [
            iterate_term_sparse_supp(activated_supp, poly, basis, NoElimination()) for
            (poly, basis) in
            zip([one(spop.objective); spop.constraints[cons_idx]], idcs_bases)
        ]
    end

    solver_config = SolverConfig(; optimizer = QUICK_SOLVER, order = d)

    result_mom =  cs_nctssos(spop, solver_config; dualize=false)
    @test isapprox(result_mom.objective, -4.0, atol = 1e-4)

    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -4.0, atol = 1e-5)
end

@testset "State Polynomial Opt 7.2.2" begin
    @ncpolyvar x[1:3] y[1:3]
    cov(a, b) = 1.0 * ς(x[a] * y[b]) - 1.0 * ς(x[a]) * ς(y[b])
    sp =
        cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) -
        cov(3, 2)

    spop = polyopt(sp*one(Monomial); is_unipotent = true, comm_gps = [x[1:3], y[1:3]])

    solver_config = SolverConfig(; optimizer = SOLVER, order = 2)

    @test cs_nctssos(spop, solver_config) ≈ -5.0 atol = 1e-2

    @ncpolyvar x[1:6]
    sp =
        -1.0 * ς(x[1] * x[4]) + 1 * ς(x[1]) * ς(x[4]) - 1 * ς(x[1] * x[5]) +
        1 * ς(x[1]) * ς(x[5]) - 1 * ς(x[1] * x[6]) + 1 * ς(x[1]) * ς(x[6]) -
        1 * ς(x[2] * x[4]) + 1 * ς(x[2]) * ς(x[4]) - 1 * ς(x[2] * x[5]) +
        1 * ς(x[2]) * ς(x[5]) +
        1 * ς(x[2] * x[6]) - 1 * ς(x[2]) * ς(x[6]) - 1 * ς(x[3] * x[4]) +
        1 * ς(x[3]) * ς(x[4]) +
        1 * ς(x[3] * x[5]) - 1 * ς(x[3]) * ς(x[5])


    spop = polyopt(sp; is_unipotent = true, comm_gps = [x[1:3], x[4:6]])

    solver_config = SolverConfig(; optimizer = SOLVER, order = 2)

    @test cs_nctssos(spop, solver_config) ≈ -5.0 atol = 1e-2

    d = 2
    cr = correlative_sparsity(spop, d, NoElimination())

    cliques_objective = [
        reduce(
            +,
            [
                issubset(variables(t[2]), clique) ? t[1] * t[2] : zero(t[2]) for
                t in terms(spop.objective)
            ],
        ) for clique in cr.cliques
    ]

    initial_activated_supp = [
        sorted_union(
            symmetric_canonicalize.(monomials(obj_part)),
            mapreduce(
                a -> monomials(a),
                vcat,
                spop.constraints[cons_idx];
                init = typeof(monomials(spop.objective)[1])[],
            ),
        ) for (obj_part, cons_idx, idcs_bases) in
        zip(cliques_objective, cr.cliques_cons, cr.cliques_idcs_bases)
    ]

    cliques_term_sparsities = map(
        zip(initial_activated_supp, cr.cliques_cons, cr.cliques_idcs_bases),
    ) do (activated_supp, cons_idx, idcs_bases)
        [
            iterate_term_sparse_supp(activated_supp, poly, basis, NoElimination()) for
            (poly, basis) in
            zip([one(spop.objective); spop.constraints[cons_idx]], idcs_bases)
        ]
    end

    mom_problem =
        moment_relax(spop, cr.cliques_cons, cr.global_cons, cliques_term_sparsities)

    sos_problem = sos_dualize(mom_problem)
    set_optimizer(
        sos_problem.model,
        optimizer_with_attributes(COSMO.Optimizer, "eps_rel" => 1e-8),
    )
    optimize!(sos_problem.model)
    @test isapprox(objective_value(sos_problem.model), -5.0, atol = 1e-3)
    result_sos = cs_nctssos(spop, solver_config)

    @test isapprox(result_sos.objective, -5.0, atol=1e-3)
end


@testset "Constrain Moment matrix" begin
    @ncpolyvar x[1:2]

    sa = SimplifyAlgorithm(comm_gps=[x], is_unipotent=false, is_projective=false)

    basis = get_state_basis(x, 1, sa)

    sp = 1.0 * ς(x[1] * x[2]) + 2.0 * ς(x[1]) + 3.0 * ς(x[2])
    nc_words = monomial.([one(x[1]), x[1], x[2]])
    ncsp =
        1.0 * ς(x[1] * x[2]) * one(Monomial) +
        2.0 * ς(x[1]) * monomial(x[1]) +
        3.0 * ς(x[2]) * monomial(x[2])
    poly = one(ncsp)

    total_basis = sort(unique([expval(neat_dot(a, b)) for a in basis for b in basis]))

    model = GenericModel{Float64}()
    @variable(model, y[1:length(total_basis)])
    wordmap = Dict(zip(total_basis, y))

    ncterms = map(
        a -> a[1] * NCStateWord(monomial.(a[2]), a[3]),
        zip([1.0, 2.0, 3.0], [[x[1] * x[2]], [x[1]], [x[2]]], nc_words),
    )

    @test map(a -> a[1] * a[2], terms(ncsp)) == ncterms
    @test substitute_variables(ncsp, wordmap) == 1.0 * y[7] + 3.0 * y[6] + 2.0 * y[4]

    true_mom_mtx = expval.([neat_dot(a, b) for a in basis, b in basis])
    mom_mtx_cons =
        constrain_moment_matrix!(model, one(ncsp), basis, wordmap, PSDCone(), identity)
    mom_mtx = constraint_object(mom_mtx_cons)
    reshape(mom_mtx.func, 5, 5)
    @test reshape(mom_mtx.func, 5, 5) == AffExpr[
        y[1] y[2] y[3] y[2] y[3];
        y[2] y[4] y[5] y[4] y[5];
        y[3] y[5] y[6] y[5] y[6];
        y[2] y[4] y[5] y[8] y[7];
        y[3] y[5] y[6] y[9] y[10]
    ]
end

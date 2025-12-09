using Test, NCTSSoS
using NCTSSoS.FastPolynomials
using JuMP

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end
using Graphs

using NCTSSoS.FastPolynomials: get_basis
using NCTSSoS: substitute_variables

@testset "Complex Polynomial Optimization" begin
    @testset "1D Trasverse Filed Ising Model" begin
        N = 3
        @ncpolyvar x[1:N] y[1:N] z[1:N]

        J = 1.0
        h = 2.0
        ham = sum(-complex(J / 4) * z[i] * z[mod1(i + 1, N)] for i in 1:N) + sum(-h / 2 * x[i] for i in 1:N)

        eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

        cpop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=1)

        order = iszero(solver_config.order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [cpop.objective; cpop.eq_constraints; cpop.ineq_constraints]]) : solver_config.order

        corr_sparsity = NCTSSoS.correlative_sparsity(cpop, order, solver_config.cs_algo)

        cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in NCTSSoS.FastPolynomials.terms(cpop.objective)]) for clique in corr_sparsity.cliques]

        initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
            NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
        end

        cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
            NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
        end

        cmp = NCTSSoS.moment_relax(cpop, corr_sparsity, cliques_term_sparsities)

        @test length(cmp.constraints) == 19
        @test length(cmp.total_basis) == 55
    end

    @testset "Naive Example" begin
        N = 1
        @ncpolyvar x[1:N] y[1:N] z[1:N]

        ham = sum(ComplexF64(1 / 2) * op[1] for op in [x, y, z])

        eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

        cpop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=1)

        order = iszero(solver_config.order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [cpop.objective; cpop.eq_constraints; cpop.ineq_constraints]]) : solver_config.order

        corr_sparsity = NCTSSoS.correlative_sparsity(cpop, order, solver_config.cs_algo)

        cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in NCTSSoS.FastPolynomials.terms(cpop.objective)]) for clique in corr_sparsity.cliques]

        initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
            NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
        end

        cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
            NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
        end

        cmp = NCTSSoS.moment_relax(cpop, corr_sparsity, cliques_term_sparsities)

        @test length(cmp.constraints) == 7

        @test cmp.constraints[1][1] == :HPSD
        @test size(cmp.constraints[1][2]) == (4, 4)
        @test length(cmp.total_basis) == 10
    end
end

@testset "Special Constraint Type " begin
    @testset "CHSH Inequality" begin
        @ncpolyvar x[1:2]
        @ncpolyvar y[1:2]

        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(f; comm_gps=[x, y], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER; order=1)

        result = cs_nctssos(pop, solver_config; dualize=false)

        @test isapprox(result.objective, -2.8284271321623193, atol=1e-6)
    end
end

@testset "CS TS Example" begin
    order = 3
    n = 10
    @ncpolyvar x[1:n]
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

    pop = polyopt(f; ineq_constraints=cons)

    solver_config =
        SolverConfig(optimizer=SOLVER; order=order, cs_algo=MF(), ts_algo=MMD())

    result = cs_nctssos(pop, solver_config; dualize=false)

    @test isapprox(result.objective, 3.011288, atol=1e-4)
end

@testset "Moment Method Heisenberg Model on Star Graph" begin
    num_sites = 10
    g = star_graph(num_sites)

    true_ans = -1.0

    vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]

    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

    @ncpolyvar pij[1:length(vec_idx2ij)]

    objective = sum(1.0 * pij[[findvaridx(ee.src, ee.dst) for ee in edges(g)]])

    gs = unique!([
        (
            pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
            pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
            pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
            pij[findvaridx(sort([i, k])...)] + 1.0
        ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
        (i != j && j != k && i != k)
    ])


    pop = polyopt(
        objective;
        eq_constraints=gs,
        is_unipotent=true,
    )
    order = 1
    cs_algo = MF()

    solver_config = SolverConfig(
        optimizer=SOLVER,
        order=order,
        cs_algo=cs_algo
    )

    result = cs_nctssos(pop, solver_config; dualize=false)
    @test isapprox(result.objective, true_ans, atol=1e-6)
end

@testset "Replace DynamicPolynomials variables with JuMP variables" begin
    @ncpolyvar x y z
    poly = 1.0 * x^2 - 2.0 * x * y - 1

    model = GenericModel{Float64}()
    @variable(model, jm[1:13])

    monomap = Dict(get_basis([x, y, z], 2) .=> jm)

    @test substitute_variables(poly, monomap) == 1.0 * jm[5] - 2.0 * jm[6] - jm[1]
end

@testset "Moment Method Example 1" begin
    order = 2
    n = 3
    @ncpolyvar x[1:n]
    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] +
        2.0x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    pop = polyopt(f)

    @testset "Dense" begin

        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            cs_algo=NoElimination(),
        )

        result = cs_nctssos(pop, solver_config; dualize=false)

        # NOTE: differs from original test case value since that one is a relaxed in terms of sparsity
        # This value here is obtained by running the master branch with no sparsity relaxation
        @test isapprox(
            result.objective,
            4.372259295498716e-10,
            atol=1e-6,
        )
    end

    @testset "Sprase" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            ts_algo=MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize=false)

        @test isapprox(result.objective, -0.0035512, atol=1e-7)
    end
end

@testset "Moment Method Example 2" begin
    order = 2
    n = 2
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    pop = polyopt(f; eq_constraints=[h1], ineq_constraints=[g])


    @testset "Dense" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order)

        result = cs_nctssos(pop, solver_config; dualize=false)
        @test isapprox(result.objective, -1.0, atol=1e-6)
    end

    @testset "Term Sparse" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            cs_algo=MF(),
            ts_algo=MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize=false)

        @test isapprox(result.objective, -1.0, atol=1e-6)
    end
end

@testset "Moment Method Correlative Sparsity" begin
    n = 3
    @ncpolyvar x[1:n]
    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] +
        2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0 / 3 for i = 1:n])
    order = 3

    pop = polyopt(f; ineq_constraints=cons)

    @testset "Correlative Sparse" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            cs_algo=MF(),
        )
        result = cs_nctssos(pop, solver_config; dualize=false)

        # FIXME: reduced accuracy
        # @test is_solved_and_feasible(moment_problem.model)
        @test isapprox(
            result.objective,
            0.9975306427277915,
            atol=1e-5,
        )
    end

    @testset "Term Sparse" begin
        solver_config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            ts_algo=MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize=false)

        result = cs_nctssos_higher(pop, result, solver_config; dualize=false)

        @test isapprox(
            result.objective,
            0.9975306427277915,
            atol=1e-5,
        )
    end
end

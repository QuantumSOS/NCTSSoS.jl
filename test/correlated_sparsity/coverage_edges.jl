# test/correlated_sparsity/coverage_edges.jl

struct _StatusStub
    term::MOI.TerminationStatusCode
    primal::MOI.ResultStatusCode
    dual::MOI.ResultStatusCode
end

JuMP.termination_status(model::_StatusStub) = model.term
JuMP.primal_status(model::_StatusStub) = model.primal
JuMP.dual_status(model::_StatusStub) = model.dual

@testset "Coverage Edges" begin
    @testset "compute_relaxation_order auto odd-degree" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        objective = 1.0 * x[1]^3 + x[2]
        pop = polyopt(objective, reg)
        @test NCTSSoS.compute_relaxation_order(pop, 0) == 2
    end

    @testset "_check_solver_status accepts iteration-limit feasible" begin
        model = _StatusStub(MOI.ITERATION_LIMIT, MOI.FEASIBLE_POINT, MOI.NO_SOLUTION)
        @test NCTSSoS._check_solver_status(model) == MOI.ITERATION_LIMIT

        model = _StatusStub(MOI.SLOW_PROGRESS, MOI.NO_SOLUTION, MOI.NEARLY_FEASIBLE_POINT)
        @test NCTSSoS._check_solver_status(model) == MOI.SLOW_PROGRESS
    end

    @testset "_check_solver_status throws with status details" begin
        model = _StatusStub(MOI.NUMERICAL_ERROR, MOI.NO_SOLUTION, MOI.NO_SOLUTION)

        @test_throws NCTSSoS.SolverStatusError NCTSSoS._check_solver_status(model)

        err = try
            NCTSSoS._check_solver_status(model)
            nothing
        catch e
            e
        end

        @test err isa NCTSSoS.SolverStatusError
        @test err.termination == MOI.NUMERICAL_ERROR
        @test err.primal == MOI.NO_SOLUTION
        @test err.dual == MOI.NO_SOLUTION
    end

    @testset "polyopt rejects integer coefficients" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])
        int_poly = NCTSSoS.Polynomial([(1, x[1])])
        @test_throws ArgumentError polyopt(int_poly, reg)
    end

    @testset "moment_relax processes global constraints branch" begin
        pop, sparsity, corr_with_global = build_global_constraint_fixture()
        @test !isempty(corr_with_global.global_cons)

        mp = NCTSSoS.moment_relax(pop, corr_with_global, sparsity.cliques_term_sparsities)
        global_poly = corr_with_global.cons[first(corr_with_global.global_cons)]
        @test any(size(mat) == (1, 1) && mat[1, 1] == global_poly for (_cone, mat) in mp.constraints)
    end

    @testset "_substitute_poly zero branch" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])
        model = Model()
        @variable(model, y)

        zero_poly = zero(typeof(1.0 * x[1]))
        monomap = Dict(1 => y)
        expr = NCTSSoS._substitute_poly(zero_poly, monomap)
        @test JuMP.coefficient(expr, y) == 0
    end

    @testset "_add_parity_constraints! branch for fermionic algebra" begin
        reg, (a, a_dag) = create_fermionic_variables(1:2)
        objective = a_dag[1] * a[1] + a_dag[2] * a[2]
        pop = polyopt(objective, reg)

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        sparsity = compute_sparsity(pop, config)
        mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        odd_poly = 1.0 * a[1]
        push!(mp.constraints, (:HPSD, reshape([odd_poly], 1, 1)))

        zero_before = count(c -> c[1] == :Zero, mp.constraints)
        NCTSSoS._add_parity_constraints!(mp)
        zero_after = count(c -> c[1] == :Zero, mp.constraints)

        @test zero_after > zero_before
        @test any(cone == :Zero && size(mat) == (1, 1) && mat[1, 1] == odd_poly for (cone, mat) in mp.constraints)
    end

    @testset "_solve_real_moment_problem rejects unknown cone" begin
        reg, (u,) = create_unipotent_variables([("u", 1:1)])
        objective = 1.0 * u[1]
        pop = polyopt(objective, reg)
        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        sparsity = compute_sparsity(pop, config)
        mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)

        bad_constraints = copy(mp.constraints)
        bad_constraints[1] = (:BadCone, bad_constraints[1][2])
        bad_mp = typeof(mp)(mp.objective, bad_constraints, mp.total_basis, mp.n_unique_moment_matrix_elements)

        err = try
            NCTSSoS._solve_real_moment_problem(bad_mp, SOLVER, true)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("Unexpected cone type", sprint(showerror, err))
    end

    @testset "show(::CorrelativeSparsity) prints global constraints" begin
        _pop, _sparsity, corr_with_global = build_global_constraint_fixture()

        io = IOBuffer()
        show(io, corr_with_global)
        shown = String(take!(io))
        @test occursin("Global Constraints", shown)
        @test occursin("1", shown)
    end
end

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

    @testset "term_sparsity_graph_supp empty-basis guard (non-state)" begin
        expected = correlated_structure_case("term_sparsity_empty_basis_nonstate")
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])
        basis0 = get_ncbasis(reg, 1)
        empty_basis = eltype(basis0)[]
        g = 1.0 * x[1]

        supp = NCTSSoS.term_sparsity_graph_supp(SimpleGraph(0), empty_basis, g)
        @test length(supp) == json_int(expected["supp_length"])
        @test isempty(supp)
    end

    @testset "term_sparsity_graph_supp empty-basis guard (state)" begin
        expected = correlated_structure_case("term_sparsity_empty_basis_state")
        reg, (x,) = create_unipotent_variables([("x", 1:1)])
        one_mono = one(typeof(x[1]))
        g = (1.0 * ς(x[1])) * one_mono
        basis0 = get_state_basis(reg, 1)
        empty_basis = eltype(basis0)[]

        supp = NCTSSoS.term_sparsity_graph_supp(SimpleGraph(0), empty_basis, g)
        @test length(supp) == json_int(expected["supp_length"])
        @test isempty(supp)
    end

    @testset "get_term_sparsity_graph edge detection (state)" begin
        expected_lr = correlated_structure_case("term_sparsity_connected_lr_state")
        expected_rl = correlated_structure_case("term_sparsity_connected_rl_state")

        # connected_lr branch
        reg_lr, (x_lr,) = create_unipotent_variables([("x", 1:1)])
        basis_lr = get_state_basis(reg_lr, 1)
        id_lr = one(eltype(basis_lr))
        w_lr = first(filter(b -> !isone(b), basis_lr))

        graph_lr = NCTSSoS.get_term_sparsity_graph([id_lr], [w_lr], [id_lr, w_lr])
        @test has_edge(graph_lr, 1, 2) == Bool(expected_lr["has_edge_1_2"])

        # connected_rl branch: activated support matches only swapped product
        reg_rl, (x_rl,) = create_unipotent_variables([("x", 1:2)])
        idx1 = x_rl[1].word[1]
        idx2 = x_rl[2].word[1]
        mono1 = typeof(x_rl[1])([idx1])
        mono12 = typeof(x_rl[1])([idx1, idx2])
        mono2 = typeof(x_rl[1])([idx2])

        basis_rl = get_state_basis(reg_rl, 2)
        id_rl = one(eltype(basis_rl))
        i_w = findfirst(b -> isone(b.sw) && b.nc_word == mono1, basis_rl)
        i_m = findfirst(b -> isone(b.sw) && b.nc_word == mono12, basis_rl)
        i_t = findfirst(b -> isone(b.sw) && b.nc_word == mono2, basis_rl)
        @test i_w !== nothing
        @test i_m !== nothing
        @test i_t !== nothing

        w = basis_rl[i_w]
        m = basis_rl[i_m]
        target = basis_rl[i_t]

        graph_rl = NCTSSoS.get_term_sparsity_graph([m], [target], [id_rl, w])
        @test has_edge(graph_rl, 1, 2) == Bool(expected_rl["has_edge_1_2"])
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

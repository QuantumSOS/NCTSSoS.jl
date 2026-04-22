using NCTSSoS: CorrelativeSparsity

function _build_chsh_trace_problem()
    reg, (vars,) = create_unipotent_variables([("v", 1:4)])
    x = vars[1:2]
    y = vars[3:4]

    p = -1.0 * tr(x[1] * y[1]) - tr(x[1] * y[2]) -
        tr(x[2] * y[1]) + tr(x[2] * y[2])
    return polyopt(p * one(typeof(x[1])), reg)
end

@testset "Trace StatePolynomial input is promoted to NCStatePolynomial" begin
    reg, (vars,) = create_unipotent_variables([("v", 1:4)])
    x = vars[1:2]
    y = vars[3:4]

    p = -1.0 * tr(x[1] * y[1]) - tr(x[1] * y[2]) - tr(x[2] * y[1]) + tr(x[2] * y[2])
    eq_p = -1.0 + tr(x[1])
    ineq_p = 2.0 + tr(y[1])
    pop_state = polyopt(p, reg)
    pop_state_cons = polyopt(p, reg; eq_constraints=[eq_p], ineq_constraints=[ineq_p])
    pop_nc = polyopt(p * one(typeof(x[1])), reg)
    pop_nc_cons = polyopt(
        p * one(typeof(x[1])),
        reg;
        eq_constraints=[eq_p * one(typeof(x[1]))],
        ineq_constraints=[ineq_p * one(typeof(x[1]))],
    )

    @test pop_state.objective isa NCTSSoS.NCStatePolynomial{Float64,MaxEntangled,UnipotentAlgebra,UInt8}
    @test pop_state.objective == pop_nc.objective
    @test pop_state_cons.eq_constraints == pop_nc_cons.eq_constraints
    @test pop_state_cons.ineq_constraints == pop_nc_cons.ineq_constraints
    @test isempty(pop_state.eq_constraints)
    @test isempty(pop_state.ineq_constraints)
    @test all(isone(ncsw.nc_word) for ncsw in monomials(pop_state.objective))

    config = SolverConfig(
        optimizer=SOLVER,
        order=1,
        cs_algo=NoElimination(),
        ts_algo=NoElimination()
    )
    result_state = cs_nctssos(pop_state, config)
    result_nc = cs_nctssos(pop_nc, config)

    @test result_state.objective ≈ result_nc.objective atol = 1e-5
    @test result_state.n_unique_moment_matrix_elements == result_nc.n_unique_moment_matrix_elements
    @test result_state.moment_matrix_sizes == result_nc.moment_matrix_sizes
end

@testset "Trace Optimization Paths" begin
    tpop = _build_chsh_trace_problem()
    dense_oracle = expectations_oracle("expectations/trace_optimization_paths.toml", "dense_auto_order")
    mmd_oracle = expectations_oracle("expectations/trace_optimization_paths.toml", "mmd_order1")

    dense_cfg = SolverConfig(
        optimizer=SOLVER,
        order=0,  # exercise auto-order path
        cs_algo=NoElimination(),
        ts_algo=NoElimination()
    )

    dense_dual = cs_nctssos(tpop, dense_cfg; dualize=true)
    @test dense_dual.objective ≈ dense_oracle.opt atol = 1e-5
    @test flatten_sizes(dense_dual.moment_matrix_sizes) == dense_oracle.sides
    @test dense_dual.n_unique_moment_matrix_elements == dense_oracle.nuniq

    dense_moment = cs_nctssos(tpop, dense_cfg; dualize=false)
    @test dense_moment.objective ≈ dense_oracle.opt atol = 1e-5
    @test dense_moment.n_unique_moment_matrix_elements == dense_oracle.nuniq

    mmd_cfg = SolverConfig(
        optimizer=SOLVER,
        order=1,
        cs_algo=NoElimination(),
        ts_algo=MMD()
    )
    sparse_dual = cs_nctssos(tpop, mmd_cfg; dualize=true)
    @test sparse_dual.objective ≈ mmd_oracle.opt atol = 1e-5
    @test flatten_sizes(sparse_dual.moment_matrix_sizes) == mmd_oracle.sides
    @test sparse_dual.n_unique_moment_matrix_elements == mmd_oracle.nuniq

    sparse_moment = cs_nctssos(tpop, mmd_cfg; dualize=false)
    @test sparse_moment.objective ≈ mmd_oracle.opt atol = 1e-5
    @test flatten_sizes(sparse_moment.moment_matrix_sizes) == mmd_oracle.sides
    @test sparse_moment.n_unique_moment_matrix_elements == mmd_oracle.nuniq

    @test occursin("State Optimization Result", sprint(show, dense_dual))
    @test occursin("Correlative Sparsity", sprint(show, dense_dual.sparsity.corr_sparsity))
    @test occursin("Number of Activated supp", sprint(show, dense_dual.sparsity.cliques_term_sparsities[1][1]))
end

@testset "Trace symmetry MVP is rejected cleanly" begin
    tpop = _build_chsh_trace_problem()
    symmetry = SymmetrySpec(SignedPermutation(1 => 2, 2 => 1))

    err = try
        cs_nctssos(
            tpop,
            SolverConfig(
                optimizer=SOLVER,
                order=1,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
                symmetry=symmetry,
            ),
        )
        nothing
    catch caught
        caught
    end

    @test err isa ArgumentError
    @test occursin(
        "state/trace polynomial optimization",
        sprint(showerror, err),
    )
end

@testset "Trace Equality Constraints" begin
    reg, (x,) = create_unipotent_variables([("x", 1:2)])
    objective = (-1.0 * tr(x[1] * x[2])) * one(typeof(x[1]))
    eq_con = (1.0 * tr(x[1])) * one(typeof(x[1]))
    tpop = polyopt(objective, reg; eq_constraints=[eq_con])

    config = SolverConfig(
        optimizer=SOLVER,
        order=1,
        cs_algo=NoElimination(),
        ts_algo=NoElimination()
    )

    result_dual = cs_nctssos(tpop, config; dualize=true)
    result_moment = cs_nctssos(tpop, config; dualize=false)

    @test isfinite(result_dual.objective)
    @test isfinite(result_moment.objective)

    mp = NCTSSoS.moment_relax(tpop, result_dual.sparsity.corr_sparsity, result_dual.sparsity.cliques_term_sparsities)
    @test any(c -> c[1] == :Zero, mp.constraints)
end

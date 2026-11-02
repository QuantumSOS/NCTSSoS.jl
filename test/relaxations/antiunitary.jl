using Test, NCTSSoS, JuMP, LinearAlgebra

@testset "Pauli antiunitary reduction" begin
    reg, (x,y,z) = create_pauli_variables(1:2)
    basis = [one(x[1]), x[1], y[1], z[1], y[2]]
    config = SolverConfig(optimizer=SOLVER, moment_basis=basis)
    pop = polyopt(1.0x[1] + 2z[1] + y[1]*y[2], reg)
    sparsity = compute_sparsity(pop, config)
    reduced = antiunitary_moment_relaxation(pop, sparsity)
    mp = reduced.moment_problem
    @test reduced.theta.x == [false,false]
    @test reduced.theta.z == [false,false]
    @test !antiunitary_parity(reduced.theta, x[1])
    @test antiunitary_parity(reduced.theta, y[1])
    @test all(!antiunitary_parity(reduced.theta,m) for m in mp.total_basis)
    @test size(mp.constraints[1][2]) == (5,5)
    @test mp.constraints[1][1] == :PSD
    @test mp.linear.psd_blocks_lin[1].meta.row_labels == reduced.original.linear.psd_blocks_lin[1].meta.row_labels
    @test all(mp.constraints[1][2][i,j] == mp.constraints[1][2][j,i] for i in 1:5, j in 1:5)
    # Cross-parity products survive: dropping these would change the feasible set.
    @test any(!iszero(mp.constraints[1][2][i,j]) for i in 1:5, j in 1:5
        if reduced.phases[1][i] != reduced.phases[1][j])
    for dualize in (false, true)
        baseline = NCTSSoS.solve_sdp(reduced.original, SOLVER; dualize)
        result = NCTSSoS.solve_sdp(mp, SOLVER; dualize)
        @test result.objective ≈ -sqrt(6) atol=2e-6
        @test result.objective ≈ baseline.objective atol=2e-6
        # Actual JuMP cone dimensions, not just symbolic metadata.
        sides = [JuMP.constraint_object(c).set.side_dimension for
            (F,S) in JuMP.list_of_constraint_types(result.model) if S <: JuMP.MOI.PositiveSemidefiniteConeTriangle
            for c in JuMP.all_constraints(result.model,F,S)]
        @test sides == [5]
    end

    # Recover a solver-produced certificate in the original, unprojected algebra.
    model = NCTSSoS.sos_dualize(mp).model
    set_optimizer(model, SOLVER)
    set_silent(model)
    optimize!(model)
    cones = [c for (F,S) in list_of_constraint_types(model) if S <: JuMP.MOI.PositiveSemidefiniteConeTriangle
        for c in all_constraints(model,F,S)]
    Z = Matrix(JuMP.reshape_vector(value.(constraint_object(only(cones)).func), JuMP.shape(constraint_object(only(cones)))))
    Q = recover_antiunitary_gram(reduced, 1, Z)
    @test eigmin(Hermitian(Q)) >= -2e-6
    original_mat = reduced.original.constraints[1][2]
    certificate = sum(Q[j,i] * original_mat[i,j] for i in 1:5, j in 1:5)
    residual = pop.objective - objective_value(model) - certificate
    @test maximum(abs, coefficients(residual); init=0.0) < 3e-6

    # Exhaust every one-qubit support, including the inconsistent XYZ case.
    for mask in 0:7
        h = sum((Float64((mask >> (i-1)) & 1) * basis[i+1] for i in 1:3))
        p = polyopt(h, reg)
        theta = detect_pauli_antiunitary(p)
        @test isnothing(theta) == (mask == 7)
        if !isnothing(theta)
            @test all(!antiunitary_parity(theta,m) for (c,m) in terms(h) if !iszero(c))
        end
    end
    bad = polyopt(1.0x[1] + z[1], reg; ineq_constraints=[1.0y[1]])
    @test isnothing(detect_pauli_antiunitary(bad))
    @test_throws ArgumentError antiunitary_moment_relaxation(bad, compute_sparsity(bad,config))
    @test_throws ArgumentError antiunitary_moment_relaxation(pop,sparsity; theta=PauliAntiunitary([1,0],[1,0]))
    @test isnothing(detect_pauli_antiunitary(polyopt(1.0x[1]+z[1],reg; eq_constraints=[1.0y[1]])))
    @test isnothing(detect_pauli_antiunitary(polyopt(1.0x[1]+z[1],reg; moment_eq_constraints=[1.0y[1]])))
    @test_throws ArgumentError detect_pauli_antiunitary(polyopt(1.0x[1],reg; eq_constraints=[x[1]+im*y[1]]))

    # Θ²=-I is valid for moment optimization (invariant states may be mixed).
    reg3, ops = create_pauli_variables(1:3)
    h = heisenberg_chain_hamiltonian(ops; periodic=false, coupling=1.0)
    p3 = polyopt(h,reg3)
    s3 = compute_sparsity(p3,SolverConfig(optimizer=SOLVER,order=1))
    minus = PauliAntiunitary(ones(Int,3),ones(Int,3))
    a3 = antiunitary_moment_relaxation(p3,s3; theta=minus)
    @test isodd(sum(a3.theta.x .* a3.theta.z))
    @test NCTSSoS.solve_sdp(a3.moment_problem,SOLVER).objective ≈ NCTSSoS.solve_sdp(a3.original,SOLVER).objective atol=3e-6
    # All weight-two axes force time reversal on every site, not just K.
    all_pairs = sum(1.0u[i]*v[j] for u in ops, v in ops, i in 1:3, j in 1:3 if i < j)
    detected_minus = detect_pauli_antiunitary(polyopt(all_pairs,reg3))
    @test detected_minus.x == trues(3)
    @test detected_minus.z == trues(3)

    # Individually fixed localizers, equalities, and sparse/term-sparse word bases.
    reg2, (a,b,c) = create_pauli_variables(1:2)
    constrained = polyopt(1.0a[1] + 0.3c[2],reg2;
        ineq_constraints=[1.0 + c[1], 2.0 + b[1]*b[2]], moment_eq_constraints=[1.0 - c[2]])
    for (cs,ts) in ((NoElimination(),NoElimination()),(MF(),MF()))
        s = compute_sparsity(constrained,SolverConfig(optimizer=SOLVER,order=2,cs_algo=cs,ts_algo=ts))
        a = antiunitary_moment_relaxation(constrained,s)
        @test length(a.phases) == length(a.original.linear.psd_blocks_lin)
        @test all(size(new.entries) == size(old.entries) for (new,old) in
            zip(a.moment_problem.linear.psd_blocks_lin,a.original.linear.psd_blocks_lin))
        @test NCTSSoS.solve_sdp(a.moment_problem,SOLVER).objective ≈ NCTSSoS.solve_sdp(a.original,SOLVER).objective atol=3e-6
        model = NCTSSoS.sos_dualize(a.moment_problem).model
        set_optimizer(model,SOLVER)
        set_silent(model)
        optimize!(model)
        cones = [q for (F,S) in list_of_constraint_types(model) if S <: JuMP.MOI.PositiveSemidefiniteConeTriangle
            for q in all_constraints(model,F,S)]
        gram_vars = Set(v for q in cones for v in constraint_object(q).func)
        equality_vars = filter(v -> !(v in gram_vars),all_variables(model))
        certificate = zero(a.original.objective)
        for (k,q) in enumerate(cones)
            Q = recover_antiunitary_gram(a,k,Matrix(JuMP.reshape_vector(value.(constraint_object(q).func),JuMP.shape(constraint_object(q)))))
            @test eigmin(Hermitian(Q)) >= -3e-6
            mat = a.original.constraints[a.original.linear.psd_block_constraint_idx[k]][2]
            certificate += sum(Q[j,i]*mat[i,j] for i in axes(mat,1),j in axes(mat,2))
        end
        @test length(equality_vars) == length(a.moment_problem.linear.zero_constraints)
        for (v,eq) in zip(equality_vars,a.moment_problem.linear.zero_constraints)
            certificate += value(v)*sum(coef*a.moment_problem.linear.key_to_monomial[key] for (key,coef) in eq.form)
        end
        residual = constrained.objective-objective_value(model)-certificate
        @test maximum(abs,coefficients(residual);init=0.0) < 5e-6
    end
end

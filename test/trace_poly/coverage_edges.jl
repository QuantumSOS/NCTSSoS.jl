using Graphs
using JuMP
const MOI = JuMP.MOI

_idx(op_id::Integer, site::Integer=1) = NCTSSoS.encode_index(UInt16, op_id, site)
_mono(ids::Integer...; site::Integer=1) =
    NormalMonomial{NonCommutativeAlgebra,UInt16}(UInt16[_idx(i, site) for i in ids])

@testset "State API Coverage Edges" begin
    m1 = _mono(1)
    m2 = _mono(2)

    sym1 = StateSymbol{Arbitrary}(m1)
    sym2 = StateSymbol{Arbitrary}(m2)
    @test one(sym1) == one(StateSymbol{Arbitrary,NonCommutativeAlgebra,UInt16})
    @test sprint(show, Arbitrary()) == "Arbitrary()"
    @test sprint(show, MaxEntangled()) == "MaxEntangled()"

    tr_sym = StateSymbol{MaxEntangled}(m1)
    @test startswith(sprint(show, tr_sym), "tr(")

    sw_from_sym = StateWord{Arbitrary}(sym1)
    sw_from_vec = StateWord{Arbitrary}([sym1, sym2])
    @test length(sw_from_vec.state_syms) == 2
    @test one(sw_from_sym) == one(StateWord{Arbitrary,NonCommutativeAlgebra,UInt16})

    @test (sym1 * sw_from_sym) isa StateWord
    @test (sw_from_sym * sym1) isa StateWord

    ncsw = NCStateWord(sw_from_sym, m2)
    @test symmetric_canon(ncsw) == symmetric_canon(expval(ncsw))
    @test one(ncsw) == one(NCStateWord{Arbitrary,NonCommutativeAlgebra,UInt16})
    @test varsigma(m1) == ς(m1)

    _, (σx, _, _) = create_pauli_variables(1:1)
    pauli_mono = first(monomials(σx[1]))
    pauli_sym = StateSymbol{Arbitrary}(pauli_mono)
    @test adjoint(pauli_sym) === pauli_sym

    _, (a, _) = create_fermionic_variables(1:1)
    fermionic_mono = first(monomials(a[1]))
    fermionic_sym = StateSymbol{Arbitrary}(fermionic_mono)
    @test_throws ErrorException adjoint(fermionic_sym)
end

@testset "StatePolynomial Operator Edges" begin
    m1 = _mono(1)
    m2 = _mono(2)
    sym1 = StateSymbol{Arbitrary}(m1)
    sym2 = StateSymbol{Arbitrary}(m2)
    sw1 = StateWord{Arbitrary}(sym1)
    sw2 = StateWord{Arbitrary}(sym2)
    sp = StatePolynomial([2.0], [sw1])

    @test coefficients(sp) == [2.0]
    @test monomials(sp) == [sw1]
    @test collect(terms(sp)) == [(2.0, sw1)]

    @test_throws ArgumentError expval(MaxEntangled, m1)
    @test iszero(zero(sp))
    @test isone(one(sp))

    t = (3.0, sw2)
    @test (sw1 + t) isa StatePolynomial
    @test (t + sw1) isa StatePolynomial
    @test (sw1 - t) isa StatePolynomial
    @test (t - sw1) isa StatePolynomial
    @test (sp + t) isa StatePolynomial
    @test (t + sp) isa StatePolynomial
    @test (sp - t) isa StatePolynomial
    @test (sp + sw2) isa StatePolynomial
    @test (sp - sw2) isa StatePolynomial
    @test (sw2 + sp) isa StatePolynomial

    @test (sp * 2.0) isa StatePolynomial
    @test (sp * m1) isa NCStatePolynomial
    @test (m1 * sp) isa NCStatePolynomial

    @test (2.0 + sw1) isa StatePolynomial
    @test (sw1 + 2.0) isa StatePolynomial
    @test (2.0 - sw1) isa StatePolynomial
    @test (2.0 * sw1) isa StatePolynomial
    @test (sw1 * 2.0) isa StatePolynomial

    @test (2.0 * sym1) isa StatePolynomial
    @test (sym1 * 2.0) isa StatePolynomial
    @test (sym1 + sym2) isa StatePolynomial
    @test (sym1 - sym2) isa StatePolynomial
    @test (-sym1) isa StatePolynomial
    @test (2.0 + sym1) isa StatePolynomial
    @test (sym1 + 2.0) isa StatePolynomial
    @test (2.0 - sym1) isa StatePolynomial
    @test (sym1 - 2.0) isa StatePolynomial
    @test (sym1 + sw1) isa StatePolynomial
    @test (sw1 + sym1) isa StatePolynomial
    @test (sym1 - sw1) isa StatePolynomial
    @test (sw1 - sym1) isa StatePolynomial
    @test (sp + sym1) isa StatePolynomial
    @test (sym1 + sp) isa StatePolynomial
    @test (sp - sym1) isa StatePolynomial
    @test_throws MethodError sym1 - sp
    @test (sym1 * sp) isa StatePolynomial
    @test (sp * sym1) isa StatePolynomial
    @test (sw1 * sp) == (sp * sw1)

    @test sprint(show, zero(sp)) == "0"
    complex_sp = StatePolynomial(ComplexF64[1.0 + 0.0im, 2.0 + 3.0im], [sw1, sw2])
    @test !isempty(sprint(show, complex_sp))
end

@testset "NCStatePolynomial Operator Edges" begin
    m1 = _mono(1)
    m2 = _mono(2)
    sw1 = tr(m1)
    sw2 = tr(m2)
    ncsw1 = NCStateWord(sw1, m1)
    ncsw2 = NCStateWord(sw2, m2)
    ncsp = NCStatePolynomial([1.0], [ncsw1])

    @test coefficients(ncsp) == [1.0]
    @test monomials(ncsp) == [ncsw1]
    @test collect(terms(ncsp)) == [(1.0, ncsw1)]
    @test degree(ncsp) == degree(ncsw1)
    @test variable_indices(ncsp) == variables(ncsp)
    @test variable_indices(ncsw1) == variables(ncsw1)
    @test iszero(zero(ncsp))
    @test isone(one(ncsp))
    @test hash(ncsp) == hash(ncsp)

    @test (ncsp * 2.0) isa NCStatePolynomial
    @test (2.0 * ncsw1) isa NCStatePolynomial
    @test (ncsw1 * 2.0) isa NCStatePolynomial
    @test (ncsw1 + ncsw2) isa NCStatePolynomial
    @test (ncsw1 - ncsw2) isa NCStatePolynomial
    @test (-ncsw1) isa NCStatePolynomial
    @test (2.0 + ncsw1) isa NCStatePolynomial
    @test (ncsw1 + 2.0) isa NCStatePolynomial
    @test (2.0 - ncsw1) isa NCStatePolynomial
    @test (ncsw1 - 2.0) isa NCStatePolynomial
    @test (ncsp + ncsw2) isa NCStatePolynomial
    @test (ncsw2 + ncsp) isa NCStatePolynomial
    @test (ncsp - ncsw2) isa NCStatePolynomial
    @test (ncsw2 - ncsp) isa NCStatePolynomial

    empty_ncsp = zero(ncsp)
    @test iszero(expval(empty_ncsp))
    @test expval(ncsp) isa StatePolynomial
    @test sprint(show, empty_ncsp) == "0"
end

@testset "Optimization Helper Coverage Edges" begin
    @test NCTSSoS._is_complex_problem(PauliAlgebra)
    @test NCTSSoS._is_complex_problem(FermionicAlgebra)
    @test NCTSSoS._is_complex_problem(BosonicAlgebra)
    @test !NCTSSoS._is_complex_problem(NonCommutativeAlgebra)
    @test !NCTSSoS._is_complex_problem(ProjectorAlgebra)
    @test !NCTSSoS._is_complex_problem(UnipotentAlgebra)
    @test Base.invokelatest(NCTSSoS._is_complex_problem, PauliAlgebra)
    @test Base.invokelatest(NCTSSoS._is_complex_problem, FermionicAlgebra)
    @test Base.invokelatest(NCTSSoS._is_complex_problem, BosonicAlgebra)
    @test !Base.invokelatest(NCTSSoS._is_complex_problem, NonCommutativeAlgebra)
    @test !Base.invokelatest(NCTSSoS._is_complex_problem, ProjectorAlgebra)
    @test !Base.invokelatest(NCTSSoS._is_complex_problem, UnipotentAlgebra)

    reg_many, (x_many,) = create_noncommutative_variables([("x", 1:11)])
    pop_many = polyopt(
        1.0 * x_many[1],
        reg_many;
        eq_constraints=[1.0 * x_many[1]],
        ineq_constraints=[1.0 * x_many[2]]
    )
    pop_str = sprint(show, pop_many)
    @test occursin("Optimization Problem (NonCommutativeAlgebra)", pop_str)
    @test occursin("Equality constraints (1)", pop_str)
    @test occursin("Inequality constraints (1)", pop_str)
    @test occursin("...", pop_str)

    G = SimpleGraph(4)
    add_edge!(G, 1, 2)
    add_edge!(G, 3, 4)
    max_cliques = clique_decomp(G, MaximalElimination())
    @test sort(length.(max_cliques)) == [2, 2]

    err = NCTSSoS.SolverStatusError(MOI.OTHER_ERROR, MOI.NO_SOLUTION, MOI.NO_SOLUTION)
    @test occursin("Solver failed: termination=", sprint(showerror, err))

    reg_small, (x_small,) = create_noncommutative_variables([("x", 1:1)])
    pop_small = polyopt(1.0 * x_small[1] * x_small[1], reg_small; ineq_constraints=[1.0 * x_small[1]])
    config = SolverConfig(
        optimizer=SOLVER,
        order=1,
        cs_algo=NoElimination(),
        ts_algo=NoElimination()
    )
    res_small = cs_nctssos(pop_small, config; dualize=true)
    @test occursin("Objective:", sprint(show, res_small))
    @test occursin("Correlative Sparsity", sprint(show, res_small.sparsity.corr_sparsity))
    @test occursin("Number of Activated supp", sprint(show, res_small.sparsity.cliques_term_sparsities[1][1]))
end

function _find_rl_only_state_case(basis)
    for a in basis, b in basis, supp in basis
        a == b && continue
        connected_lr = simplify(NCTSSoS._neat_dot3(a, supp, b))
        connected_rl = simplify(NCTSSoS._neat_dot3(b, supp, a))
        lr_sws = Set(symmetric_canon(expval(ncsw)) for ncsw in monomials(connected_lr))
        for ncsw in monomials(connected_rl)
            sw = symmetric_canon(expval(ncsw))
            if !(sw in lr_sws)
                return (a, b, supp, ncsw)
            end
        end
    end
    return nothing
end

@testset "Sparsity and Moment Edge Branches" begin
    reg_nc, (x_nc,) = create_noncommutative_variables([("x", 1:2)])
    pop_nc = polyopt(1.0 * x_nc[1], reg_nc; eq_constraints=[1.0 * x_nc[1] * x_nc[1] * x_nc[2]])
    corr_nc = correlative_sparsity(pop_nc, 1, NoElimination())
    @test isempty(corr_nc.clq_localizing_mtx_bases[1][1])

    objective_state = (1.0 * tr(x_nc[1] * x_nc[2])) * one(typeof(x_nc[1]))
    eq_state = (1.0 * tr(x_nc[1] * x_nc[1] * x_nc[2])) * one(typeof(x_nc[1]))
    pop_state = polyopt(objective_state, reg_nc; eq_constraints=[eq_state])
    corr_state = correlative_sparsity(pop_state, 1, NoElimination())
    @test isempty(corr_state.clq_localizing_mtx_bases[1][1])

    compound_obj = (1.0 * (tr(x_nc[1]) * tr(x_nc[2]))) * one(typeof(x_nc[1]))
    compound_words = NCTSSoS._extract_compound_state_words(compound_obj, [compound_obj])
    @test any(sw -> length(sw.state_syms) > 1, compound_words)

    graph_state, sorted_indices, idx_to_node = NCTSSoS.get_state_correlative_graph(reg_nc, compound_obj, [eq_state])
    @test nv(graph_state) == length(sorted_indices) == length(idx_to_node)

    idx1 = only(variable_indices(x_nc[1]))
    idx2 = only(variable_indices(x_nc[2]))
    mixed_constraint = (1.0 * tr(x_nc[1] * x_nc[2])) * one(typeof(x_nc[1]))
    clq_cons, global_cons = NCTSSoS.assign_state_constraint([[idx1], [idx2]], [mixed_constraint], reg_nc)
    @test clq_cons == [Int[], Int[]]
    @test global_cons == [1]

    init_state = NCTSSoS.init_activated_supp(pop_state.objective, [eq_state], corr_state.clq_mom_mtx_bases[1])
    ts_state = NCTSSoS.term_sparsities(
        init_state,
        [eq_state],
        corr_state.clq_mom_mtx_bases[1],
        corr_state.clq_localizing_mtx_bases[1],
        NoElimination()
    )
    @test length(ts_state) == 1 + length(corr_state.clq_localizing_mtx_bases[1])

    basis_state = get_state_basis(reg_nc, 2; state_type=Arbitrary)
    rl_case = _find_rl_only_state_case(basis_state)
    @test rl_case !== nothing
    a_case, b_case, supp_case, activated_case = rl_case
    graph_rl = NCTSSoS.get_term_sparsity_graph([supp_case], [activated_case], [a_case, b_case])
    @test has_edge(graph_rl, 1, 2)

    @test NCTSSoS._conj_coef(PauliAlgebra, UInt8(1)) == UInt8(3)
    @test !NCTSSoS._has_odd_parity_only(1.0 * x_nc[1])
    @test !Base.invokelatest(NCTSSoS._has_odd_parity_only, 1.0 * x_nc[1])

    pauli_model = GenericModel{Float64}()
    @variable(pauli_model, y_re[1:1], set_string_name=false)
    @variable(pauli_model, y_im[1:1], set_string_name=false)
    pauli_zero_poly = zero(Polynomial{PauliAlgebra,UInt8,ComplexF64})
    id_pauli = symmetric_canon(expval(one(NormalMonomial{PauliAlgebra,UInt8})))
    re0, im0 = NCTSSoS._substitute_complex_poly(pauli_zero_poly, Dict(id_pauli => 1), y_re, y_im)
    @test iszero(re0)
    @test iszero(im0)

    reg_pauli, (σx, _, _) = create_pauli_variables(1:1)
    pauli_nonbasis_poly = (1.0 + 1.0im) * σx[1]
    re_skip, im_skip = NCTSSoS._substitute_complex_poly(pauli_nonbasis_poly, Dict(id_pauli => 1), y_re, y_im)
    @test iszero(re_skip)
    @test iszero(im_skip)
    pop_pauli = polyopt(1.0 * σx[1], reg_pauli; eq_constraints=[1.0 * σx[1]])
    sparsity_pauli = compute_sparsity(pop_pauli, SolverConfig(optimizer=SOLVER, order=1))
    mp_pauli = NCTSSoS.moment_relax(pop_pauli, sparsity_pauli.corr_sparsity, sparsity_pauli.cliques_term_sparsities)
    @test any(c -> c[1] == :Zero, mp_pauli.constraints)
    pauli_result = NCTSSoS.solve_moment_problem(mp_pauli, SOLVER)
    @test isfinite(real(pauli_result.objective))

    sos_pauli = NCTSSoS.sos_dualize(mp_pauli)
    @test sos_pauli.n_unique_elements > 0

    bad_mp_pauli = NCTSSoS.MomentProblem(
        mp_pauli.objective,
        [(:PSD, mp_pauli.constraints[1][2])],
        mp_pauli.total_basis,
        mp_pauli.n_unique_moment_matrix_elements
    )
    @test_throws ErrorException NCTSSoS._solve_complex_moment_problem(bad_mp_pauli, SOLVER, true)

    mp_state_manual = NCTSSoS.moment_relax(pop_state, corr_state, [ts_state])
    @test !isempty(mp_state_manual.constraints)
    pop_state_len0 = polyopt(
        (1.0 * tr(x_nc[1])) * one(typeof(x_nc[1])),
        reg_nc;
        eq_constraints=[(1.0 * tr(x_nc[1] * x_nc[1] * x_nc[1])) * one(typeof(x_nc[1]))]
    )
    corr_state_len0 = correlative_sparsity(pop_state_len0, 1, NoElimination())
    @test isempty(corr_state_len0.clq_localizing_mtx_bases[1][1])

    M_state = eltype(corr_state.clq_mom_mtx_bases[1])
    P_state = typeof(eq_state)
    manual_corr_state = CorrelativeSparsity{NonCommutativeAlgebra,UInt8,P_state,M_state,MaxEntangled}(
        [[idx1]],
        reg_nc,
        [eq_state],
        [Int[]],
        [1],
        [[one(M_state)]],
        [Vector{M_state}[]]
    )
    manual_ts = [[NCTSSoS.TermSparsity([one(M_state)], [[one(M_state)]])]]
    mp_state_global = NCTSSoS.moment_relax(pop_state, manual_corr_state, manual_ts)
    @test any(c -> c[1] == :Zero, mp_state_global.constraints)

    sos_state = NCTSSoS.sos_dualize(mp_state_manual)
    @test sos_state.n_unique_elements > 0

    bad_mp_state = NCTSSoS.StateMomentProblem(
        mp_state_manual.objective,
        [(:HPSD, mp_state_manual.constraints[1][2], mp_state_manual.constraints[1][3])],
        mp_state_manual.total_basis,
        mp_state_manual.n_unique_moment_matrix_elements
    )
    @test_throws ErrorException NCTSSoS.solve_moment_problem(bad_mp_state, SOLVER)
    bad_identity_mp = NCTSSoS.StateMomentProblem(
        mp_state_manual.objective,
        mp_state_manual.constraints,
        eltype(mp_state_manual.total_basis)[],
        mp_state_manual.n_unique_moment_matrix_elements
    )
    @test_throws ErrorException NCTSSoS.solve_moment_problem(bad_identity_mp, SOLVER)

    sw_identity = one(StateWord{MaxEntangled,NonCommutativeAlgebra,UInt16})
    zero_expr = NCTSSoS._substitute_state_poly(
        zero(NCStatePolynomial{Float64,MaxEntangled,NonCommutativeAlgebra,UInt16}),
        Dict(sw_identity => 1),
        [1.0]
    )
    @test zero_expr == 0.0

    non_basis_ncsp = (1.0 * tr(_mono(1) * _mono(2))) * one(_mono(1))
    skipped_expr = NCTSSoS._substitute_state_poly(non_basis_ncsp, Dict(sw_identity => 1), [1.0])
    @test skipped_expr == 0.0
end

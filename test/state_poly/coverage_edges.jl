# test/state_poly/coverage_edges.jl

using Test, NCTSSoS, JuMP
using NCTSSoS: NCStateWord, NCStatePolynomial, Arbitrary, MaxEntangled
using Graphs

@testset "Coverage Edges" begin
    @testset "states/types.jl show" begin
        @test sprint(show, Arbitrary()) == "Arbitrary()"
        @test sprint(show, MaxEntangled()) == "MaxEntangled()"
    end

    @testset "states/word.jl misc branches" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))

        sym = StateSymbol{Arbitrary}(m)
        @test one(sym) == one(StateSymbol{Arbitrary,NonCommutativeAlgebra,UInt16})

        sym_tr = StateSymbol{MaxEntangled}(m)
        @test occursin("tr(", sprint(show, sym_tr))

        # Adjoint: TwistedGroupAlgebra returns self, PBWAlgebra throws
        _, (σx, _, _) = create_pauli_variables(1:1)
        sym_pauli = StateSymbol{Arbitrary}(σx[1])
        @test adjoint(sym_pauli) === sym_pauli

        _, (a, _) = create_fermionic_variables(1:1)
        sym_ferm = StateSymbol{Arbitrary}(a[1])
        @test_throws ErrorException adjoint(sym_ferm)

        sw_from_sym = StateWord{Arbitrary}(sym)
        sw_from_syms = StateWord{Arbitrary}([sym])
        @test sw_from_sym == sw_from_syms
        @test one(sw_from_sym) == one(StateWord{Arbitrary,NonCommutativeAlgebra,UInt16})

        # StateSymbol <-> StateWord multiplication overloads
        sw_double = StateWord{Arbitrary}([sym, sym])
        @test sym * sw_from_sym == sw_double
        @test sw_from_sym * sym == sw_double

        # NCStateWord canonicalization + one(ncsw)
        ncsw = NCStateWord(one(StateWord{Arbitrary,NonCommutativeAlgebra,UInt16}), m)
        @test symmetric_canon(ncsw) isa StateWord
        @test one(ncsw) == one(NCStateWord{Arbitrary,NonCommutativeAlgebra,UInt16})

        # varsigma ASCII alias
        @test varsigma(m) == ς(m)
    end

    @testset "states/polynomial.jl accessors + lifted arithmetic" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))
        sw1 = ς(m1)
        sw2 = ς(m2)

        # Accessors
        sp = StatePolynomial([1.0, 2.0], [sw1, sw2])
        @test coefficients(sp) == sp.coeffs
        @test monomials(sp) == sp.state_words
        @test collect(terms(sp)) == collect(zip(sp.coeffs, sp.state_words))

        # expval(::Type{ST}, ::NormalMonomial) is rejected (use expval(m) / tr(m))
        @test_throws ArgumentError NCTSSoS.expval(MaxEntangled, m1)

        # instance zero/one
        @test zero(sp) == zero(typeof(sp))
        @test one(sp) == one(typeof(sp))

        # StateWord +/- Tuple and Tuple +/- StateWord
        t = (2.0, sw2)
        @test sw1 + t == StatePolynomial([1.0, 2.0], [sw1, sw2])
        @test t + sw1 == StatePolynomial([1.0, 2.0], [sw1, sw2])
        @test sw1 - t == StatePolynomial([1.0, -2.0], [sw1, sw2])
        @test t - sw1 == StatePolynomial([-1.0, 2.0], [sw1, sw2])

        # StatePolynomial +/- Tuple and Tuple + StatePolynomial
        @test sp + t == StatePolynomial([1.0, 4.0], [sw1, sw2])
        @test t + sp == StatePolynomial([1.0, 4.0], [sw1, sw2])
        @test sp - t == StatePolynomial([1.0], [sw1])

        # StateWord + StatePolynomial
        @test sw1 + sp == StatePolynomial([2.0, 2.0], [sw1, sw2])

        # Monomial * StatePolynomial (lifts to NCStatePolynomial)
        @test m1 * sp == NCStatePolynomial([1.0, 2.0], [NCStateWord(sw1, m1), NCStateWord(sw2, m1)])

        # Scalar +/-/* StateWord
        id_sw = one(StateWord{Arbitrary,NonCommutativeAlgebra,UInt16})
        @test 2 + sw1 == StatePolynomial([2.0, 1.0], [id_sw, sw1])
        @test sw1 + 2 == StatePolynomial([2.0, 1.0], [id_sw, sw1])
        @test 2 - sw1 == StatePolynomial([2.0, -1.0], [id_sw, sw1])
        @test 2 * sw1 == StatePolynomial([2.0], [sw1])
        @test sw1 * 2 == StatePolynomial([2.0], [sw1])

        # StateSymbol lifting helpers
        sym1 = StateSymbol{Arbitrary}(m1)
        sym2 = StateSymbol{Arbitrary}(m2)
        @test 2 * sym1 == StatePolynomial([2.0], [sw1])
        @test sym1 * 2 == StatePolynomial([2.0], [sw1])
        @test sym1 + sym2 == StatePolynomial([1.0, 1.0], [sw1, sw2])
        @test sym1 - sym2 == StatePolynomial([1.0, -1.0], [sw1, sw2])
        @test -sym1 == StatePolynomial([-1.0], [sw1])
        @test 2 + sym1 == StatePolynomial([2.0, 1.0], [id_sw, sw1])
        @test sym1 + 2 == StatePolynomial([2.0, 1.0], [id_sw, sw1])
        @test 2 - sym1 == StatePolynomial([2.0, -1.0], [id_sw, sw1])
        @test sym1 - 2 == StatePolynomial([-2.0, 1.0], [id_sw, sw1])
        @test sym1 + sw1 == StatePolynomial([2.0], [sw1])
        @test sw1 + sym1 == StatePolynomial([2.0], [sw1])
        @test iszero(sym1 - sw1)
        @test iszero(sw1 - sym1)
        @test sp + sym1 == StatePolynomial([2.0, 2.0], [sw1, sw2])
        @test sym1 + sp == StatePolynomial([2.0, 2.0], [sw1, sw2])
        @test sp - sym1 == StatePolynomial([2.0], [sw2])
        @test sym1 - sp == StatePolynomial([-2.0], [sw2])
        sw1sw1 = sw1 * sw1
        sw1sw2 = sw1 * sw2
        @test sym1 * sp == StatePolynomial([1.0, 2.0], [sw1sw1, sw1sw2])
        @test sp * sym1 == StatePolynomial([1.0, 2.0], [sw1sw1, sw1sw2])

        # show: cover non-first, non-unit, positive coefficient branch
        str = sprint(show, sp)
        @test occursin(" + 2", str)

        # show: empty polynomial prints 0
        sp0 = zero(typeof(sp))
        @test sprint(show, sp0) == "0"

        # NCStatePolynomial accessors/ops not hit elsewhere
        ncsw = NCStateWord(sw1, m2)
        ncsp = NCStatePolynomial([1.0], [ncsw])
        @test collect(terms(ncsp)) == collect(zip(coefficients(ncsp), monomials(ncsp)))
        @test zero(ncsp) == zero(typeof(ncsp))
        @test hash(ncsp) == hash(ncsp)

        # NCStateWord scalar/arith helpers
        id_ncsw = one(typeof(ncsw))
        @test ncsw * 2 == NCStatePolynomial([2.0], [ncsw])
        @test ncsw + ncsw == NCStatePolynomial([2.0], [ncsw])
        @test iszero(ncsw - ncsw)
        @test -ncsw == NCStatePolynomial([-1.0], [ncsw])
        @test 2 + ncsw == NCStatePolynomial([2.0, 1.0], [id_ncsw, ncsw])
        @test ncsw + 2 == NCStatePolynomial([2.0, 1.0], [id_ncsw, ncsw])
        @test 2 - ncsw == NCStatePolynomial([2.0, -1.0], [id_ncsw, ncsw])
        @test ncsw - 2 == NCStatePolynomial([-2.0, 1.0], [id_ncsw, ncsw])
        @test ncsp + ncsw == NCStatePolynomial([2.0], [ncsw])
        @test ncsw + ncsp == NCStatePolynomial([2.0], [ncsw])
        @test iszero(ncsp - ncsw)
        @test iszero(ncsw - ncsp)
        @test NCTSSoS.expval(ncsp) == StatePolynomial([1.0], [NCTSSoS.expval(ncsw)])

        # show: empty NCStatePolynomial prints 0
        @test sprint(show, zero(typeof(ncsp))) == "0"

        # Scalar on right: ensure `*(ncsp::NCStatePolynomial, c::Number)` is exercised
        @test Base.invokelatest(*, ncsp, 2.0) == NCStatePolynomial([2.0], [ncsw])
    end

    @testset "optimization/problem.jl show + traits" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        pop = polyopt(1.0 * x[1], reg; eq_constraints=[1.0 * x[2]], ineq_constraints=[1.0 * x[3]])
        @test occursin("Optimization Problem", sprint(show, pop))

        reg_large, (xl,) = create_noncommutative_variables([("x", 1:11)])
        pop_large = polyopt(1.0 * xl[1], reg_large)
        @test occursin("Variables (11)", sprint(show, pop_large))

        @test NCTSSoS._is_complex_problem(PauliAlgebra) == true
        @test NCTSSoS._is_complex_problem(FermionicAlgebra) == true
        @test NCTSSoS._is_complex_problem(BosonicAlgebra) == true
        @test NCTSSoS._is_complex_problem(NonCommutativeAlgebra) == false
        @test NCTSSoS._is_complex_problem(ProjectorAlgebra) == false
        @test NCTSSoS._is_complex_problem(UnipotentAlgebra) == false
    end

    @testset "optimization/interface.jl show + showerror" begin
        err = NCTSSoS.SolverStatusError(MOI.NUMERICAL_ERROR, MOI.NO_SOLUTION, MOI.NO_SOLUTION)
        @test occursin("Solver failed:", sprint(showerror, err))
    end

    @testset "optimization/moment.jl and glue coverage" begin
        # _conj_coef special-case for Pauli phases
        @test NCTSSoS._conj_coef(PauliAlgebra, 0x01) == 0x03

        # _has_odd_parity_only fallback for non-fermionic problems
        reg_nc, (xn,) = create_noncommutative_variables([("x", 1:1)])
        poly_any = 1.0 * xn[1]
        @test NCTSSoS._has_odd_parity_only(poly_any) == false

        # _substitute_complex_poly: zero + missing-basis skip
        reg_p, (σx, _, _) = create_pauli_variables(1:1)
        model = Model()
        @variable(model, y_re[1:1])
        @variable(model, y_im[1:1])
        poly_zero = zero(typeof(1.0 * σx[1]))
        z_re, z_im = NCTSSoS._substitute_complex_poly(poly_zero, Dict{Any,Int}(), y_re, y_im)
        @test JuMP.coefficient(z_re, y_re[1]) == 0
        @test JuMP.coefficient(z_im, y_im[1]) == 0

        _re, _im = NCTSSoS._substitute_complex_poly(1.0im * σx[1], Dict{Any,Int}(), y_re, y_im)
        @test JuMP.coefficient(_re, y_re[1]) == 0
        @test JuMP.coefficient(_im, y_im[1]) == 0

        # State: build a tiny feasible problem; exercise term/correlative sparsity with constraints
        reg_u, (u,) = create_unipotent_variables([("u", 1:1)])
        TU = eltype(u[1].word)
        base = (1.0 * ς(u[1])) * one(typeof(u[1]))
        obj0 = 0.0 * base
        eq = 1.0 * base
        spop = polyopt(obj0, reg_u; eq_constraints=[eq])

        config = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )

        sparsity = compute_sparsity(spop, config)
        mp_state = NCTSSoS.moment_relax(spop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)

        # Force the state-moment global-constraint branch (normally empty due to clique-connecting)
        push!(sparsity.corr_sparsity.global_cons, 1)
        mp_state_global = NCTSSoS.moment_relax(spop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        pop!(sparsity.corr_sparsity.global_cons)

        global_poly = sparsity.corr_sparsity.cons[1]
        @test any(size(mat) == (1, 1) && mat[1, 1] == global_poly for (_cone, mat, _basis) in mp_state_global.constraints)

        # Direct solve (dualize=false) to cover solve_moment_problem(::StateMomentProblem)
        sol_state = NCTSSoS.solve_sdp(mp_state, config.optimizer; dualize=false)
        res_state = NCTSSoS.PolyOptResult(sol_state.objective, sparsity, sol_state.model, sol_state.n_unique_elements)
        @test isfinite(res_state.objective)
        @test sprint(show, res_state) isa String

        # State solve: missing identity branch
        one_sw = one(StateWord{Arbitrary,UnipotentAlgebra,TU})
        ncsw_no_id = NCStateWord(one_sw, u[1])
        obj_bad = NCStatePolynomial([1.0], [ncsw_no_id])
        cons_bad = Tuple{Symbol, Matrix{typeof(obj_bad)}, Vector{typeof(ncsw_no_id)}}[]
        mp_no_id = NCTSSoS.StateMomentProblem{UnipotentAlgebra,Arbitrary,TU,typeof(ncsw_no_id),typeof(obj_bad)}(
            obj_bad, cons_bad, [ncsw_no_id], 1
        )
        err_no_id = try
            NCTSSoS.solve_moment_problem(mp_no_id, SOLVER)
            nothing
        catch e
            e
        end
        @test err_no_id isa ErrorException
        @test occursin("Identity StateWord not found", sprint(showerror, err_no_id))

        # State solve: unexpected cone branch
        bad_constraints = copy(mp_state.constraints)
        bad_constraints[1] = (:BadCone, bad_constraints[1][2], bad_constraints[1][3])
        bad_mp_state = typeof(mp_state)(mp_state.objective, bad_constraints, mp_state.total_basis, mp_state.n_unique_moment_matrix_elements)
        err_bad_cone = try
            NCTSSoS.solve_moment_problem(bad_mp_state, SOLVER)
            nothing
        catch e
            e
        end
        @test err_bad_cone isa ErrorException
        @test occursin("Unexpected cone type", sprint(showerror, err_bad_cone))

        # Complex direct solve with a Zero cone constraint to hit :Zero handling
        pop_p = polyopt(0.0 * σx[1], reg_p; eq_constraints=[1.0 * σx[1]])
        sparsity_p = compute_sparsity(pop_p, config)
        mp_p = NCTSSoS.moment_relax(pop_p, sparsity_p.corr_sparsity, sparsity_p.cliques_term_sparsities)
        sol_p = NCTSSoS.solve_sdp(mp_p, config.optimizer; dualize=false)
        res_p = NCTSSoS.PolyOptResult(sol_p.objective, sparsity_p, sol_p.model, sol_p.n_unique_elements)
        @test sprint(show, res_p) isa String

        # Complex solve: unexpected cone branch
        bad_constraints_p = copy(mp_p.constraints)
        bad_constraints_p[1] = (:BadCone, bad_constraints_p[1][2])
        bad_mp_p = typeof(mp_p)(mp_p.objective, bad_constraints_p, mp_p.total_basis, mp_p.n_unique_moment_matrix_elements)
        err_bad_p = try
            NCTSSoS._solve_complex_moment_problem(bad_mp_p, SOLVER, true)
            nothing
        catch e
            e
        end
        @test err_bad_p isa ErrorException
        @test occursin("Unexpected cone type", sprint(showerror, err_bad_p))
    end

    @testset "optimization/sos.jl :Zero branches" begin
        # State SOS dualization: need at least one :Zero constraint
        reg_u, (u,) = create_unipotent_variables([("u", 1:1)])
        base = (1.0 * ς(u[1])) * one(typeof(u[1]))
        obj0 = 0.0 * base
        eq = 1.0 * base
        spop = polyopt(obj0, reg_u; eq_constraints=[eq])
        config = SolverConfig(optimizer=SOLVER, order=1, cs_algo=NoElimination(), ts_algo=NoElimination())
        sparsity = compute_sparsity(spop, config)
        mp_state = NCTSSoS.moment_relax(spop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        @test NCTSSoS.sos_dualize(mp_state) isa NCTSSoS.SOSProblem

        # Hermitian SOS dualization (Pauli) with :Zero constraint
        reg_p, (σx, _, _) = create_pauli_variables(1:1)
        pop_p = polyopt(0.0 * σx[1], reg_p; eq_constraints=[1.0 * σx[1]])
        sparsity_p = compute_sparsity(pop_p, config)
        mp_p = NCTSSoS.moment_relax(pop_p, sparsity_p.corr_sparsity, sparsity_p.cliques_term_sparsities)
        @test NCTSSoS.sos_dualize(mp_p) isa NCTSSoS.SOSProblem
    end

    @testset "optimization/sparsity.jl state-only helpers" begin
        reg_u, (u,) = create_unipotent_variables([("u", 1:2)])
        sw_compound = ς(u[1]) * ς(u[2])
        obj = NCStatePolynomial([1.0], [NCStateWord(sw_compound, one(typeof(u[1])))])

        compound = NCTSSoS._extract_compound_state_words(obj, [obj])
        @test sw_compound in compound

        # Localizing-basis truncation branch (len == 0) for *regular* correlative sparsity
        reg_nc, (x,) = create_noncommutative_variables([("x", 1:3)])
        obj_nc = 0.0 * x[1]
        con_hi_nc = 1.0 * x[1] * x[2] * x[3]  # degree 3 -> reduced order negative at order=1
        pop_hi_nc = polyopt(obj_nc, reg_nc; eq_constraints=[con_hi_nc])
        cs_nc = NCTSSoS.correlative_sparsity(pop_hi_nc, 1, NoElimination())
        @test any(isempty(b) for clique_bases in cs_nc.clq_localizing_mtx_bases for b in clique_bases)

        # Localizing-basis truncation branch (len == 0) for state correlative sparsity
        reg3, (v,) = create_unipotent_variables([("v", 1:3)])
        obj1 = (1.0 * ς(v[1])) * one(typeof(v[1]))
        # degree 3 operator part -> reduced order negative for order=1
        v123 = NormalMonomial{UnipotentAlgebra,eltype(v[1].word)}(vcat(v[1].word, v[2].word, v[3].word))
        con_hi = NCStatePolynomial([1.0], [NCStateWord(ς(one(typeof(v[1]))), v123)])
        pop_hi = polyopt(obj1, reg3; eq_constraints=[con_hi])
        cs = NCTSSoS.correlative_sparsity(pop_hi, 1, NoElimination())
        @test any(isempty(b) for clique_bases in cs.clq_localizing_mtx_bases for b in clique_bases)

        # get_term_sparsity_graph: force connected_rl path (lines 825-826)
        sw_id = one(StateWord{Arbitrary,NonCommutativeAlgebra,UInt16})
        b1 = NCStateWord(sw_id, NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1)))
        b2 = NCStateWord(sw_id, NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2)))
        supp = NCStateWord(sw_id, NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 4)))

        target = NCStateWord(sw_id, NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 3, 4, 1)))
        G = NCTSSoS.get_term_sparsity_graph([supp], [target], [b1, b2])
        @test has_edge(G, 1, 2)
    end
end

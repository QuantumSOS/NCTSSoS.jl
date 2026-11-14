# Tests of the redesign work area `NCTSSoS.Next`.
using Test, NCTSSoS, LinearAlgebra, JuMP, Random

isdefined(Main, :SOLVER) || include("../TestUtils.jl")
include("contract_checks.jl")
include("reference.jl")
include("fixtures.jl")

function solve_reference(program)
    model, y = Nx.to_jump(program)
    set_optimizer(model, SOLVER)
    optimize!(model)
    return objective_value(model), Nx.moment_values(program, y)
end

@testset "Next: contracts and references" begin
    @testset "Phase" begin
        @test [ComplexF64(1, 2) * im^k for k in 0:3] == [Nx.scale(Nx.Phase(k), ComplexF64(1, 2)) for k in 0:3]
        @test Nx.scale(Nx.Phase(2), 3.0) === -3.0
        @test_throws InexactError Nx.scale(Nx.Phase(1), 3.0)
        @test Nx.Phase(3) * Nx.Phase(3) == Nx.Phase(2)
        @test isone(Nx.Phase(1) * conj(Nx.Phase(1)))
    end

    @testset "contract 1 on the toy Pauli word" begin
        rng = MersenneTwister(1)
        samples = unique!([toy_words(3, 2)[1:12]; [ToyPauli(rand(rng, UInt64), rand(rng, UInt64)) for _ in 1:12]])
        check_word_contract(ToyPauli, samples)
        @test Nx.act!(TermDict{ToyPauli,ComplexF64}(), ToyShift(4), toy_letter(4, 1), 2.0 + 0im).terms == Dict(toy_letter(1, 1) => 2.0 + 0im)
        @test !Nx.is_antiunitary(ToyShift(4))
    end

    @testset "toy products equal the current package" begin
        reg, _ = create_pauli_variables(1:6)
        T = eltype(indices(reg))
        words = toy_words(6, 2)
        rng = MersenneTwister(2)
        for _ in 1:500
            a, b = rand(rng, words), rand(rng, words)
            @test all(w == Nx.word_from_legacy(ToyPauli, Nx.legacy_letters(w)) for w in (a, b))
            @test Nx.legacy_terms(a) == [(1.0 + 0im, Nx.legacy_letters(a))]
            @test Nx.terms_from_legacy!(TermDict{ToyPauli,ComplexF64}(), ToyPauli, Nx.legacy_letters(a), 2.0 + 0im).terms == Dict(a => 2.0 + 0im)
            new = Dict(Nx.legacy_letters(w) => c for (w, c) in product_terms(a, b, ComplexF64).terms)
            @test new == legacy_product_terms(PauliAlgebra, T, Nx.legacy_letters(a), Nx.legacy_letters(b))
        end
    end

    @testset "contract 2: format checks" begin
        p = Nx.MomentProgram{Float64}(2)
        Nx.push_term!(p.forms, 1.0, 1); Nx.close_form!(p.forms)
        Nx.push_term!(p.forms, 1.0, 2); Nx.close_form!(p.forms)
        Nx.push_term!(p.forms, 1.0, 1); Nx.close_form!(p.forms)
        @test_throws ArgumentError Nx.validate(p)                       # forms without a block
        @test_throws ArgumentError Nx.push_block!(p, Nx.PSD, 3)         # 6 forms needed, 3 closed
        Nx.push_block!(p, Nx.PSD, 2)
        @test Nx.validate(p) === p
        @test Nx.tri_offset(3, 3, 1) == 2 && Nx.tri_offset(3, 2, 2) == 3 && Nx.tri_offset(3, 3, 3) == 5
        Nx.push_term!(p.forms, 1.0, 3); Nx.close_form!(p.forms); Nx.push_block!(p, Nx.ZERO, 1)
        @test_throws ArgumentError Nx.validate(p)                       # id 3 with nvars == 2
        bad = Nx.MomentProgram{Float64}(1)
        Nx.push_term!(bad.forms, 1.0, 1); Nx.close_form!(bad.forms); Nx.push_block!(bad, Nx.HPSD, 1)
        @test_throws ArgumentError Nx.validate(bad)                     # HPSD with real coefficients
    end

    @testset "contract 2: small programs with a known optimum" begin
        # minimise y2 with [1 y2; y2 1] PSD: optimum -1
        p = Nx.MomentProgram{Float64}(2)
        for id in (1, 2, 1)
            Nx.push_term!(p.forms, 1.0, id); Nx.close_form!(p.forms)
        end
        Nx.push_block!(p, Nx.PSD, 2)
        push!(p.objective_ids, 2); push!(p.objective_coefs, 1.0)
        @test solve_reference(p)[1] ≈ -1 atol = 1e-5
        # add y2 + 0.5 >= 0 to the same program: optimum -0.5
        Nx.push_term!(p.forms, 1.0, 2); Nx.push_term!(p.forms, 0.5, 1); Nx.close_form!(p.forms)
        Nx.push_block!(p, Nx.NONNEGATIVE, 1)
        @test solve_reference(p)[1] ≈ -0.5 atol = 1e-5
        # minimise y2 + y3 with norm((y2, y3)) <= 1 and y3 == 0.6: optimum -0.8 + 0.6
        q = Nx.MomentProgram{Float64}(3)
        for id in (1, 2, 3)
            Nx.push_term!(q.forms, 1.0, id); Nx.close_form!(q.forms)
        end
        Nx.push_block!(q, Nx.SOC, 3)
        Nx.push_term!(q.forms, 1.0, 3); Nx.push_term!(q.forms, -0.6, 1); Nx.close_form!(q.forms)
        Nx.push_block!(q, Nx.ZERO, 1)
        append!(q.objective_ids, (2, 3)); append!(q.objective_coefs, (1.0, 1.0))
        value, y = solve_reference(q)
        @test value ≈ -0.2 atol = 1e-5
        @test y ≈ [1, -0.8, 0.6] atol = 1e-4
        # Hermitian 2x2 block [1 y2+im*y3; conj 1]: minimise y2 + y3, optimum -sqrt(2)
        h = Nx.MomentProgram{ComplexF64}(3)
        Nx.push_term!(h.forms, 1.0, 1); Nx.close_form!(h.forms)
        Nx.push_term!(h.forms, 1.0, 2); Nx.push_term!(h.forms, -1.0im, 3); Nx.close_form!(h.forms)   # entry (2,1) = conj of (1,2)
        Nx.push_term!(h.forms, 1.0, 1); Nx.close_form!(h.forms)
        Nx.push_block!(h, Nx.HPSD, 2)
        append!(h.objective_ids, (2, 3)); append!(h.objective_coefs, (1.0, 1.0))
        @test solve_reference(h)[1] ≈ -sqrt(2) atol = 1e-5
    end

    @testset "dense reference equals the current package" begin
        for (n, order) in ((4, 1), (3, 2))
            reg, (x, y, z) = create_pauli_variables(1:n)
            pop = polyopt(sum(x[i] * x[i + 1] + y[i] * y[i + 1] + z[i] * z[i + 1] for i in 1:(n - 1)), reg)
            config = SolverConfig(optimizer=SOLVER, order=order)
            sparsity = compute_sparsity(pop, config)
            old = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
            objective = [(1.0, Nx.mul(toy_letter(i, a), toy_letter(i + 1, a))[2]) for i in 1:(n - 1) for a in 0:2]
            program, table = dense_program(toy_words(n, order), objective)
            @test Nx.nvars(program) == old.n_unique_moment_matrix_elements
            @test Set(Nx.legacy_letters(l.word) for l in table.labels) == Set(Int.(m.word) for m in values(old.linear.key_to_monomial))
            @test all(l.part == 0x01 for l in table.labels[2:end])       # each Pauli moment is real
            @test solve_reference(program)[1] ≈ cs_nctssos(pop, config; dualize=false).objective atol = 1e-5
        end
    end
end

include("core/runtests.jl")
include("conic/runtests.jl")
include("pipeline/runtests.jl")
include("model/runtests.jl")

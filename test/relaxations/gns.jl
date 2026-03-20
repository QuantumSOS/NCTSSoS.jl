using Test, NCTSSoS, LinearAlgebra

using NCTSSoS: get_ncbasis, monomials

function _moment_lookup(monomap, mono)
    key = symmetric_canon(NCTSSoS.expval(mono))
    return get(monomap, key, 0.0)
end

function _evaluate_monomial(matrices::Dict, mono::NormalMonomial)
    sample_matrix = first(values(matrices))
    T = eltype(sample_matrix)
    dim = size(sample_matrix, 1)
    result = Matrix{T}(I, dim, dim)
    for idx in mono.word
        result *= matrices[idx]
    end
    return result
end

function _expectation(matrices::Dict, xi::AbstractVector, mono::NormalMonomial)
    op = _evaluate_monomial(matrices, mono)
    return dot(xi, op * xi)
end

function _moment_key(mono::NormalMonomial)
    return symmetric_canon(NCTSSoS.expval(mono))
end

@testset "GNS Construction" begin
    @testset "Hankel dictionary and localizing matrix" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])
        basis_full = get_ncbasis(reg, 2)
        basis = get_ncbasis(reg, 1)

        # Moments m_k = L(x^k)
        moments = Dict(
            0 => 1.0,
            1 => 0.2,
            2 => 0.7,
            3 => 0.5,
            4 => 0.4,
        )

        H = [moments[i + j - 2] for i in 1:3, j in 1:3]
        dict = NCTSSoS.hankel_entries_dict(H, basis_full)

        @test dict[only(monomials(x[1]))] == moments[1]
        @test dict[only(monomials(x[1]^2))] == moments[2]

        x_idx = first(indices(reg))
        K = NCTSSoS.construct_localizing_matrix(dict, x_idx, basis)

        @test size(K) == (2, 2)
        @test K ≈ [moments[1] moments[2]; moments[2] moments[3]] atol = 1e-12

        gns = @test_logs (:warn, r"Flatness condition violated") gns_reconstruct(
            H,
            reg,
            2;
            hankel_deg=1,
            atol=1e-12,
        )
        Xmat = gns.matrices[x_idx]
        @test size(Xmat) == (2, 2)
        @test real(dot(gns.xi, Xmat * gns.xi)) ≈ moments[1] atol = 1e-10
        @test real(dot(gns.xi, Xmat^2 * gns.xi)) ≈ moments[2] atol = 1e-10
    end

    @testset "Real moment solve returns numeric moments" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])
        pop = polyopt(1.0 * x[1]^2 + 1.0, reg)
        cfg = SolverConfig(optimizer=SOLVER, order=1)

        sparsity = compute_sparsity(pop, cfg)
        mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        result = NCTSSoS.solve_moment_problem(mp, SOLVER)

        @test !isempty(result.monomap)
        @test all(v -> v isa Real, values(result.monomap))
        @test _moment_lookup(result.monomap, one(x[1])) ≈ 1.0 atol = 1e-8
    end

    @testset "Fermionic reconstruction preserves adjoint pairs" begin
        reg, (a, a_dag) = create_fermionic_variables(1:1)

        monomap = Dict(
            _moment_key(one(a[1])) => 1.0,
            _moment_key(only(monomials(a[1]))) => 0.0,
            _moment_key(only(monomials(a_dag[1]))) => 0.0,
            _moment_key(only(monomials(a_dag[1] * a[1]))) => 0.0,
        )

        gns = gns_reconstruct(monomap, reg, 2; hankel_deg=1, atol=1e-10)
        a_idx = only(only(monomials(a[1])).word)
        adag_idx = only(only(monomials(a_dag[1])).word)
        A = gns.matrices[a_idx]
        A_dag = gns.matrices[adag_idx]

        @test A_dag ≈ A' atol = 1e-12
        @test !isapprox(A_dag, A; atol=1e-12)
        @test A^2 ≈ zeros(eltype(A), size(A)) atol = 1e-12
        @test A_dag^2 ≈ zeros(eltype(A_dag), size(A_dag)) atol = 1e-12
    end

    @testset "Monomap reconstruction rejects missing moments" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])

        monomap = Dict(
            _moment_key(one(x[1])) => 1.0,
            _moment_key(only(monomials(x[1]))) => 0.2,
        )

        err = try
            gns_reconstruct(monomap, reg, 2; hankel_deg=1)
            nothing
        catch caught
            caught
        end

        @test err isa ArgumentError
        @test occursin("missing", lowercase(sprint(showerror, err)))
    end

    @testset "Sparse GNS reconstructs disjoint cliques" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        pop = polyopt(1.0 * x[1]^2 + 1.0 * x[2]^2, reg)
        cfg = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MaximalElimination(),
            ts_algo=NoElimination(),
        )

        sparsity = compute_sparsity(pop, cfg)
        @test length(sparsity.corr_sparsity.cliques) == 2
        @test all(length(clique) == 1 for clique in sparsity.corr_sparsity.cliques)

        monomap = Dict(
            _moment_key(one(x[1])) => 1.0,
            _moment_key(only(monomials(x[1]))) => 0.0,
            _moment_key(only(monomials(x[1]^2))) => 1.0,
            _moment_key(only(monomials(x[1]^3))) => 0.0,
            _moment_key(only(monomials(x[1]^4))) => 1.0,
            _moment_key(only(monomials(x[2]))) => 0.0,
            _moment_key(only(monomials(x[2]^2))) => 1.0,
            _moment_key(only(monomials(x[2]^3))) => 0.0,
            _moment_key(only(monomials(x[2]^4))) => 1.0,
        )

        gns = @test_logs (:warn, r"Flatness condition violated") gns_reconstruct(
            monomap,
            sparsity,
            2;
            hankel_deg=1,
            atol=1e-10,
        )
        x1_idx = reg[:x₁]
        x2_idx = reg[:x₂]
        X1 = gns.matrices[x1_idx]
        X2 = gns.matrices[x2_idx]

        @test size(X1) == size(X2)
        @test norm(X1 - X1') ≤ 1e-10
        @test norm(X2 - X2') ≤ 1e-10
        @test norm(gns.xi) ≈ 1.0 atol = 1e-10

        expected_moments = [
            (one(x[1]), 1.0),
            (only(monomials(x[1])), 0.0),
            (only(monomials(x[2])), 0.0),
            (only(monomials(x[1]^2)), 1.0),
            (only(monomials(x[2]^2)), 1.0),
            (only(monomials(x[1] * x[2])), 0.0),
            (only(monomials(x[2] * x[1])), 0.0),
        ]

        for (mono, expected) in expected_moments
            observed = real(_expectation(gns.matrices, gns.xi, mono))
            @test observed ≈ expected atol = 1e-8
        end
    end

    @testset "Sparse GNS rejects overlapping cliques" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        pop = polyopt(1.0 * x[1] * x[2] + 1.0 * x[2] * x[3], reg)
        cfg = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=AsIsElimination(),
            ts_algo=NoElimination(),
        )

        sparsity = compute_sparsity(pop, cfg)
        @test length(sparsity.corr_sparsity.cliques) > 1
        @test !all(isempty, [intersect(Set(Int.(a)), Set(Int.(b))) for (i, a) in enumerate(sparsity.corr_sparsity.cliques) for b in sparsity.corr_sparsity.cliques[i+1:end]])

        err = try
            gns_reconstruct(Dict{Any,Float64}(), sparsity, 2; hankel_deg=1)
            nothing
        catch caught
            caught
        end

        @test err isa ArgumentError
        @test occursin("pairwise-disjoint", sprint(showerror, err))
    end

    @testset "Example 5.3 reproduction on the nc unit ball" begin
        reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
        X, Y = vars

        f = 2.0 - X^2 + X * Y^2 * X - Y^2
        g = 1.0 - X^2 - Y^2
        pop = polyopt(
            f,
            reg;
            ineq_constraints=[g],
        )

        cfg = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )

        sparsity = compute_sparsity(pop, cfg)
        mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        result = NCTSSoS.solve_moment_problem(mp, SOLVER)

        @test result.objective ≈ 1.0 atol = 1e-5
        @test all(v -> v isa Real, values(result.monomap))

        gns = @test_logs (:warn, r"Flatness condition violated") gns_reconstruct(
            result.monomap,
            reg,
            3;
            hankel_deg=2,
            atol=1e-6,
        )
        A1 = gns.matrices[reg[:X₁]]
        A2 = gns.matrices[reg[:X₂]]

        @test gns.rank == 5
        @test gns.full_rank >= gns.rank
        @test size(A1) == (5, 5)
        @test size(A2) == (5, 5)
        @test norm(A1 - A1') ≤ 1e-8
        @test norm(A2 - A2') ≤ 1e-8
        @test norm(gns.xi) ≈ 1.0 atol = 1e-8

        F = Matrix(2I - A1^2 + A1 * A2^2 * A1 - A2^2)
        F = Symmetric((F + F') / 2)
        λmin = minimum(eigvals(F))

        @test λmin ≈ 1.0 atol = 1e-5
        @test real(dot(gns.xi, F * gns.xi)) ≈ 1.0 atol = 1e-5

        for mono in get_ncbasis(reg, 3)
            observed = _expectation(gns.matrices, gns.xi, mono)
            expected = _moment_lookup(result.monomap, mono)
            @test observed ≈ expected atol = 5e-5
        end
    end
end

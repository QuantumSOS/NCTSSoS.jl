using Test, NCTSSoS, JuMP, LinearAlgebra
import Graphs, TOML

@testset "Local Pauli invariant SOCs" begin
    # Matrix oracle is independent of the separator's binary encoding.
    matrices = Dict('I' => ComplexF64[1 0; 0 1], 'X' => ComplexF64[0 1; 1 0],
        'Y' => ComplexF64[0 -im; im 0], 'Z' => ComplexF64[1 0; 0 -1])
    pauli(label) = reduce(kron, [matrices[c] for c in label])

    reg, (x, y, z) = create_pauli_variables(1:3)
    h = -(x[1]*y[2] + x[1]*y[3] + x[2]*x[3] + y[2]*y[3] + z[2]*z[3])
    pop = polyopt(h, reg)
    config = SolverConfig(optimizer=SOLVER, order=1)
    bare = pauli_invariant_soc_nctssos(pop, config; max_rounds=0)
    strengthened = pauli_invariant_soc_nctssos(pop, config; max_rounds=1, batch_size=1)
    @test bare.objective ≈ -5 atol=2e-5
    @test strengthened.objective ≈ -sqrt(15) atol=2e-5
    # COSMO 0.8.11's reported dual objective need not be tight for these cones.
    # The benchmark separately certifies the bound using exact dual repair.
    @test strengthened.bound <= strengthened.objective + 2e-5
    @test keys(bare.monomap) == keys(strengthened.monomap)
    @test num_variables(bare.model) == num_variables(strengthened.model)
    @test length(strengthened.records) == 1
    record = only(strengthened.records)
    @test isodd(count(==('Y'), record.r))
    @test record.constant_squared == 3
    @test all(label -> label in record.labels, ["XYI", "XIY", "IXX", "IYY", "IZZ"])
    H = -sum(pauli(s) for s in ["XYI", "XIY", "IXX", "IYY", "IZZ"])
    @test eigmin(Hermitian(H)) ≈ -3

    # The advertised arbitrary subset, and validation before mutating the model.
    model, _ = build_jump_model(bare.moment_problem)
    subset = add_pauli_invariant_soc!(model, bare.moment_problem, [1,2,3], "YII";
        labels=["XYI", "XIY", "IXX", "IYY", "IZZ"])
    @test length(subset.labels) == 5
    @test subset.r == "YII"
    @test_throws ArgumentError add_pauli_invariant_soc!(model, bare.moment_problem, [1,2,3], "III")
    for labels in (["III"], ["XII"], ["XXX"], ["IXX", "IXX"])
        @test_throws ArgumentError add_pauli_invariant_soc!(model, bare.moment_problem, [1,2,3], "YII"; labels)
    end
    @test_throws ArgumentError pauli_invariant_soc_nctssos(pop, config; clusters=[[3,2,1]])
    @test_throws ArgumentError pauli_invariant_soc_nctssos(pop, config; batch_size=0)

    # Check the selected maximal score by direct matrix conjugation, including
    # asymmetry and unavailable higher-weight words. This catches a missing
    # x/z swap or Y-transpose sign in the Walsh transform through the public API.
    for k in (3, 4)
        reg, (x, y, z) = create_pauli_variables(1:k)
        h = -(x[1]*y[2] + x[1]*y[3] + x[2]*x[3] + y[2]*y[3] + z[2]*z[3]) +
            0.23z[1] + 0.17x[2] - 0.31z[1]*y[3]
        k == 4 && (h -= 0.7*(x[1]*y[4] + x[2]*x[4] + y[2]*y[4] + z[2]*z[4]))
        pop = polyopt(h, reg)
        config = SolverConfig(optimizer=SOLVER, order=1)
        bare = pauli_invariant_soc_nctssos(pop, config; max_rounds=0)
        strengthened = pauli_invariant_soc_nctssos(pop, config;
            clusters=[collect(1:k)], max_rounds=1, batch_size=1)
        @test length(strengthened.records) == 1
        rec = only(strengthened.records)
        @test rec.cluster == collect(1:k)
        @test rec.constant_squared == 2^(k-1)-1
        @test keys(bare.monomap) == keys(strengthened.monomap)
        @test strengthened.objective > bare.objective + 1e-3

        labeled = Dict{String,Float64}()
        for (key, value) in bare.monomap
            letters = fill('I', k)
            for letter in bare.moment_problem.linear.key_to_monomial[key].word
                letters[(letter-1) ÷ 3 + 1] = "XYZ"[(letter-1) % 3 + 1]
            end
            all(==('I'), letters) && continue
            labeled[String(letters)] = real(value)
        end
        scores = Float64[]
        for letters in Iterators.product(ntuple(_ -> "IXYZ", k)...)
            r = join(letters)
            isodd(count(==('Y'), r)) || continue
            R = pauli(r)
            score = sum(value^2 for (label, value) in labeled if
                R * conj(pauli(label)) * R ≈ pauli(label))
            push!(scores, score)
        end
        @test rec.violation ≈ maximum(scores) - rec.constant_squared atol=2e-4
        R = pauli(rec.r)
        @test R * conj(R) ≈ -Matrix{ComplexF64}(I, 2^k, 2^k)
        for label in rec.labels
            @test R * conj(pauli(label)) * R ≈ pauli(label)
        end
        # Independent physical state: selected expectations obey every cut.
        ψ = normalize([sin(j) + im*cos(2j) for j in 1:2^k])
        @test sum(abs2(dot(ψ, pauli(label)*ψ)) for label in rec.labels) <= rec.constant_squared + 1e-12
    end

    @testset "Composition with exhaustive anticommutation constraints" begin
        fixture = only(TOML.parsefile(joinpath(@__DIR__, "..", "data", "expectations", "pauli_invariant_soc.toml"))["cases"])
        @test fixture["id"] == "dev382_rich3_combined"
        reg, _ = create_pauli_variables(1:3)
        pairs = filter(m -> NCTSSoS.degree(m) == 2, get_ncbasis(reg, 2))
        h = sum((c/5)*m for (c,m) in zip(fixture["weights_times_five"], pairs))
        bare = pauli_invariant_soc_nctssos(polyopt(h, reg), SolverConfig(optimizer=SOLVER, order=1); max_rounds=0)
        model, mp = bare.model, bare.moment_problem
        monos = [mp.linear.key_to_monomial[key] for key in mp.linear.moments]
        # Use explicit matrices to define the comparison graph independently.
        operators = map(monos) do mono
            letters = fill('I', 3)
            for letter in mono.word
                letters[(letter-1) ÷ 3 + 1] = "XYZ"[(letter-1) % 3 + 1]
            end
            pauli(String(letters))
        end
        graph = Graphs.SimpleGraph(length(monos))
        for i in eachindex(monos), j in 1:(i-1)
            isapprox(operators[i]*operators[j], -operators[j]*operators[i]; atol=1e-12) &&
                Graphs.add_edge!(graph, i, j)
        end
        cliques = filter(c -> length(c) > 1, Graphs.maximal_cliques(graph))
        @test length(cliques) == 363
        for clique in cliques
            @constraint(model, [1.0; model[:y_re][clique]] in SecondOrderCone())
        end
        before_variables = num_variables(model)
        result = separate_pauli_invariant_socs!(model, mp, [[1,2,3]]; max_rounds=3, batch_size=2)
        @test result.objectives[1] ≈ fixture["expected"]["anticommutation"] atol=2e-4
        @test result.objective ≈ fixture["expected"]["objective"] atol=2e-4
        @test result.objective > result.objectives[1] + 0.15
        @test num_variables(model) == before_variables
        @test length(mp.linear.moments) == fixture["expected"]["nuniq"]
        @test length(result.records) == 1
        @test only(result.records).r == "YYY"
        @test length(only(result.records).labels) == 27
        @test length(result.objectives) == 2 # Stops when no new cut is violated.
    end
end

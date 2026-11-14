# Action fixtures live with the pipeline tests, not in the shared contracts.
struct TestPauliMap <: Nx.AbstractAction
    images::Vector{Tuple{Nx.Phase,ToyPauli}}
    anti::Bool
end
Nx.is_antiunitary(g::TestPauliMap) = g.anti
function Nx.act(g::TestPauliMap, w::ToyPauli)
    p, u = one(Nx.Phase), one(ToyPauli)
    for site in Nx.sites(w)
        bit = UInt64(1) << (site - 1)
        x, z = w.x & bit != 0, w.z & bit != 0
        l = 3 * (site - 1) + (x && z ? 2 : x ? 1 : 3)
        q, v = g.images[l]
        r, u = Nx.mul(u, v)
        p = p * q * r
    end
    return p, u
end
Base.:(==)(g::TestPauliMap, h::TestPauliMap) = g.anti == h.anti && g.images == h.images
Base.hash(g::TestPauliMap, h::UInt) = hash(g.anti, hash(g.images, h))
function Base.one(g::TestPauliMap)
    n = length(g.images) ÷ 3
    TestPauliMap([(one(Nx.Phase), toy_letter(i, a)) for i in 1:n for a in 0:2], false)
end
Base.isone(g::TestPauliMap) = g == one(g)
function Nx.compose(g::TestPauliMap, h::TestPauliMap)
    TestPauliMap(map(h.images) do (p, w)
        q, v = Nx.act(g, w)
        ((g.anti ? conj(p) : p) * q, v)
    end, xor(g.anti, h.anti))
end
function Base.inv(g::TestPauliMap)
    h = one(g)
    for _ in 1:10000
        Nx.compose(g, h) == one(g) && return h
        h = Nx.compose(g, h)
    end
    error("test action order exceeds 10000")
end
function test_map(n; sites=collect(1:n), axes=(0, 1, 2), signs=(1, 1, 1), anti=false)
    TestPauliMap([(Nx.Phase(signs[a + 1] == 1 ? 0 : 2), toy_letter(sites[i], axes[a + 1])) for i in 1:n for a in 0:2], anti)
end
test_shift(n) = test_map(n; sites=[mod1(i + 1, n) for i in 1:n])
test_reflect(n) = test_map(n; sites=[mod1(2 - i, n) for i in 1:n])
test_flip(n) = test_map(n; signs=(-1, -1, 1))
test_conjugate(n) = test_map(n; signs=(1, -1, 1), anti=true)
function ring_objective(n, model)
    if model == :heisenberg
        return [(1.0, Nx.mul(toy_letter(i, a), toy_letter(mod1(i + 1, n), a))[2]) for i in 1:n for a in 0:2]
    end
    return vcat([(1.0, Nx.mul(toy_letter(i, 2), toy_letter(mod1(i + 1, n), 2))[2]) for i in 1:n],
        [(0.7, toy_letter(i, 0)) for i in 1:n])
end
function evaluated_blocks(p, y)
    out = Matrix{ComplexF64}[]
    for b in p.blocks
        A = zeros(ComplexF64, b.side, b.side)
        for j in 1:b.side, i in j:b.side
            f = b.first_form + Nx.tri_offset(b.side, i, j)
            A[i,j] = sum((p.forms.coefs[t] * y[p.forms.ids[t]] for t in Nx.term_range(p.forms, f)); init=0im)
            A[j,i] = conj(A[i,j])
        end
        push!(out, A)
    end
    out
end

@testset "Direct cyclic orbit programs" begin
    # Unequal orbit sizes and signed stabilizers exercise both normalization
    # and sector selection; compare all coefficients through arbitrary moments.
    for n in (3, 4), signed in (false, true)
        rows = toy_words(n, 2)
        g = signed ? Nx.compose(test_shift(n), test_flip(n)) : test_shift(n)
        order = signed && isodd(n) ? 2n : n
        reduced, table, sectors = Nx.cyclic_program(rows, [], g, order)
        dense, dt = dense_program(rows, [])
        rng = MersenneTwister(71)
        y = randn(rng, reduced.nvars); y[1] = 1
        yd = map(dt.labels) do l
            id, phase = Nx.orbit_moment!(table, l.word)
            id == 0 ? 0.0 : Nx.scale(phase, y[id])
        end
        D = only(evaluated_blocks(dense, yd))
        orbits = Nx.CyclicRows(rows, g, order)
        index = Dict(w => i for (i, w) in enumerate(rows))
        for ((k, _), A) in zip(sectors, evaluated_blocks(reduced, y))
            vectors = Vector{ComplexF64}[]
            for (ws, ps, stabilizer) in zip(orbits.words, orbits.phases, orbits.stabilizers)
                mod(4k * length(ws) - Int(stabilizer.k) * order, 4order) == 0 || continue
                v = zeros(ComplexF64, length(rows))
                for (t, (w, phase)) in enumerate(zip(ws, ps))
                    v[index[w]] = Nx.scale(phase, cispi(-2k * (t - 1) / order)) / sqrt(length(ws))
                end
                push!(vectors, v)
            end
            V = hcat(vectors...)
            @test V'V ≈ I atol=1e-12
            @test V'D*V ≈ A atol=1e-11
        end
        @test sum(b.side for b in reduced.blocks) == length(rows)
    end
    @test_throws ArgumentError Nx.cyclic_program(toy_words(3, 1), [(1., toy_letter(1, 0))], test_shift(3), 3)
    @test_throws ArgumentError Nx.CyclicRows([one(ToyPauli), toy_letter(1, 0)], test_shift(3), 3)
    @test_throws ArgumentError Nx.CyclicRows(toy_words(3, 1), test_shift(3), 2)
    @testset "Bounds match unreduced reference" begin
        for model in (:heisenberg, :ising), n in (3, 4)
            rows = toy_words(n, 1)
            objective = ring_objective(n, model)
            dense, _ = dense_program(rows, objective)
            reference = solve_reference(dense)[1]
            for extra in ((), (test_reflect(n),), (test_conjugate(n),))
                p, _, _ = Nx.cyclic_program(rows, objective, test_shift(n), n; moment_generators=extra)
                @test solve_reference(p)[1] ≈ reference atol=1e-6
            end
            if model == :heisenberg
                p, _, _ = Nx.cyclic_program(rows, objective, test_shift(n), n; labels=(test_flip(n),))
                @test solve_reference(p)[1] ≈ reference atol=1e-6
            end
        end
    end
end

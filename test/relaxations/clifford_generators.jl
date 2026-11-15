using Test, NCTSSoS, LinearAlgebra

@testset "Complete Clifford detection and generator-only reduction" begin
    N = NCTSSoS
    reg, (x,y,z) = create_pauli_variables(1:2)
    rows = get_ncbasis(reg, 2)
    lookup = Dict(w => i for (i,w) in enumerate(rows))
    signature(g) = Tuple(begin
        s, w = N._act_monomial(g, row)
        s * lookup[w]
    end for row in rows)
    compose(a,b) = Tuple(sign(i) * a[abs(i)] for i in b)
    function closure(generators)
        queue = [Tuple(1:length(first(generators)))]
        seen = Set(queue)
        pos = 1
        while pos <= length(queue)
            current = queue[pos]
            pos += 1
            for g in generators
                image = compose(g, current)
                image in seen && continue
                push!(seen, image)
                push!(queue, image)
            end
        end
        return seen
    end

    # Independent oracle: actual 4 x 4 unitary conjugation, not tableau lifts.
    X = ComplexF64[0 1;1 0]
    Y = ComplexF64[0 -im;im 0]
    Z = ComplexF64[1 0;0 -1]
    E = Matrix{ComplexF64}(I,2,2)
    paulis = [E,X,Y,Z]
    word_matrices = [kron(paulis[1 + sum(N._pauli_type(k)+1 for k in w.word if N._pauli_site(k)==1;init=0)],
                          paulis[1 + sum(N._pauli_type(k)+1 for k in w.word if N._pauli_site(k)==2;init=0)]) for w in rows]
    h = (X+Z)/sqrt(2)
    s = ComplexF64[1 0;0 im]
    cnot = ComplexF64[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]
    function matrix_signature(U)
        Tuple(begin
            overlaps = [real(LinearAlgebra.tr(P' * U * Q * U'))/4 for P in word_matrices]
            index = argmax(abs.(overlaps))
            round(Int, overlaps[index]) * index
        end for Q in word_matrices)
    end
    full = closure(matrix_signature.([kron(h,E),kron(E,h),kron(s,E),kron(E,s),cnot]))
    @test length(full) == 11520

    heisenberg = x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
    for H in (heisenberg, Polynomial(z[2]), x[1]-z[1]+2z[2],
              x[1]+z[2]+x[1]*z[2], x[1]+z[1]+y[1]+2z[2], 3+heisenberg)
        coefficients = zeros(ComplexF64, length(rows))
        for (c,w) in terms(H)
            coefficients[lookup[w]] = c
        end
        expected = Set(g for g in full if all(
            coefficients[abs(g[i])] == sign(g[i])*coefficients[i] for i in eachindex(g)))
        detected = sympleq_symmetry_spec(H)
        actual = closure(signature.(detected.clifford_generators))
        @test actual == expected
        @test all(g -> N._act_polynomial(g,H) == H, detected.clifford_generators)
        if H == heisenberg
            @test length(actual) == 192
            @test signature(CliffordSymmetry(:SWAP,1,2)) in actual
        end
    end
    @test_throws ArgumentError sympleq_symmetry_spec(heisenberg;max_search_nodes=1)

    # The ring translation changes a chosen nullspace basis but preserves its
    # span. Check membership in the generated term action, not generator names.
    _, (rx,ry,rz) = create_pauli_variables(1:4)
    ring = sum(v[i]*v[mod1(i+1,4)] for v in (rx,ry,rz) for i in 1:4)
    open_chain = sum(v[i]*v[i+1] for v in (rx,ry,rz) for i in 1:3)
    translation = pauli_site_permutation([2,3,4,1])
    term_rows = monomials(ring)
    term_lookup = Dict(w=>i for (i,w) in enumerate(term_rows))
    term_signature(g) = Tuple(term_lookup[last(N._act_monomial(g,w))] for w in term_rows)
    ring_spec = sympleq_symmetry_spec(ring)
    @test term_signature(translation) in closure(term_signature.(ring_spec.clifford_generators))
    @test N._act_polynomial(translation,ring) == ring
    @test N._act_polynomial(translation,open_chain) != open_chain

    # A real cyclic action with complex irreps, an entangling action, and the
    # full detected subgroup. Compare actual SDP bounds against unreduced ones.
    global_s = CliffordSymmetry(Dict(x[1]=>(1,y[1]), y[1]=>(-1,x[1]),
                                    x[2]=>(1,y[2]), y[2]=>(-1,x[2])); nqubits=2)
    cases = [
        (x[1]*x[2]+y[1]*y[2]+2z[1]+3z[2], [global_s]),
        (z[2]+z[1]*z[2]+0.3x[1]+0.3x[1]*x[2]+0.2(y[1]*y[2]-x[1]*z[2]),
         [CliffordSymmetry(:CNOT,1,2)]),
        (heisenberg, sympleq_symmetry_spec(heisenberg).clifford_generators),
    ]
    for (H, gs) in cases
        pop = polyopt(H,reg)
        baseline = cs_nctssos(pop,SolverConfig(optimizer=SOLVER,moment_basis=rows,
            cs_algo=NoElimination(),ts_algo=NoElimination()))
        spec = SymmetrySpec(gs;clifford_backend=:generators,offblock_check=:full)
        result = cs_nctssos(pop,SolverConfig(optimizer=SOLVER,moment_basis=rows,
            cs_algo=NoElimination(),ts_algo=NoElimination(),symmetry=spec))
        @test result.objective ≈ baseline.objective atol=2e-6
        @test result.symmetry.group_order === nothing
        @test all(==(:generators), result.symmetry.block_provenance)
        @test maximum(result.symmetry.psd_block_sizes) < length(rows)
        H == heisenberg && @test result.objective ≈ -3 atol=2e-6
    end

    # Repeated irreducible copies must give multiplicity-size blocks [4,3],
    # not irreducible-size blocks or one block per eigenvector of a sample.
    r3,(a,b,c) = create_pauli_variables(1:3)
    H3 = sum(v[i]*v[j] for v in (a,b,c) for i in 1:3 for j in i+1:3)
    gs = [pauli_site_permutation([2,1,3]),pauli_site_permutation([2,3,1])]
    spec = SymmetrySpec(gs;clifford_backend=:generators,offblock_check=:full)
    result = cs_nctssos(polyopt(H3,r3),SolverConfig(optimizer=SOLVER,order=1,
        cs_algo=NoElimination(),ts_algo=NoElimination(),symmetry=spec))
    @test sort(result.symmetry.psd_block_sizes) == [3,4]
    # Each axis has a 3 x 3 site Gram matrix with diagonal 1. Its PSD
    # constraint gives sum_{i<j} c_ij >= -3/2, attained at c_ij=-1/2.
    @test result.objective ≈ -4.5 atol=2e-6

    cnot_spec = SymmetrySpec(CliffordSymmetry(:CNOT,1,2);clifford_backend=:generators)
    @test_throws ArgumentError cs_nctssos(polyopt(cases[2][1],reg),
        SolverConfig(optimizer=SOLVER,moment_basis=[one(x[1]);x;y;z],
            cs_algo=NoElimination(),ts_algo=NoElimination(),symmetry=cnot_spec))
end

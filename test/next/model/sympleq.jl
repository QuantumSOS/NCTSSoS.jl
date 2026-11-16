# Retained detector laws from test/relaxations/symmetry.jl, now on packed words.
module SympleQChecks
using Test
import GroupsCore
import NCTSSoS: Next as N

word(p) = last(only(N.terms(p)))

function mapped(g, p)
    W = eltype(p.data.words)
    data = N.PolynomialBuffer{W,ComplexF64}(length(p.data))
    for (c,w) in N.terms(p)
        phase, image = N.act(g,w)
        N.push_term!(data,N.scale(phase,c),image)
    end
    N._polynomial(p.registry,data)
end

# Include all letters and their images: the signed action distinguishes maps
# that agree on Hamiltonian terms. This avoids enumerating Clifford objects.
function signed_word_group(generators, letters)
    PG = N.PermutationGroups
    words = unique([letters; [N.act(g,w)[2] for g in generators for w in letters]])
    lookup = Dict(w => i for (i,w) in enumerate(words))
    d = length(words)
    function permutation(g)
        images = zeros(Int,2d)
        for (i,w) in enumerate(words)
            phase, image = N.act(g,w)
            sign = N.scale(phase,1)
            @assert sign in (-1,1)
            j = lookup[image] # Fails if the proposed domain is not closed.
            images[i] = sign == 1 ? j : j+d
            images[i+d] = sign == 1 ? j+d : j
        end
        PG.Perm{Int}(images)
    end
    PG.PermGroup(permutation.(generators)), permutation
end

function check()
    @testset "SympleQ preserves the relation space, not a chosen basis" begin
        for (model,sizes) in ((:heisenberg,(4,6,8,16,32)),(:ising,(16,32,64)))
            for n in sizes
                _, ops = N.pauli_variables(n)
                x,y,z = ops
                H = model == :heisenberg ?
                    sum(o[i]*o[mod1(i+1,n)] for o in ops for i in 1:n) :
                    sum(z[i]*z[mod1(i+1,n)]+0.7x[i] for i in 1:n)
                # Test the retained SympleQ reduction, before public detection
                # appends lattice translations and classifies moment-only maps.
                generators = N._sympleq_clifford_candidates(H)
                @test all(g -> mapped(g,H) == H,generators)
                @test length(generators) <= 10
                group, permutation = signed_word_group(generators,word.([x;y;z]))
                @test permutation(N.SitePermutation([2:n;1])) in group
                @test permutation(N.SitePermutation(collect(n:-1:1))) in group
                @test GroupsCore.order(Int,group) == (model == :heisenberg ? 96n : 4n)
                if model == :heisenberg
                    W = typeof(word(x[1]))
                    rotations = [N.CliffordMap([N.pauli_letter(W,i,a) for a in axes for i in 1:n])
                        for axes in ((2,0),(1,2))]
                    @test all(g -> permutation(g) in group,rotations)
                end
            end
        end

        # Four commuting terms have one relation; the fifth is independent.
        # Of 5! graph permutations, exactly 4! preserve it (identity omitted).
        _, (x,_,_) = N.pauli_variables(4)
        H = 1.0*(x[1]+x[2]+x[3]+x[1]*x[2]*x[3]+x[4])
        tab = N.SymplecticTableau(H)
        @test only(length.(N.pauli_cycle_basis(tab))) == 4
        generators = N.sympleq_generators(tab;max_automorphisms=25)
        @test length(generators) == 23
        independent = findfirst(==(word(x[4])),tab.monomials)
        @test all(g -> g.permutation[independent] == independent,generators)
        @test all(g -> mapped(g,H) == H,N._sympleq_clifford_candidates(H))
    end

    @testset "SympleQ deterministic internals" begin
        _, (x,y,z) = N.pauli_variables(3)
        hamiltonian = 2.0x[1]+2.0z[1]+3.0x[2]*z[3]
        tab = N.SymplecticTableau(hamiltonian)
        W = typeof(word(x[1]))
        @test length(tab) == 3
        @test size(tab) == (3,6)
        @test sprint(show,tab) == "SymplecticTableau(3 terms, 3 qubits)"
        @test tab.coeffs == [2.0+0.0im,2.0+0.0im,3.0+0.0im]
        @test tab.nqubits == 3

        # Word-to-row conversion now belongs to the packed tableau constructor.
        ytab = N.SymplecticTableau(y[2]+z[3])
        i = findfirst(==(word(y[2])),ytab.monomials)
        y_row, y_eta = ytab.paulis[i,:], ytab.eta[i]
        @test y_eta == 1
        @test y_row == UInt8[0,1,0,0,1,0]
        y_word, eta_back = N._symplectic_row_to_pauli_word(y_row,W)
        @test eta_back == 1
        @test y_word == word(y[2])

        @test N.symplectic_product_matrix(tab) == UInt8[0 1 0;1 0 0;0 0 0]
        graph = N.anticommutation_graph(tab)
        @test sprint(show,graph) == "SympleQGraph(3 vertices, 1 edges, 3 Pauli vertices, 0 cycle vertices)"
        @test graph.colors == [(:pauli,2.0),(:pauli,2.0),(:pauli,3.0)]
        cycle_tab = N.SymplecticTableau(1.0*(x[1]+x[2]+x[1]*x[2]))
        @test N.pauli_cycle_basis(cycle_tab) == [[1,2,3]]
        augmented = N.cycle_augmented_graph(cycle_tab)
        @test isempty(augmented.cycles)
        @test isempty(augmented.auxiliary_vertices)
        @test_throws ArgumentError N.pauli_cycle_basis(cycle_tab;cycle_strategy=:minimal_circuits)

        @test_throws ArgumentError N.TermPermutation([1,1])
        swap_perm = N.TermPermutation([2,1,3])
        @test length(swap_perm) == 3
        @test swap_perm[1] == 2
        @test !isone(swap_perm)
        @test sprint(show,swap_perm) == "TermPermutation([2, 1, 3])"
        automorphisms = N.automorphism_generators(graph;backend=:backtracking)
        @test any(perm -> perm.images == [2,1,3],automorphisms)
        @test N.automorphism_generators(graph;backend=:bliss) == automorphisms
        @test_throws ArgumentError N.automorphism_generators(graph;backend=:nauty)

        S = N.symplectic_matrix_from_permutation(tab,swap_perm)
        @test sprint(show,S) == "SymplecticMatrix(6×6)"
        @test size(S) == (6,6)
        @test Matrix(S) == S.data
        @test N.is_symplectic_matrix(S)
        @test N._gf2_matmul(tab.paulis,S.data) == tab.paulis[swap_perm.images,:]
        @test_throws ArgumentError N.symplectic_matrix_from_permutation(tab,N.TermPermutation([1,2]))
        @test_throws ArgumentError N.SymplecticMatrix(ones(Int,2,3))
        @test_throws ArgumentError N.SymplecticMatrix(ones(Int,3,3))
        phase = N.PhaseVector(zeros(Int8,6),true)
        @test length(phase) == 6
        @test phase[1] == 0
        generator = N.SympleQGenerator(swap_perm,S,phase,true)
        clifford = N._sympleq_clifford(generator.S,generator.phase,W)
        @test mapped(clifford,hamiltonian) == hamiltonian
        p,image = N.act(clifford,word(x[1]))
        @test (N.scale(p,1),image) == (1,word(z[1]))
        p,image = N.act(clifford,word(z[1]))
        @test (N.scale(p,1),image) == (1,word(x[1]))
        @test N.recover_phase_vector(tab,swap_perm,S).verified
        @test N.sympleq_generators(tab) isa Vector{N.SympleQGenerator}
        # The low-level search still rejects constants. The public front door
        # now deliberately returns an empty symmetry for a constant objective.
        @test_throws ArgumentError N._sympleq_clifford_candidates(1.0one(x[1]))
        constant = N.detect_symmetry(1.0one(x[1]))
        @test constant.detected
        @test isempty(constant.generators) && isempty(constant.moment_only)
    end

    @testset "SympleQ Witt extension scales past 10 qubits" begin
        A = UInt8[1 0 1;0 1 1]
        b = UInt8[1,0]
        particular, kernel = N._gf2_affine_solution_space(A,b)
        @test mod.(Int.(A)*Int.(particular),2) == Int.(b)
        @test length(kernel) == 1
        @test mod.(Int.(A)*Int.(only(kernel)),2) == [0,0]
        @test N._gf2_affine_solution_space(UInt8[1 0;1 0],UInt8[0,1]) === nothing
        D = UInt8[1 0 0 0 0 0;0 1 0 0 0 0]
        R = UInt8[0 1 0 0 0 0;1 0 0 0 0 0]
        d = N._choose_domain_extension(D)
        r = N._choose_target_extension(d,D,R)
        @test !N._gf2_row_in_span(r,R)
        for i in axes(D,1)
            @test N._gf2_pairing(d,view(D,i,:)) == N._gf2_pairing(r,view(R,i,:))
        end

        # Dimension 24, rank deficiency 2: exercises the linear Witt extension
        # beyond the old exhaustive solver's dimension-20 limit.
        n = 12
        _, ops = N.pauli_variables(n)
        ring_h = sum(ComplexF64(1/4)*op[i]*op[mod1(i+1,n)] for op in ops for i in 1:n)
        generators = N._sympleq_clifford_candidates(ring_h)
        @test length(generators) >= 3
        @test all(g -> mapped(g,ring_h) == ring_h,generators)
    end
end
end

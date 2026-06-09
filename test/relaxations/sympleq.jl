# test/relaxations/sympleq.jl
# Tests: SympleQ find-side bridge into CliffordSymmetry

using Test, NCTSSoS

@testset "SympleQ bridge" begin
    @testset "tableau encodes Pauli words" begin
        _, (σx, σy, σz) = create_pauli_variables(1:2)
        H = 2.0 * σx[1] + 3.0 * σy[2] + 5.0 * σz[1] * σx[2]
        tab = NCTSSoS.SymplecticTableau(H)

        @test length(tab) == 3
        @test tab.nqubits == 2
        @test size(tab.paulis) == (3, 4)
        @test Set(tab.coeffs) == Set([2.0, 3.0, 5.0])
        @test any(row -> collect(row) == UInt8[1, 0, 0, 0], eachrow(tab.paulis))
        @test any(row -> collect(row) == UInt8[0, 1, 0, 1], eachrow(tab.paulis))
        @test any(row -> collect(row) == UInt8[0, 1, 1, 0], eachrow(tab.paulis))
    end

    @testset "binary Clifford conversion preserves Pauli SWAP action" begin
        _, (σx, σy, σz) = create_pauli_variables(1:2)
        mono(poly) = only(monomials(poly))
        act(g, mono) = NCTSSoS._act_monomial(g, mono)

        swap_rows = UInt8[
            0 1 0 0;
            1 0 0 0;
            0 0 0 1;
            0 0 1 0;
        ]
        g = clifford_to_clifford_symmetry(
            NCTSSoS.SymplecticMatrix(swap_rows),
            NCTSSoS.PhaseVector(zeros(Int8, 4), true);
            nqubits=2,
            integer_type=UInt8,
        )

        @test act(g, σx[1]) == (1, σx[2])
        @test act(g, σy[1]) == (1, σy[2])
        @test act(g, σz[2]) == (1, σz[1])
        @test act(g, mono(σx[1] * σz[2])) == (1, mono(σz[1] * σx[2]))

        cnot_rows = UInt8[
            1 1 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 1 1;
        ]
        cnot = clifford_to_clifford_symmetry(
            NCTSSoS.SymplecticMatrix(cnot_rows),
            NCTSSoS.PhaseVector(zeros(Int8, 4), true);
            nqubits=2,
            integer_type=UInt8,
        )

        @test act(cnot, σx[1]) == (1, mono(σx[1] * σx[2]))
        @test act(cnot, σz[2]) == (1, mono(σz[1] * σz[2]))
    end

    @testset "SympleQ discovers a Clifford symmetry spec" begin
        _, (_, _, σz) = create_pauli_variables(1:2)
        H = 1.0 * σz[1] + 1.0 * σz[2]
        spec = sympleq_symmetry_spec(H)

        @test isempty(spec.generators)
        @test isempty(spec.fermionic_generators)
        @test length(spec.clifford_generators) >= 1
        @test any(g -> NCTSSoS._act_polynomial(g, H) == H, spec.clifford_generators)
    end
end

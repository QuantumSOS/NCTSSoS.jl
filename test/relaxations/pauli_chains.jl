# Tests for sparse Pauli spin-chain bases and charge-compatible sign symmetry.

using Test, NCTSSoS

function _pauli_support(mono)
    return sort!(Int[NCTSSoS._pauli_site(idx) for idx in mono.word])
end

function _is_periodic_contiguous(sites::Vector{Int}, n::Int)
    isempty(sites) && return true
    target = Set(sites)
    width = length(sites)
    return any(Set(mod1(start + offset, n) for offset in 0:(width - 1)) == target for start in 1:n)
end

@testset "Pauli contiguous chain basis" begin
    reg4, (σx4, _, _) = create_pauli_variables(1:4)

    @test length(pauli_contiguous_chain_basis(reg4, 0)) == 1
    @test length(pauli_contiguous_chain_basis(reg4, 1)) == 13
    @test length(pauli_contiguous_chain_basis(reg4, 2)) == 49
    @test length(pauli_contiguous_chain_basis(reg4, 3)) == 157
    @test length(pauli_contiguous_chain_basis(reg4, 4)) == 238
    @test length(pauli_contiguous_chain_basis(reg4, 2; periodic=false)) == 40

    basis = pauli_contiguous_chain_basis(reg4, 3)
    @test one(σx4[1]) in basis
    @test all(mono -> degree(mono) <= 3, basis)
    @test all(mono -> _is_periodic_contiguous(_pauli_support(mono), 4), basis)

    reg10, _ = create_pauli_variables(1:10)
    @test length(pauli_contiguous_chain_basis(reg10, 2)) == 1 + 10 * (3 + 9)

    reg100, _ = create_pauli_variables(1:100)
    @test length(pauli_contiguous_chain_basis(reg100, 4)) == 1 + 100 * (3 + 9 + 27 + 81)

    @test_throws ArgumentError pauli_contiguous_chain_basis(reg4, -1)
end

@testset "Pauli chain basis closure under spatial and sign symmetries" begin
    n = 6
    reg, (σx, σy, σz) = create_pauli_variables(1:n)
    basis = pauli_contiguous_chain_basis(reg, 4)
    lookup = Set(basis)

    translation = pauli_site_permutation([2:n; 1])
    reflection = pauli_site_permutation(reverse(1:n))
    sign = pauli_sign_symmetry(n; integer_type=eltype(σx[1].word))

    for g in (translation, reflection, sign), mono in basis
        _, image = NCTSSoS._act_monomial(g, mono)
        @test image in lookup
    end

    @test NCTSSoS._act_monomial(sign, σx[1]) == (-1, σx[1])
    @test NCTSSoS._act_monomial(sign, σy[1]) == (-1, σy[1])
    @test NCTSSoS._act_monomial(sign, σz[1]) == (1, σz[1])
end

@testset "Sparse Pauli charge words follow the supplied chain basis" begin
    n = 6
    reg, _ = create_pauli_variables(1:n)
    basis = pauli_contiguous_chain_basis(reg, 4)

    charge_groups = NCTSSoS._pauli_charge_transform_groups(
        basis,
        PauliChargeSectorSpec(nqubits=n, max_degree=4),
        nothing,
    )
    blocks = collect(Iterators.flatten(charge_groups))

    @test sum(size(block.row_basis, 1) for block in blocks) == length(basis)
    @test Set(block.label.charge for block in blocks) == Set(-4:4)
    @test all(block -> block.provenance == :charge_sector, blocks)

    sign_group = CliffordSymmetryGroup(
        pauli_sign_symmetry(n; integer_type=eltype(basis[1].word));
        nqubits=n,
        integer_type=eltype(basis[1].word),
    )
    signed_groups = NCTSSoS._pauli_charge_transform_groups(
        basis,
        PauliChargeSectorSpec(nqubits=n, max_degree=4),
        sign_group,
    )
    signed_blocks = collect(Iterators.flatten(signed_groups))

    @test sum(size(block.row_basis, 1) for block in signed_blocks) == length(basis)
    @test all(block -> block.label.group_order == 2, signed_blocks)
    @test all(block -> block.provenance == :charge_sector, signed_blocks)
end

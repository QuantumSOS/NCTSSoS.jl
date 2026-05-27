using Test
using NCTSSoS

function _toy_periodic_integrals(; nk::Int=2, norb::Int=1)
    h1e = Dict{Int,Matrix{ComplexF64}}(
        k => fill(ComplexF64(1.0 + 0.25k), norb, norb) for k in 0:(nk - 1)
    )

    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()
    for k1 in 0:(nk - 1), k2 in 0:(nk - 1), k3 in 0:(nk - 1)
        k4 = mod(k1 + k2 - k3, nk)
        block = zeros(ComplexF64, norb, norb, norb, norb)
        block[1, 1, 1, 1] = 0.1 * (1 + k1 + 2k2 + 3k3 + 4k4)
        eri[(k1, k2, k3, k4)] = block
    end
    return h1e, eri
end

@testset "structured PQG V2RDM assembler" begin
    @testset "canonical moment helper normal-orders small fermionic words" begin
        @test NCTSSoS.canonical_moment_id(:identity, 0, 0, ()) == Int[]
        @test NCTSSoS.canonical_moment_id(:one, 0, 0, (1, 1)) == Int[-1, 1]
        @test NCTSSoS.canonical_moment_id(:raw, 0, 0, [-1, 2]) == Int[-1, 2]
        @test_throws ArgumentError NCTSSoS.canonical_moment_id(:raw, 0, 0, [1, -1])
        @test_throws ArgumentError NCTSSoS.canonical_moment_id(:bogus, 0, 0, ())
    end

    @testset "toy periodic PQG data is self-contained and deterministic" begin
        h1e, eri = _toy_periodic_integrals()
        linear = build_pqg_moment_data(
            h1e,
            eri;
            nk=2,
            norb=1,
            nelec_per_cell=1,
            blocking=:momentum,
            spin_resolved_trace=false,
            singlet_s2=false,
            include_one_d=false,
        )

        @test NCTSSoS.assert_moment_linear_data_invariants(linear) === nothing
        @test isempty(linear.identity)
        @test !isempty(linear.moments)
        @test !isempty(linear.objective_lin)
        @test !isempty(linear.zero_constraints)
        @test length(linear.psd_block_constraint_idx) == length(linear.psd_blocks_lin)
        @test issorted(linear.moments; lt=NCTSSoS.key_lt)

        families = [block.meta.origin.family for block in linear.psd_blocks_lin]
        @test count(==(:D), families) == 2
        @test count(==(:Q), families) == 2
        @test count(==(:G), families) == 2
        @test all(block.meta.cone == :HPSD for block in linear.psd_blocks_lin)
        @test all(length(block.meta.origin.block_key) == 1 for block in linear.psd_blocks_lin)
        @test isempty(linear.free_keys)
    end

    @testset "layout options change block provenance without external assets" begin
        h1e, eri = _toy_periodic_integrals()

        unblocked = build_pqg_moment_data(
            h1e,
            eri;
            nk=2,
            norb=1,
            nelec_per_cell=1,
            blocking=:none,
            spin_resolved_trace=false,
            singlet_s2=false,
            include_one_d=false,
        )
        @test [block.meta.origin.family for block in unblocked.psd_blocks_lin] == [:D, :Q, :G]

        spin_blocked = build_pqg_moment_data(
            h1e,
            eri;
            nk=2,
            norb=1,
            nelec_per_cell=1,
            blocking=:spin,
            spin_resolved_trace=false,
            singlet_s2=false,
            include_one_d=true,
        )
        @test any(block.meta.origin.family == :oneD for block in spin_blocked.psd_blocks_lin)
        @test all(length(block.meta.origin.block_key) == 2 for block in spin_blocked.psd_blocks_lin)
        @test NCTSSoS.assert_moment_linear_data_invariants(spin_blocked) === nothing

        @test_throws ArgumentError build_pqg_moment_data(
            h1e,
            eri;
            nk=2,
            norb=1,
            nelec_per_cell=1,
            blocking=:bad,
        )
    end
end

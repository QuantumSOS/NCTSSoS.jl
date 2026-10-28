using Test
using NCTSSoS
using JuMP, LinearAlgebra

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
    end

    @testset "real moments preserve an off-diagonal one-body ground energy" begin
        h = ComplexF64[-1 0.3; 0.3 2]
        eri = Dict((0, 0, 0, 0) => zeros(ComplexF64, 2, 2, 2, 2))
        options = (; nk=1, norb=2, nelec_per_cell=2, blocking=:spin,
                    singlet_s2=false, include_one_d=true)
        for real_moments in (false, true)
            linear = build_pqg_moment_data(Dict(0 => h), eri; options..., real_moments)
            @test NCTSSoS.assert_moment_linear_data_invariants(linear) === nothing
            @test all(b.meta.cone == (real_moments ? :PSD : :HPSD) for b in linear.psd_blocks_lin)
            M = typeof(first(values(linear.key_to_monomial)))
            T = eltype(linear.identity)
            P = Polynomial{FermionicAlgebra,T,ComplexF64}
            mp = NCTSSoS._moment_problem_with_linear(
                FermionicAlgebra, T, M, P, zero(P), Tuple{Symbol,Matrix{P}}[], M[], 0, linear,
            )
            sos = NCTSSoS.sos_dualize(mp)
            set_optimizer(sos.model, SOLVER)
            optimize!(sos.model)
            @test termination_status(sos.model) == MOI.OPTIMAL
            # Two opposite-spin electrons occupy the lowest spatial orbital.
            @test objective_value(sos.model) ≈ 2eigmin(Hermitian(h)) atol=1e-6
        end
        h[1, 2] += 0.2im
        h[2, 1] -= 0.2im
        @test_throws ArgumentError build_pqg_moment_data(Dict(0 => h), eri; options..., real_moments=true)
        @test !isempty(build_pqg_moment_data(Dict(0 => h), eri; options...).objective_lin)
        h = real.(h)
        eri[(0, 0, 0, 0)][1, 2, 1, 2] = 0.1im
        @test_throws ArgumentError build_pqg_moment_data(Dict(0 => h), eri; options..., real_moments=true)
    end

    @testset "fixed spin populations have explicit zero-direction equations" begin
        for (nk, norb, nelec) in ((1, 4, 4), (2, 1, 1)),
            blocking in (:momentum, :spin, :none), real_moments in (false, true)
            @testset "nk=$nk blocking=$blocking real=$real_moments" begin
                h = Dict(k => zeros(norb, norb) for k in 0:(nk - 1))
                eri = Dict((k1, k2, k3, mod(k1 + k2 - k3, nk)) =>
                    zeros(norb, norb, norb, norb)
                    for k1 in 0:(nk - 1), k2 in 0:(nk - 1), k3 in 0:(nk - 1))
                linear = build_pqg_moment_data(h, eri; nk, norb, nelec_per_cell=nelec,
                    blocking, singlet_s2=false, real_moments)
                @test NCTSSoS.assert_moment_linear_data_invariants(linear) === nothing
                @test all(b.meta.cone == (real_moments ? :PSD : :HPSD) for b in linear.psd_blocks_lin)
                # Physical ordering: (α, β) within each momentum sector.
                # For nk=2,norb=1 the signs are +,-,+,-, not +,+,-,-.
                signs = repeat(vcat(ones(Int, norb), -ones(Int, norb)), nk)
                coefficients = Dict((-p, p) => signs[p] for p in eachindex(signs))
                block = only(filter(b -> b.meta.origin.family == :G &&
                    any(row -> haskey(coefficients, Tuple(row.word)), b.meta.row_labels),
                    linear.psd_blocks_lin))
                s = [get(coefficients, Tuple(row.word), 0) for row in block.meta.row_labels]
                @test count(!iszero, s) == 2nk*norb
                F = typeof(linear.objective_lin)
                for i in 1:block.size
                    form = F(key => s[j] * coefficient for j in 1:block.size
                        for (key, coefficient) in block.entries[i, j])
                    # Split f=0 into (f+f†)/2=0 and (f-f†)/(2im)=0,
                    # independently of the assembler's splitting helper.
                    parts = if real_moments
                        (form,)
                    else
                        adj = F(linear.adjoint_key[key] => conj(c) for (key, c) in form)
                        (F(vcat([key => c/2 for (key, c) in form],
                                [key => c/2 for (key, c) in adj])),
                         F(vcat([key => c/(2im) for (key, c) in form],
                                [key => -c/(2im) for (key, c) in adj])))
                    end
                    for part in parts
                        @test isempty(part) || any(row -> row.form.terms == part.terms, linear.zero_constraints)
                    end
                end
            end
        end
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

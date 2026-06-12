using Test, NCTSSoS
import LinearAlgebra

function _test_pauli_matrix(op::Int)
    op == 0 && return ComplexF64[1 0; 0 1]
    op == 1 && return ComplexF64[0 1; 1 0]
    op == 2 && return ComplexF64[0 -im; im 0]
    op == 3 && return ComplexF64[1 0; 0 -1]
    error("bad Pauli op $op")
end

function _test_kron_all(mats)
    out = mats[1]
    for mat in mats[2:end]
        out = LinearAlgebra.kron(out, mat)
    end
    return out
end

function _test_pauli_word_matrix(mono, nsites::Int)
    mats = [_test_pauli_matrix(0) for _ in 1:nsites]
    for idx in mono.word
        site = NCTSSoS._pauli_site(idx)
        pauli_type = NCTSSoS._pauli_type(idx) + 1
        mats[site] = _test_pauli_matrix(pauli_type)
    end
    return _test_kron_all(mats)
end

function _test_expectation(poly, rho, nsites::Int)
    value = 0.0 + 0.0im
    for (coef, mono) in terms(poly)
        value += coef * LinearAlgebra.tr(rho * _test_pauli_word_matrix(mono, nsites))
    end
    return value
end

@testset "Ground-state scalar and PSD physical constraints" begin
    reg, (x, y, z) = create_pauli_variables(1:2)
    H = sum(0.25 * op[1] * op[2] for op in (x, y, z))
    basis = [one(x[1]), x[1], y[1]]

    comms = commutator_constraints(H, basis)
    @test length(comms) == 2
    @test all(!iszero, comms)
    @test iszero(comms[1] - (H * x[1] - x[1] * H))

    pop = polyopt(H, reg; scalar_moment_eq_constraints=comms)
    sparsity = compute_sparsity(pop, SolverConfig(optimizer=nothing, order=2))
    mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    @test !isempty(mp.linear.zero_constraints)
    @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing

    curv = curvature_block(H, basis)
    @test curv.cone == :HPSD
    @test size(curv) == (3, 3)
    for i in 1:3, j in 1:3
        @test iszero(curv.matrix[i, j] - adjoint(curv.matrix[j, i]))
    end

    rdm1 = rdm_block(reg, [1])
    mp_extra = NCTSSoS.moment_relax(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities;
        extra_psd_blocks=[rdm1],
    )
    @test any(block -> block.meta.origin isa NCTSSoS.PhysicalPSDOrigin, mp_extra.linear.psd_blocks_lin)
    @test NCTSSoS.assert_moment_linear_data_invariants(mp_extra.linear, mp_extra.constraints) === nothing
end

@testset "Spin RDM block reconstructs a two-qubit density matrix" begin
    reg, _ = create_pauli_variables(1:2)
    block = rdm_block(reg, [1, 2])
    @test block.cone == :HPSD
    @test size(block) == (4, 4)

    ψ = ComplexF64[1, 2im, -1, 0.5]
    ψ ./= LinearAlgebra.norm(ψ)
    ρ = ψ * ψ'

    assembled = Matrix{ComplexF64}(undef, 4, 4)
    for i in 1:4, j in 1:4
        assembled[i, j] = _test_expectation(block.matrix[i, j], ρ, 2)
    end

    @test assembled ≈ ρ atol = 1e-12
end

if @isdefined(SOLVER)
    @testset "Physical constraints solve a two-spin Heisenberg instance" begin
        reg, (x, y, z) = create_pauli_variables(1:2)
        H = sum(0.25 * op[1] * op[2] for op in (x, y, z))
        basis = [one(x[1]), x[1], y[1]]
        cfg = SolverConfig(optimizer=SOLVER, order=1)

        plain = cs_nctssos(polyopt(H, reg), cfg)
        strengthened = cs_nctssos(
            polyopt(H, reg; scalar_moment_eq_constraints=commutator_constraints(H, basis)),
            cfg;
            extra_psd_blocks=[rdm_block(reg, [1])],
        )

        @test strengthened.objective >= plain.objective - 1e-5
        @test strengthened.objective ≈ -0.75 atol = 1e-4
    end
end

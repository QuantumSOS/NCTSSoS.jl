# Public API regression for Clifford symmetry reduction on a 2-site Heisenberg model.

using Test, NCTSSoS, JuMP

if !@isdefined(SOLVER)
    using COSMO
    const SOLVER = optimizer_with_attributes(
        COSMO.Optimizer,
        "verbose" => false,
        "eps_abs" => 1e-8,
        "eps_rel" => 1e-8,
        "eps_prim_inf" => 1e-6,
        "eps_dual_inf" => 1e-6,
        "max_iter" => 50_000,
        "rho" => 1.0,
        "adaptive_rho" => true,
        "alpha" => 1.0,
        "scaling" => 10,
    )
end

if !@isdefined(flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

@testset "2-site Heisenberg model with SWAP Clifford symmetry" begin
    registry, (σx, σy, σz) = create_pauli_variables(1:2)
    heisenberg = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))
    pop = polyopt(heisenberg, registry)
    basis = [one(σx[1]); σx; σy; σz]

    plain_config = SolverConfig(
        optimizer=SOLVER,
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
    )
    symmetric_config = SolverConfig(
        optimizer=SOLVER,
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=SymmetrySpec(CliffordSymmetry(:SWAP, 1, 2)),
    )

    plain = cs_nctssos(pop, plain_config)
    symmetric = cs_nctssos(pop, symmetric_config)

    @test plain.objective ≈ -0.75 atol = 1e-6
    @test symmetric.objective ≈ plain.objective atol = 1e-6
    @test !isnothing(symmetric.symmetry)
    @test symmetric.symmetry.group_order == 2
    @test symmetric.symmetry.psd_block_sizes == [4, 3]
    @test maximum(symmetric.symmetry.psd_block_sizes) < only(flatten_sizes(plain.moment_matrix_sizes))
end

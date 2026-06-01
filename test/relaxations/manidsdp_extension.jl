# Optional extension smoke tests.  ManiDSDP is a weak dependency, so this file is
# deliberately silent unless the test environment explicitly provides it.

using LinearAlgebra

if Base.find_package("ManiDSDP") !== nothing
    @eval using ManiDSDP

    @testset "ManiDSDP Pauli extension" begin
        MOI = NCTSSoS.MOI
        silent_clarabel = MOI.OptimizerWithAttributes(
            NCTSSoS.Clarabel.Optimizer,
            MOI.Silent() => true,
        )
        manidsdp_fast = MOI.OptimizerWithAttributes(
            ManiDSDP.Optimizer,
            MOI.RawOptimizerAttribute("rank") => 2,
            MOI.RawOptimizerAttribute("max_iter") => 200,
            MOI.RawOptimizerAttribute("tol") => 1e-7,
            MOI.RawOptimizerAttribute("inner_max_iter") => 1_000,
            MOI.RawOptimizerAttribute("inner_tol") => 1e-8,
            MOI.Silent() => true,
        )

        @test Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt) !== nothing

        registry1, (_, _, σz1) = create_pauli_variables(1:1)
        pop1 = polyopt(1.0 * σz1[1], registry1)
        mani1 = cs_nctssos(pop1, SolverConfig(optimizer=manidsdp_fast, order=1))
        clar1 = cs_nctssos(pop1, SolverConfig(optimizer=silent_clarabel, order=1))
        @test mani1.objective ≈ -1.0 atol=1e-8
        @test mani1.objective ≈ clar1.objective atol=1e-6
        @test isempty(NCTSSoS.all_variables(mani1.model))

        mani1_type = cs_nctssos(pop1, SolverConfig(optimizer=ManiDSDP.Optimizer, order=1))
        @test mani1_type.objective ≈ -1.0 atol=1e-8

        @test_throws ArgumentError cs_nctssos(pop1, SolverConfig(optimizer=manidsdp_fast, order=1); dualize=false)

        constrained = polyopt(
            1.0 * σz1[1],
            registry1;
            ineq_constraints=[1.0 * one(σz1[1]) + 1.0 * σz1[1]],
        )
        @test_throws ArgumentError cs_nctssos(constrained, SolverConfig(optimizer=manidsdp_fast, order=1))

        registry2, (σx, σy, σz) = create_pauli_variables(1:2)
        heisenberg = sum(ComplexF64(0.25) * op[1] * op[2] for op in (σx, σy, σz))
        pop2 = polyopt(heisenberg, registry2)
        mani2 = cs_nctssos(pop2, SolverConfig(optimizer=manidsdp_fast, order=1))
        clar2 = cs_nctssos(pop2, SolverConfig(optimizer=silent_clarabel, order=1))
        @test mani2.objective ≈ -0.75 atol=1e-8
        @test mani2.objective ≈ clar2.objective atol=1e-6

        ext = Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt)
        sparsity2 = compute_sparsity(pop2, SolverConfig(optimizer=manidsdp_fast, order=1))
        mp2 = NCTSSoS.moment_relax(pop2, sparsity2.corr_sparsity, sparsity2.cliques_term_sparsities)
        opt2 = ext._manidsdp_optimizer_from_attributes(manidsdp_fast)
        native_problem2, native_constant2 = ext._pauli_moment_problem_to_manidsdp(mp2, opt2)
        @test native_constant2 == 0.0
        @test native_problem2.gram isa Diagonal
        dense_problem2 = ManiDSDP.DualSDP(native_problem2.C, native_problem2.A, native_problem2.b; rank=native_problem2.rank)
        @test Matrix(native_problem2.gram) ≈ dense_problem2.gram atol=1e-12
        @test native_problem2.D ≈ dense_problem2.D atol=1e-12
    end
end

using Test, NCTSSoS, JuMP, SparseArrays
const MOI = JuMP.MOI

const RUN_H2_NK2_BRIDGE_SMOKE = get(ENV, "NCTSSOS_RUN_H2_NK2_BRIDGE_SMOKE", "0") == "1"
const H2_NK2_DEMO_PATH = joinpath(@__DIR__, "..", "..", "demos", "h2_periodic_nk2_moment_sos.jl")

module H2Nk2BridgeDemo end

if !RUN_H2_NK2_BRIDGE_SMOKE
    @testset "H2/Nk=2 :psd_blocks BPSDP bridge smoke" begin
        @test_skip "set NCTSSOS_RUN_H2_NK2_BRIDGE_SMOKE=1 to run the HAI-only BPSDP bridge smoke"
    end
else
    import BPSDP

    Base.include(H2Nk2BridgeDemo, H2_NK2_DEMO_PATH)

    function _h2_nk2_moment_problem()
        options = H2Nk2BridgeDemo.parse_options(String[
            "--bpsdp-print-level=0",
            "--output-dir=$(mktempdir())",
        ])
        data = H2Nk2BridgeDemo.build_h2_pqg_moment_problem(options)
        return data.moment_problem
    end

    function _h2_bpsdp_copy_optimizer()
        return BPSDP.Optimizer(
            max_iter = 0,
            cg_max_iter = 1,
            mu_update_frequency = 25,
            penalty_parameter = 0.1,
            cg_convergence = 1e-12,
            dynamic_cg_convergence = false,
            sdp_objective_convergence = 1e-8,
            sdp_error_convergence = 1e-8,
            guess_type = :zero,
            print_level = 0,
            dependent_rows = :keep,
        )
    end

    function _num_constraints(backend, F::Type, S::Type)
        return MOI.get(backend, MOI.NumberOfConstraints{F,S}())
    end

    function _count_1x1_constraints(backend, F::Type, S::Type)
        total = 0
        for ci in MOI.get(backend, MOI.ListOfConstraintIndices{F,S}())
            set = MOI.get(backend, MOI.ConstraintSet(), ci)
            total += Int(getfield(set, :side_dimension) == 1)
        end
        return total
    end

    function _moi_scalar_psd_1x1_count(backend)
        Fs = (MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64})
        Ss = (MOI.PositiveSemidefiniteConeTriangle, MOI.HermitianPositiveSemidefiniteConeTriangle)
        return sum(_count_1x1_constraints(backend, F, S) for F in Fs for S in Ss; init = 0)
    end

    function _raw_scalar_1x1_psd_count(raw)
        dims = getfield(raw, :block_dims)
        kinds = getfield(raw, :block_kinds)
        return count(idx -> dims[idx] == 1 && kinds[idx] == :real, eachindex(dims))
    end

    function _copy_to_bpsdp!(model)
        set_optimizer(model, _h2_bpsdp_copy_optimizer)
        set_silent(model)
        optimize!(model)
        return H2Nk2BridgeDemo.raw_optimizer(model)
    end

    @testset "H2/Nk=2 :psd_blocks delivers BPSDP-favorable cones" begin
        mp = _h2_nk2_moment_problem()
        free_keys = mp.linear.free_keys
        n_free = length(free_keys)
        n_symbolic_hpsd = count(block -> block.meta.cone == :HPSD, mp.linear.psd_blocks_lin)
        n_aux_hpsd = cld(n_free, NCTSSoS.AUX_ORPHANS_PER_BLOCK)

        # The explicit N̂ - N constraint is intentionally absent here. TrD/TrG
        # plus the identity-augmented ²G block imply particle number without
        # producing orphan moment keys.
        @test n_free == 0
        @test isempty(NCTSSoS.orphan_keys(mp))

        model, _extract = build_jump_model(mp;
            formulation = :psd_blocks,
            representation = :complex,
            orphan_policy = :aux_psd_free,
        )

        backend = JuMP.backend(model)
        @test _num_constraints(backend, MOI.VectorOfVariables, MOI.HermitianPositiveSemidefiniteConeTriangle) ==
              n_symbolic_hpsd + n_aux_hpsd
        @test _moi_scalar_psd_1x1_count(backend) < 100

        raw = _copy_to_bpsdp!(model)
        dims = getfield(raw, :block_dims)
        kinds = getfield(raw, :block_kinds)

        @test count(==(:hermitian), kinds) == n_symbolic_hpsd + n_aux_hpsd
        @test _raw_scalar_1x1_psd_count(raw) == 0
        @test _raw_scalar_1x1_psd_count(raw) < 100  # hard ceiling vs old 406,706
        @test length(dims) == n_symbolic_hpsd + n_aux_hpsd
        @test length(getfield(raw, :b)) > 0
        @test nnz(getfield(raw, :A)) > 0
    end
end

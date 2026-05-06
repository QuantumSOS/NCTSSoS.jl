using Test, NCTSSoS, JuMP
const MOI = JuMP.MOI

const RUN_H2_NK2_BRIDGE_SMOKE = get(ENV, "NCTSSOS_RUN_H2_NK2_BRIDGE_SMOKE", "0") == "1"

if !RUN_H2_NK2_BRIDGE_SMOKE
    @testset "H2/Nk=2 :psd_blocks BPSDP bridge smoke" begin
        @test_skip false
    end
else
    import BPSDP

    include(joinpath(@__DIR__, "..", "..", "demos", "h2_periodic_nk2_moment_sos.jl"))

    function _h2_nk2_moment_problem()
        options = parse_options(String[
            "--bpsdp-print-level=0",
            "--output-dir=$(mktempdir())",
        ])
        data = build_h2_pqg_moment_problem(options)
        return data.moment_problem
    end

    function _h2_bpsdp_optimizer()
        return BPSDP.Optimizer(
            max_iter = 1,
            cg_max_iter = 5000,
            mu_update_frequency = 25,
            penalty_parameter = 0.1,
            cg_convergence = 1e-12,
            dynamic_cg_convergence = false,
            sdp_objective_convergence = 1e-8,
            sdp_error_convergence = 1e-8,
            guess_type = :zero,
            print_level = 0,
            dependent_rows = :drop,
        )
    end

    @testset "H2/Nk=2 :psd_blocks uses free orphan variables and survives first BPSDP step" begin
        mp = _h2_nk2_moment_problem()
        n_orphans = length(NCTSSoS.orphan_keys(mp))
        n_symbolic_hpsd = count(c -> c[1] == :HPSD, mp.constraints)
        @test n_orphans > 0

        model, _extract = build_jump_model(mp;
            formulation = :psd_blocks,
            representation = :complex,
            orphan_policy = :free_variables,
        )
        set_optimizer(model, _h2_bpsdp_optimizer)
        optimize!(model)

        raw = raw_optimizer(model)
        dims = getfield(raw, :block_dims)
        kinds = getfield(raw, :block_kinds)
        A = getfield(raw, :A)
        c = getfield(raw, :c)
        state = getfield(raw, :state)

        @test count(==(:hermitian), kinds) == n_symbolic_hpsd
        @test count(==(1), dims) >= 2 * n_orphans
        @test count(j -> diff(A.colptr)[j] == 0 && iszero(c[j]), eachindex(c)) == 0
        @test size(A, 1) < 42336

        @test termination_status(model) != MOI.NUMERICAL_ERROR
        @test MOI.get(model, MOI.RawStatusString()) != "cg_failure"
        @test getfield(state, :status) != :cg_failure
        @test getfield(state, :iterations) >= 1 || getfield(state, :status) == :optimal
        @test getfield(state, :inner_iterations) > 0
    end
end

using Test

if !isdefined(@__MODULE__, :TestExpectations)
    include("Expectations.jl")
end
using .TestExpectations: expectations_oracle

module QMBReferenceHarnessTests
include("../perf/qmbcertify_reference_runs.jl")
end

module QMBProfileProbeHarnessTests
include("../perf/qmbcertify_profile_probe.jl")
end

module PauliTranslationProfileHarnessTests
include("../perf/pauli_translation_profile.jl")
end

module PauliSparseChainD4HarnessTests
include("../perf/pauli_sparse_chain_d4_blocks.jl")
end

module PauliChargeSingletPrepHarnessTests
include("../perf/pauli_charge_singlet_prep.jl")
end

module PauliTranslationCompareHarnessTests
include("../perf/pauli_translation_compare.jl")
end

module PauliTranslationSolverProbeHarnessTests
include("../perf/pauli_translation_solver_probe.jl")
end

module FakeMosekVersionHarnessTests
getversion() = (11, 2, 0)
end

function _test_structural_target_only(target)
    @test target["status"] == "analytic_target"
    @test target["solve_supported"] == false
    @test target["solve_blocker"] == "structural_target_only"
    @test occursin("Structural target only", target["solve_blocker_reason"])
    @test target["estimated_model_size_gate_status"] ==
        "blocked_missing_scalar_equality_estimate"
    @test occursin(
        "Structural target only",
        target["estimated_model_size_gate_reason"],
    )
    @test target["solve_unsupported_block_features"] == []
    @test target["solve_unsupported_zero_features"] == []
    @test target["requires_construction"] == false
end

function _test_block_domain_histograms(target, coefficient_histogram, exact_histogram)
    @test target["block_coefficient_domain_histogram"] == coefficient_histogram
    @test target["block_exact_coefficient_domain_histogram"] == exact_histogram
end

function _check_source_base_solved_reference_contract(probe)
    get(probe, "qmbcertify_base_construct", false) == true || return nothing
    get(probe, "solved", false) == true || return nothing

    required_fields = (
        "qmbcertify_probe_id",
        "qmbcertify_objective_total_estimate",
        "objective_minus_qmbcertify_total_estimate",
        "objective",
    )
    missing_fields = [field for field in required_fields if !haskey(probe, field)]
    if !isempty(missing_fields)
        throw(
            ArgumentError(
                "solved source-base probe $(get(probe, "id", "<unknown>")) " *
                "is missing QMBCertify reference fields: " *
                join(missing_fields, ", "),
            ),
        )
    end

    qmbcertify_probe_id = probe["qmbcertify_probe_id"]
    if !(qmbcertify_probe_id isa AbstractString) || isempty(strip(qmbcertify_probe_id))
        throw(
            ArgumentError(
                "solved source-base probe $(get(probe, "id", "<unknown>")) " *
                "has an empty qmbcertify_probe_id",
            ),
        )
    end

    expected_delta = probe["objective"] - probe["qmbcertify_objective_total_estimate"]
    actual_delta = probe["objective_minus_qmbcertify_total_estimate"]
    if !isapprox(actual_delta, expected_delta; atol = 1e-12, rtol = 1e-12)
        throw(
            ArgumentError(
                "solved source-base probe $(get(probe, "id", "<unknown>")) " *
                "has inconsistent objective_minus_qmbcertify_total_estimate",
            ),
        )
    end

    return nothing
end

function _check_deleted_primal_su2_solved_fixture_contract(probe)
    get(probe, "solved", false) == true || return nothing
    get(probe, "model_mode", "") == "primal" || return nothing
    get(probe, "formulation", "") == "moment_variables" || return nothing
    get(probe, "su2_symmetry", false) == true || return nothing
    get(probe, "su2_moment_quotient", false) == true && return nothing
    get(probe, "contiguous_rdm_decomposition", "") == "su2" || return nothing
    get(probe, "linear_state_opt_width", nothing) == 7 || return nothing
    get(probe, "psd_state_opt_width", nothing) == 3 || return nothing

    haskey(probe, "length") || throw(
        ArgumentError(
            "solved primal SU(2) probe $(get(probe, "id", "<unknown>")) " *
            "is missing length, so the deleted-run gate cannot be evaluated",
        ),
    )
    n = Int(probe["length"])
    n < 14 && return nothing

    throw(
        ArgumentError(
            "solved primal moment-variable SU(2) RDM/LSO/PSO fixture row " *
            "$(get(probe, "id", "<unknown>")) at length $n is deleted under " *
            "the current formulation; keep no-solve sizing rows or smaller " *
            "smoke solves until a formulation change updates the gate",
        ),
    )
    return nothing
end

function _check_su2_moment_quotient_probe_contract(probe)
    get(probe, "su2_moment_quotient", false) == true || return nothing
    required = (
        "su2_moment_raw_count",
        "su2_moment_quotient_count",
        "su2_moment_quotient_reduction_ratio",
        "su2_moment_support_orbit_count",
        "su2_moment_max_pivot_residual",
        "su2_moment_max_invariant_residual",
        "su2_moment_max_reconstruction_residual",
        "su2_moment_max_condition",
        "su2_moment_eliminated_zero_row_count",
        "estimated_scalar_equalities_upper_bound",
        "estimated_dense_schur_bytes",
    )
    missing = [field for field in required if !haskey(probe, field)]
    isempty(missing) || throw(ArgumentError(
        "SU(2) moment quotient probe $(get(probe, "id", "<unknown>")) is " *
        "missing fields: $(join(missing, ", "))",
    ))

    raw_count = Int(probe["su2_moment_raw_count"])
    quotient_count = Int(probe["su2_moment_quotient_count"])
    0 < quotient_count < raw_count || throw(ArgumentError(
        "SU(2) moment quotient probe has invalid raw/quotient counts",
    ))
    haskey(probe, "linear_moments") &&
        Int(probe["linear_moments"]) != quotient_count && throw(ArgumentError(
            "SU(2) moment quotient count does not match linear_moments",
        ))
    ratio = Float64(probe["su2_moment_quotient_reduction_ratio"])
    isfinite(ratio) && 0 < ratio < 1 &&
        isapprox(ratio, quotient_count / raw_count; atol=1e-12, rtol=1e-12) ||
        throw(ArgumentError("SU(2) moment quotient reduction ratio is inconsistent"))
    Int(probe["su2_moment_support_orbit_count"]) > 0 || throw(ArgumentError(
        "SU(2) moment quotient probe has no support orbits",
    ))
    for field in (
        "su2_moment_max_pivot_residual",
        "su2_moment_max_invariant_residual",
        "su2_moment_max_reconstruction_residual",
    )
        residual = Float64(probe[field])
        isfinite(residual) && 0 <= residual <= 1e-8 || throw(ArgumentError(
            "SU(2) moment quotient probe has invalid $field",
        ))
    end
    condition = Float64(probe["su2_moment_max_condition"])
    isfinite(condition) && condition > 0 || throw(ArgumentError(
        "SU(2) moment quotient probe has invalid coordinate condition",
    ))
    Int(probe["su2_moment_eliminated_zero_row_count"]) >= 0 || throw(ArgumentError(
        "SU(2) moment quotient probe has a negative eliminated-row count",
    ))

    estimated_equalities = Int(probe["estimated_scalar_equalities_upper_bound"])
    estimated_equalities > 0 || throw(ArgumentError(
        "SU(2) moment quotient probe has no scalar-equality estimate",
    ))
    Int(probe["estimated_dense_schur_bytes"]) == 8 * estimated_equalities^2 ||
        throw(ArgumentError("SU(2) moment quotient dense-Schur estimate is inconsistent"))
    if get(probe, "model_built", false)
        haskey(probe, "scalar_equalities") && haskey(probe, "dense_schur_bytes") ||
            throw(ArgumentError(
                "model-built SU(2) moment quotient probe lacks actual equation telemetry",
            ))
        actual_equalities = Int(probe["scalar_equalities"])
        0 < actual_equalities <= estimated_equalities || throw(ArgumentError(
            "SU(2) moment quotient actual scalar-equality count is inconsistent",
        ))
        Int(probe["dense_schur_bytes"]) == 8 * actual_equalities^2 ||
            throw(ArgumentError("SU(2) moment quotient actual dense-Schur proxy is inconsistent"))
    end
    if get(probe, "solved", false)
        for field in (
            "primal_status",
            "dual_status",
            "raw_status",
            "result_count",
        )
            haskey(probe, field) || throw(ArgumentError(
                "solved SU(2) moment quotient probe lacks $field telemetry",
            ))
        end
        Int(probe["result_count"]) > 0 || throw(ArgumentError(
            "solved SU(2) moment quotient probe has no solver result",
        ))
        if haskey(probe, "relative_gap")
            gap = Float64(probe["relative_gap"])
            isfinite(gap) && gap >= 0 || throw(ArgumentError(
                "solved SU(2) moment quotient probe has an invalid relative gap",
            ))
        end
    end
    return nothing
end

function _test_probe_execution_state_contract(probe)
    @test haskey(probe, "construction_only")
    @test haskey(probe, "model_built")
    @test haskey(probe, "solved")
    @test haskey(probe, "termination_status")

    construction_only = probe["construction_only"]
    model_built = probe["model_built"]
    solved = probe["solved"]
    @test construction_only isa Bool
    @test model_built isa Bool
    @test solved isa Bool

    @test !(construction_only && model_built)
    @test !(construction_only && solved)
    @test !(solved && !model_built)

    termination_status = probe["termination_status"]
    @test termination_status in ("OPTIMAL", "not_solved", "timeout")

    if construction_only
        @test termination_status == "not_solved"
        @test !haskey(probe, "dualization_time_seconds")
        @test !haskey(probe, "lowering_time_seconds")
        @test !haskey(probe, "solve_time_seconds")
        @test !haskey(probe, "objective")
        @test !haskey(probe, "sos_certificate_moment_count")
    elseif model_built
        @test haskey(probe, "model_mode")
    end

    if solved
        @test termination_status == "OPTIMAL"
        @test haskey(probe, "solve_time_seconds")
        @test haskey(probe, "objective")
    else
        @test termination_status in ("not_solved", "timeout")
        @test !haskey(probe, "objective")
        @test !haskey(probe, "sos_certificate_moment_count")
    end

    @test _check_source_base_solved_reference_contract(probe) === nothing
    @test _check_deleted_primal_su2_solved_fixture_contract(probe) === nothing
    @test _check_su2_moment_quotient_probe_contract(probe) === nothing

    if termination_status == "timeout"
        @test model_built
        @test !construction_only
        @test !solved
        @test probe["timeout_exit_status"] == 124
        @test probe["run_queue_status"] == "deleted_until_evidence_gate"
        @test occursin("formulation_change", probe["restore_gate"])
        @test occursin("wall_rss_load_estimate", probe["restore_gate"])
        @test get(probe, "larger_solved_probes_active", true) == false
        @test haskey(probe, "deleted_solved_probe_lengths")
        @test !isempty(probe["deleted_solved_probe_lengths"])
    end
end

function _test_fixture_probe_execution_states(qmb)
    for field in (
        "nctssos_source_like_solver_probes",
        "nctssos_large_construction_probes",
        "nctssos_qmbcertify_source_base_model_probes",
        "nctssos_su2_formulation_probes",
    )
        haskey(qmb, field) || continue
        for probe in qmb[field]
            _test_probe_execution_state_contract(probe)
        end
    end
end

function _with_env(f, pairs)
    old_values = Dict{String,Union{Nothing,String}}()
    for (key, value) in pairs
        name = String(key)
        old_values[name] = get(ENV, name, nothing)
        value === nothing ? delete!(ENV, name) : (ENV[name] = String(value))
    end
    try
        return f()
    finally
        for (name, value) in old_values
            value === nothing ? delete!(ENV, name) : (ENV[name] = value)
        end
    end
end

@testset "Source-base Solved Probe QMBCertify Reference Contract" begin
    solved_source_base_probe = Dict{String,Any}(
        "qmbcertify_base_construct" => true,
        "construction_only" => false,
        "model_built" => true,
        "solved" => true,
        "termination_status" => "OPTIMAL",
        "solve_time_seconds" => 1.0,
        "objective" => -8.0,
    )

    @test_throws ArgumentError _check_source_base_solved_reference_contract(
        solved_source_base_probe,
    )

    solved_source_base_probe["qmbcertify_probe_id"] = "A1_L20"
    solved_source_base_probe["qmbcertify_objective_total_estimate"] = -7.5
    solved_source_base_probe["objective_minus_qmbcertify_total_estimate"] = -0.5
    @test _check_source_base_solved_reference_contract(solved_source_base_probe) ===
        nothing

    inconsistent_delta_probe = copy(solved_source_base_probe)
    inconsistent_delta_probe["objective_minus_qmbcertify_total_estimate"] = 1.0
    @test_throws ArgumentError _check_source_base_solved_reference_contract(
        inconsistent_delta_probe,
    )

    construction_only_source_base_probe = copy(solved_source_base_probe)
    construction_only_source_base_probe["construction_only"] = true
    construction_only_source_base_probe["model_built"] = false
    construction_only_source_base_probe["solved"] = false
    construction_only_source_base_probe["termination_status"] = "not_solved"
    delete!(construction_only_source_base_probe, "qmbcertify_probe_id")
    delete!(
        construction_only_source_base_probe,
        "qmbcertify_objective_total_estimate",
    )
    delete!(
        construction_only_source_base_probe,
        "objective_minus_qmbcertify_total_estimate",
    )
    @test _check_source_base_solved_reference_contract(
        construction_only_source_base_probe,
    ) === nothing
end

@testset "Deleted Primal SU(2) Solved Fixture Contract" begin
    solved_large_primal_probe = Dict{String,Any}(
        "id" => "synthetic_deleted_primal_su2_solve",
        "length" => 14,
        "model_mode" => "primal",
        "formulation" => "moment_variables",
        "su2_symmetry" => true,
        "contiguous_rdm_decomposition" => "su2",
        "linear_state_opt_width" => 7,
        "psd_state_opt_width" => 3,
        "construction_only" => false,
        "model_built" => true,
        "solved" => true,
        "termination_status" => "OPTIMAL",
        "solve_time_seconds" => 1.0,
        "objective" => -1.0,
    )

    @test_throws ArgumentError _check_deleted_primal_su2_solved_fixture_contract(
        solved_large_primal_probe,
    )

    no_solve_large_primal_probe = copy(solved_large_primal_probe)
    no_solve_large_primal_probe["solved"] = false
    no_solve_large_primal_probe["termination_status"] = "not_solved"
    delete!(no_solve_large_primal_probe, "solve_time_seconds")
    delete!(no_solve_large_primal_probe, "objective")
    @test _check_deleted_primal_su2_solved_fixture_contract(
        no_solve_large_primal_probe,
    ) === nothing

    solved_small_primal_probe = copy(solved_large_primal_probe)
    solved_small_primal_probe["length"] = 12
    @test _check_deleted_primal_su2_solved_fixture_contract(
        solved_small_primal_probe,
    ) === nothing
end

@testset "SU(2) Moment Quotient Probe Contract" begin
    probe = Dict{String,Any}(
        "id" => "synthetic_su2_moment_quotient",
        "su2_moment_quotient" => true,
        "su2_moment_raw_count" => 100,
        "su2_moment_quotient_count" => 10,
        "su2_moment_quotient_reduction_ratio" => 0.1,
        "su2_moment_support_orbit_count" => 4,
        "su2_moment_max_pivot_residual" => 1.0e-13,
        "su2_moment_max_invariant_residual" => 2.0e-13,
        "su2_moment_max_reconstruction_residual" => 1.0e-12,
        "su2_moment_max_condition" => 5.0,
        "su2_moment_eliminated_zero_row_count" => 90,
        "linear_moments" => 10,
        "estimated_scalar_equalities_upper_bound" => 9,
        "estimated_dense_schur_bytes" => 8 * 9^2,
        "scalar_equalities" => 8,
        "dense_schur_bytes" => 8 * 8^2,
        "construction_only" => false,
        "model_built" => true,
        "solved" => false,
        "termination_status" => "not_solved",
        "model_mode" => "sos_dual",
    )
    @test _check_su2_moment_quotient_probe_contract(probe) === nothing
    @test _test_probe_execution_state_contract(probe) === nothing

    bad_count = copy(probe)
    bad_count["su2_moment_quotient_count"] = 100
    @test_throws ArgumentError _check_su2_moment_quotient_probe_contract(bad_count)

    bad_residual = copy(probe)
    bad_residual["su2_moment_max_reconstruction_residual"] = Inf
    @test_throws ArgumentError _check_su2_moment_quotient_probe_contract(bad_residual)

    missing_actual = copy(probe)
    delete!(missing_actual, "scalar_equalities")
    @test_throws ArgumentError _check_su2_moment_quotient_probe_contract(missing_actual)

    solved_probe = copy(probe)
    solved_probe["solved"] = true
    solved_probe["termination_status"] = "OPTIMAL"
    solved_probe["solve_time_seconds"] = 1.0
    solved_probe["objective"] = -1.0
    @test_throws ArgumentError _check_su2_moment_quotient_probe_contract(solved_probe)
    solved_probe["primal_status"] = "FEASIBLE_POINT"
    solved_probe["dual_status"] = "FEASIBLE_POINT"
    solved_probe["raw_status"] = "solved"
    solved_probe["result_count"] = 1
    solved_probe["relative_gap"] = 1.0e-10
    @test _check_su2_moment_quotient_probe_contract(solved_probe) === nothing
end

@testset "QMBCertify Profile Probe Lengths" begin
    _with_env(("NCTS_QMB_NS" => "8,10",)) do
        @test QMBProfileProbeHarnessTests._profile_probe_lengths() == [8, 10]
    end
end

@testset "QMBCertify Harness Guards and Probe Overrides" begin
    _with_env((
        "NCTS_QMB_ALLOW_LARGE" => nothing,
        "NCTS_QMB_MAX_L" => "13",
    )) do
        @test QMBReferenceHarnessTests._check_size_guard([13]) === nothing
        @test_throws ArgumentError QMBReferenceHarnessTests._check_size_guard([14])
    end

    _with_env((
        "NCTS_QMB_ALLOW_LARGE" => nothing,
        "NCTS_QMB_MAX_L" => "13",
        "NCTS_QMB_NS" => "14",
        "NCTS_QMB_PROFILE" => "A1",
        "NCTS_QMB_BOOTSTRAP_ONLY" => "true",
    )) do
        @test_throws ArgumentError QMBProfileProbeHarnessTests.probe_main()
    end

    _with_env((
        "NCTS_QMB_ALLOW_LARGE" => "true",
        "NCTS_QMB_MAX_L" => "13",
        "NCTS_QMB_ESTIMATED_WALL_SECONDS" => nothing,
        "NCTS_QMB_ESTIMATED_RSS_GIB" => nothing,
        "NCTS_QMB_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test_throws ArgumentError QMBReferenceHarnessTests._check_size_guard([14])
    end

    _with_env((
        "NCTS_QMB_ALLOW_LARGE" => "true",
        "NCTS_QMB_MAX_L" => "13",
        "NCTS_QMB_NS" => "14",
        "NCTS_QMB_PROFILE" => "A1",
        "NCTS_QMB_BOOTSTRAP_ONLY" => "false",
        "NCTS_QMB_ESTIMATED_WALL_SECONDS" => nothing,
        "NCTS_QMB_ESTIMATED_RSS_GIB" => nothing,
        "NCTS_QMB_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test_throws ArgumentError QMBProfileProbeHarnessTests.probe_main()
    end

    _with_env((
        "NCTS_QMB_ALLOW_LARGE" => "true",
        "NCTS_QMB_MAX_L" => "13",
        "NCTS_QMB_ESTIMATED_WALL_SECONDS" => "1",
        "NCTS_QMB_ESTIMATED_RSS_GIB" => "0.001",
        "NCTS_QMB_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test QMBReferenceHarnessTests._check_size_guard([14]) === nothing
    end

    _with_env((
        "NCTS_QMB_ARG_OVERRIDES" => "lso = false\npso = 0",
        "NCTS_QMB_ARG_OVERRIDES_FILE" => nothing,
    )) do
        overrides = QMBProfileProbeHarnessTests._probe_overrides()
        @test overrides["lso"] == false
        @test overrides["pso"] == 0
        args, parsed_overrides = QMBProfileProbeHarnessTests._probe_args("A1")
        @test parsed_overrides == overrides
        @test args["lso"] == false
        @test args["pso"] == 0
        @test args["rdm"] == 8
    end

    _with_env((
        "NCTS_QMB_ARG_OVERRIDES" => "lso = false",
        "NCTS_QMB_ARG_OVERRIDES_FILE" => "unused.toml",
    )) do
        @test_throws ArgumentError QMBProfileProbeHarnessTests._probe_overrides()
    end

    @test QMBProfileProbeHarnessTests._probe_accepted_statuses("A2") ==
        ("OPTIMAL", "ALMOST_OPTIMAL")
    @test QMBProfileProbeHarnessTests._toml_string_array(
        QMBProfileProbeHarnessTests._probe_accepted_statuses("A2"),
    ) == "[\"OPTIMAL\", \"ALMOST_OPTIMAL\"]"
    _with_env((
        "NCTS_QMB_ESTIMATED_WALL_SECONDS" => "12",
        "NCTS_QMB_ESTIMATED_RSS_GIB" => "3.5",
        "NCTS_LOAD_GUARD_ESTIMATED_WALL_SECONDS" => nothing,
        "NCTS_LOAD_GUARD_ESTIMATED_RSS_GIB" => nothing,
    )) do
        @test QMBReferenceHarnessTests._qmb_estimated_wall_seconds() == 12.0
        @test QMBReferenceHarnessTests._qmb_estimated_rss_gib() == 3.5
    end
    _with_env((
        "NCTS_QMB_ESTIMATED_WALL_SECONDS" => nothing,
        "NCTS_QMB_ESTIMATED_RSS_GIB" => nothing,
        "NCTS_LOAD_GUARD_ESTIMATED_WALL_SECONDS" => "21",
        "NCTS_LOAD_GUARD_ESTIMATED_RSS_GIB" => "4.25",
    )) do
        @test QMBReferenceHarnessTests._qmb_estimated_wall_seconds() == 21.0
        @test QMBReferenceHarnessTests._qmb_estimated_rss_gib() == 4.25
    end
    @test QMBProfileProbeHarnessTests._check_large_profile_probe_variant(
        "A1",
        [20, 30],
        Dict{String,Any}(),
    ) === nothing
    @test QMBProfileProbeHarnessTests._check_large_profile_probe_variant(
        "A2",
        [12],
        Dict{String,Any}(),
    ) === nothing
    @test_throws ArgumentError QMBProfileProbeHarnessTests._check_large_profile_probe_variant(
        "A2",
        [20],
        Dict{String,Any}(),
    )
    @test QMBProfileProbeHarnessTests._check_large_profile_probe_variant(
        "A2",
        [20],
        Dict{String,Any}("pso" => 2),
    ) === nothing

    profile_probe_source = read(
        normpath(joinpath(@__DIR__, "..", "perf", "qmbcertify_profile_probe.jl")),
        String,
    )
    bootstrap_branch = findfirst("if bootstrap_only", profile_probe_source)
    load_call = findfirst("qmb = _load_qmbcertify", profile_probe_source)
    @test bootstrap_branch !== nothing
    @test load_call !== nothing
    @test first(bootstrap_branch) < first(load_call)
end

@testset "Pauli Perf Harness Guards" begin
    _with_env((
        "NCTS_TRANSLATION_ALLOW_LARGE" => nothing,
        "NCTS_TRANSLATION_MAX_N" => "13",
        "NCTS_TRANSLATION_NS" => "14",
        "NCTS_TRANSLATION_TARGET_ONLY" => "true",
    )) do
        @test PauliTranslationProfileHarnessTests._check_size_guard([13]) === nothing
        @test_throws ArgumentError PauliTranslationProfileHarnessTests._check_size_guard([14])
        @test_throws ArgumentError PauliTranslationProfileHarnessTests.main()
    end

    _with_env((
        "NCTS_TRANSLATION_ALLOW_LARGE" => "true",
        "NCTS_TRANSLATION_MAX_N" => "13",
        "NCTS_TRANSLATION_ESTIMATED_WALL_SECONDS" => nothing,
        "NCTS_TRANSLATION_ESTIMATED_RSS_GIB" => nothing,
        "NCTS_TRANSLATION_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test_throws ArgumentError PauliTranslationProfileHarnessTests._check_size_guard([14])
    end

    _with_env((
        "NCTS_TRANSLATION_ALLOW_LARGE" => "true",
        "NCTS_TRANSLATION_MAX_N" => "13",
        "NCTS_TRANSLATION_ESTIMATED_WALL_SECONDS" => "1",
        "NCTS_TRANSLATION_ESTIMATED_RSS_GIB" => "0.001",
        "NCTS_TRANSLATION_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test PauliTranslationProfileHarnessTests._check_size_guard([14]) === nothing
    end

    _with_env((
        "NCTS_TRANSLATION_WARMUP" => "false",
    )) do
        @test !PauliTranslationProfileHarnessTests._parse_bool_env(
            "NCTS_TRANSLATION_WARMUP",
            true,
        )
    end
    _with_env((
        "NCTS_TRANSLATION_REPEATS" => "2",
        "NCTS_TRANSLATION_REPEAT_STAGE_LIMIT" => "12",
    )) do
        @test PauliTranslationProfileHarnessTests._parse_positive_int_env(
            "NCTS_TRANSLATION_REPEATS",
            1,
        ) == 2
        @test PauliTranslationProfileHarnessTests._parse_positive_int_env(
            "NCTS_TRANSLATION_REPEAT_STAGE_LIMIT",
            6,
        ) == 12
    end
    _with_env((
        "NCTS_TRANSLATION_REPEATS" => "0",
    )) do
        @test_throws ArgumentError PauliTranslationProfileHarnessTests._parse_positive_int_env(
            "NCTS_TRANSLATION_REPEATS",
            1,
        )
    end
    repeat_result_1 = (
        n=4,
        run_index=1,
        construction_time_seconds=1.0,
        outer_relaxation_wall_time_seconds=1.25,
        allocated_bytes=1024,
        gc_time_seconds=0.0,
        linear_moments=7,
        zero_constraints=11,
        psd_blocks=3,
        solver_facing_max_block=5,
        contiguous_rdm_blocks=2,
        psd_state_opt_blocks=0,
        construction_stage_time_seconds=Dict(
            :su2_base_wigner_rows => 3.0,
            :setup => 0.5,
        ),
    )
    repeat_result_2 = merge(
        repeat_result_1,
        (
            run_index=2,
            construction_time_seconds=1.5,
            outer_relaxation_wall_time_seconds=1.75,
            allocated_bytes=2048,
            construction_stage_time_seconds=Dict(
                :su2_base_wigner_rows => 4.0,
                :setup => 0.25,
                :other => 0.1,
            ),
        ),
    )
    @test PauliTranslationProfileHarnessTests._profile_repeat_shape(repeat_result_1) ==
        PauliTranslationProfileHarnessTests._profile_repeat_shape(repeat_result_2)
    @test PauliTranslationProfileHarnessTests._format_seconds_range([
        repeat_result_1.construction_time_seconds,
        repeat_result_2.construction_time_seconds,
    ]) == "1.000000000 .. 1.500000000"
    @test PauliTranslationProfileHarnessTests._format_bytes_range([
        repeat_result_1.allocated_bytes,
        repeat_result_2.allocated_bytes,
    ]) == "1.00 KiB .. 2.00 KiB"
    repeat_io = IOBuffer()
    PauliTranslationProfileHarnessTests._print_repeat_summary(
        [
            repeat_result_1,
            repeat_result_2,
        ],
        repeat_io,
    )
    repeat_summary = String(take!(repeat_io))
    @test occursin("#### Repeat Summary", repeat_summary)
    @test occursin("| 4 | 2 |", repeat_summary)
    @test occursin("linear=7, zero=11, psd=3, max=5, rdm=2, pso=0", repeat_summary)
    @test occursin("#### Repeat Stage Summary", repeat_summary)
    @test occursin("| 4 | su2_base_wigner_rows | 2 | 3.000000000 .. 4.000000000 |", repeat_summary)
    @test occursin("| 4 | setup | 2 | 0.250000000 .. 0.500000000 |", repeat_summary)

    repeat_limited_io = IOBuffer()
    PauliTranslationProfileHarnessTests._print_repeat_summary(
        [
            repeat_result_1,
            repeat_result_2,
        ],
        repeat_limited_io;
        stage_limit=1,
    )
    repeat_limited_summary = String(take!(repeat_limited_io))
    @test occursin("| 4 | su2_base_wigner_rows | 2 | 3.000000000 .. 4.000000000 |", repeat_limited_summary)
    @test !occursin("| 4 | setup | 2 | 0.250000000 .. 0.500000000 |", repeat_limited_summary)

    @test PauliTranslationProfileHarnessTests._check_base_su2_extend_rdm_profile_guard(
        target_only=false,
        direct_linear=true,
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=true,
        reflection_symmetry=true,
    ) === nothing
    @test PauliTranslationProfileHarnessTests._check_base_su2_extend_rdm_profile_guard(
        target_only=true,
        direct_linear=false,
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=true,
        reflection_symmetry=true,
    ) === nothing
    @test PauliTranslationProfileHarnessTests._check_base_su2_extend_rdm_profile_guard(
        target_only=false,
        direct_linear=true,
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=false,
        reflection_symmetry=false,
    ) === nothing
    @test_throws ArgumentError PauliTranslationProfileHarnessTests._check_base_su2_extend_rdm_profile_guard(
        target_only=false,
        direct_linear=false,
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=true,
        reflection_symmetry=true,
    )
    @test_throws ArgumentError PauliTranslationProfileHarnessTests._check_base_su2_extend_rdm_profile_guard(
        target_only=false,
        direct_linear=true,
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:closed,
        real_moment_matrix=true,
        reflection_symmetry=true,
    )
    @test_throws ArgumentError PauliTranslationProfileHarnessTests._check_base_su2_extend_rdm_profile_guard(
        target_only=true,
        direct_linear=false,
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=false,
        reflection_symmetry=true,
    )
    @test_throws ArgumentError PauliTranslationProfileHarnessTests._check_base_su2_extend_rdm_profile_guard(
        target_only=false,
        direct_linear=true,
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=false,
        reflection_symmetry=true,
    )
    _with_env((
        "NCTS_TRANSLATION_NS" => "8",
        "NCTS_TRANSLATION_DIRECT_LINEAR" => "true",
        "NCTS_TRANSLATION_SU2" => "true",
        "NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM" => "true",
        "NCTS_TRANSLATION_RDM_DECOMPOSITION" => "su2",
        "NCTS_TRANSLATION_RDM_SUPPORT" => "closed",
        "NCTS_TRANSLATION_WARMUP" => "false",
    )) do
        @test_throws ArgumentError PauliTranslationProfileHarnessTests.main()
    end

    @test PauliTranslationProfileHarnessTests._check_singlet_channel_profile_guard(
        singlet_channel_equalities=false,
        su2_symmetry=false,
    ) === nothing
    @test PauliTranslationProfileHarnessTests._check_singlet_channel_profile_guard(
        singlet_channel_equalities=true,
        su2_symmetry=true,
    ) === nothing
    @test_throws ArgumentError PauliTranslationProfileHarnessTests._check_singlet_channel_profile_guard(
        singlet_channel_equalities=true,
        su2_symmetry=false,
    )
    _with_env((
        "NCTS_TRANSLATION_NS" => "8",
        "NCTS_TRANSLATION_SINGLET_EQUALITIES" => "true",
        "NCTS_TRANSLATION_SU2" => "false",
        "NCTS_TRANSLATION_WARMUP" => "false",
    )) do
        @test_throws ArgumentError PauliTranslationProfileHarnessTests.main()
    end

    _with_env((
        "NCTS_PERF_ALLOW_LARGE" => nothing,
        "NCTS_PERF_MAX_N" => "13",
        "NCTS_PERF_NS" => "14",
        "NCTS_PERF_TARGET_ONLY" => "true",
    )) do
        @test PauliSparseChainD4HarnessTests._check_size_guard([13]) === nothing
        @test_throws ArgumentError PauliSparseChainD4HarnessTests._check_size_guard([14])
        @test_throws ArgumentError PauliSparseChainD4HarnessTests.main()
        @test PauliChargeSingletPrepHarnessTests._check_size_guard([13]) === nothing
        @test_throws ArgumentError PauliChargeSingletPrepHarnessTests._check_size_guard([14])
        @test_throws ArgumentError PauliChargeSingletPrepHarnessTests.main()
    end

    _with_env((
        "NCTS_PERF_ALLOW_LARGE" => "true",
        "NCTS_PERF_MAX_N" => "13",
        "NCTS_PERF_ESTIMATED_WALL_SECONDS" => nothing,
        "NCTS_PERF_ESTIMATED_RSS_GIB" => nothing,
        "NCTS_PERF_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test_throws ArgumentError PauliSparseChainD4HarnessTests._check_size_guard([14])
        @test_throws ArgumentError PauliChargeSingletPrepHarnessTests._check_size_guard([14])
    end

    _with_env((
        "NCTS_PERF_ALLOW_LARGE" => "true",
        "NCTS_PERF_MAX_N" => "13",
        "NCTS_PERF_ESTIMATED_WALL_SECONDS" => "1",
        "NCTS_PERF_ESTIMATED_RSS_GIB" => "0.001",
        "NCTS_PERF_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test PauliSparseChainD4HarnessTests._check_size_guard([14]) === nothing
        @test PauliChargeSingletPrepHarnessTests._check_size_guard([14]) === nothing
    end

    _with_env((
        "NCTS_COMPARE_ALLOW_LARGE" => nothing,
        "NCTS_COMPARE_MAX_N" => "13",
        "NCTS_COMPARE_NS" => "14",
    )) do
        @test PauliTranslationCompareHarnessTests._check_size_guard([13]) === nothing
        @test_throws ArgumentError PauliTranslationCompareHarnessTests._check_size_guard([14])
        @test_throws ArgumentError PauliTranslationCompareHarnessTests.main()
    end

    _with_env((
        "NCTS_COMPARE_ALLOW_LARGE" => "true",
        "NCTS_COMPARE_MAX_N" => "13",
        "NCTS_COMPARE_ESTIMATED_WALL_SECONDS" => nothing,
        "NCTS_COMPARE_ESTIMATED_RSS_GIB" => nothing,
        "NCTS_COMPARE_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test_throws ArgumentError PauliTranslationCompareHarnessTests._check_size_guard([14])
    end

    _with_env((
        "NCTS_COMPARE_ALLOW_LARGE" => "true",
        "NCTS_COMPARE_MAX_N" => "13",
        "NCTS_COMPARE_ESTIMATED_WALL_SECONDS" => "1",
        "NCTS_COMPARE_ESTIMATED_RSS_GIB" => "0.001",
        "NCTS_COMPARE_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test PauliTranslationCompareHarnessTests._check_size_guard([14]) === nothing
    end

    _with_env((
        "NCTS_SOLVER_PROBE_ALLOW_LARGE" => nothing,
        "NCTS_SOLVER_PROBE_MAX_N" => "13",
        "NCTS_SOLVER_PROBE_N" => "14",
    )) do
        @test PauliTranslationSolverProbeHarnessTests._check_size_guard(13) === nothing
        @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._check_size_guard(14)
        @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests.main()
    end

    _with_env((
        "NCTS_SOLVER_PROBE_ALLOW_LARGE" => "true",
        "NCTS_SOLVER_PROBE_MAX_N" => "13",
        "NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS" => nothing,
        "NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB" => nothing,
        "NCTS_SOLVER_PROBE_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._check_size_guard(14)
    end

    _with_env((
        "NCTS_SOLVER_PROBE_ALLOW_LARGE" => "true",
        "NCTS_SOLVER_PROBE_MAX_N" => "13",
        "NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS" => "1",
        "NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB" => "0.001",
        "NCTS_SOLVER_PROBE_MAX_LOAD_PER_CPU" => "1000",
    )) do
        @test PauliTranslationSolverProbeHarnessTests._check_size_guard(14) === nothing
    end

    _with_env((
        "NCTS_SOLVER_PROBE_ALLOW_DELETED_PRIMAL_SU2_SOLVE" => nothing,
        "NCTS_SOLVER_PROBE_DELETED_PRIMAL_SU2_SOLVE_MIN_N" => nothing,
    )) do
        @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._check_deleted_primal_su2_solved_probe_guard(
            n=14,
            solve=true,
            dualize=false,
            formulation=:moment_variables,
            su2_symmetry=true,
            contiguous_rdm_decomposition=:su2,
            linear_state_opt_width=7,
            psd_state_opt_width=3,
        )
        @test PauliTranslationSolverProbeHarnessTests._check_deleted_primal_su2_solved_probe_guard(
            n=12,
            solve=true,
            dualize=false,
            formulation=:moment_variables,
            su2_symmetry=true,
            contiguous_rdm_decomposition=:su2,
            linear_state_opt_width=7,
            psd_state_opt_width=3,
        ) === nothing
        @test PauliTranslationSolverProbeHarnessTests._check_deleted_primal_su2_solved_probe_guard(
            n=14,
            solve=false,
            dualize=false,
            formulation=:moment_variables,
            su2_symmetry=true,
            contiguous_rdm_decomposition=:su2,
            linear_state_opt_width=7,
            psd_state_opt_width=3,
        ) === nothing
        @test PauliTranslationSolverProbeHarnessTests._check_deleted_primal_su2_solved_probe_guard(
            n=14,
            solve=true,
            dualize=true,
            formulation=:moment_variables,
            su2_symmetry=true,
            contiguous_rdm_decomposition=:su2,
            linear_state_opt_width=7,
            psd_state_opt_width=3,
        ) === nothing
    end
    _with_env((
        "NCTS_SOLVER_PROBE_ALLOW_DELETED_PRIMAL_SU2_SOLVE" => "true",
        "NCTS_SOLVER_PROBE_DELETED_PRIMAL_SU2_SOLVE_MIN_N" => nothing,
    )) do
        @test PauliTranslationSolverProbeHarnessTests._check_deleted_primal_su2_solved_probe_guard(
            n=14,
            solve=true,
            dualize=false,
            formulation=:moment_variables,
            su2_symmetry=true,
            contiguous_rdm_decomposition=:su2,
            linear_state_opt_width=7,
            psd_state_opt_width=3,
        ) === nothing
    end

    @test PauliTranslationSolverProbeHarnessTests._check_base_su2_extend_rdm_solver_probe_guard(
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=true,
        reflection_symmetry=true,
    ) === nothing
    @test PauliTranslationSolverProbeHarnessTests._check_base_su2_extend_rdm_solver_probe_guard(
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=false,
        reflection_symmetry=false,
    ) === nothing
    @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._check_base_su2_extend_rdm_solver_probe_guard(
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:closed,
        real_moment_matrix=true,
        reflection_symmetry=true,
    )
    @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._check_base_su2_extend_rdm_solver_probe_guard(
        base_su2_extend_rdm=true,
        su2_symmetry=true,
        contiguous_rdm_decomposition=:su2,
        contiguous_rdm_support=:extend,
        real_moment_matrix=false,
        reflection_symmetry=true,
    )
    _with_env((
        "NCTS_SOLVER_PROBE_N" => "8",
        "NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM" => "true",
        "NCTS_SOLVER_PROBE_SU2" => "true",
        "NCTS_SOLVER_PROBE_RDM_DECOMPOSITION" => "su2",
        "NCTS_SOLVER_PROBE_RDM_SUPPORT" => "closed",
        "NCTS_SOLVER_PROBE_LOWER_MODEL" => "false",
    )) do
        @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests.main()
    end

    _with_env((
        "NCTS_SOLVER_PROBE_LOWER_MODEL" => "false",
        "NCTS_SOLVER_PROBE_LOWER" => nothing,
    )) do
        @test !PauliTranslationSolverProbeHarnessTests._parse_lower_model_env()
    end

    _with_env((
        "NCTS_SOLVER_PROBE_LOWER_MODEL" => nothing,
        "NCTS_SOLVER_PROBE_LOWER" => "false",
    )) do
        @test !PauliTranslationSolverProbeHarnessTests._parse_lower_model_env()
    end

    _with_env((
        "NCTS_SOLVER_PROBE_LOWER_MODEL" => "true",
        "NCTS_SOLVER_PROBE_LOWER" => "false",
    )) do
        @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._parse_lower_model_env()
    end

    psd_diagnostics = [
        (native_hermitian_residual=2.0e-9, native_min_eigenvalue=-1.0e-8),
        (native_hermitian_residual=4.0e-10, native_min_eigenvalue=3.0e-7),
    ]
    zero_diagnostics = [(abs_value=0.25,), (abs_value=0.5,)]
    @test PauliTranslationSolverProbeHarnessTests._max_native_hermitian_residual(
        psd_diagnostics,
    ) == 2.0e-9
    @test PauliTranslationSolverProbeHarnessTests._min_native_eigenvalue(
        psd_diagnostics,
    ) == -1.0e-8
    @test PauliTranslationSolverProbeHarnessTests._max_zero_dual_abs_value(
        zero_diagnostics,
    ) == 0.5
    telemetry_model = PauliTranslationSolverProbeHarnessTests.JuMP.Model(
        PauliTranslationSolverProbeHarnessTests.Clarabel.Optimizer,
    )
    PauliTranslationSolverProbeHarnessTests.JuMP.set_silent(telemetry_model)
    PauliTranslationSolverProbeHarnessTests.JuMP.@variable(telemetry_model, x <= 1)
    PauliTranslationSolverProbeHarnessTests.JuMP.@objective(telemetry_model, Max, x)
    PauliTranslationSolverProbeHarnessTests.JuMP.optimize!(telemetry_model)
    telemetry =
        PauliTranslationSolverProbeHarnessTests._solver_probe_solve_result(
            telemetry_model,
        )
    @test telemetry.status == PauliTranslationSolverProbeHarnessTests.JuMP.MOI.OPTIMAL
    @test telemetry.primal_status ==
        PauliTranslationSolverProbeHarnessTests.JuMP.MOI.FEASIBLE_POINT
    @test telemetry.result_count == 1
    @test telemetry.raw_status isa String
    @test telemetry.objective ≈ 1.0 atol = 1e-8
    @test PauliTranslationSolverProbeHarnessTests._probe_execution_state(
        nothing,
        nothing,
        nothing,
    ) == (construction_only=true, model_built=false, solved=false)
    @test PauliTranslationSolverProbeHarnessTests._probe_execution_state(
        (time=1.0,),
        nothing,
        nothing,
    ) == (construction_only=false, model_built=true, solved=false)
    @test PauliTranslationSolverProbeHarnessTests._probe_execution_state(
        nothing,
        (time=1.0,),
        nothing,
    ) == (construction_only=false, model_built=true, solved=false)
    @test PauliTranslationSolverProbeHarnessTests._probe_execution_state(
        nothing,
        (time=1.0,),
        (status=:OPTIMAL, objective=-1.0),
    ) == (construction_only=false, model_built=true, solved=true)
    real_primal_estimate =
        PauliTranslationSolverProbeHarnessTests._solver_probe_model_size_estimate(
            block_sizes=[3, 4],
            block_cones=[:PSD, :PSD],
            n_moments=10,
            n_free_keys=2,
            n_zero_constraints=5,
            formulation=:psd_blocks,
            representation=:real,
            orphan_policy=:free_variables,
            dualize=false,
            sos_hermitian_representation=:real_lift,
        )
    @test real_primal_estimate.model_mode == "primal"
    @test real_primal_estimate.model_variables == 18
    @test real_primal_estimate.psd_cone_scalar_variables == 16
    @test real_primal_estimate.scalar_equalities_upper_bound == 22
    @test real_primal_estimate.dense_schur_bytes == 8 * 22^2
    @test PauliTranslationSolverProbeHarnessTests._solver_probe_model_size_gate_status(
        real_primal_estimate;
        mem_available_bytes=10_000,
        max_rss_fraction=0.5,
    ) == "ok"
    @test PauliTranslationSolverProbeHarnessTests._solver_probe_model_size_gate_status(
        real_primal_estimate;
        mem_available_bytes=3_000,
        max_rss_fraction=1.0,
    ) == "blocked_insufficient_memory"
    @test PauliTranslationSolverProbeHarnessTests._solver_probe_model_size_gate_status(
        real_primal_estimate;
        mem_available_bytes=4_000,
        max_rss_fraction=0.8,
    ) == "blocked_insufficient_memory_headroom"
    real_primal_gate =
        PauliTranslationSolverProbeHarnessTests._solver_probe_model_size_gate(
            real_primal_estimate;
            mem_available_bytes=4_000,
            max_rss_fraction=0.8,
        )
    @test real_primal_gate.status == "blocked_insufficient_memory_headroom"
    @test real_primal_gate.estimated_rss_bytes == real_primal_estimate.dense_schur_bytes
    @test real_primal_gate.mem_available_bytes == 4_000
    @test real_primal_gate.max_rss_fraction == 0.8
    @test occursin(
        "would consume more than 80.0% of available memory",
        real_primal_gate.reason,
    )
    @test real_primal_estimate.zero_dual_variables == 0

    hermitian_lift_estimate =
        PauliTranslationSolverProbeHarnessTests._solver_probe_model_size_estimate(
            block_sizes=[3],
            block_cones=[:HPSD],
            n_moments=10,
            n_free_keys=0,
            n_zero_constraints=5,
            formulation=:psd_blocks,
            representation=:complex,
            orphan_policy=:error,
            dualize=true,
            sos_hermitian_representation=:real_lift,
        )
    @test hermitian_lift_estimate.model_mode == "sos_dual"
    @test hermitian_lift_estimate.model_variables == 26
    @test hermitian_lift_estimate.psd_cone_scalar_variables == 21
    @test hermitian_lift_estimate.scalar_equalities_upper_bound == 19
    @test hermitian_lift_estimate.dense_schur_bytes == 8 * 19^2
    @test hermitian_lift_estimate.zero_dual_variables == 5

    native_hermitian_estimate =
        PauliTranslationSolverProbeHarnessTests._solver_probe_model_size_estimate(
            block_sizes=[3],
            block_cones=[:HPSD],
            n_moments=10,
            n_free_keys=0,
            n_zero_constraints=5,
            formulation=:psd_blocks,
            representation=:complex,
            orphan_policy=:error,
            dualize=true,
            sos_hermitian_representation=:native,
        )
    @test native_hermitian_estimate.model_variables == 14
    @test native_hermitian_estimate.psd_cone_scalar_variables == 9
    @test native_hermitian_estimate.dense_schur_bytes == 8 * 19^2
    @test PauliTranslationSolverProbeHarnessTests._normalize_solver_probe_coefficient_scaling(
        :none,
    ) == :none
    @test PauliTranslationSolverProbeHarnessTests._normalize_solver_probe_coefficient_scaling(
        :max_abs,
    ) == :max_abs
    @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._normalize_solver_probe_coefficient_scaling(
        :mystery,
    )
    @test PauliTranslationSolverProbeHarnessTests._normalize_solver_probe_coefficient_scaling_floor(
        0.1,
    ) == 0.1
    @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._normalize_solver_probe_coefficient_scaling_floor(
        -0.1,
    )
    @test PauliTranslationSolverProbeHarnessTests._construction_stage_timing_fields(
        Dict(
            :su2_base_wigner_rows => 3.0,
            :su2_wigner_spin_form_build => 2.0,
        ),
    ) == [
        "construction_stage_su2_base_wigner_rows_seconds" => 3.0,
        "construction_stage_su2_wigner_spin_form_build_seconds" => 2.0,
    ]
    @test PauliTranslationSolverProbeHarnessTests._mosek_library_version(
        FakeMosekVersionHarnessTests,
    ) == "11.2.0"
    @test PauliTranslationSolverProbeHarnessTests._max_native_hermitian_residual([]) ==
        0.0
    @test PauliTranslationSolverProbeHarnessTests._min_native_eigenvalue([]) === nothing
    @test PauliTranslationSolverProbeHarnessTests._max_zero_dual_abs_value([]) == 0.0
    _with_env((
        "NCTS_SOLVER_PROBE_QMBCERTIFY_PROBE_ID" => "A1_L20",
        "NCTS_SOLVER_PROBE_QMBCERTIFY_OBJECTIVE_TOTAL_ESTIMATE" => "-8.904954518917137",
    )) do
        reference = PauliTranslationSolverProbeHarnessTests._solver_probe_qmbcertify_reference()
        @test reference.probe_id == "A1_L20"
        @test reference.objective_total_estimate ≈ -8.904954518917137
    end
    _with_env((
        "NCTS_SOLVER_PROBE_QMBCERTIFY_PROBE_ID" => "A1_L20",
        "NCTS_SOLVER_PROBE_QMBCERTIFY_OBJECTIVE_TOTAL_ESTIMATE" => nothing,
    )) do
        @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._solver_probe_qmbcertify_reference()
    end
    qmb_ref_missing = (probe_id=nothing, objective_total_estimate=nothing)
    qmb_ref_present = (probe_id="A1_L20", objective_total_estimate=-8.904954518917137)
    solved_probe = (status=:OPTIMAL, objective=-8.904954795042105)
    @test PauliTranslationSolverProbeHarnessTests._check_source_base_solved_probe_qmbcertify_reference(
        false,
        solved_probe,
        qmb_ref_missing,
    ) === nothing
    @test PauliTranslationSolverProbeHarnessTests._check_source_base_solved_probe_qmbcertify_reference(
        true,
        nothing,
        qmb_ref_missing,
    ) === nothing
    @test_throws ArgumentError PauliTranslationSolverProbeHarnessTests._check_source_base_solved_probe_qmbcertify_reference(
        true,
        solved_probe,
        qmb_ref_missing,
    )
    @test PauliTranslationSolverProbeHarnessTests._check_source_base_solved_probe_qmbcertify_reference(
        true,
        solved_probe,
        qmb_ref_present,
    ) === nothing
end

@testset "Perf Harness Entry Guards" begin
    perf_scripts = (
        "heisenberg_mosek_scaling.jl",
        "pauli_charge_singlet_prep.jl",
        "pauli_sparse_chain_d4_blocks.jl",
        "pauli_translation_compare.jl",
        "pauli_translation_profile.jl",
        "pauli_translation_solver_probe.jl",
        "qmbcertify_formulation_audit.jl",
        "qmbcertify_profile_probe.jl",
        "qmbcertify_reference_runs.jl",
    )
    for script in perf_scripts
        source = read(normpath(joinpath(@__DIR__, "..", "perf", script)), String)
        @test occursin("if abspath(PROGRAM_FILE) == @__FILE__", source)
    end

    profile_probe_source = read(
        normpath(joinpath(@__DIR__, "..", "perf", "qmbcertify_profile_probe.jl")),
        String,
    )
    normalized_profile_probe_source = replace(profile_probe_source, r"\n#\s*" => " ")
    @test occursin("status_policy_profile=profile", profile_probe_source)
    @test occursin("NCTS_QMB_ALLOW_LARGE pressure gate", normalized_profile_probe_source)
    @test occursin("explicit wall/RSS estimates", normalized_profile_probe_source)
    @test occursin("safe load/memory telemetry", normalized_profile_probe_source)

    translation_profile_source = read(
        normpath(joinpath(@__DIR__, "..", "perf", "pauli_translation_profile.jl")),
        String,
    )
    @test !occursin(
        "NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true` requires `NCTS_TRANSLATION_TARGET_ONLY=false`",
        translation_profile_source,
    )
    @test occursin(
        "qmbcertify_base_construct ? :qmbcertify : :full",
        translation_profile_source,
    )
    @test occursin(
        "qmbcertify_base_construct ? :extend : :closed",
        translation_profile_source,
    )
    @test occursin("model-size gate status", translation_profile_source)
    @test occursin("estimated_model_size_gate_status", translation_profile_source)
    @test occursin("SOS-dual equality upper bound", translation_profile_source)
    @test occursin("SOS-dual dense-Schur bytes", translation_profile_source)
    @test occursin("SOS-dual model-size gate status", translation_profile_source)
    @test occursin("SOS-dual model-size gate reason", translation_profile_source)
    @test occursin("SOS-dual model-size gate estimated RSS bytes", translation_profile_source)
    @test occursin("SOS-dual model-size gate MemAvailable bytes", translation_profile_source)
    @test occursin("SOS-dual model-size gate max RSS fraction", translation_profile_source)

    formulation_audit_source = read(
        normpath(joinpath(@__DIR__, "..", "perf", "qmbcertify_formulation_audit.jl")),
        String,
    )
    @test occursin("_try_nctssos_qmb_source_base_formulation", formulation_audit_source)
    @test occursin("NCTSSoS.simplify(PauliAlgebra", formulation_audit_source)
    @test occursin("nctssos_source_base_linear_moments", formulation_audit_source)
    @test occursin("nctssos_source_base_lso_row_count", formulation_audit_source)
    @test occursin("nctssos_source_base_rdm_blocks", formulation_audit_source)
    @test occursin("nctssos_source_base_block_histogram", formulation_audit_source)
    @test occursin("_nctssos_qmb_objective_metrics", formulation_audit_source)
    @test occursin("nctssos_qmb_objective_total_match", formulation_audit_source)
    @test occursin("nctssos_qmb_objective_max_abs_delta", formulation_audit_source)
    @test occursin("mode=:qmbcertify", formulation_audit_source)
    @test occursin("nctssos_qmbcertify_lso_distinct_row_count", formulation_audit_source)
    @test occursin("qmbcertify_lso_row_qmb_only_count", formulation_audit_source)
    @test occursin("qmbcertify_lso_row_nctssos_only_count", formulation_audit_source)
    @test occursin("NCTS_FORMULATION_AUDIT_DIRECT_COMPARE", formulation_audit_source)
    @test occursin("nctssos_direct_compare", formulation_audit_source)
    @test occursin("nctssos_linear_moments=skipped", formulation_audit_source)
    @test occursin("_qmbcertify_rdm_source_entry_metrics", formulation_audit_source)
    @test occursin("NCTS_FORMULATION_AUDIT_RDM_ENTRY_MAX_K", formulation_audit_source)
    @test occursin("qmb_rdm_source_row_block_exact_match", formulation_audit_source)
    @test occursin("qmb_rdm_source_entry_match", formulation_audit_source)
    @test occursin("qmb_rdm_source_entry_max_abs_delta", formulation_audit_source)
    @test occursin("_qmbcertify_base_source_entry_metrics", formulation_audit_source)
    @test occursin("_qmb_source_form_from_pairs!", formulation_audit_source)
    @test occursin("abs(coef) <= atol", formulation_audit_source)
    @test occursin("nctssos_source_base_entry_block_label_match", formulation_audit_source)
    @test occursin("nctssos_source_base_entry_compared_count", formulation_audit_source)
    @test occursin("nctssos_source_base_entry_match", formulation_audit_source)
    @test occursin("nctssos_source_base_entry_max_abs_delta", formulation_audit_source)

    pauli_chains_source = read(
        normpath(joinpath(@__DIR__, "..", "src", "optimization", "pauli_chains.jl")),
        String,
    )
    @test occursin("estimated_model_size_gate_status", pauli_chains_source)
    @test occursin("blocked_missing_scalar_equality_estimate", pauli_chains_source)
    @test occursin("estimated_model_size_gate_reason", pauli_chains_source)
    @test occursin("estimated_sos_dual_scalar_equalities_upper_bound", pauli_chains_source)
    @test occursin("estimated_sos_dual_dense_schur_bytes", pauli_chains_source)
    @test occursin("_qmbcertify_linear_state_opt_row_key", pauli_chains_source)
    @test occursin("_qmbcertify_realified_phase_number", pauli_chains_source)
    @test occursin("_qmbcertify_chain_support_term", pauli_chains_source)
    @test occursin("realify=true", pauli_chains_source)
    @test occursin(
        "all(mono -> mono in moment_set, monomials(reduced)) || continue\n" *
        "        row_key = _qmbcertify_linear_state_opt_row_key(reduced)",
        pauli_chains_source,
    )
    @test occursin(
        "covered = all(mono -> mono in moment_set, monomials(reduced))\n" *
        "        if !covered\n" *
        "            mode == :qmbcertify && continue",
        pauli_chains_source,
    )
    @test occursin(
        "cross_scale = xor(isone(row), isone(col)) ? C(inv(sqrt(R(n)))) : one(C)",
        pauli_chains_source,
    )
    @test occursin(
        "_qmbcertify_real_fixed_linear_block(complex_entries, LC; atol)",
        pauli_chains_source,
    )
end

@testset "Solver Probe Certificate Diagnostics Fixture Fields" begin
    source = read(
        normpath(joinpath(@__DIR__, "..", "perf", "pauli_translation_solver_probe.jl")),
        String,
    )
    @test occursin("sos_dual_certificate_diagnostics", source)
    @test occursin("sos_certificate_diagnostics_source", source)
    @test occursin("sos_certificate_psd_block_count", source)
    @test occursin("sos_certificate_zero_dual_count", source)
    @test occursin("sos_certificate_max_native_hermitian_residual", source)
    @test occursin("sos_certificate_min_native_eigenvalue", source)
    @test occursin("sos_certificate_max_zero_dual_abs_value", source)
    @test occursin("construction_only", source)
    @test occursin("model_built", source)
    @test occursin("solved", source)
    @test occursin("estimated_model_variables", source)
    @test occursin("estimated_risk_gate_status", source)
    @test occursin("estimated_model_size_gate_status", source)
    @test occursin("estimated_model_size_gate_estimated_rss_bytes", source)
    @test occursin("estimated_model_size_gate_mem_available_bytes", source)
    @test occursin("estimated_model_size_gate_max_rss_fraction", source)
    @test occursin("estimated_model_size_gate_reason", source)
    @test occursin(r"_print_toml_field\(\s*\"formulation\"", source)
    @test occursin(r"_print_toml_field\(\s*\"representation\"", source)
    @test occursin("blocked_lowering_orphan_policy", source)
    @test occursin("estimated_scalar_equalities_upper_bound", source)
    @test occursin("estimated_dense_schur_bytes", source)
    @test occursin("estimated_psd_cone_scalar_variables", source)
    @test occursin("NCTS_SOLVER_PROBE_QMBCERTIFY_PROBE_ID", source)
    @test occursin("NCTS_SOLVER_PROBE_QMBCERTIFY_OBJECTIVE_TOTAL_ESTIMATE", source)
    @test occursin("NCTS_SOLVER_PROBE_MOSEK_TOL_PFEAS", source)
    @test occursin("NCTS_SOLVER_PROBE_MOSEK_TOL_DFEAS", source)
    @test occursin("NCTS_SOLVER_PROBE_MOSEK_TOL_REL_GAP", source)
    @test occursin("_check_deleted_primal_su2_solved_probe_guard", source)
    @test occursin("NCTS_SOLVER_PROBE_ALLOW_DELETED_PRIMAL_SU2_SOLVE", source)
    @test occursin("NCTS_SOLVER_PROBE_DELETED_PRIMAL_SU2_SOLVE_MIN_N", source)
    @test occursin("The guarded N=12 solve timed out", source)
    @test occursin("MSK_DPAR_INTPNT_CO_TOL_PFEAS", source)
    @test occursin("MSK_DPAR_INTPNT_CO_TOL_DFEAS", source)
    @test occursin("MSK_DPAR_INTPNT_CO_TOL_REL_GAP", source)
    @test occursin("qmbcertify_probe_id", source)
    @test occursin("qmbcertify_objective_total_estimate", source)
    @test occursin("objective_minus_qmbcertify_total_estimate", source)
    @test occursin(r"_print_toml_field\(\s*\"su2_base_zero_row_count\"", source)
    @test occursin(r"_print_toml_field\(\s*\"su2_base_spin_offblock_row_count\"", source)
    @test occursin(r"_print_toml_field\(\s*\"su2_base_magnetic_offdiag_row_count\"", source)
    @test occursin(r"_print_toml_field\(\s*\"su2_base_magnetic_copy_row_count\"", source)
end

@testset "QMBCertify Parity Plan Run Gates" begin
    plan_source = read(
        normpath(joinpath(@__DIR__, "..", "plan", "qmbcertify_parity.md")),
        String,
    )
    @test occursin("If `estimated_risk_gate_status` is", plan_source)
    @test occursin("If `estimated_model_size_gate_status` is", plan_source)
    @test occursin("blocked_missing_scalar_equality_estimate", plan_source)
    normalized_plan_source = replace(plan_source, r"\s+" => " ")
    @test occursin(
        "These entries are deleted from the active run plan until fresh evidence proves both scientific necessity and remote safety.",
        normalized_plan_source,
    )
    @test occursin(
        "These templates are historical command recipes, not a run queue.",
        normalized_plan_source,
    )
    @test occursin(
        "The large-run guard now prints independent load, estimate, and memory sub-statuses",
        normalized_plan_source,
    )
    @test occursin(
        "The estimate sub-status is now split further into wall-time and RSS estimate statuses",
        normalized_plan_source,
    )
    @test occursin(
        "The guard now prints per-substatus reason fields",
        normalized_plan_source,
    )
    @test occursin(
        "Future source-base solve probes must set `NCTS_SOLVER_PROBE_QMBCERTIFY_PROBE_ID`",
        normalized_plan_source,
    )
    @test occursin(
        "The solver probe harness now enforces that deletion",
        normalized_plan_source,
    )
    @test !occursin(
        "Deferred larger probes, to re-add only after",
        normalized_plan_source,
    )
    @test !occursin(
        "base half-basis/block mismatch before treating any larger solve",
        normalized_plan_source,
    )
    @test !occursin(
        "remaining dominant formulation mismatch is the base half-basis/block construction",
        normalized_plan_source,
    )
    @test !occursin(
        "finite-axis NCTSSoS base still has 18 solver-facing blocks of size 120",
        normalized_plan_source,
    )
    @test !occursin(
        "A2_L30 solve is removed from the active run queue until the base formulation mismatch is resolved",
        normalized_plan_source,
    )
end

@testset "Perf Harness Large-Run Pressure Guards" begin
    repo_root = normpath(joinpath(@__DIR__, ".."))
    guard_path = joinpath(repo_root, "perf", "shared_load_guard.jl")
    @test isfile(guard_path)

    if isfile(guard_path)
        guard_source = read(guard_path, String)
        @test occursin("_ncts_check_large_run_pressure_guard", guard_source)
        @test occursin("blocked_overloaded_remote", guard_source)
        @test occursin("_ncts_loadavg_triplet", guard_source)
        @test occursin("large_run_pressure_gate_load5", guard_source)
        @test occursin("large_run_pressure_gate_load15", guard_source)
        @test occursin("large_run_pressure_gate_reason", guard_source)
        @test occursin("blocked_missing_loadavg", guard_source)
        @test occursin("_ncts_mem_available_bytes", guard_source)
        @test occursin("blocked_insufficient_memory", guard_source)
        @test occursin("blocked_insufficient_memory_headroom", guard_source)
        @test occursin("blocked_missing_wall_estimate", guard_source)
        @test occursin("blocked_excessive_wall_estimate", guard_source)
        @test occursin("blocked_missing_rss_estimate", guard_source)
        @test occursin("_ncts_large_run_estimate_status", guard_source)
        @test occursin("_ncts_large_run_pressure_failure_detail", guard_source)
        @test occursin("large_run_pressure_gate_estimated_wall_seconds", guard_source)
        @test occursin("large_run_pressure_gate_max_wall_seconds", guard_source)
        @test occursin("large_run_pressure_gate_max_rss_fraction", guard_source)
        @test occursin("NCTS_LOAD_GUARD_ESTIMATED_WALL_SECONDS", guard_source)
        @test occursin("NCTS_LOAD_GUARD_MAX_WALL_SECONDS", guard_source)
        @test occursin("NCTS_LOAD_GUARD_MAX_RSS_FRACTION", guard_source)
        @test occursin("NCTS_LOAD_GUARD_ESTIMATED_RSS_GIB", guard_source)
        @test occursin("_ncts_cpu_max_quota_count", guard_source)
        @test occursin("_ncts_effective_cpu_count", guard_source)
        @test occursin("cpu.max", guard_source)
        @test occursin("override_overcommitted_remote", guard_source)
        @test occursin("large_run_pressure_gate_override", guard_source)
        @test occursin("large_run_pressure_gate_blocked_status", guard_source)
        @test occursin("large_run_pressure_gate_load_status", guard_source)
        @test occursin("large_run_pressure_gate_estimate_status", guard_source)
        @test occursin("large_run_pressure_gate_wall_estimate_status", guard_source)
        @test occursin("large_run_pressure_gate_rss_estimate_status", guard_source)
        @test occursin("large_run_pressure_gate_memory_status", guard_source)
        @test occursin("large_run_pressure_gate_load_reason", guard_source)
        @test occursin("large_run_pressure_gate_estimate_reason", guard_source)
        @test occursin("large_run_pressure_gate_wall_estimate_reason", guard_source)
        @test occursin("large_run_pressure_gate_rss_estimate_reason", guard_source)
        @test occursin("large_run_pressure_gate_memory_reason", guard_source)
        @test occursin("NCTS_LOAD_GUARD_ALLOW_OVERCOMMITTED", guard_source)
    end

    guard_module = Module(:LargeRunPressureGuardReasonProbe)
    Base.include(guard_module, guard_path)
    output_path = tempname()
    open(output_path, "w") do io
        redirect_stdout(io) do
            guard_module._ncts_print_large_run_pressure_gate(
                status="override_overcommitted_remote",
                reason="override accepted despite blocked_overloaded_remote",
                blocked_status="blocked_overloaded_remote",
                load_status="blocked_overloaded_remote",
                load_reason="load averages 42.0, 41.0, 40.0 exceed 25.0 * 1.0",
                estimate_status="blocked_missing_wall_estimate",
                estimate_reason="no wall-time estimate was provided",
                wall_estimate_status="blocked_missing_wall_estimate",
                wall_estimate_reason="no wall-time estimate was provided",
                rss_estimate_status="blocked_missing_rss_estimate",
                rss_estimate_reason="no RSS estimate was provided",
                memory_status="ok",
                memory_reason="",
                label="synthetic large-run reason probe",
                load1=42.0,
                load5=41.0,
                load15=40.0,
                cpu_count=25.0,
                max_load_per_cpu=1.0,
                estimated_wall_seconds=nothing,
                max_wall_seconds=1800.0,
                mem_available_bytes=700 * 1024^3,
                estimated_rss_bytes=nothing,
                max_rss_fraction=0.8,
                override=true,
            )
        end
    end
    printed_guard = read(output_path, String)
    rm(output_path; force=true)
    @test occursin("large_run_pressure_gate_load_reason=load averages", printed_guard)
    @test occursin(
        "large_run_pressure_gate_estimate_reason=no wall-time estimate was provided",
        printed_guard,
    )
    @test occursin(
        "large_run_pressure_gate_wall_estimate_reason=no wall-time estimate was provided",
        printed_guard,
    )
    @test occursin(
        "large_run_pressure_gate_rss_estimate_reason=no RSS estimate was provided",
        printed_guard,
    )
    @test occursin("large_run_pressure_gate_memory_reason=", printed_guard)

    harnesses = [
        (
            "pauli_translation_solver_probe.jl",
            "NCTS_SOLVER_PROBE",
            "_solver_probe_preimport_large_run_pressure_guard",
            "using NCTSSoS",
        ),
        (
            "pauli_translation_profile.jl",
            "NCTS_TRANSLATION",
            "_translation_profile_preimport_large_run_pressure_guard",
            "using NCTSSoS",
        ),
        (
            "pauli_translation_compare.jl",
            "NCTS_COMPARE",
            "_translation_compare_preimport_large_run_pressure_guard",
            "using NCTSSoS",
        ),
        (
            "heisenberg_mosek_scaling.jl",
            "NCTS_PERF",
            "_heisenberg_scaling_preimport_large_run_pressure_guard",
            "using JuMP",
        ),
        (
            "pauli_charge_singlet_prep.jl",
            "NCTS_PERF",
            "_charge_singlet_preimport_large_run_pressure_guard",
            "using JuMP",
        ),
        (
            "pauli_sparse_chain_d4_blocks.jl",
            "NCTS_PERF",
            "_sparse_chain_preimport_large_run_pressure_guard",
            "using JuMP",
        ),
        (
            "qmbcertify_reference_runs.jl",
            "NCTS_QMB",
            "_qmb_reference_preimport_large_run_pressure_guard",
            "using Pkg",
        ),
        (
            "qmbcertify_formulation_audit.jl",
            "NCTS_QMB",
            "_formulation_audit_preimport_large_run_pressure_guard",
            "using NCTSSoS",
        ),
    ]

    for (filename, env_prefix, preimport_function, first_heavy_import) in harnesses
        source = read(joinpath(repo_root, "perf", filename), String)
        normalized_source = replace(source, r"\n#\s*" => " ")
        @test occursin("shared_load_guard.jl", source)
        @test occursin("_ncts_check_large_run_pressure_guard", source)
        @test occursin("env_prefix=\"$env_prefix\"", source)
        @test occursin(preimport_function, source)
        @test occursin("explicit wall/RSS", normalized_source) ||
            occursin("ESTIMATED_WALL_SECONDS", normalized_source)
        @test occursin("safe load/memory telemetry", normalized_source) ||
            occursin("pressure gate", normalized_source)
        shared_idx = findfirst("shared_load_guard.jl", source)
        import_idx = findfirst(first_heavy_import, source)
        @test shared_idx !== nothing
        @test import_idx !== nothing
        @test first(shared_idx) < first(import_idx)
    end
end

@testset "Pauli Perf Harness Axis Diagnostics" begin
    profile_source = read(
        normpath(joinpath(@__DIR__, "..", "perf", "pauli_translation_profile.jl")),
        String,
    )
    compare_source = read(
        normpath(joinpath(@__DIR__, "..", "perf", "pauli_translation_compare.jl")),
        String,
    )

    @test occursin("axis-orbit reduction ratio", profile_source)
    @test occursin("axis-orbit reduction ratio", compare_source)
    @test !occursin("missing axis-reduction ratio", profile_source)
    @test !occursin("axis-reduction ratio still missing", compare_source)
end

@testset "QMBCertify Reference Harness Provenance" begin
    @test QMBReferenceHarnessTests._reviewed_reference_commit("abc123"; allow_missing=false) == "abc123"
    @test QMBReferenceHarnessTests._reviewed_reference_commit(nothing; allow_missing=true) === nothing
    @test QMBReferenceHarnessTests._reviewed_solve_status_policy() ==
        ("OPTIMAL", "ALMOST_OPTIMAL")
    @test QMBReferenceHarnessTests._reviewed_solve_status_policy("A0") ==
        ("OPTIMAL", "ALMOST_OPTIMAL")
    @test QMBReferenceHarnessTests._reviewed_solve_status_policy("A2") ==
        ("OPTIMAL", "ALMOST_OPTIMAL")
    _with_env(("NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID" => nothing,)) do
        @test QMBReferenceHarnessTests._check_deleted_reference_case_restore_gate(
            "A1",
            [20, 30],
        ) === nothing
        @test_throws ArgumentError QMBReferenceHarnessTests._check_deleted_reference_case_restore_gate(
            "A0",
            [20],
        )
        @test_throws ArgumentError QMBReferenceHarnessTests._check_deleted_reference_case_restore_gate(
            "A2",
            [20],
        )
    end
    _with_env(("NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID" => "missing_probe",)) do
        @test_throws ArgumentError QMBReferenceHarnessTests._check_deleted_reference_case_restore_gate(
            "A2",
            [20],
        )
    end
    _with_env(("NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID" => "A1_L12_profile_probe_L12_2026_06_30",)) do
        @test_throws ArgumentError QMBReferenceHarnessTests._check_deleted_reference_case_restore_gate(
            "A2",
            [20],
        )
    end
    _with_env(("NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID" => "A2_L12_profile_probe_L12_2026_06_30",)) do
        @test_throws ArgumentError QMBReferenceHarnessTests._check_deleted_reference_case_restore_gate(
            "A2",
            [20, 30],
        )
    end
    a2_fixture = QMBReferenceHarnessTests._profile_fixture("A2")
    @test_throws ArgumentError QMBReferenceHarnessTests._check_accepted_profile_variant_restore_evidence(
        a2_fixture,
        "A2",
        "A2_L12_profile_probe_L12_2026_06_30",
    )
    a2_variant_fixture = deepcopy(a2_fixture)
    a2_variant_probe = deepcopy(QMBReferenceHarnessTests._profile_probe_by_id(
        a2_fixture,
        "A2_L12_profile_probe_L12_2026_06_30",
    ))
    a2_variant_probe["id"] = "A2_L12_pso2_variant_probe"
    a2_variant_probe["overrides_summary"] = "pso=2"
    a2_variant_probe["estimated_wall_seconds"] = 120
    a2_variant_probe["estimated_rss_gib"] = 2.5
    push!(a2_variant_fixture["profile_probes"], a2_variant_probe)
    @test QMBReferenceHarnessTests._check_accepted_profile_variant_restore_evidence(
        a2_variant_fixture,
        "A2",
        "A2_L12_pso2_variant_probe",
    ) === nothing
    a2_variant_missing_rss_fixture = deepcopy(a2_variant_fixture)
    delete!(last(a2_variant_missing_rss_fixture["profile_probes"]), "peak_rss_bytes")
    @test_throws ArgumentError QMBReferenceHarnessTests._check_accepted_profile_variant_restore_evidence(
        a2_variant_missing_rss_fixture,
        "A2",
        "A2_L12_pso2_variant_probe",
    )
    a2_variant_missing_estimate_fixture = deepcopy(a2_variant_fixture)
    delete!(
        last(a2_variant_missing_estimate_fixture["profile_probes"]),
        "estimated_wall_seconds",
    )
    @test_throws ArgumentError QMBReferenceHarnessTests._check_accepted_profile_variant_restore_evidence(
        a2_variant_missing_estimate_fixture,
        "A2",
        "A2_L12_pso2_variant_probe",
    )
    @test_throws QMBReferenceHarnessTests.NonOptimalReferenceSolve QMBReferenceHarnessTests._check_reviewed_solve_status!(
        (termination_status=missing, solution_status=missing),
        -0.4515446349992884,
    )
    @test QMBReferenceHarnessTests._check_reviewed_solve_status!(
        (termination_status="OPTIMAL", solution_status=missing),
        -0.4515446349992884,
    ) === nothing
    successful_log = """
    SDP size: n = 114, m = 1044
    SDP solving time: 5.550101014 seconds.
    optimum = -0.4489492424398082
    """
    successful_timings = QMBReferenceHarnessTests._qmb_log_timings(successful_log)
    @test successful_timings.termination_status == "OPTIMAL"
    @test ismissing(successful_timings.solution_status)
    @test_throws QMBReferenceHarnessTests.NonOptimalReferenceSolve QMBReferenceHarnessTests._check_reviewed_solve_status!(
        (
            termination_status="DUAL_INFEASIBLE",
            solution_status="INFEASIBILITY_CERTIFICATE",
        ),
        0.4633637063586052,
    )
    @test_throws QMBReferenceHarnessTests.NonOptimalReferenceSolve QMBReferenceHarnessTests._check_reviewed_solve_status!(
        (
            termination_status="SLOW_PROGRESS",
            solution_status="FEASIBLE_POINT",
        ),
        -0.4452212294966909,
    )

    err = try
        QMBReferenceHarnessTests._reviewed_reference_commit(nothing; allow_missing=false)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("NCTS_NCTSSOS_COMMIT", sprint(showerror, err))
    @test_throws ArgumentError QMBReferenceHarnessTests._check_size_guard([20])
end

@testset "QMBCertify Reference Harness Include Side Effects" begin
    old_ns = get(ENV, "NCTS_QMB_NS", nothing)
    old_allow_large = get(ENV, "NCTS_QMB_ALLOW_LARGE", nothing)
    try
        ENV["NCTS_QMB_NS"] = "20"
        delete!(ENV, "NCTS_QMB_ALLOW_LARGE")
        harness = normpath(joinpath(@__DIR__, "..", "perf", "qmbcertify_reference_runs.jl"))
        mod = Module(:QMBReferenceHarnessLargeEnvInclude)
        include_err = try
            Base.include(mod, harness)
            nothing
        catch caught
            caught
        end

        @test include_err === nothing
        @test isdefined(mod, :_requested_lengths)
    finally
        old_ns === nothing ? delete!(ENV, "NCTS_QMB_NS") : (ENV["NCTS_QMB_NS"] = old_ns)
        old_allow_large === nothing ?
            delete!(ENV, "NCTS_QMB_ALLOW_LARGE") :
            (ENV["NCTS_QMB_ALLOW_LARGE"] = old_allow_large)
    end
end

@testset "Expectations TOML Fixtures" begin
    data = TestExpectations.expectations_load("expectations/chsh_simple.toml")
    @test data["schema_version"] == 1

    case = TestExpectations.expectations_case(data, "Dense_d1")
    @test haskey(case, "expected")
    @test haskey(case["expected"], "objective")

    oracle = expectations_oracle("expectations/chsh_simple.toml", "Dense_d1")
    @test oracle.opt isa Float64
    @test oracle.nuniq isa Int
    @test oracle.sides isa Vector{Int}

    for relpath in (
        "expectations/heisenberg_qmbcertify_base.toml",
        "expectations/heisenberg_qmbcertify_rdm.toml",
    )
        qmb = TestExpectations.expectations_load(relpath)
        _test_fixture_probe_execution_states(qmb)
        @test qmb["schema_version"] == 1
        @test qmb["status"] == "pending_qmbcertify_runs"
        @test qmb["source"]["commit"] == "b18830a9460de4daa03013e389808d522c7823cf"
        @test qmb["source"]["remote"] == "git@github.com:wangjie212/QMBCertify.git"
        @test qmb["reviewed_reference_status_policy"]["accepted_termination_statuses"] ==
            ["OPTIMAL", "ALMOST_OPTIMAL"]
        @test qmb["reviewed_reference_status_policy"]["rejected_examples"] ==
            ["SLOW_PROGRESS/FEASIBLE_POINT"]
        @test occursin(
            "not reviewed",
            qmb["reviewed_reference_status_policy"]["rejection_note"],
        )
        @test qmb["qmbcertify_api"]["entrypoint"] == "GSB"
        @test qmb["qmbcertify_api"]["source_file"] == "src/sdp.jl"
        @test qmb["qmbcertify_api"]["defaults"]["lso"] == true
        @test qmb["qmbcertify_api"]["defaults"]["pso"] == 3
        @test qmb["qmbcertify_api"]["defaults"]["rdm"] == false
        @test qmb["qmbcertify_api"]["defaults"]["lattice"] == "chain"
        @test qmb["qmbcertify_api"]["rdm"]["supported_chain_values"] == [8, 9, 10]
        @test qmb["qmbcertify_api"]["rdm"]["local_word_filter"] == "even_axis_counts"
        @test qmb["qmbcertify_api"]["rdm"]["unsupported_behavior"] == "print_only"
        @test qmb["qmbcertify_api"]["linear_state_opt"]["candidate_generator"] ==
            "generate_mons(L, lol, rdm-1)"
        @test qmb["qmbcertify_api"]["linear_state_opt"]["closure_filter"] ==
            "filter_mons"
        @test qmb["qmbcertify_api"]["linear_state_opt"]["deduplication"] ==
            "random_weight_row_hash"
        @test qmb["qmbcertify_api"]["psd_state_opt"]["basis_filter"] ==
            "length.(basis[i]) .<= pso"
        @test qmb["qmbcertify_api"]["psd_state_opt"]["entry_function"] ==
            "PSDstate_entry"
        @test qmb["qmbcertify_api"]["psd_state_opt"]["endpoint_rule"] ==
            "exactly_one_endpoint_changes"
        @test qmb["required_runs"]["lengths"] == [10, 20, 30]
        @test haskey(qmb["required_runs"], "pending_cases")
        if haskey(qmb["required_runs"], "pending_cases")
            pending_cases = qmb["required_runs"]["pending_cases"]
            reviewed_cases = get(qmb, "cases", [])
            reviewed_ids = Set(case["id"] for case in reviewed_cases)
            @test length(pending_cases) == length(qmb["required_runs"]["profiles"]) *
                length(qmb["required_runs"]["lengths"])
            @test Set(case["profile"] for case in pending_cases) ==
                Set(qmb["required_runs"]["profiles"])
            @test Set(case["length"] for case in pending_cases) ==
                Set(qmb["required_runs"]["lengths"])
            @test all(
                case["status"] == (case["id"] in reviewed_ids ? "reviewed_run" : "missing_reviewed_run")
                for case in pending_cases
            )
            deleted_queue_cases = [
                case for case in pending_cases if
                get(case, "run_queue_status", "queued") == "deleted_until_evidence_gate"
            ]
            @test all(case["status"] == "missing_reviewed_run" for case in deleted_queue_cases)
            @test all(haskey(case, "restore_gate") for case in deleted_queue_cases)
            @test all(
                occursin("wall_rss_load_estimate", case["restore_gate"]) for
                case in deleted_queue_cases
            )
        end
        @test haskey(qmb, "profile_definitions")
        if haskey(qmb, "profile_definitions")
            @test all(
                haskey(qmb["profile_definitions"], profile)
                for profile in qmb["required_runs"]["profiles"]
            )
            @test all(
                haskey(qmb["profile_definitions"][profile], "command_args")
                for profile in qmb["required_runs"]["profiles"]
            )
        end
        @test all(
            field -> field in qmb["required_runs"]["required_fields"],
            [
                "objective",
                "psd_block_histogram",
                "moment_count",
                "build_time_seconds",
                "solve_time_seconds",
                "peak_rss_bytes",
            ],
        )
        @test all(
            field -> field in qmb["required_runs"]["environment_required_fields"],
            [
                "nctssos_commit",
                "julia_version",
                "mosek_version",
                "cpu_model",
                "ram_bytes",
                "thread_count",
                "blas_vendor",
            ],
        )
        if relpath == "expectations/heisenberg_qmbcertify_base.toml"
            @test length(qmb["cases"]) == 1
            reviewed = only(qmb["cases"])
            @test reviewed["id"] == "A0_L10"
            @test reviewed["profile"] == "A0"
            @test reviewed["length"] == 10
            @test reviewed["status"] == "reviewed_run"
            @test reviewed["environment"] == "autodl_2026_06_28"
            @test reviewed["objective"] ≈ -0.4515446349992884
            @test reviewed["objective_per_site"] ≈ reviewed["objective"]
            @test reviewed["objective_total_estimate"] ≈ 10 * reviewed["objective"]
            @test reviewed["moment_count"] == 527
            @test reviewed["psd_block_histogram"] ==
                [[11, 2], [15, 2], [22, 4], [30, 4], [38, 2], [42, 1], [43, 1], [76, 4], [84, 4]]
            @test reviewed["psd_max_block"] == 84
            @test reviewed["psd_scalar_variables"] == 32_559
            @test reviewed["build_time_seconds"] > 0
            @test reviewed["solve_time_seconds"] > 0
            @test reviewed["peak_rss_bytes"] > 0
            env = qmb["reference_environments"]["autodl_2026_06_28"]
            @test env["nctssos_commit"] == "1033964d36e273ed1e3f2cc457e47db134cf0a4e"
            @test env["julia_version"] == "1.12.6"
            @test env["mosek_version"] == "11.2.0"
            @test env["thread_count"] == 1
            @test env["qmbcertify_path"] == "/root/autodl-tmp/QMBCertify"
            pending_by_id = Dict(case["id"] => case for case in qmb["required_runs"]["pending_cases"])
            @test pending_by_id["A0_L20"]["run_queue_status"] ==
                "deleted_until_evidence_gate"
            @test pending_by_id["A0_L20"]["restore_gate"] ==
                "accepted_profile_variant_and_fresh_wall_rss_load_estimate"
            @test pending_by_id["A0_L20"]["failed_attempt_id"] ==
                "A0_L20_2026_07_06"
            @test pending_by_id["A0_L30"]["run_queue_status"] ==
                "deleted_until_evidence_gate"
            @test pending_by_id["A0_L30"]["restore_gate"] ==
                "accepted_profile_variant_and_fresh_wall_rss_load_estimate"
        else
            @test length(qmb["cases"]) == 2
            reviewed_by_id = Dict(case["id"] => case for case in qmb["cases"])
            @test Set(keys(reviewed_by_id)) == Set(["A1_L20", "A1_L30"])
            a1_l20 = reviewed_by_id["A1_L20"]
            a1_l30 = reviewed_by_id["A1_L30"]
            @test a1_l20["profile"] == "A1"
            @test a1_l20["length"] == 20
            @test a1_l20["status"] == "reviewed_run"
            @test a1_l20["objective"] ≈ -0.4452477259458568
            @test a1_l20["objective_total_estimate"] ≈
                20 * a1_l20["objective_per_site"]
            @test a1_l20["moment_count"] == 3534
            @test a1_l20["psd_block_histogram"] ==
                [[48, 2], [56, 1], [57, 1], [58, 1], [64, 1], [72, 1], [96, 9], [114, 9]]
            @test a1_l20["psd_max_block"] == 114
            @test a1_l20["psd_scalar_variables"] == 112_919
            @test a1_l20["termination_status"] == "OPTIMAL"
            @test a1_l20["peak_rss_bytes"] > 0
            @test a1_l30["profile"] == "A1"
            @test a1_l30["length"] == 30
            @test a1_l30["status"] == "reviewed_run"
            @test a1_l30["objective"] ≈ -0.4441690264512694
            @test a1_l30["objective_total_estimate"] ≈
                30 * a1_l30["objective_per_site"]
            @test a1_l30["moment_count"] == 7405
            @test a1_l30["psd_block_histogram"] ==
                [[48, 2], [56, 1], [57, 1], [58, 1], [64, 1], [72, 1], [96, 14], [114, 14]]
            @test a1_l30["psd_max_block"] == 114
            @test a1_l30["psd_scalar_variables"] == 168_974
            @test a1_l30["termination_status"] == "OPTIMAL"
            @test a1_l30["peak_rss_bytes"] > a1_l20["peak_rss_bytes"]
            reviewed_env = qmb["reference_environments"][a1_l20["environment"]]
            @test reviewed_env["nctssos_commit"] == "1033964d36e273ed1e3f2cc457e47db134cf0a4e"
            @test reviewed_env["julia_version"] == "1.12.6"
            @test reviewed_env["mosek_version"] == "11.2.0"
            @test reviewed_env["thread_count"] == 1
            @test reviewed_env["qmbcertify_path"] == "/root/autodl-tmp/QMBCertify"
            @test haskey(qmb, "failed_attempts")
            pending_by_id = Dict(case["id"] => case for case in qmb["required_runs"]["pending_cases"])
            @test pending_by_id["A1_L10"]["source_mapping"] == "extrapolated_from_example_L20"
            @test pending_by_id["A1_L20"]["status"] == "reviewed_run"
            @test pending_by_id["A1_L20"]["source_mapping"] == "examples/example.jl_line_10"
            @test pending_by_id["A1_L30"]["status"] == "reviewed_run"
            @test pending_by_id["A1_L30"]["source_mapping"] == "extrapolated_from_example_L20"
            @test pending_by_id["A2_L10"]["source_mapping"] == "candidate_from_gsb_defaults"
            @test pending_by_id["A2_L20"]["run_queue_status"] ==
                "deleted_until_evidence_gate"
            @test pending_by_id["A2_L20"]["restore_gate"] ==
                "accepted_profile_variant_and_fresh_wall_rss_load_estimate"
            @test pending_by_id["A2_L20"]["failed_attempt_id"] ==
                "A2_L20_2026_07_06"
            @test pending_by_id["A2_L30"]["run_queue_status"] ==
                "deleted_until_evidence_gate"
            @test pending_by_id["A2_L30"]["restore_gate"] ==
                "accepted_profile_variant_and_fresh_wall_rss_load_estimate"
            @test qmb["profile_definitions"]["A1"]["source_length"] == 20
            @test qmb["profile_definitions"]["A2"]["source_kind"] == "candidate_from_gsb_defaults"
            @test length(qmb["failed_attempts"]) == 4
            failed_by_id = Dict(attempt["id"] => attempt for attempt in qmb["failed_attempts"])
            a1_failed = failed_by_id["A1_L10_2026_06_29"]
            a1_rdm_only_probe = failed_by_id["A1_L10_rdm_only_probe_2026_06_29"]
            a2_failed = failed_by_id["A2_L10_2026_06_29"]
            a2_l20_failed = failed_by_id["A2_L20_2026_07_06"]

            for failed in (a1_failed, a1_rdm_only_probe, a2_failed)
                @test failed["length"] == 10
                @test failed["status"] == "rejected_nonoptimal_solve"
                @test failed["termination_status"] == "DUAL_INFEASIBLE"
                @test failed["solution_status"] == "INFEASIBILITY_CERTIFICATE"
                @test failed["reported_sdp_n"] == 114
                @test failed["reported_sdp_m"] == 641
                @test failed["basis_time_seconds"] > 0
                @test failed["rdm_time_seconds"] > failed["basis_time_seconds"]
                @test failed["solve_time_seconds"] > 0
                failed_env = qmb["failed_attempt_environments"][failed["environment"]]
                @test failed_env["julia_version"] == "1.12.6"
                @test failed_env["mosek_version"] == "11.2.0"
                @test failed_env["thread_count"] == 1
                @test failed_env["qmbcertify_path"] == "/root/autodl-tmp/QMBCertify"
            end

            @test a1_failed["case_id"] == "A1_L10"
            @test a1_failed["profile"] == "A1"
            @test a1_failed["objective_reported"] ≈ 0.4633637063586052
            @test !haskey(a1_failed, "psd_state_opt_time_seconds")

            @test a1_rdm_only_probe["case_id"] == "A1_L10"
            @test a1_rdm_only_probe["profile"] == "A1_rdm_only_probe"
            @test a1_rdm_only_probe["parent_profile"] == "A1"
            @test occursin("lso=false", a1_rdm_only_probe["command_args_summary"])
            @test a1_rdm_only_probe["overrides_summary"] == "lso=false"
            @test a1_rdm_only_probe["objective_reported"] ≈ 0.463363705028278
            @test !haskey(a1_rdm_only_probe, "linear_state_opt_time_seconds")
            @test !haskey(a1_rdm_only_probe, "psd_state_opt_time_seconds")

            @test a2_failed["case_id"] == "A2_L10"
            @test a2_failed["profile"] == "A2"
            @test a2_failed["objective_reported"] ≈ 0.6273534182534728
            @test a2_failed["psd_state_opt_time_seconds"] > 0

            @test a2_l20_failed["case_id"] == "A2_L20"
            @test a2_l20_failed["profile"] == "A2"
            @test a2_l20_failed["length"] == 20
            @test a2_l20_failed["status"] == "rejected_nonoptimal_solve"
            @test a2_l20_failed["termination_status"] == "SLOW_PROGRESS"
            @test a2_l20_failed["solution_status"] == "FEASIBLE_POINT"
            @test a2_l20_failed["objective_reported"] ≈ -0.4452212294966909
            @test a2_l20_failed["reported_sdp_n"] == 114
            @test a2_l20_failed["reported_sdp_m"] == 5034
            @test a2_l20_failed["rdm_time_seconds"] > 0
            @test a2_l20_failed["psd_state_opt_time_seconds"] > 0
            @test a2_l20_failed["solve_time_seconds"] > 0

            @test haskey(qmb, "profile_probes")
            @test length(qmb["profile_probes"]) == 2
            probes_by_id = Dict(probe["id"] => probe for probe in qmb["profile_probes"])
            a1_l12_probe = probes_by_id["A1_L12_profile_probe_L12_2026_06_30"]
            a2_l12_probe = probes_by_id["A2_L12_profile_probe_L12_2026_06_30"]
            @test a1_l12_probe["id"] == "A1_L12_profile_probe_L12_2026_06_30"
            @test a1_l12_probe["case_id"] == "A1_L12_probe"
            @test a1_l12_probe["profile"] == "A1_L12_profile_probe"
            @test a1_l12_probe["parent_profile"] == "A1"
            @test a1_l12_probe["length"] == 12
            @test a1_l12_probe["status"] == "probe_run"
            @test a1_l12_probe["objective"] ≈ -0.4489492424398082
            @test a1_l12_probe["objective_total_estimate"] ≈
                12 * a1_l12_probe["objective_per_site"]
            @test a1_l12_probe["moment_count"] == 1044
            @test a1_l12_probe["psd_block_histogram"] ==
                [[48, 2], [56, 1], [57, 1], [58, 1], [64, 1], [72, 1], [96, 5], [114, 5]]
            @test a1_l12_probe["psd_max_block"] == 114
            @test a1_l12_probe["psd_scalar_variables"] == 68_075
            @test a1_l12_probe["build_time_seconds"] > 0
            @test a1_l12_probe["solve_time_seconds"] > 0
            @test a1_l12_probe["termination_status"] == "OPTIMAL"
            @test a1_l12_probe["solution_status"] == "missing"
            @test occursin(
                "optimal termination status is inferred",
                a1_l12_probe["notes"],
            )
            @test a1_l12_probe["peak_rss_bytes"] > 0
            @test occursin("rdm=8", a1_l12_probe["command_args_summary"])
            probe_env = qmb["profile_probe_environments"][a1_l12_probe["environment"]]
            @test probe_env["julia_version"] == "1.12.6"
            @test probe_env["mosek_version"] == "11.2.0"
            @test probe_env["thread_count"] == 1
            @test probe_env["qmbcertify_path"] == "/root/autodl-tmp/QMBCertify"

            @test a2_l12_probe["case_id"] == "A2_L12_probe"
            @test a2_l12_probe["profile"] == "A2_L12_profile_probe"
            @test a2_l12_probe["parent_profile"] == "A2"
            @test a2_l12_probe["length"] == 12
            @test a2_l12_probe["status"] == "probe_run"
            @test a2_l12_probe["objective"] ≈ -0.4489491908788246
            @test a2_l12_probe["objective_total_estimate"] ≈
                12 * a2_l12_probe["objective_per_site"]
            @test a2_l12_probe["moment_count"] == 1184
            @test a2_l12_probe["psd_block_histogram"] ==
                [[28, 2], [36, 2], [48, 2], [56, 6], [57, 1], [58, 1], [64, 1], [72, 6], [96, 5], [114, 5]]
            @test a2_l12_probe["psd_max_block"] == 114
            @test a2_l12_probe["psd_scalar_variables"] == 91_339
            @test a2_l12_probe["build_time_seconds"] > a1_l12_probe["build_time_seconds"]
            @test a2_l12_probe["solve_time_seconds"] > a1_l12_probe["solve_time_seconds"]
            @test a2_l12_probe["termination_status"] == "OPTIMAL"
            @test a2_l12_probe["solution_status"] == "missing"
            @test a2_l12_probe["peak_rss_bytes"] > a1_l12_probe["peak_rss_bytes"]
            @test occursin("pso=3", a2_l12_probe["command_args_summary"])
            @test a2_l12_probe["environment"] == a1_l12_probe["environment"]

            @test haskey(qmb, "nctssos_source_like_solver_probes")
            @test length(qmb["nctssos_source_like_solver_probes"]) == 16
            nctssos_by_id =
                Dict(probe["id"] => probe for probe in qmb["nctssos_source_like_solver_probes"])
            @test Set(keys(nctssos_by_id)) ==
                Set([
                    "NCTSSOS_A0_L12_source_like_sos_dual_mosek_2026_07_02",
                    "NCTSSOS_A1_L12_source_like_sos_dual_mosek_2026_07_02",
                    "NCTSSOS_A2_L12_source_like_sos_dual_mosek_2026_07_02",
                    "NCTSSOS_A2_SU2_RDM_PSO_L8_source_like_sos_dual_mosek_2026_07_02",
                    "NCTSSOS_A2_SU2_RDM_LSO_PSO_L8_source_like_sos_dual_mosek_2026_07_02",
                    "NCTSSOS_A2_SU2_RDM_LSO_PSO_REALLIFT_COMPLEX_L8_source_like_sos_dual_mosek_2026_07_02",
                    "NCTSSOS_A2_SU2_RDM_LSO_PSO_NATIVE_L8_source_like_sos_dual_mosek_2026_07_02",
                    "NCTSSOS_A2_L20_source_like_sos_dual_clarabel_2026_07_06",
                    "NCTSSOS_A2_L30_source_like_sos_dual_clarabel_2026_07_06",
                    "NCTSSOS_A1_L20_source_like_sos_dual_mosek_2026_07_06",
                    "NCTSSOS_A1_L20_source_like_sos_dual_mosek_2026_07_08",
                    "NCTSSOS_A1_L30_source_like_sos_dual_clarabel_2026_07_08",
                    "NCTSSOS_A1_L30_source_like_sos_dual_mosek_2026_07_08",
                    "NCTSSOS_A2_L20_source_like_sos_dual_mosek_2026_07_06",
                    "NCTSSOS_A1_L30_source_like_sos_dual_mosek_2026_07_06",
                    "NCTSSOS_A2_SU2_MOMENT_QUOTIENT_L30_source_like_sos_dual_mosek_2026_07_11",
                ])
            nctssos_a0 =
                nctssos_by_id["NCTSSOS_A0_L12_source_like_sos_dual_mosek_2026_07_02"]
            nctssos_a1 =
                nctssos_by_id["NCTSSOS_A1_L12_source_like_sos_dual_mosek_2026_07_02"]
            nctssos_a2 =
                nctssos_by_id["NCTSSOS_A2_L12_source_like_sos_dual_mosek_2026_07_02"]
            nctssos_su2_pso =
                nctssos_by_id["NCTSSOS_A2_SU2_RDM_PSO_L8_source_like_sos_dual_mosek_2026_07_02"]
            nctssos_su2_lso_pso =
                nctssos_by_id["NCTSSOS_A2_SU2_RDM_LSO_PSO_L8_source_like_sos_dual_mosek_2026_07_02"]
            nctssos_su2_complex_lift =
                nctssos_by_id["NCTSSOS_A2_SU2_RDM_LSO_PSO_REALLIFT_COMPLEX_L8_source_like_sos_dual_mosek_2026_07_02"]
            nctssos_su2_native =
                nctssos_by_id["NCTSSOS_A2_SU2_RDM_LSO_PSO_NATIVE_L8_source_like_sos_dual_mosek_2026_07_02"]
            nctssos_a2_l20_model =
                nctssos_by_id["NCTSSOS_A2_L20_source_like_sos_dual_clarabel_2026_07_06"]
            nctssos_a2_l30_model =
                nctssos_by_id["NCTSSOS_A2_L30_source_like_sos_dual_clarabel_2026_07_06"]
            nctssos_a1_l20_solve =
                nctssos_by_id["NCTSSOS_A1_L20_source_like_sos_dual_mosek_2026_07_06"]
            nctssos_source_base_a1_l20_solve =
                nctssos_by_id["NCTSSOS_A1_L20_source_like_sos_dual_mosek_2026_07_08"]
            nctssos_source_base_a1_l30_model =
                nctssos_by_id["NCTSSOS_A1_L30_source_like_sos_dual_clarabel_2026_07_08"]
            nctssos_source_base_a1_l30_solve =
                nctssos_by_id["NCTSSOS_A1_L30_source_like_sos_dual_mosek_2026_07_08"]
            nctssos_a2_l20_solve =
                nctssos_by_id["NCTSSOS_A2_L20_source_like_sos_dual_mosek_2026_07_06"]
            nctssos_a1_l30_solve =
                nctssos_by_id["NCTSSOS_A1_L30_source_like_sos_dual_mosek_2026_07_06"]
            nctssos_su2_quotient_l30_solve = nctssos_by_id[
                "NCTSSOS_A2_SU2_MOMENT_QUOTIENT_L30_source_like_sos_dual_mosek_2026_07_11"
            ]

            for probe in values(nctssos_by_id)
                @test probe["status"] == "nctssos_solver_probe"
                @test probe["model_mode"] == "sos_dual"
                @test probe["optimizer"] in ["Clarabel", "Mosek"]
                @test probe["probe_harness"] == "perf/pauli_translation_solver_probe.jl"
                if probe["id"] in (
                    "NCTSSOS_A1_L20_source_like_sos_dual_mosek_2026_07_08",
                    "NCTSSOS_A1_L30_source_like_sos_dual_mosek_2026_07_08",
                )
                    @test occursin("bound quality parity is met", probe["comparison_scope"])
                else
                    @test occursin("not a reviewed QMBCertify parity run", probe["comparison_scope"]) ||
                        occursin("bound quality parity is not met", probe["comparison_scope"])
                end
                @test probe["construction_only"] == false
                @test probe["model_built"] == true
                @test probe["su2_base_zero_row_count"] ==
                    probe["su2_base_spin_offblock_row_count"] +
                    probe["su2_base_magnetic_offdiag_row_count"] +
                    probe["su2_base_magnetic_copy_row_count"]
                @test probe["su2_base_zero_row_count"] >= 0
                @test probe["su2_base_spin_offblock_row_count"] >= 0
                @test probe["su2_base_magnetic_offdiag_row_count"] >= 0
                @test probe["su2_base_magnetic_copy_row_count"] >= 0
                @test probe["order"] == 4
                @test probe["sos_hermitian_representation"] in ["real_lift", "native"]
                @test probe["construction_time_seconds"] > 0
                @test probe["dualization_time_seconds"] > 0
                @test probe["script_wall_seconds"] > 0
                @test probe["remote_wall_seconds"] > 0
                @test probe["peak_rss_bytes"] > 0
                if probe["solved"]
                    @test probe["termination_status"] == "OPTIMAL"
                    @test probe["solve_time_seconds"] > 0
                    @test haskey(probe, "objective")
                    @test probe["sos_certificate_moment_count"] == probe["linear_moments"]
                    @test probe["sos_certificate_residual_source"] ==
                        "sos_dual_certificate_residual"
                    @test probe["sos_certificate_diagnostics_source"] ==
                        "sos_dual_certificate_diagnostics"
                    @test probe["sos_certificate_max_abs_residual"] <= 1e-6
                    @test probe["sos_certificate_identity_abs_residual"] <= 1e-6
                    @test probe["sos_certificate_worst_abs_residual"] ==
                        probe["sos_certificate_max_abs_residual"]
                    @test probe["sos_certificate_psd_block_count"] == probe["psd_cones"]
                    @test probe["sos_certificate_zero_dual_count"] ==
                        probe["zero_dual_variables"]
                    @test probe["sos_certificate_max_native_hermitian_residual"] <= 1e-10
                    @test probe["sos_certificate_min_native_eigenvalue"] >= -1e-8
                    @test probe["sos_certificate_max_zero_dual_abs_value"] >= 0
                    if probe["zero_dual_variables"] == 0
                        @test probe["sos_certificate_max_zero_dual_abs_value"] == 0
                    end
                else
                    @test probe["termination_status"] == "not_solved"
                    @test !haskey(probe, "solve_time_seconds")
                    @test !haskey(probe, "objective")
                    @test !haskey(probe, "sos_certificate_moment_count")
                end
                env = qmb["nctssos_solver_probe_environments"][probe["environment"]]
                @test env["julia_version"] == "1.12.6"
                @test env["thread_count"] == 1
                if probe["optimizer"] == "Mosek"
                    @test env["mosek_version"] == "11.2.2"
                    @test get(env, "mosek_thread_count", 1) == 1
                else
                    @test !haskey(env, "mosek_version")
                end
            end

            for probe in (nctssos_a0, nctssos_a1, nctssos_a2)
                @test probe["length"] == 12
                @test probe["reflection_symmetry"] == true
                @test probe["sign_symmetry"] == true
                @test probe["real_moment_matrix"] == true
                @test probe["u1_symmetry"] == false
                @test probe["su2_symmetry"] == false
                @test probe["base_su2_extend_rdm"] == false
                @test probe["axis_rotation_symmetry"] == true
                @test probe["axis_rotation_equalities"] == true
                @test probe["axis_rotation_quotient"] == true
                @test probe["report_su2_moment_symmetry"] == false
                @test probe["report_su2_rdm_symmetry"] == false
                @test probe["report_u1_rdm_symmetry"] == false
                @test probe["report_axis_rotation_symmetry"] == true
                @test probe["report_axis_quotient"] == true
                @test probe["solve_supported"] == true
                @test probe["solve_blocker"] == "none"
                @test probe["solve_support_reason"] == ""
                @test probe["solve_unsupported_block_features"] == []
                @test probe["solve_unsupported_zero_features"] == []
                @test probe["linear_state_opt_width"] in [0, 7]
                @test probe["psd_state_opt_width"] in [0, 3]
            end

            for probe in (nctssos_su2_pso, nctssos_su2_lso_pso)
                @test probe["length"] == 8
                @test probe["reflection_symmetry"] == true
                @test probe["sign_symmetry"] == true
                @test probe["real_moment_matrix"] == true
                @test probe["u1_symmetry"] == false
                @test probe["su2_symmetry"] == true
                @test probe["base_su2_extend_rdm"] == false
                @test probe["axis_rotation_symmetry"] == false
                @test probe["axis_rotation_equalities"] == false
                @test probe["axis_rotation_quotient"] == false
                @test probe["report_su2_moment_symmetry"] == false
                @test probe["report_su2_rdm_symmetry"] == true
                @test probe["report_u1_rdm_symmetry"] == false
                @test probe["report_axis_rotation_symmetry"] == false
                @test probe["report_axis_quotient"] == false
                @test probe["solve_supported"] == true
                @test probe["solve_blocker"] == "none"
                @test probe["solve_support_reason"] == ""
                @test probe["solve_unsupported_block_features"] == []
                @test probe["solve_unsupported_zero_features"] == []
                @test probe["linear_moments"] == 4_098
                @test probe["free_keys"] == 4_097
                @test probe["contiguous_rdm_k"] == 8
                @test probe["contiguous_rdm_decomposition"] == "su2"
                @test probe["contiguous_rdm_support"] == "extend"
                @test probe["qmbcertify_rdm_block_sizes"] == [28, 56, 40, 14, 1]
                @test probe["psd_state_opt_width"] == 3
                @test probe["psd_state_opt_block_sizes"] == fill(18, 5)
                @test probe["psd_cones"] == 50
                @test probe["max_block"] == 56
            end

            for probe in (nctssos_su2_complex_lift, nctssos_su2_native)
                @test probe["length"] == 8
                @test probe["reflection_symmetry"] == false
                @test probe["sign_symmetry"] == true
                @test probe["real_moment_matrix"] == false
                @test probe["u1_symmetry"] == false
                @test probe["su2_symmetry"] == true
                @test probe["base_su2_extend_rdm"] == false
                @test probe["axis_rotation_symmetry"] == false
                @test probe["axis_rotation_equalities"] == false
                @test probe["axis_rotation_quotient"] == false
                @test probe["report_su2_moment_symmetry"] == false
                @test probe["report_su2_rdm_symmetry"] == true
                @test probe["report_u1_rdm_symmetry"] == false
                @test probe["report_axis_rotation_symmetry"] == false
                @test probe["report_axis_quotient"] == false
                @test probe["solve_supported"] == true
                @test probe["solve_blocker"] == "none"
                @test probe["solve_support_reason"] == ""
                @test probe["solve_unsupported_block_features"] == []
                @test probe["solve_unsupported_zero_features"] == []
                @test probe["linear_moments"] == 4_098
                @test probe["free_keys"] == 4_097
                @test probe["linear_state_opt_width"] == 7
                @test probe["linear_state_opt_row_count"] == 819
                @test probe["contiguous_rdm_k"] == 8
                @test probe["contiguous_rdm_decomposition"] == "su2"
                @test probe["contiguous_rdm_support"] == "extend"
                @test probe["qmbcertify_rdm_block_sizes"] == [14, 28, 20, 7, 1]
                @test probe["psd_state_opt_width"] == 3
                @test probe["psd_state_opt_block_sizes"] == fill(9, 8)
                @test probe["psd_cones"] == 45
                @test probe["zero_dual_variables"] == 819
                @test probe["max_block"] == 31
            end

            @test nctssos_su2_complex_lift["sos_hermitian_representation"] == "real_lift"
            @test nctssos_su2_complex_lift["dual_variables"] == 63_800
            @test nctssos_su2_complex_lift["objective"] ≈ -3.6510934102555996
            @test nctssos_su2_native["sos_hermitian_representation"] == "native"
            @test nctssos_su2_native["dual_variables"] == 31_758
            @test nctssos_su2_native["objective"] ≈ -3.6510934109179667
            @test nctssos_su2_native["dual_variables"] <
                nctssos_su2_complex_lift["dual_variables"]
            @test nctssos_su2_native["dual_variables"] <=
                nctssos_su2_complex_lift["dual_variables"] / 2
            @test abs(
                nctssos_su2_native["objective"] -
                nctssos_su2_complex_lift["objective"],
            ) <= 1e-8

            @test nctssos_su2_quotient_l30_solve["profile"] ==
                "A2_SU2_MOMENT_QUOTIENT"
            @test nctssos_su2_quotient_l30_solve["length"] == 30
            @test nctssos_su2_quotient_l30_solve["comparison_scope"] ==
                "large-N NCTSSoS backend evidence; not a reviewed QMBCertify parity run"
            @test nctssos_su2_quotient_l30_solve["representation"] == "complex"
            @test nctssos_su2_quotient_l30_solve["sos_hermitian_representation"] ==
                "native"
            @test nctssos_su2_quotient_l30_solve["sos_coefficient_scaling"] ==
                "max_abs"
            @test nctssos_su2_quotient_l30_solve["sos_coefficient_scaling_floor"] ==
                0.1
            @test nctssos_su2_quotient_l30_solve["su2_moment_quotient"] == true
            @test nctssos_su2_quotient_l30_solve["su2_moment_raw_count"] == 94_129
            @test nctssos_su2_quotient_l30_solve["su2_moment_quotient_count"] ==
                3_529
            @test nctssos_su2_quotient_l30_solve["termination_status"] == "OPTIMAL"
            @test nctssos_su2_quotient_l30_solve["primal_status"] == "FEASIBLE_POINT"
            @test nctssos_su2_quotient_l30_solve["dual_status"] == "FEASIBLE_POINT"
            @test nctssos_su2_quotient_l30_solve["relative_gap"] <= 2e-8
            @test nctssos_su2_quotient_l30_solve["sos_certificate_max_abs_residual"] <=
                1e-6
            @test nctssos_su2_quotient_l30_solve["sos_certificate_min_native_eigenvalue"] >=
                -1e-8
            @test nctssos_su2_quotient_l30_solve["dual_variables"] == 12_001
            @test nctssos_su2_quotient_l30_solve["scalar_equalities"] == 4_528
            @test nctssos_su2_quotient_l30_solve["max_block"] == 28
            @test nctssos_su2_quotient_l30_solve["remote_wall_seconds"] <= 960
            @test nctssos_su2_quotient_l30_solve["peak_rss_bytes"] <= 15 * 2^30

            @test nctssos_a0["objective"] ≈ -5.47924827890285
            @test nctssos_a0["linear_moments"] == 1621
            @test nctssos_a0["free_keys"] == 1620
            @test nctssos_a0["psd_cones"] == 28
            @test nctssos_a0["dual_variables"] == 86_379
            @test nctssos_a0["zero_dual_variables"] == 0
            @test nctssos_a0["contiguous_rdm_k"] == 0
            @test nctssos_a0["contiguous_rdm_decomposition"] == "none"
            @test nctssos_a0["qmbcertify_rdm_block_sizes"] == []
            @test nctssos_a0["psd_state_opt_block_sizes"] == []

            @test nctssos_a1["objective"] ≈ -5.392564604362093
            @test nctssos_a1["qmbcertify_probe_id"] ==
                a1_l12_probe["id"]
            @test nctssos_a1["qmbcertify_objective_total_estimate"] ≈
                a1_l12_probe["objective_total_estimate"]
            @test nctssos_a1["objective_minus_qmbcertify_total_estimate"] ≈
                nctssos_a1["objective"] - a1_l12_probe["objective_total_estimate"]
            @test nctssos_a1["objective_minus_qmbcertify_total_estimate"] < 0
            @test nctssos_a1["linear_moments"] == 2143
            @test nctssos_a1["free_keys"] == 2142
            @test nctssos_a1["linear_state_opt_row_count"] == 810
            @test nctssos_a1["psd_cones"] == 31
            @test nctssos_a1["dual_variables"] == 93_493
            @test nctssos_a1["zero_dual_variables"] == 810
            @test nctssos_a1["contiguous_rdm_k"] == 8
            @test nctssos_a1["contiguous_rdm_decomposition"] == "qmbcertify"
            @test nctssos_a1["contiguous_rdm_support"] == "extend"
            @test nctssos_a1["qmbcertify_rdm_block_sizes"] == [72, 64, 56]
            @test nctssos_a1["psd_state_opt_block_sizes"] == []

            @test nctssos_a2["objective"] ≈ -5.392562452506382
            @test nctssos_a2["qmbcertify_probe_id"] ==
                a2_l12_probe["id"]
            @test nctssos_a2["qmbcertify_objective_total_estimate"] ≈
                a2_l12_probe["objective_total_estimate"]
            @test nctssos_a2["objective_minus_qmbcertify_total_estimate"] ≈
                nctssos_a2["objective"] - a2_l12_probe["objective_total_estimate"]
            @test nctssos_a2["objective_minus_qmbcertify_total_estimate"] < 0
            @test nctssos_a2["linear_moments"] == 2158
            @test nctssos_a2["free_keys"] == 2157
            @test nctssos_a2["linear_state_opt_row_count"] == 810
            @test nctssos_a2["psd_cones"] == 38
            @test nctssos_a2["dual_variables"] == 94_690
            @test nctssos_a2["zero_dual_variables"] == 810
            @test nctssos_a2["contiguous_rdm_k"] == 8
            @test nctssos_a2["contiguous_rdm_decomposition"] == "qmbcertify"
            @test nctssos_a2["contiguous_rdm_support"] == "extend"
            @test nctssos_a2["qmbcertify_rdm_block_sizes"] == [72, 64, 56]
            @test nctssos_a2["psd_state_opt_block_sizes"] == fill(18, 7)
            @test nctssos_a1["objective"] - nctssos_a0["objective"] ≈ 0.08668367454075705
            @test nctssos_a2["objective"] - nctssos_a1["objective"] ≈ 2.151855711174119e-6

            for probe in (
                nctssos_a2_l20_model,
                nctssos_a2_l30_model,
                nctssos_a1_l20_solve,
                nctssos_a2_l20_solve,
                nctssos_a1_l30_solve,
            )
                @test probe["length"] in [20, 30]
                @test probe["reflection_symmetry"] == true
                @test probe["sign_symmetry"] == true
                @test probe["real_moment_matrix"] == true
                @test probe["u1_symmetry"] == false
                @test probe["su2_symmetry"] == false
                @test probe["base_su2_extend_rdm"] == false
                @test probe["axis_rotation_symmetry"] == true
                @test probe["axis_rotation_equalities"] == true
                @test probe["axis_rotation_quotient"] == true
                @test probe["report_su2_moment_symmetry"] == false
                @test probe["report_su2_rdm_symmetry"] == false
                @test probe["report_u1_rdm_symmetry"] == false
                @test probe["report_axis_rotation_symmetry"] == true
                @test probe["report_axis_quotient"] == true
                @test probe["solve_supported"] == true
                @test probe["solve_blocker"] == "none"
                @test probe["solve_support_reason"] == ""
                @test probe["solve_unsupported_block_features"] == []
                @test probe["solve_unsupported_zero_features"] == []
                @test probe["linear_state_opt_width"] == 7
                @test probe["linear_state_opt_row_count"] == 810
                @test probe["contiguous_rdm_k"] == 8
                @test probe["contiguous_rdm_decomposition"] == "qmbcertify"
                @test probe["contiguous_rdm_support"] == "extend"
                @test probe["qmbcertify_rdm_block_sizes"] == [72, 64, 56]
                @test probe["zero_dual_variables"] == 810
                @test probe["max_block"] == 120
                @test probe["peak_rss_bytes"] < 8 * 1024^3
            end

            for probe in (nctssos_a2_l20_model, nctssos_a2_l30_model)
                @test probe["optimizer"] == "Clarabel"
                @test probe["termination_status"] == "not_solved"
                @test probe["solved"] == false
                @test probe["model_built"] == true
                @test probe["psd_state_opt_width"] == 3
            end
            @test nctssos_a2_l20_model["length"] == 20
            @test nctssos_a2_l20_model["linear_moments"] == 4614
            @test nctssos_a2_l20_model["free_keys"] == 4613
            @test nctssos_a2_l20_model["psd_state_opt_block_sizes"] == fill(18, 11)
            @test nctssos_a2_l20_model["psd_cones"] == 50
            @test nctssos_a2_l20_model["dual_variables"] == 153_454
            @test nctssos_a2_l20_model["construction_time_seconds"] ≈ 160.60223414
            @test nctssos_a2_l20_model["dualization_time_seconds"] ≈ 3.136656289
            @test nctssos_a2_l20_model["remote_wall_seconds"] ≈ 242.63

            @test nctssos_a2_l30_model["length"] == 30
            @test nctssos_a2_l30_model["linear_moments"] == 7684
            @test nctssos_a2_l30_model["free_keys"] == 7683
            @test nctssos_a2_l30_model["psd_state_opt_block_sizes"] == fill(18, 16)
            @test nctssos_a2_l30_model["psd_cones"] == 65
            @test nctssos_a2_l30_model["dual_variables"] == 226_909
            @test nctssos_a2_l30_model["construction_time_seconds"] ≈ 430.666762024
            @test nctssos_a2_l30_model["dualization_time_seconds"] ≈ 7.840913287
            @test nctssos_a2_l30_model["remote_wall_seconds"] ≈ 517.50

            @test nctssos_a1_l20_solve["optimizer"] == "Mosek"
            @test nctssos_a1_l20_solve["profile"] == "A1"
            @test nctssos_a1_l20_solve["length"] == 20
            @test nctssos_a1_l20_solve["objective"] ≈ -9.092707010191038
            @test nctssos_a1_l20_solve["qmbcertify_probe_id"] == "A1_L20"
            @test nctssos_a1_l20_solve["qmbcertify_objective_total_estimate"] ≈
                a1_l20["objective_total_estimate"]
            @test nctssos_a1_l20_solve["objective_minus_qmbcertify_total_estimate"] ≈
                nctssos_a1_l20_solve["objective"] -
                a1_l20["objective_total_estimate"]
            @test nctssos_a1_l20_solve["objective_minus_qmbcertify_total_estimate"] < -0.18
            @test nctssos_a1_l20_solve["linear_moments"] == 4559
            @test nctssos_a1_l20_solve["free_keys"] == 4558
            @test nctssos_a1_l20_solve["psd_state_opt_width"] == 0
            @test nctssos_a1_l20_solve["psd_state_opt_block_sizes"] == []
            @test nctssos_a1_l20_solve["psd_cones"] == 39
            @test nctssos_a1_l20_solve["dual_variables"] == 151_573
            @test nctssos_a1_l20_solve["solve_time_seconds"] ≈ 118.795514017
            @test nctssos_a1_l20_solve["remote_wall_seconds"] ≈ 377.39
            @test nctssos_a1_l20_solve["sos_certificate_max_abs_residual"] <= 4e-9

            @test nctssos_source_base_a1_l20_solve["optimizer"] == "Mosek"
            @test nctssos_source_base_a1_l20_solve["profile"] == "A1"
            @test nctssos_source_base_a1_l20_solve["length"] == 20
            @test occursin(
                "source-base",
                nctssos_source_base_a1_l20_solve["comparison_scope"],
            )
            @test nctssos_source_base_a1_l20_solve["qmbcertify_base_construct"] == true
            @test nctssos_source_base_a1_l20_solve["qmbcertify_base_extra"] == 9
            @test nctssos_source_base_a1_l20_solve["reflection_symmetry"] == false
            @test nctssos_source_base_a1_l20_solve["axis_rotation_symmetry"] == false
            @test nctssos_source_base_a1_l20_solve["objective"] ≈ -8.904954795042105
            @test nctssos_source_base_a1_l20_solve["qmbcertify_probe_id"] == "A1_L20"
            @test nctssos_source_base_a1_l20_solve["qmbcertify_objective_total_estimate"] ≈
                a1_l20["objective_total_estimate"]
            @test nctssos_source_base_a1_l20_solve[
                "objective_minus_qmbcertify_total_estimate"
            ] ≈ nctssos_source_base_a1_l20_solve["objective"] -
                 a1_l20["objective_total_estimate"]
            @test nctssos_source_base_a1_l20_solve[
                "objective_minus_qmbcertify_total_estimate"
            ] < 0
            @test abs(
                nctssos_source_base_a1_l20_solve[
                    "objective_minus_qmbcertify_total_estimate"
                ],
            ) < 1e-6
            @test nctssos_source_base_a1_l20_solve["linear_moments"] == 3534
            @test nctssos_source_base_a1_l20_solve["free_keys"] == 3505
            @test nctssos_source_base_a1_l20_solve["linear_state_opt_width"] == 7
            @test nctssos_source_base_a1_l20_solve["linear_state_opt_mode"] == "qmbcertify"
            @test nctssos_source_base_a1_l20_solve["linear_state_opt_row_count"] == 198
            @test nctssos_source_base_a1_l20_solve["contiguous_rdm_k"] == 8
            @test nctssos_source_base_a1_l20_solve["contiguous_rdm_decomposition"] ==
                "qmbcertify"
            @test nctssos_source_base_a1_l20_solve["contiguous_rdm_support"] == "extend"
            @test nctssos_source_base_a1_l20_solve["qmbcertify_rdm_block_sizes"] ==
                [72, 64, 56]
            @test nctssos_source_base_a1_l20_solve["psd_state_opt_width"] == 0
            @test nctssos_source_base_a1_l20_solve["psd_cones"] == 25
            @test nctssos_source_base_a1_l20_solve["dual_variables"] == 113_117
            @test nctssos_source_base_a1_l20_solve["zero_dual_variables"] == 198
            @test nctssos_source_base_a1_l20_solve["max_block"] == 114
            @test nctssos_source_base_a1_l20_solve["solve_time_seconds"] ≈ 35.763360136
            @test nctssos_source_base_a1_l20_solve["remote_wall_seconds"] ≈ 142.57
            @test nctssos_source_base_a1_l20_solve["peak_rss_bytes"] < 3 * 1024^3
            @test nctssos_source_base_a1_l20_solve[
                "sos_certificate_max_abs_residual"
            ] <= 2e-8

            @test nctssos_source_base_a1_l30_model["optimizer"] == "Clarabel"
            @test nctssos_source_base_a1_l30_model["profile"] == "A1"
            @test nctssos_source_base_a1_l30_model["length"] == 30
            @test nctssos_source_base_a1_l30_model["qmbcertify_base_construct"] == true
            @test nctssos_source_base_a1_l30_model["qmbcertify_base_extra"] == 9
            @test nctssos_source_base_a1_l30_model["reflection_symmetry"] == false
            @test nctssos_source_base_a1_l30_model["axis_rotation_symmetry"] == false
            @test nctssos_source_base_a1_l30_model["termination_status"] == "not_solved"
            @test nctssos_source_base_a1_l30_model["model_built"] == true
            @test nctssos_source_base_a1_l30_model["solved"] == false
            @test nctssos_source_base_a1_l30_model["estimated_model_size_gate_status"] ==
                "ok"
            @test nctssos_source_base_a1_l30_model["estimated_model_variables"] == 169_600
            @test nctssos_source_base_a1_l30_model[
                "estimated_scalar_equalities_upper_bound"
            ] == 7404
            @test nctssos_source_base_a1_l30_model["estimated_dense_schur_bytes"] ==
                438_553_728
            @test nctssos_source_base_a1_l30_model["linear_moments"] == 7405
            @test nctssos_source_base_a1_l30_model["free_keys"] == 7404
            @test nctssos_source_base_a1_l30_model["linear_state_opt_width"] == 7
            @test nctssos_source_base_a1_l30_model["linear_state_opt_mode"] == "qmbcertify"
            @test nctssos_source_base_a1_l30_model["linear_state_opt_row_count"] == 626
            @test nctssos_source_base_a1_l30_model["contiguous_rdm_k"] == 8
            @test nctssos_source_base_a1_l30_model["contiguous_rdm_decomposition"] ==
                "qmbcertify"
            @test nctssos_source_base_a1_l30_model["contiguous_rdm_support"] == "extend"
            @test nctssos_source_base_a1_l30_model["qmbcertify_rdm_block_sizes"] ==
                [72, 64, 56]
            @test nctssos_source_base_a1_l30_model["psd_state_opt_width"] == 0
            @test nctssos_source_base_a1_l30_model["psd_cones"] == 35
            @test nctssos_source_base_a1_l30_model["dual_variables"] == 169_600
            @test nctssos_source_base_a1_l30_model["zero_dual_variables"] == 626
            @test nctssos_source_base_a1_l30_model["max_block"] == 114
            @test nctssos_source_base_a1_l30_model["construction_time_seconds"] ≈
                16.723815959
            @test nctssos_source_base_a1_l30_model["dualization_time_seconds"] ≈
                2.731044266
            @test nctssos_source_base_a1_l30_model["remote_wall_seconds"] ≈ 97.16
            @test nctssos_source_base_a1_l30_model["peak_rss_bytes"] < 3 * 1024^3

            @test nctssos_source_base_a1_l30_solve["optimizer"] == "Mosek"
            @test nctssos_source_base_a1_l30_solve["profile"] == "A1"
            @test nctssos_source_base_a1_l30_solve["length"] == 30
            @test occursin(
                "source-base",
                nctssos_source_base_a1_l30_solve["comparison_scope"],
            )
            @test nctssos_source_base_a1_l30_solve["qmbcertify_base_construct"] == true
            @test nctssos_source_base_a1_l30_solve["qmbcertify_base_extra"] == 9
            @test nctssos_source_base_a1_l30_solve["reflection_symmetry"] == false
            @test nctssos_source_base_a1_l30_solve["axis_rotation_symmetry"] == false
            @test nctssos_source_base_a1_l30_solve["objective"] ≈ -13.325070927620173
            @test nctssos_source_base_a1_l30_solve["qmbcertify_probe_id"] == "A1_L30"
            @test nctssos_source_base_a1_l30_solve["qmbcertify_objective_total_estimate"] ≈
                a1_l30["objective_total_estimate"]
            @test nctssos_source_base_a1_l30_solve[
                "objective_minus_qmbcertify_total_estimate"
            ] ≈ nctssos_source_base_a1_l30_solve["objective"] -
                 a1_l30["objective_total_estimate"]
            @test nctssos_source_base_a1_l30_solve[
                "objective_minus_qmbcertify_total_estimate"
            ] < 0
            @test abs(
                nctssos_source_base_a1_l30_solve[
                    "objective_minus_qmbcertify_total_estimate"
                ],
            ) < 1e-6
            @test nctssos_source_base_a1_l30_solve["linear_moments"] == 7405
            @test nctssos_source_base_a1_l30_solve["free_keys"] == 7404
            @test nctssos_source_base_a1_l30_solve["linear_state_opt_width"] == 7
            @test nctssos_source_base_a1_l30_solve["linear_state_opt_mode"] == "qmbcertify"
            @test nctssos_source_base_a1_l30_solve["linear_state_opt_row_count"] == 313
            @test nctssos_source_base_a1_l30_solve["contiguous_rdm_k"] == 8
            @test nctssos_source_base_a1_l30_solve["contiguous_rdm_decomposition"] ==
                "qmbcertify"
            @test nctssos_source_base_a1_l30_solve["contiguous_rdm_support"] == "extend"
            @test nctssos_source_base_a1_l30_solve["qmbcertify_rdm_block_sizes"] ==
                [72, 64, 56]
            @test nctssos_source_base_a1_l30_solve["psd_state_opt_width"] == 0
            @test nctssos_source_base_a1_l30_solve["psd_cones"] == 35
            @test nctssos_source_base_a1_l30_solve["dual_variables"] == 169_287
            @test nctssos_source_base_a1_l30_solve["zero_dual_variables"] == 313
            @test nctssos_source_base_a1_l30_solve["max_block"] == 114
            @test nctssos_source_base_a1_l30_solve["solve_time_seconds"] ≈ 139.54671609
            @test nctssos_source_base_a1_l30_solve["remote_wall_seconds"] ≈ 249.77
            @test nctssos_source_base_a1_l30_solve["peak_rss_bytes"] < 5 * 1024^3
            @test nctssos_source_base_a1_l30_solve[
                "sos_certificate_max_abs_residual"
            ] <= 2e-8

            @test nctssos_a2_l20_solve["optimizer"] == "Mosek"
            @test nctssos_a2_l20_solve["length"] == 20
            @test nctssos_a2_l20_solve["objective"] ≈ -9.092707017660949
            @test nctssos_a2_l20_solve["linear_moments"] == 4614
            @test nctssos_a2_l20_solve["free_keys"] == 4613
            @test nctssos_a2_l20_solve["psd_state_opt_width"] == 3
            @test nctssos_a2_l20_solve["psd_state_opt_block_sizes"] == fill(18, 11)
            @test nctssos_a2_l20_solve["psd_cones"] == 50
            @test nctssos_a2_l20_solve["dual_variables"] == 153_454
            @test nctssos_a2_l20_solve["solve_time_seconds"] ≈ 127.170103655
            @test nctssos_a2_l20_solve["remote_wall_seconds"] ≈ 384.37
            @test nctssos_a2_l20_solve["sos_certificate_max_abs_residual"] <= 2e-9

            @test nctssos_a1_l30_solve["optimizer"] == "Mosek"
            @test nctssos_a1_l30_solve["profile"] == "A1"
            @test nctssos_a1_l30_solve["length"] == 30
            @test nctssos_a1_l30_solve["objective"] ≈ -13.664267269723961
            @test nctssos_a1_l30_solve["qmbcertify_probe_id"] == "A1_L30"
            @test nctssos_a1_l30_solve["qmbcertify_objective_total_estimate"] ≈
                a1_l30["objective_total_estimate"]
            @test nctssos_a1_l30_solve["objective_minus_qmbcertify_total_estimate"] ≈
                nctssos_a1_l30_solve["objective"] -
                a1_l30["objective_total_estimate"]
            @test nctssos_a1_l30_solve["objective_minus_qmbcertify_total_estimate"] < -0.3
            @test nctssos_a1_l30_solve["linear_moments"] == 7579
            @test nctssos_a1_l30_solve["free_keys"] == 7578
            @test nctssos_a1_l30_solve["psd_state_opt_width"] == 0
            @test nctssos_a1_l30_solve["psd_state_opt_block_sizes"] == []
            @test nctssos_a1_l30_solve["psd_cones"] == 49
            @test nctssos_a1_l30_solve["dual_variables"] == 224_173
            @test nctssos_a1_l30_solve["solve_time_seconds"] ≈ 542.494610598
            @test nctssos_a1_l30_solve["remote_wall_seconds"] ≈ 1092.76
            @test nctssos_a1_l30_solve["sos_certificate_max_abs_residual"] <= 5e-10

            @test nctssos_su2_pso["objective"] ≈ -3.651093409167063
            @test nctssos_su2_pso["linear_state_opt_width"] == 0
            @test nctssos_su2_pso["linear_state_opt_row_count"] == 0
            @test nctssos_su2_pso["dual_variables"] == 22_854
            @test nctssos_su2_pso["zero_dual_variables"] == 0

            @test nctssos_su2_lso_pso["objective"] ≈ -3.6510934100794663
            @test nctssos_su2_lso_pso["linear_state_opt_width"] == 7
            @test nctssos_su2_lso_pso["linear_state_opt_row_count"] == 819
            @test nctssos_su2_lso_pso["dual_variables"] ==
                nctssos_su2_pso["dual_variables"] +
                nctssos_su2_lso_pso["zero_dual_variables"]
            @test nctssos_su2_lso_pso["zero_dual_variables"] == 819
            @test nctssos_su2_lso_pso["objective"] - nctssos_su2_pso["objective"] ≈
                -9.12403486034362e-10

            @test haskey(qmb, "nctssos_large_construction_probes")
            large_probes = qmb["nctssos_large_construction_probes"]
            @test length(large_probes) == 4
            large_by_id = Dict(probe["id"] => probe for probe in large_probes)
            @test Set(keys(large_by_id)) == Set([
                "NCTSSOS_A1_L20_source_like_construct_only_2026_07_06",
                "NCTSSOS_A2_L20_source_like_construct_only_2026_07_06",
                "NCTSSOS_A1_L30_source_like_construct_only_2026_07_06",
                "NCTSSOS_A2_L30_source_like_construct_only_2026_07_06",
            ])
            large_a1_l20 =
                large_by_id["NCTSSOS_A1_L20_source_like_construct_only_2026_07_06"]
            large_a2_l20 =
                large_by_id["NCTSSOS_A2_L20_source_like_construct_only_2026_07_06"]
            large_a1_l30 =
                large_by_id["NCTSSOS_A1_L30_source_like_construct_only_2026_07_06"]
            large_a2_l30 =
                large_by_id["NCTSSOS_A2_L30_source_like_construct_only_2026_07_06"]

            for probe in values(large_by_id)
                @test probe["status"] == "nctssos_construction_probe"
                @test occursin("not a solved parity row", probe["comparison_scope"])
                @test probe["model_mode"] == "primal"
                @test probe["optimizer"] == "Clarabel"
                @test probe["probe_harness"] == "perf/pauli_translation_solver_probe.jl"
                @test probe["termination_status"] == "not_solved"
                @test probe["construction_only"] == true
                @test probe["model_built"] == false
                @test probe["solved"] == false
                @test !haskey(probe, "objective")
                @test !haskey(probe, "dualization_time_seconds")
                @test !haskey(probe, "solve_time_seconds")
                @test probe["reflection_symmetry"] == true
                @test probe["sign_symmetry"] == true
                @test probe["real_moment_matrix"] == true
                @test probe["u1_symmetry"] == false
                @test probe["su2_symmetry"] == false
                @test probe["base_su2_extend_rdm"] == false
                @test probe["axis_rotation_symmetry"] == true
                @test probe["axis_rotation_equalities"] == true
                @test probe["axis_rotation_quotient"] == true
                @test probe["report_axis_rotation_symmetry"] == true
                @test probe["report_axis_quotient"] == true
                @test probe["solve_supported"] == true
                @test probe["solve_blocker"] == "none"
                @test probe["solve_support_reason"] == ""
                @test probe["solve_unsupported_block_features"] == []
                @test probe["solve_unsupported_zero_features"] == []
                @test probe["order"] == 4
                @test probe["linear_state_opt_width"] == 7
                @test probe["linear_state_opt_row_count"] == 810
                @test probe["contiguous_rdm_k"] == 8
                @test probe["contiguous_rdm_decomposition"] == "qmbcertify"
                @test probe["contiguous_rdm_support"] == "extend"
                @test probe["qmbcertify_rdm_block_sizes"] == [72, 64, 56]
                @test probe["max_block"] == 120
                @test probe["construction_time_seconds"] > 0
                @test probe["script_wall_seconds"] > probe["construction_time_seconds"]
                @test probe["remote_wall_seconds"] > probe["script_wall_seconds"]
                @test probe["peak_rss_bytes"] > 0
                @test probe["peak_rss_bytes"] < 10 * 1024^3
                @test probe["construction_stage_block_assembly_seconds"] >
                    0.95 * probe["construction_time_seconds"]
                env = qmb["nctssos_large_probe_environments"][probe["environment"]]
                @test env["julia_version"] == "1.12.6"
                @test env["thread_count"] == 1
                @test env["host"] == "autodl"
                @test env["timeout_seconds"] == 1500
            end

            @test large_a1_l20["profile"] == "A1"
            @test large_a1_l20["length"] == 20
            @test large_a1_l20["linear_moments"] == 4559
            @test large_a1_l20["free_keys"] == 4558
            @test large_a1_l20["psd_state_opt_width"] == 0
            @test large_a1_l20["psd_state_opt_block_sizes"] == []
            @test large_a1_l20["psd_cones"] == 39
            @test large_a1_l20["product_cache_hit_rate"] ≈ 0.9077155824508321

            @test large_a2_l20["profile"] == "A2"
            @test large_a2_l20["length"] == 20
            @test large_a2_l20["linear_moments"] == 4614
            @test large_a2_l20["free_keys"] == 4613
            @test large_a2_l20["psd_state_opt_width"] == 3
            @test large_a2_l20["psd_state_opt_block_sizes"] == fill(18, 11)
            @test large_a2_l20["psd_cones"] == 50
            @test large_a2_l20["product_cache_hit_rate"] ≈ 0.9077155824508321

            @test large_a1_l30["profile"] == "A1"
            @test large_a1_l30["length"] == 30
            @test large_a1_l30["linear_moments"] == 7579
            @test large_a1_l30["free_keys"] == 7578
            @test large_a1_l30["psd_state_opt_width"] == 0
            @test large_a1_l30["psd_state_opt_block_sizes"] == []
            @test large_a1_l30["psd_cones"] == 49
            @test large_a1_l30["product_cache_hit_rate"] ≈ 0.9365244536940687

            @test large_a2_l30["profile"] == "A2"
            @test large_a2_l30["length"] == 30
            @test large_a2_l30["linear_moments"] == 7684
            @test large_a2_l30["free_keys"] == 7683
            @test large_a2_l30["psd_state_opt_width"] == 3
            @test large_a2_l30["psd_state_opt_block_sizes"] == fill(18, 16)
            @test large_a2_l30["psd_cones"] == 65
            @test large_a2_l30["product_cache_hit_rate"] ≈ 0.9365244536940687

            @test large_a1_l30["remote_wall_seconds"] > large_a1_l20["remote_wall_seconds"]
            @test large_a2_l30["remote_wall_seconds"] > large_a2_l20["remote_wall_seconds"]
            @test large_a1_l30["peak_rss_bytes"] > large_a1_l20["peak_rss_bytes"]
            @test large_a2_l30["peak_rss_bytes"] > large_a2_l20["peak_rss_bytes"]

            @test haskey(qmb, "nctssos_qmbcertify_source_base_model_probes")
            source_base_probes = qmb["nctssos_qmbcertify_source_base_model_probes"]
            @test length(source_base_probes) == 3
            source_base_by_id = Dict(probe["id"] => probe for probe in source_base_probes)
            @test Set(keys(source_base_by_id)) == Set([
                "NCTSSOS_QMBCERTIFY_SOURCE_BASE_A0_L20_sos_dual_clarabel_2026_07_07",
                "NCTSSOS_QMBCERTIFY_SOURCE_BASE_A1_L20_sos_dual_clarabel_2026_07_07",
                "NCTSSOS_QMBCERTIFY_SOURCE_BASE_A2_L20_sos_dual_clarabel_2026_07_07",
            ])
            source_base_a0 = source_base_by_id[
                "NCTSSOS_QMBCERTIFY_SOURCE_BASE_A0_L20_sos_dual_clarabel_2026_07_07"
            ]
            source_base_a1 = source_base_by_id[
                "NCTSSOS_QMBCERTIFY_SOURCE_BASE_A1_L20_sos_dual_clarabel_2026_07_07"
            ]
            source_base_a2 = source_base_by_id[
                "NCTSSOS_QMBCERTIFY_SOURCE_BASE_A2_L20_sos_dual_clarabel_2026_07_07"
            ]

            for probe in values(source_base_by_id)
                @test probe["status"] == "nctssos_source_base_model_probe"
                @test occursin("not a solved parity row", probe["comparison_scope"])
                @test probe["length"] == 20
                @test probe["order"] == 4
                @test probe["model_mode"] == "sos_dual"
                @test probe["sos_hermitian_representation"] == "real_lift"
                @test probe["optimizer"] == "Clarabel"
                @test probe["qmbcertify_base_construct"] == true
                @test probe["qmbcertify_base_extra"] == 9
                @test probe["qmbcertify_base_three_type"] == [1, 1]
                @test probe["reflection_symmetry"] == false
                @test probe["sign_symmetry"] == true
                @test probe["real_moment_matrix"] == true
                @test probe["u1_symmetry"] == false
                @test probe["su2_symmetry"] == false
                @test probe["axis_rotation_symmetry"] == false
                @test probe["axis_rotation_equalities"] == false
                @test probe["axis_rotation_quotient"] == false
                @test probe["termination_status"] == "not_solved"
                @test probe["construction_only"] == false
                @test probe["model_built"] == true
                @test probe["solved"] == false
                @test !haskey(probe, "objective")
                @test !haskey(probe, "solve_time_seconds")
                @test probe["solve_supported"] == true
                @test probe["solve_blocker"] == "none"
                @test probe["dualization_time_seconds"] > 0
                @test probe["remote_wall_seconds"] < 120
                @test probe["peak_rss_gib"] < 3
                @test probe["estimated_model_mode"] == "sos_dual"
                @test probe["estimated_model_variables"] == probe["dual_variables"]
                @test probe["estimated_zero_dual_variables"] == probe["zero_dual_variables"]
                @test probe["estimated_lowering_would_error"] == false
                @test probe["estimated_risk_gate_status"] == "ok"
                @test probe["estimated_model_size_gate_status"] == "ok"
                @test probe["estimated_model_size_gate_estimated_rss_bytes"] ==
                    probe["estimated_dense_schur_bytes"]
                @test probe["estimated_model_size_gate_mem_available_bytes"] > 0
                @test probe["estimated_model_size_gate_max_rss_fraction"] == 0.8
                @test probe["estimated_model_size_gate_reason"] == ""
                @test probe["estimated_scalar_equalities_upper_bound"] ==
                    probe["linear_moments"] - 1
                @test probe["estimated_dense_schur_bytes"] ==
                    8 * probe["estimated_scalar_equalities_upper_bound"]^2
                @test probe["dual_psd_block_count"] == probe["psd_cones"]
                @test probe["max_block"] == 114
                env = qmb["nctssos_solver_probe_environments"][probe["environment"]]
                @test env["thread_count"] == 1
                @test env["host"] == "autodl"
                @test env["timeout_seconds"] == 240
            end

            @test source_base_a0["profile"] == "A0_source_base"
            @test source_base_a0["linear_moments"] == 3420
            @test source_base_a0["free_keys"] == 3391
            @test source_base_a0["linear_zero_constraints"] == 0
            @test source_base_a0["linear_state_opt_width"] == 0
            @test source_base_a0["contiguous_rdm_k"] == 0
            @test source_base_a0["qmbcertify_rdm_block_sizes"] == []
            @test source_base_a0["psd_state_opt_width"] == 0
            @test source_base_a0["psd_state_opt_block_sizes"] == []
            @test source_base_a0["psd_cones"] == 22
            @test source_base_a0["estimated_psd_cone_scalar_variables"] == 106_615
            @test source_base_a0["zero_dual_variables"] == 0

            @test source_base_a1["profile"] == "A1_source_base"
            @test source_base_a1["linear_moments"] == 3534
            @test source_base_a1["free_keys"] == 3505
            @test source_base_a1["linear_zero_constraints"] == 396
            @test source_base_a1["linear_state_opt_width"] == 7
            @test source_base_a1["linear_state_opt_mode"] == "qmbcertify"
            @test source_base_a1["linear_state_opt_row_count"] == 396
            @test source_base_a1["contiguous_rdm_k"] == 8
            @test source_base_a1["contiguous_rdm_decomposition"] == "qmbcertify"
            @test source_base_a1["contiguous_rdm_support"] == "extend"
            @test source_base_a1["qmbcertify_rdm_block_sizes"] == [72, 64, 56]
            @test source_base_a1["psd_state_opt_width"] == 0
            @test source_base_a1["psd_state_opt_block_sizes"] == []
            @test source_base_a1["psd_cones"] == 25
            @test source_base_a1["dual_variables"] == 113_315
            @test source_base_a1["zero_dual_variables"] == 396
            @test source_base_a1["preflight_load_average_1m"] > 25

            @test source_base_a2["profile"] == "A2_source_base"
            @test source_base_a2["linear_moments"] == 5034
            @test source_base_a2["free_keys"] == 5005
            @test source_base_a2["linear_zero_constraints"] == 672
            @test source_base_a2["linear_state_opt_width"] == 7
            @test source_base_a2["linear_state_opt_row_count"] == 672
            @test source_base_a2["contiguous_rdm_k"] == 8
            @test source_base_a2["contiguous_rdm_decomposition"] == "qmbcertify"
            @test source_base_a2["contiguous_rdm_support"] == "extend"
            @test source_base_a2["qmbcertify_rdm_block_sizes"] == [72, 64, 56]
            @test source_base_a2["psd_state_opt_width"] == 3
            @test source_base_a2["psd_state_opt_block_sizes"] == [
                36, 72, 72, 72, 72, 72, 72, 72, 72, 72, 36,
                28, 56, 56, 56, 56, 56, 56, 56, 56, 56, 28,
            ]
            @test source_base_a2["psd_cones"] == 47
            @test source_base_a2["dual_variables"] == 153_751
            @test source_base_a2["zero_dual_variables"] == 672
            @test source_base_a0["dual_variables"] <
                source_base_a1["dual_variables"] <
                source_base_a2["dual_variables"]

            @test haskey(qmb, "nctssos_su2_formulation_probes")
            su2_formulation_probes = qmb["nctssos_su2_formulation_probes"]
            @test length(su2_formulation_probes) == 6
            su2_formulation_by_id =
                Dict(probe["id"] => probe for probe in su2_formulation_probes)
            @test Set(keys(su2_formulation_by_id)) == Set([
                "NCTSSOS_SU2_A2_L8_native_sos_dual_nosolve_zero_row_skip_2026_07_08",
                "NCTSSOS_SU2_A2_L8_primal_moment_variables_nosolve_2026_07_08",
                "NCTSSOS_SU2_A2_L30_primal_moment_variables_nosolve_2026_07_08",
                "NCTSSOS_SU2_A2_L8_primal_moment_variables_solve_2026_07_08",
                "NCTSSOS_SU2_A2_L12_primal_moment_variables_nosolve_2026_07_08",
                "NCTSSOS_SU2_A2_L12_primal_moment_variables_timeout_2026_07_08",
            ])

            for probe in values(su2_formulation_by_id)
                @test probe["status"] == "nctssos_su2_formulation_probe"
                @test probe["profile"] == "A2_SU2_RDM_LSO_PSO"
                @test probe["order"] == 4
                @test probe["optimizer"] == "Mosek"
                @test probe["probe_harness"] == "perf/pauli_translation_solver_probe.jl"
                @test probe["sign_symmetry"] == true
                @test probe["su2_symmetry"] == true
                @test probe["u1_symmetry"] == false
                @test probe["contiguous_rdm_k"] == 8
                @test probe["contiguous_rdm_decomposition"] == "su2"
                @test probe["contiguous_rdm_support"] == "extend"
                @test probe["linear_state_opt_width"] == 7
                @test probe["psd_state_opt_width"] == 3
                @test probe["construction_only"] == false
                @test probe["model_built"] == true
                @test probe["remote_wall_seconds"] > 0
                @test probe["peak_rss_bytes"] > 0
                @test probe["peak_rss_bytes"] < 32 * 1024^3
                @test probe["estimated_scalar_equalities_upper_bound"] >=
                    probe["scalar_equalities"]
                env = qmb["nctssos_solver_probe_environments"][probe["environment"]]
                @test env["host"] == "autodl"
                @test env["thread_count"] == 1
                @test env["mosek_thread_count"] == 1
                @test env["mosek_version"] == "11.2.2"
            end

            su2_native_l8 = su2_formulation_by_id[
                "NCTSSOS_SU2_A2_L8_native_sos_dual_nosolve_zero_row_skip_2026_07_08"
            ]
            @test su2_native_l8["model_mode"] == "sos_dual"
            @test su2_native_l8["sos_hermitian_representation"] == "native"
            @test su2_native_l8["termination_status"] == "not_solved"
            @test su2_native_l8["solved"] == false
            @test su2_native_l8["zero_row_skip_removed_equalities"] ==
                su2_native_l8["estimated_scalar_equalities_upper_bound"] -
                su2_native_l8["scalar_equalities"]
            @test su2_native_l8["zero_row_skip_removed_equalities"] == 56

            su2_primal_l8_nosolve = su2_formulation_by_id[
                "NCTSSOS_SU2_A2_L8_primal_moment_variables_nosolve_2026_07_08"
            ]
            su2_primal_l30_nosolve = su2_formulation_by_id[
                "NCTSSOS_SU2_A2_L30_primal_moment_variables_nosolve_2026_07_08"
            ]
            su2_primal_l8_solve = su2_formulation_by_id[
                "NCTSSOS_SU2_A2_L8_primal_moment_variables_solve_2026_07_08"
            ]
            su2_primal_l12_nosolve = su2_formulation_by_id[
                "NCTSSOS_SU2_A2_L12_primal_moment_variables_nosolve_2026_07_08"
            ]
            su2_primal_l12_timeout = su2_formulation_by_id[
                "NCTSSOS_SU2_A2_L12_primal_moment_variables_timeout_2026_07_08"
            ]

            for probe in (
                su2_primal_l8_nosolve,
                su2_primal_l30_nosolve,
                su2_primal_l8_solve,
                su2_primal_l12_nosolve,
                su2_primal_l12_timeout,
            )
                @test probe["model_mode"] == "primal"
                @test probe["formulation"] == "moment_variables"
                @test probe["representation"] == "real"
                @test probe["scalar_equalities"] == 821
                @test probe["estimated_scalar_equalities_upper_bound"] == 821
                @test probe["estimated_dense_schur_bytes"] == 8 * 821^2
            end

            @test su2_primal_l30_nosolve["length"] == 30
            @test su2_primal_l30_nosolve["model_variables"] == 122_386
            @test su2_primal_l30_nosolve["native_sos_dual_scalar_equalities_reference"] ==
                122_385
            @test su2_primal_l30_nosolve["native_sos_dual_dense_schur_bytes_reference"] >
                10_000 * su2_primal_l30_nosolve["estimated_dense_schur_bytes"]

            @test su2_primal_l8_solve["termination_status"] == "OPTIMAL"
            @test su2_primal_l8_solve["solved"] == true
            @test abs(
                su2_primal_l8_solve["objective"] -
                su2_primal_l8_solve["native_dual_objective_reference"],
            ) <= 1e-6
            @test abs(su2_primal_l8_solve["objective_minus_native_dual_reference"]) <=
                1e-6

            @test su2_primal_l12_nosolve["model_variables"] ==
                su2_primal_l12_timeout["model_variables"]
            @test su2_primal_l12_nosolve["max_block"] == 31
            @test su2_primal_l12_timeout["max_block"] == 31
            @test su2_primal_l12_timeout["termination_status"] == "timeout"
            @test su2_primal_l12_timeout["timeout_exit_status"] == 124
            @test su2_primal_l12_timeout["run_queue_status"] ==
                "deleted_until_evidence_gate"
            @test su2_primal_l12_timeout["restore_gate"] ==
                "post_formulation_change_and_fresh_wall_rss_load_estimate"
            @test su2_primal_l12_timeout["remote_wall_seconds"] >=
                su2_primal_l12_timeout["timeout_seconds"]
            @test su2_primal_l12_timeout["swaps"] == 0
            @test su2_primal_l12_timeout["larger_solved_probes_active"] == false
            @test su2_primal_l12_timeout["deleted_solved_probe_lengths"] ==
                [14, 20, 30]
        end
    end

    rdm_fixture = TestExpectations.expectations_load("expectations/heisenberg_qmbcertify_rdm.toml")
    @test rdm_fixture["qmbcertify_rdm_blocks"]["8"] == [72, 64, 56]
    @test rdm_fixture["qmbcertify_rdm_blocks"]["9"] == [136, 120]
    @test rdm_fixture["qmbcertify_rdm_blocks"]["10"] == [256, 272, 240]
    qmb_rdm8_metrics = rdm_fixture["qmbcertify_rdm_block_metrics"]["rdm8"]
    @test qmb_rdm8_metrics["status"] == "source_derived"
    @test qmb_rdm8_metrics["block_sizes"] == [72, 64, 56]
    @test qmb_rdm8_metrics["n_blocks"] == 3
    @test qmb_rdm8_metrics["max_block"] == 72
    @test qmb_rdm8_metrics["total_block_side"] == 192
    @test qmb_rdm8_metrics["dense_entries"] == 12_416
    @test qmb_rdm8_metrics["symmetric_entries"] == 6_304
    @test qmb_rdm8_metrics["dense_bytes"] == 8 * qmb_rdm8_metrics["dense_entries"]
    @test qmb_rdm8_metrics["symmetric_bytes"] ==
        8 * qmb_rdm8_metrics["symmetric_entries"]
    @test qmb_rdm8_metrics["requires_construction"] == false
    qmb_rdm9_metrics = rdm_fixture["qmbcertify_rdm_block_metrics"]["rdm9"]
    @test qmb_rdm9_metrics["status"] == "source_derived"
    @test qmb_rdm9_metrics["block_sizes"] == [136, 120]
    @test qmb_rdm9_metrics["n_blocks"] == 2
    @test qmb_rdm9_metrics["max_block"] == 136
    @test qmb_rdm9_metrics["total_block_side"] == 256
    @test qmb_rdm9_metrics["dense_entries"] == 32_896
    @test qmb_rdm9_metrics["symmetric_entries"] == 16_576
    @test qmb_rdm9_metrics["dense_bytes"] == 8 * qmb_rdm9_metrics["dense_entries"]
    @test qmb_rdm9_metrics["symmetric_bytes"] ==
        8 * qmb_rdm9_metrics["symmetric_entries"]
    @test qmb_rdm9_metrics["requires_construction"] == false
    qmb_rdm10_metrics = rdm_fixture["qmbcertify_rdm_block_metrics"]["rdm10"]
    @test qmb_rdm10_metrics["status"] == "source_derived"
    @test qmb_rdm10_metrics["block_sizes"] == [256, 272, 240]
    @test qmb_rdm10_metrics["n_blocks"] == 3
    @test qmb_rdm10_metrics["max_block"] == 272
    @test qmb_rdm10_metrics["total_block_side"] == 768
    @test qmb_rdm10_metrics["dense_entries"] == 197_120
    @test qmb_rdm10_metrics["symmetric_entries"] == 98_944
    @test qmb_rdm10_metrics["dense_bytes"] == 8 * qmb_rdm10_metrics["dense_entries"]
    @test qmb_rdm10_metrics["symmetric_bytes"] ==
        8 * qmb_rdm10_metrics["symmetric_entries"]
    @test qmb_rdm10_metrics["requires_construction"] == false
    rdm_structural = rdm_fixture["paper_structural_targets"]["n100_d4_rdm8_full_pso3"]
    _test_structural_target_only(rdm_structural)
    @test rdm_structural["n_sites"] == 100
    @test rdm_structural["order"] == 4
    @test rdm_structural["contiguous_rdm_k"] == [8]
    @test rdm_structural["contiguous_rdm_support"] == "extend"
    @test rdm_structural["psd_state_opt_width"] == 3
    @test rdm_structural["basis_size"] == 12_001
    @test rdm_structural["orbit_basis_size"] == 121
    @test rdm_structural["momentum_sector_count"] == 51
    @test rdm_structural["n_blocks"] == 256
    @test rdm_structural["logical_max_block"] == 256
    @test rdm_structural["psd_max_block"] == 512
    @test rdm_structural["logical_block_histogram"] ==
        [[9, 51], [30, 203], [31, 1], [256, 1]]
    @test rdm_structural["psd_block_histogram"] ==
        [[9, 1], [18, 50], [60, 203], [62, 1], [512, 1]]
    _test_block_domain_histograms(
        rdm_structural,
        [["cyclotomic_float64", 255], ["nothing", 1]],
        [["cyclotomic", 255], ["nothing", 1]],
    )
    @test rdm_structural["psd_symmetric_entries"] == 513_366
    @test rdm_structural["psd_dense_bytes"] ==
        8 * rdm_structural["psd_dense_entries"]
    @test rdm_structural["psd_symmetric_bytes"] ==
        8 * rdm_structural["psd_symmetric_entries"]
    @test rdm_structural["rdm_logical_block_sizes"] == [256]
    @test rdm_structural["rdm_psd_block_sizes"] == [512]
    @test rdm_structural["qmbcertify_rdm_reference_block_sizes"] == [72, 64, 56]
    @test rdm_structural["qmbcertify_rdm_reference_n_blocks"] == 3
    @test rdm_structural["qmbcertify_rdm_reference_max_block"] == 72
    @test rdm_structural["qmbcertify_rdm_reference_total_block_side"] == 192
    @test rdm_structural["qmbcertify_rdm_reference_dense_entries"] == 12_416
    @test rdm_structural["qmbcertify_rdm_reference_symmetric_entries"] == 6_304
    @test rdm_structural["qmbcertify_rdm_reference_dense_bytes"] == 99_328
    @test rdm_structural["qmbcertify_rdm_reference_symmetric_bytes"] == 50_432
    @test !rdm_structural["qmbcertify_rdm_reference_requires_construction"]
    @test rdm_structural["psd_state_opt_candidate_count"] == 9
    @test rdm_structural["psd_state_opt_block_count"] == 51
    @test rdm_structural["psd_state_opt_psd_block_histogram"] == [[9, 1], [18, 50]]
    @test rdm_structural["contiguous_rdm_zero_row_count"] == 0
    @test rdm_structural["add_on_zero_row_count"] == 0
    @test rdm_structural["known_zero_constraint_feature_histogram"] == []

    qmb_structural = rdm_fixture["paper_structural_targets"]["n100_d4_rdm8_qmbcertify_pso3"]
    _test_structural_target_only(qmb_structural)
    @test qmb_structural["contiguous_rdm_decomposition"] == "qmbcertify"
    @test qmb_structural["contiguous_rdm_support"] == "extend"
    @test qmb_structural["n_blocks"] == 258
    @test qmb_structural["logical_max_block"] == 72
    @test qmb_structural["psd_max_block"] == 72
    @test qmb_structural["psd_dense_entries"] == 763_341
    @test qmb_structural["psd_symmetric_entries"] == 388_342
    @test qmb_structural["psd_dense_bytes"] ==
        8 * qmb_structural["psd_dense_entries"]
    @test qmb_structural["psd_symmetric_bytes"] ==
        8 * qmb_structural["psd_symmetric_entries"]
    _test_block_domain_histograms(
        qmb_structural,
        [["cyclotomic_float64", 255], ["nothing", 3]],
        [["cyclotomic", 255], ["nothing", 3]],
    )
    @test qmb_structural["rdm_logical_block_sizes"] == [72, 64, 56]
    @test qmb_structural["rdm_psd_block_sizes"] == [72, 64, 56]
    @test qmb_structural["qmbcertify_rdm_reference_block_sizes"] == [72, 64, 56]
    @test qmb_structural["qmbcertify_rdm_reference_n_blocks"] == 3
    @test qmb_structural["qmbcertify_rdm_reference_max_block"] == 72
    @test qmb_structural["qmbcertify_rdm_reference_total_block_side"] == 192
    @test qmb_structural["qmbcertify_rdm_reference_dense_entries"] == 12_416
    @test qmb_structural["qmbcertify_rdm_reference_symmetric_entries"] == 6_304
    @test qmb_structural["qmbcertify_rdm_reference_dense_bytes"] == 99_328
    @test qmb_structural["qmbcertify_rdm_reference_symmetric_bytes"] == 50_432
    @test !qmb_structural["qmbcertify_rdm_reference_requires_construction"]
    @test qmb_structural["contiguous_rdm_zero_row_count"] == 0
    @test qmb_structural["add_on_zero_row_count"] == 0
    @test qmb_structural["known_zero_constraint_feature_histogram"] == []

    base_fixture = TestExpectations.expectations_load("expectations/heisenberg_qmbcertify_base.toml")
    structural = base_fixture["paper_structural_targets"]["n100_d4"]
    _test_structural_target_only(structural)
    @test structural["n_sites"] == 100
    @test structural["order"] == 4
    @test structural["basis_size"] == 12_001
    @test structural["orbit_basis_size"] == 121
    @test structural["momentum_sector_count"] == 51
    @test structural["n_blocks"] == 204
    @test structural["logical_max_block"] == 31
    @test structural["psd_max_block"] == 62
    _test_block_domain_histograms(
        structural,
        [["cyclotomic_float64", 204]],
        [["cyclotomic", 204]],
    )
    @test structural["product_cache_hit_rate"] > 0.98
    @test structural["max_contiguous_rdm_k"] == 8
    @test structural["max_linear_state_opt_width"] == 7
    @test structural["max_psd_state_opt_width"] == 3
    @test structural["axis_rotation_equality_row_count"] == 0
    @test structural["moment_equality_row_count"] == 0
    @test structural["linear_state_opt_row_count"] == 0
    @test structural["add_on_zero_row_count"] == 0
    @test structural["known_zero_constraint_feature_histogram"] == []

    reflection_structural = base_fixture["paper_structural_targets"]["n100_d4_reflection"]
    _test_structural_target_only(reflection_structural)
    @test reflection_structural["reflection_symmetry"] == true
    @test reflection_structural["n_sites"] == 100
    @test reflection_structural["order"] == 4
    @test reflection_structural["basis_size"] == 12_001
    @test reflection_structural["orbit_basis_size"] == 121
    @test reflection_structural["momentum_sector_count"] == 51
    @test reflection_structural["n_blocks"] == 408
    @test reflection_structural["logical_max_block"] == 30
    @test reflection_structural["psd_max_block"] == 44
    _test_block_domain_histograms(
        reflection_structural,
        [["cyclotomic_float64", 408]],
        [["cyclotomic_sqrt_rational", 408]],
    )
    @test reflection_structural["psd_symmetric_entries"] == 190_191
    @test reflection_structural["add_on_zero_row_count"] == 0
    @test reflection_structural["known_zero_constraint_feature_histogram"] == []

    reflection_complex = base_fixture["paper_structural_targets"]["n100_d4_reflection_complex_k0"]
    _test_structural_target_only(reflection_complex)
    @test reflection_complex["reflection_symmetry"] == true
    @test reflection_complex["real_moment_matrix"] == false
    @test reflection_complex["momentum_sectors"] == [0]
    @test reflection_complex["n_sites"] == 100
    @test reflection_complex["order"] == 4
    @test reflection_complex["n_blocks"] == 8
    @test reflection_complex["logical_max_block"] == 22
    @test reflection_complex["psd_max_block"] == 22
    _test_block_domain_histograms(
        reflection_complex,
        [["cyclotomic_float64", 8]],
        [["cyclotomic_sqrt_rational", 8]],
    )
    @test reflection_complex["logical_block_histogram"] ==
        reflection_complex["psd_block_histogram"]
    @test reflection_complex["psd_dense_entries"] == 1_939
    @test reflection_complex["psd_symmetric_entries"] == 1_939
    @test reflection_complex["psd_dense_bytes"] == 15_512
    @test reflection_complex["psd_symmetric_bytes"] == 15_512
    @test reflection_complex["product_cache_hit_rate"] == 0.0

    su2_full = base_fixture["paper_structural_targets"]["n100_d4_su2_full_basis"]
    _test_structural_target_only(su2_full)
    @test su2_full["n_sites"] == 100
    @test su2_full["order"] == 4
    @test su2_full["basis_size"] == 12_001
    @test su2_full["support_counts"] == [[0, 1], [1, 100], [2, 100], [3, 100], [4, 100]]
    @test su2_full["singlet_channel_count"] == 501
    @test su2_full["singlet_channel_support_counts"] ==
        [[0, 1], [2, 100], [3, 100], [4, 300]]
    @test su2_full["singlet_channel_equality_row_count"] ==
        su2_full["basis_size"] - su2_full["singlet_channel_count"]
    @test su2_full["reduced_block_sizes"] == [501, 1_100, 900, 400, 100]
    @test su2_full["transformed_block_sizes"] == [501, 3_300, 4_500, 2_800, 900]
    @test su2_full["real_psd_block_sizes"] == su2_full["reduced_block_sizes"]
    @test su2_full["reduced_max_block"] == 1_100
    @test su2_full["transformed_max_block"] == 4_500
    @test su2_full["real_psd_max_block"] == su2_full["reduced_max_block"]
    @test su2_full["transformed_total_block_side"] == su2_full["basis_size"]
    @test su2_full["real_psd_total_block_side"] ==
        su2_full["reduced_total_block_side"]
    _test_block_domain_histograms(
        su2_full,
        [["complex_algebraic_float64", 5]],
        [["complex_sqrt_rational", 5]],
    )
    @test su2_full["reduced_dense_entries"] == 2_441_001
    @test su2_full["reduced_symmetric_entries"] == 1_222_001
    @test su2_full["real_psd_dense_entries"] == su2_full["reduced_dense_entries"]
    @test su2_full["real_psd_symmetric_entries"] ==
        su2_full["reduced_symmetric_entries"]
    @test su2_full["full_dense_bytes"] == 8 * su2_full["full_dense_entries"]
    @test su2_full["full_symmetric_bytes"] == 8 * su2_full["full_symmetric_entries"]
    @test su2_full["reduced_dense_bytes"] == 8 * su2_full["reduced_dense_entries"]
    @test su2_full["reduced_symmetric_bytes"] ==
        8 * su2_full["reduced_symmetric_entries"]
    @test su2_full["real_psd_dense_bytes"] ==
        8 * su2_full["real_psd_dense_entries"]
    @test su2_full["real_psd_symmetric_bytes"] ==
        8 * su2_full["real_psd_symmetric_entries"]
    @test su2_full["offblock_entry_count"] == 134_883_000
    @test su2_full["copy_entry_count"] == 6_700_000
    @test su2_full["accounted_entry_count"] ==
        su2_full["full_dense_entries"]
    @test su2_full["translation_combined"] == false

    su2_translation = base_fixture["paper_structural_targets"]["n100_d4_su2_translation_orbit"]
    _test_structural_target_only(su2_translation)
    @test su2_translation["orbit_basis_size"] == 121
    @test su2_translation["singlet_channel_count"] == 6
    @test su2_translation["singlet_channel_support_counts"] ==
        [[0, 1], [2, 1], [3, 1], [4, 3]]
    @test su2_translation["singlet_channel_equality_row_count"] ==
        su2_translation["orbit_basis_size"] - su2_translation["singlet_channel_count"]
    @test su2_translation["momentum_sector_count"] == 51
    @test su2_translation["zero_momentum_blocks"] ==
        [[0, 6], [2, 11], [4, 9], [6, 4], [8, 1]]
    @test su2_translation["nonzero_momentum_blocks"] ==
        [[0, 5], [2, 11], [4, 9], [6, 4], [8, 1]]
    @test su2_translation["n_blocks"] == 255
    @test su2_translation["logical_max_block"] == 11
    @test su2_translation["psd_max_block"] == 22
    @test su2_translation["logical_block_histogram"] ==
        [[1, 51], [4, 51], [5, 50], [6, 1], [9, 51], [11, 51]]
    @test su2_translation["psd_block_histogram"] ==
        [[2, 51], [8, 51], [10, 50], [12, 1], [18, 51], [22, 51]]
    _test_block_domain_histograms(
        su2_translation,
        [["complex_algebraic_float64", 255]],
        [["complex_sqrt_rational", 255]],
    )
    @test su2_translation["psd_symmetric_entries"] == 26_441
    @test su2_translation["psd_dense_bytes"] ==
        8 * su2_translation["psd_dense_entries"]
    @test su2_translation["psd_symmetric_bytes"] ==
        8 * su2_translation["psd_symmetric_entries"]
    @test su2_translation["su2_full_dense_entries"] == 734_641
    @test su2_translation["su2_active_dense_entries"] == 46_625
    @test su2_translation["su2_reduced_dense_entries"] == 12_455
    @test su2_translation["offblock_entry_count"] == 688_016
    @test su2_translation["copy_entry_count"] == 34_170
    @test su2_translation["wigner_offblock_zero_row_budget"] ==
        2 * su2_translation["offblock_entry_count"]
    @test su2_translation["wigner_magnetic_copy_zero_row_budget"] ==
        su2_translation["copy_entry_count"]
    @test su2_translation["wigner_zero_row_budget"] ==
        su2_translation["wigner_offblock_zero_row_budget"] +
        su2_translation["wigner_magnetic_copy_zero_row_budget"]
    @test su2_translation["accounted_entry_count"] ==
        su2_translation["su2_full_dense_entries"]
    @test su2_translation["translation_combined"] == true
    @test su2_translation["reflection_combined"] == false

    su2_reflection = base_fixture["paper_structural_targets"]["n100_d4_su2_translation_reflection"]
    _test_structural_target_only(su2_reflection)
    @test su2_reflection["orbit_basis_size"] == 121
    @test su2_reflection["singlet_channel_count"] == 6
    @test su2_reflection["singlet_channel_support_counts"] ==
        [[0, 1], [2, 1], [3, 1], [4, 3]]
    @test su2_reflection["singlet_channel_equality_row_count"] ==
        su2_reflection["orbit_basis_size"] - su2_reflection["singlet_channel_count"]
    @test su2_reflection["momentum_sector_count"] == 51
    @test su2_reflection["n_blocks"] == 507
    @test su2_reflection["logical_max_block"] == 11
    @test su2_reflection["psd_max_block"] == 16
    @test su2_reflection["logical_block_histogram"] == [
        [1, 102], [2, 2], [3, 4], [4, 98], [5, 101],
        [6, 3], [8, 1], [9, 98], [11, 98],
    ]
    @test su2_reflection["psd_block_histogram"] == [
        [1, 98], [2, 4], [4, 100], [5, 98], [6, 4],
        [9, 98], [10, 3], [11, 98], [12, 3], [16, 1],
    ]
    _test_block_domain_histograms(
        su2_reflection,
        [["complex_algebraic_float64", 507]],
        [["complex_sqrt_rational", 507]],
    )
    @test su2_reflection["psd_dense_entries"] == 25_092
    @test su2_reflection["psd_symmetric_entries"] == 14_077
    @test su2_reflection["psd_dense_bytes"] ==
        8 * su2_reflection["psd_dense_entries"]
    @test su2_reflection["psd_symmetric_bytes"] ==
        8 * su2_reflection["psd_symmetric_entries"]
    @test su2_reflection["su2_full_dense_entries"] == 1_426_033
    @test su2_reflection["su2_active_dense_entries"] == 90_619
    @test su2_reflection["su2_reduced_dense_entries"] == 24_207
    @test su2_reflection["offblock_entry_count"] == 1_335_414
    @test su2_reflection["copy_entry_count"] == 66_412
    @test su2_reflection["wigner_offblock_zero_row_budget"] ==
        2 * su2_reflection["offblock_entry_count"]
    @test su2_reflection["wigner_magnetic_copy_zero_row_budget"] ==
        su2_reflection["copy_entry_count"]
    @test su2_reflection["wigner_zero_row_budget"] ==
        su2_reflection["wigner_offblock_zero_row_budget"] +
        su2_reflection["wigner_magnetic_copy_zero_row_budget"]
    @test su2_reflection["accounted_entry_count"] ==
        su2_reflection["su2_full_dense_entries"]
    @test su2_reflection["translation_combined"] == true
    @test su2_reflection["reflection_combined"] == true
    @test su2_reflection["solve_blocker"] == su2_translation["solve_blocker"]
    @test su2_reflection["solve_blocker_reason"] == su2_translation["solve_blocker_reason"]
end

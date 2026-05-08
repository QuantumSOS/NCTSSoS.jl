using Test
using Random
using NCTSSoS

const V2RDM_STRUCTURED_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const H2_NK2_INTEGRALS = joinpath(V2RDM_STRUCTURED_ROOT, "test", "data", "assets", "h2_chain_nk2_active_2e4o_integrals.txt")
const H4_NK2_TEST_INTEGRALS = joinpath(V2RDM_STRUCTURED_ROOT, "test", "data", "assets", "h4_chain_nk2_integrals.txt")
const H4_NK2_PHYSICAL_INTEGRALS = joinpath(V2RDM_STRUCTURED_ROOT, "output", "h4_chain_nk2_figure1_integrals_drop1e-12.txt")
const H4_NK2_INTEGRALS = isfile(H4_NK2_PHYSICAL_INTEGRALS) ? H4_NK2_PHYSICAL_INTEGRALS : H4_NK2_TEST_INTEGRALS
const H4_NK2_REFERENCE_DIR = joinpath(V2RDM_STRUCTURED_ROOT, "output", "h4_nk_sweep_physical_cpu_1h", "nk2")

module V2RDMStructuredSymbolicOracle end
Base.include(V2RDMStructuredSymbolicOracle,
    joinpath(@__DIR__, "..", "..", "probes", "h4_nk_periodic_moment_sos_export_libsdp.jl"))

module V2RDMStructuredDriver end
Base.include(V2RDMStructuredDriver,
    joinpath(@__DIR__, "..", "..", "probes", "h4_nk_periodic_structured_export.jl"))

function _symbolic_h2_linear()
    opts = V2RDMStructuredSymbolicOracle.Options(
        2,      # nk
        4,      # norb per k
        2,      # nelec per cell
        H2_NK2_INTEGRALS,
        :momentum,
        false,  # include_one_d
        false,  # spin_resolved_trace
        false,  # singlet_s2
        mktempdir(),
        "h2_nk2_symbolic_oracle",
        :min,
        32,
    )
    data = V2RDMStructuredSymbolicOracle.build_h4_pqg_moment_problem(opts)
    return data.moment_problem.linear
end

function _structured_h2_linear()
    h1e, eri = V2RDMStructuredDriver.load_integrals_txt(H2_NK2_INTEGRALS; nk = 2, norb = 4)
    return build_pqg_moment_data(
        h1e,
        eri;
        nk = 2,
        norb = 4,
        nelec_per_cell = 2,
        blocking = :momentum,
        spin_resolved_trace = false,
        singlet_s2 = false,
        include_one_d = false,
    )
end

_key_string(key) = join(Int.(key), ",")

function _form_map(form)
    out = Dict{String,ComplexF64}()
    for (key, coef) in form
        out[_key_string(key)] = ComplexF64(coef)
    end
    return out
end

function _test_form_equal(a, b; atol = 1e-12)
    ma = _form_map(a)
    mb = _form_map(b)
    @test keys(ma) == keys(mb)
    for key in keys(ma)
        @test isapprox(ma[key], mb[key]; atol, rtol = 0)
    end
end

function _parse_datc_counts(path::AbstractString)
    data_lines = String[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "*") || startswith(line, "\"")) && continue
        push!(data_lines, line)
    end
    length(data_lines) >= 4 || error("not enough data lines in $path")
    n_constraints = parse(Int, data_lines[1])
    n_blocks = parse(Int, data_lines[2])
    block_dims = parse.(Int, split(data_lines[3]))
    entries = data_lines[5:end]
    n_objective_entries = count(line -> first(split(line)) == "0", entries)
    n_constraint_entries = length(entries) - n_objective_entries
    return (; n_blocks, block_dims, n_constraints, n_objective_entries, n_constraint_entries)
end

function _find_reference_datc()
    isdir(H4_NK2_REFERENCE_DIR) || return nothing
    candidates = String[]
    for (root, _, files) in walkdir(H4_NK2_REFERENCE_DIR)
        for file in files
            endswith(file, ".dat-c") && push!(candidates, joinpath(root, file))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates)
    return first(candidates)
end

@testset "structured PQG V2RDM assembler" begin
    @testset "H2/Nk=2 canonical moment enumeration" begin
        symbolic = _symbolic_h2_linear()
        structured = _structured_h2_linear()

        @test length(structured.moments) == length(symbolic.moments)
        @test length(structured.psd_blocks_lin) == length(symbolic.psd_blocks_lin)
        @test [b.size for b in structured.psd_blocks_lin] == [b.size for b in symbolic.psd_blocks_lin]
        @test length(structured.free_keys) == length(symbolic.free_keys)
        _test_form_equal(structured.objective_lin, symbolic.objective_lin; atol = 1e-10)
        @test length(structured.zero_constraints) == length(symbolic.zero_constraints)
        for idx in eachindex(symbolic.zero_constraints)
            _test_form_equal(structured.zero_constraints[idx].form, symbolic.zero_constraints[idx].form; atol = 1e-10)
        end
    end

    @testset "H2/Nk=2 random PQG block entries match symbolic path" begin
        symbolic = _symbolic_h2_linear()
        structured = _structured_h2_linear()
        rng = MersenneTwister(0x5eed)

        positions = Tuple{Int,Int,Int}[]
        for (b, block) in enumerate(symbolic.psd_blocks_lin)
            for i in 1:block.size, j in 1:block.size
                push!(positions, (b, i, j))
            end
        end

        for _ in 1:200
            b, i, j = positions[rand(rng, eachindex(positions))]
            _test_form_equal(structured.psd_blocks_lin[b].entries[i, j], symbolic.psd_blocks_lin[b].entries[i, j])
        end
    end

    @testset "H2/Nk=2 direct MomentLinearData dat-c export smoke" begin
        linear = _structured_h2_linear()
        outdir = mktempdir()
        summary = V2RDMStructuredDriver.export_libsdp(linear; outdir, basename = "h2_nk2_structured", sense = :min)
        counts = _parse_datc_counts(joinpath(outdir, "h2_nk2_structured.dat-c"))

        @test counts.n_blocks == summary.n_blocks
        @test counts.block_dims == summary.block_dims
        @test counts.n_constraints == summary.n_constraints
        @test counts.n_objective_entries == summary.n_objective_entries
        @test counts.n_constraint_entries == summary.n_constraint_entries
        @test counts.n_blocks > 0
        @test counts.n_constraints > 0
        @test counts.n_objective_entries > 0
    end

    @testset "H4/Nk=2 structured dat-c count diff against archived symbolic export" begin
        ref = _find_reference_datc()
        if ref === nothing
            @test_skip "archived H4/Nk=2 reference .dat-c not present under $H4_NK2_REFERENCE_DIR"
        else
            h1e, eri = V2RDMStructuredDriver.load_integrals_txt(H4_NK2_INTEGRALS; nk = 2, norb = 8)
            linear = build_pqg_moment_data(
                h1e,
                eri;
                nk = 2,
                norb = 8,
                nelec_per_cell = 4,
                blocking = :momentum,
                spin_resolved_trace = true,
                singlet_s2 = true,
                include_one_d = false,
            )
            outdir = mktempdir()
            summary = V2RDMStructuredDriver.export_libsdp(linear; outdir, basename = "h4_nk2_structured", sense = :min, legacy_zero_split = true)
            got = _parse_datc_counts(joinpath(outdir, "h4_nk2_structured.dat-c"))
            want = _parse_datc_counts(ref)

            @test got.n_blocks == want.n_blocks
            @test got.block_dims == want.block_dims
            @test got.n_constraints == want.n_constraints
            # The direct MomentLinearData exporter reuses cached pivots, while
            # the archived symbolic exporter rediscovered pivots during export.
            # That can change only the sparsity of F₀, not the SDP layout or
            # equality system. Keep this as a tight sanity check, not byte-count
            # cosplay.
            @test got.n_objective_entries > 0
            @test isapprox(got.n_objective_entries, want.n_objective_entries; rtol = 0.01)
            @test got.n_constraint_entries == want.n_constraint_entries
            @test summary.n_blocks == want.n_blocks
        end
    end

    @testset "H4/Nk=2 structured BPSDP objective smoke" begin
        if get(ENV, "NCTSSOS_RUN_V2RDM_STRUCTURED_BPSDP", "0") != "1"
            @test_skip "set NCTSSOS_RUN_V2RDM_STRUCTURED_BPSDP=1 on HAI to solve the structured H4/Nk=2 .dat-c with BPSDP"
        else
            integrals = get(ENV, "NCTSSOS_V2RDM_H4_NK2_INTEGRALS", H4_NK2_INTEGRALS)
            expected_active = parse(Float64, get(ENV, "NCTSSOS_V2RDM_H4_NK2_ACTIVE_OBJECTIVE", "-39.300912225773"))
            expected_total = parse(Float64, get(ENV, "NCTSSOS_V2RDM_H4_NK2_TOTAL_OBJECTIVE", "-2.188110787316"))
            energy_shift = expected_total - expected_active

            h1e, eri = V2RDMStructuredDriver.load_integrals_txt(integrals; nk = 2, norb = 8)
            linear = build_pqg_moment_data(
                h1e,
                eri;
                nk = 2,
                norb = 8,
                nelec_per_cell = 4,
                blocking = :momentum,
                spin_resolved_trace = true,
                singlet_s2 = true,
                include_one_d = false,
            )
            outdir = mktempdir()
            V2RDMStructuredDriver.export_libsdp(linear; outdir, basename = "h4_nk2_structured", sense = :min, legacy_zero_split = true)
            datc = joinpath(outdir, "h4_nk2_structured.dat-c")
            summary_path = joinpath(outdir, "bpsdp_summary.txt")
            solver = joinpath(V2RDM_STRUCTURED_ROOT, "probes", "solve_datc_bpsdp_primitive.jl")
            max_seconds = get(ENV, "NCTSSOS_V2RDM_BPSDP_MAX_SECONDS", "2400")

            project_dir = Base.active_project() === nothing ? V2RDM_STRUCTURED_ROOT : dirname(Base.active_project())
            cmd = `$(Base.julia_cmd()) --project=$project_dir $solver --problem=$datc --backend=cpu --max-seconds=$max_seconds --summary=$summary_path --energy-shift=$energy_shift`
            run(cmd)

            summary = Dict{String,String}()
            for line in eachline(summary_path)
                key, value = split(line, " = ", limit = 2)
                summary[key] = value
            end
            @test isapprox(parse(Float64, summary["primal_objective"]), expected_active; atol = 5e-4, rtol = 0)
            @test isapprox(parse(Float64, summary["recovered_total_primal"]), expected_total; atol = 5e-4, rtol = 0)
        end
    end
end

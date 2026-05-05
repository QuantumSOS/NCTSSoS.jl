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

    function _bridge_inspection_backend(model)
        for candidate in (backend(model), unsafe_backend(model))
            try
                MOI.get(candidate, MOI.ListOfConstraintTypesPresent())
                return candidate
            catch err
                err isa MOI.GetAttributeNotAllowed || rethrow()
            end
        end
        error("No backend layer exposes MOI.ListOfConstraintTypesPresent()")
    end

    function _constraint_type_counts(backend)
        types = MOI.get(backend, MOI.ListOfConstraintTypesPresent())
        counts = Dict{Tuple{DataType,DataType},Int}()
        for (F, S) in types
            counts[(F, S)] = MOI.get(backend, MOI.NumberOfConstraints{F,S}())
        end
        return types, counts
    end

    function _count_psd_1x1(backend, types)
        total = 0
        for (F, S) in types
            S <: MOI.PositiveSemidefiniteConeTriangle || continue
            for ci in MOI.get(backend, MOI.ListOfConstraintIndices{F,S}())
                set = MOI.get(backend, MOI.ConstraintSet(), ci)
                MOI.dimension(set) == 1 && (total += 1)
            end
        end
        return total
    end

    @testset "H2/Nk=2 :psd_blocks delivers BPSDP-favorable cones" begin
        mp = _h2_nk2_moment_problem()
        model, _extract = build_jump_model(mp;
            formulation=:psd_blocks,
            representation=:complex,
            orphan_policy=:aux_psd_free,
        )
        set_optimizer(model, BPSDP.Optimizer)
        MOI.Utilities.attach_optimizer(model)

        inspected_backend = _bridge_inspection_backend(model)
        types, counts = _constraint_type_counts(inspected_backend)

        n_hermitian = 0
        for ((_, S), count) in counts
            S <: MOI.HermitianPositiveSemidefiniteConeTriangle && (n_hermitian += count)
        end
        n_psd_1x1 = _count_psd_1x1(inspected_backend, types)
        n_symbolic_hpsd = count(c -> c[1] == :HPSD, mp.constraints)
        n_orphans = length(NCTSSoS.orphan_keys(mp))
        n_aux_blocks = cld(n_orphans, NCTSSoS.AUX_ORPHANS_PER_BLOCK)

        @test n_hermitian == n_symbolic_hpsd + n_aux_blocks
        @test n_psd_1x1 <= n_orphans
        @test n_psd_1x1 < 100
    end
end

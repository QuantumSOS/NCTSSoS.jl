using Test
using LinearAlgebra
using Random
using NCTSSoS
using JuMP
import Clarabel

function _right_total_generator(rows, local_generator, support_size)
    result = zeros(eltype(rows), size(rows))
    radix = 3
    for site in 1:support_size
        stride = radix^(site - 1)
        for code0 in 0:(radix^support_size - 1)
            input_digit = (code0 ÷ stride) % radix + 1
            for output_digit in 1:radix
                coefficient = local_generator[input_digit, output_digit]
                iszero(coefficient) && continue
                output_code0 = code0 + (output_digit - input_digit) * stride
                @views result[:, output_code0 + 1] .+=
                    coefficient .* rows[:, code0 + 1]
            end
        end
    end
    return result
end

function _form_value(form, values)
    return sum(
        coefficient * values[key]
        for (key, coefficient) in form;
        init=0.0 + 0.0im,
    )
end

@testset "SU(2)-invariant moment quotient" begin
    @testset "thin singlet rows" begin
        expected = [1, 0, 1, 1, 3, 6, 15, 36, 91]
        for support_size in 0:8
            rows = NCTSSoS._pauli_su2_word_singlet_rows(support_size)
            @test size(rows) == (expected[support_size + 1], 3^support_size)
            dense = NCTSSoS._dense_sparse_transform_rows(rows)
            @test dense * dense' ≈ I atol = 1e-11
            if support_size <= 4
                reference = pauli_su2_word_singlet_channels(support_size).transform
                @test dense' * dense ≈ reference' * reference atol = 1e-10
            end
            n_sites = max(1, 2 * support_size + 1)
            support = Tuple(1:support_size)
            patterns = NCTSSoS._pauli_su2_stabilizer_axis_patterns(support, n_sites)
            @test NCTSSoS._pauli_su2_exact_projected_singlet_rank(
                support,
                n_sites,
                patterns,
            ) == expected[support_size + 1]
        end

        local_transform = NCTSSoS._pauli_su2_word_local_spherical_transform()
        jz_spherical = Matrix(Diagonal(ComplexF64[-1, 0, 1]))
        jplus_spherical = zeros(ComplexF64, 3, 3)
        jplus_spherical[2, 1] = sqrt(2)
        jplus_spherical[3, 2] = sqrt(2)
        jz_cartesian = local_transform' * jz_spherical * local_transform
        jplus_cartesian = local_transform' * jplus_spherical * local_transform
        for support_size in 1:8
            rows = NCTSSoS._pauli_su2_word_singlet_rows(support_size)
            dense = NCTSSoS._dense_sparse_transform_rows(rows)
            @test maximum(abs, _right_total_generator(
                dense,
                jz_cartesian,
                support_size,
            ); init=0.0) <= 1e-10
            @test maximum(abs, _right_total_generator(
                dense,
                jplus_cartesian,
                support_size,
            ); init=0.0) <= 1e-10
        end

        support8 = NCTSSoS._pauli_su2_word_singlet_rows(8)
        @test sum(length, support8.rows) <= 91 * 3^8
        @test_throws DomainError NCTSSoS._pauli_su2_word_singlet_rows(-1)
        @test_throws DomainError NCTSSoS._pauli_su2_word_singlet_rows(2; atol=-1)
        @test_throws ArgumentError NCTSSoS._pauli_su2_word_singlet_rows(9)
        @test isdefined(NCTSSoS, :_PAULI_SU2_SPIN1_CG_CACHE_LOCK)
        @test NCTSSoS._PAULI_SU2_SPIN1_CG_CACHE_LOCK isa ReentrantLock
    end

    @testset "deterministic singlet pivot charts" begin
        dependent_rows = ComplexF64[
            1 0 1
            2 0 2
            0 1 1
        ]
        independent_rows = NCTSSoS._pauli_su2_select_expected_singlet_rows(
            dependent_rows,
            2;
            atol=1e-12,
        )
        @test size(independent_rows) == (2, 3)
        row_coefficients = adjoint(independent_rows) \ adjoint(dependent_rows)
        @test adjoint(row_coefficients) * independent_rows ≈ dependent_rows atol = 1e-12

        for support_size in (0, 2, 3, 4, 8)
            singlet = NCTSSoS._dense_sparse_transform_rows(
                NCTSSoS._pauli_su2_word_singlet_rows(support_size),
            )
            pivots, reconstruction, chart_condition =
                NCTSSoS._pauli_su2_select_pivot_columns(
                    singlet;
                    condition_limit=1e10,
                )
            channel_count = size(singlet, 1)
            @test length(pivots) == channel_count
            @test issorted(pivots)
            @test size(reconstruction) == (3^support_size, channel_count)
            @test isfinite(chart_condition)
            if !isempty(pivots)
                chart = adjoint(singlet)[pivots, :]
                @test reconstruction[pivots, :] ≈ I atol = 1e-11
                @test reconstruction * chart ≈ adjoint(singlet) atol = 1e-10
                @test reconstruction ≈
                    adjoint(singlet) * (singlet * reconstruction) atol = 1e-10
            end
            repeated = NCTSSoS._pauli_su2_select_pivot_columns(
                singlet;
                condition_limit=1e10,
            )
            @test repeated[1] == pivots
            @test repeated[2] ≈ reconstruction atol = 0 rtol = 0
        end

        singlet2 = NCTSSoS._dense_sparse_transform_rows(
            NCTSSoS._pauli_su2_word_singlet_rows(2),
        )
        @test_throws ArgumentError NCTSSoS._pauli_su2_select_pivot_columns(
            singlet2;
            condition_limit=0.5,
        )
        @test_throws DomainError NCTSSoS._pauli_su2_select_pivot_columns(
            singlet2;
            condition_limit=0,
        )
        singular_chart = ComplexF64[
            1 0
            2 0
        ]
        @test_throws ArgumentError NCTSSoS._pauli_su2_select_pivot_columns(
            singular_chart;
            condition_limit=1e10,
        )

        swapped_patterns = NCTSSoS._pauli_su2_stabilizer_axis_patterns((1, 4), 6)
        expected_patterns = Set(
            (UInt8(left), UInt8(right))
            for left in 0:2 for right in left:2
        )
        @test swapped_patterns == expected_patterns
    end

    @testset "translation moment quotient descriptor" begin
        registry, ops = create_pauli_variables(1:6)
        hamiltonian = heisenberg_chain_hamiltonian(ops)
        pop = polyopt(hamiltonian, registry)
        linear, raw_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            real_moment_matrix=true,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
        )
        descriptor = NCTSSoS._pauli_su2_moment_quotient_descriptor(
            linear,
            6;
            atol=1e-11,
            condition_limit=1e10,
        )
        @test descriptor.raw_moment_count == length(linear.moments)
        @test descriptor.coordinate_provenance == :strict_input
        @test descriptor.quotient_moment_count ==
            sum(length(orbit.pivot_keys) for orbit in descriptor.orbits; init=0)
        @test 0 < descriptor.quotient_moment_count < descriptor.raw_moment_count
        @test descriptor.support_orbit_count == length(descriptor.orbits)
        @test length(descriptor.source_to_orbit_row) == descriptor.raw_moment_count
        @test descriptor.max_pivot_residual <= 1e-11
        @test descriptor.max_invariant_residual <= 1e-10
        @test descriptor.max_reconstruction_residual <= 1e-10
        @test isfinite(descriptor.max_condition)
        @test all(orbit -> issorted(orbit.pivot_keys; lt=NCTSSoS.key_lt), descriptor.orbits)
        repeated = NCTSSoS._pauli_su2_moment_quotient_descriptor(
            linear,
            6;
            atol=1e-11,
            condition_limit=1e10,
        )
        @test [orbit.pivot_keys for orbit in repeated.orbits] ==
            [orbit.pivot_keys for orbit in descriptor.orbits]
        @test [orbit.reconstruction for orbit in repeated.orbits] ==
            [orbit.reconstruction for orbit in descriptor.orbits]

        permuted = deepcopy(linear)
        reverse!(permuted.moments)
        empty!(permuted.moment_index)
        for (idx, key) in pairs(permuted.moments)
            permuted.moment_index[copy(key)] = idx
        end
        permuted_descriptor = NCTSSoS._pauli_su2_moment_quotient_descriptor(
            permuted,
            6;
            atol=1e-11,
            condition_limit=1e10,
        )
        @test [orbit.source_keys for orbit in permuted_descriptor.orbits] ==
            [orbit.source_keys for orbit in descriptor.orbits]
        @test [orbit.pivot_keys for orbit in permuted_descriptor.orbits] ==
            [orbit.pivot_keys for orbit in descriptor.orbits]
        @test [orbit.reconstruction for orbit in permuted_descriptor.orbits] ==
            [orbit.reconstruction for orbit in descriptor.orbits]

        incomplete = deepcopy(linear)
        incomplete_orbit = first(filter(
            orbit -> !isempty(orbit.pivot_keys) && length(orbit.source_keys) > 1,
            descriptor.orbits,
        ))
        missing_key = first(incomplete_orbit.source_keys)
        missing_idx = findfirst(
            key -> NCTSSoS.key_isequal(key, missing_key),
            incomplete.moments,
        )
        @test missing_idx !== nothing
        deleteat!(incomplete.moments, something(missing_idx))
        delete!(incomplete.moment_index, missing_key)
        delete!(incomplete.key_to_monomial, missing_key)
        @test_throws ArgumentError NCTSSoS._pauli_su2_moment_quotient_descriptor(
            incomplete,
            6;
            atol=1e-11,
            condition_limit=1e10,
        )
        projected_descriptor = NCTSSoS._pauli_su2_moment_quotient_descriptor(
            incomplete,
            6;
            allow_registered_projection=true,
            atol=1e-11,
            condition_limit=1e10,
        )
        @test projected_descriptor.coordinate_provenance ==
            :generated_translation_relaxation
        @test any(
            orbit.source_coordinate_mask == :generated_projection
            for orbit in projected_descriptor.orbits
        )

        rng = Random.Xoshiro(0x5a17)
        K = typeof(linear.identity)
        raw_values = Dict{K,ComplexF64}()
        quotient_values = Dict{K,ComplexF64}()
        for orbit in descriptor.orbits
            coordinates = randn(rng, length(orbit.pivot_keys))
            for (pivot_idx, pivot_key) in pairs(orbit.pivot_keys)
                quotient_values[copy(pivot_key)] = coordinates[pivot_idx]
            end
            for (source_idx, source_key) in pairs(orbit.source_keys)
                raw_values[copy(source_key)] = sum(
                    orbit.reconstruction[source_idx, pivot_idx] *
                        coordinates[pivot_idx]
                    for pivot_idx in eachindex(coordinates);
                    init=0.0 + 0.0im,
                )
            end
        end
        raw_forms = Any[linear.objective_lin]
        append!(raw_forms, form for block in linear.psd_blocks_lin for form in block.entries)
        append!(raw_forms, constraint.form for constraint in linear.zero_constraints)
        rewrite_plan = NCTSSoS._pauli_su2_moment_rewrite_plan(descriptor)
        rewrite_workspace = NCTSSoS._pauli_su2_moment_rewrite_workspace(
            rewrite_plan,
            Float64,
        )
        for raw_form in raw_forms
            rewritten = NCTSSoS._pauli_su2_rewrite_form(
                raw_form,
                descriptor,
                Float64;
                atol=1e-11,
            )
            @test _form_value(raw_form, raw_values) ≈
                _form_value(rewritten, quotient_values) atol = 1e-10
            fast_rewritten = NCTSSoS._pauli_su2_rewrite_form_with_plan!(
                rewrite_workspace,
                raw_form,
                rewrite_plan,
                Float64;
                atol=1e-11,
            )
            @test first.(fast_rewritten.terms) == first.(rewritten.terms)
            @test last.(fast_rewritten.terms) ≈ last.(rewritten.terms) atol = 1e-12
        end

        quotient_result = NCTSSoS._pauli_su2_quotient_linear_data(
            linear,
            6;
            atol=1e-11,
            condition_limit=1e10,
        )
        quotient = quotient_result.linear
        @test NCTSSoS.assert_moment_linear_data_invariants(quotient) === nothing
        @test length(quotient.moments) == descriptor.quotient_moment_count
        @test length(quotient.moments) < length(linear.moments)
        @test length(quotient.psd_blocks_lin) == length(linear.psd_blocks_lin)
        @test quotient.psd_block_constraint_idx == linear.psd_block_constraint_idx
        for (raw_block, quotient_block) in zip(
            linear.psd_blocks_lin,
            quotient.psd_blocks_lin,
        )
            @test quotient_block.size == raw_block.size
            @test quotient_block.meta.cone == raw_block.meta.cone
            @test quotient_block.meta.row_labels == raw_block.meta.row_labels
            @test quotient_block.meta.origin isa NCTSSoS.PauliSU2QuotientBlockOrigin
            @test quotient_block.meta.origin.source_origin === raw_block.meta.origin
            @test quotient_block.meta.origin.label ==
                (hasproperty(raw_block.meta.origin, :label) ?
                    getproperty(raw_block.meta.origin, :label) : nothing)
            @test quotient_block.meta.origin.logical_row_labels ==
                (hasproperty(raw_block.meta.origin, :logical_row_labels) ?
                    getproperty(raw_block.meta.origin, :logical_row_labels) : Any[])
        end
        @test quotient_result.eliminated_zero_row_count > 0
        @test quotient_result.eliminated_zero_row_count ==
            length(linear.zero_constraints) - length(quotient.zero_constraints)
        @test sum(last, quotient_result.eliminated_zero_feature_histogram; init=0) ==
            quotient_result.eliminated_zero_row_count
        @test all(
            constraint -> isconcretetype(typeof(constraint.origin)),
            quotient.zero_constraints,
        )
        @test all(
            constraint -> constraint.trusted_self_adjoint,
            filter(constraint -> constraint.trusted_self_adjoint, quotient.zero_constraints),
        )

        redundant_indices = findall(
            constraint -> NCTSSoS._pauli_su2_quotient_redundant_zero_origin(
                constraint.origin,
            ),
            linear.zero_constraints,
        )
        @test !isempty(redundant_indices)
        for constraint_idx in redundant_indices
            rewritten = NCTSSoS._pauli_su2_rewrite_form(
                linear.zero_constraints[constraint_idx].form,
                descriptor,
                Float64;
                atol=1e-11,
            )
            @test isempty(rewritten)
        end

        mislabeled = deepcopy(linear)
        mislabeled_idx = first(redundant_indices)
        mislabeled_constraint = mislabeled.zero_constraints[mislabeled_idx]
        mislabeled.zero_constraints[mislabeled_idx] = NCTSSoS.ScalarLinearConstraint(
            linear.objective_lin,
            :zero,
            mislabeled_constraint.origin;
            trusted_self_adjoint=mislabeled_constraint.trusted_self_adjoint,
        )
        @test_throws ArgumentError NCTSSoS._pauli_su2_quotient_linear_data(
            mislabeled,
            6;
            atol=1e-11,
            condition_limit=1e10,
        )

        pivot_orbit = only(filter(orbit -> !isempty(orbit.pivot_keys), descriptor.orbits)[1:1])
        caller_key = copy(first(pivot_orbit.pivot_keys))
        caller_form = NCTSSoS.LinearMomentForm{K,Float64}([caller_key => 1.0])
        owned_rewrite = NCTSSoS._pauli_su2_rewrite_form(
            caller_form,
            descriptor,
            Float64;
            atol=1e-11,
        )
        saved_terms = deepcopy(owned_rewrite.terms)
        isempty(caller_key) || fill!(caller_key, typemax(eltype(caller_key)))
        @test owned_rewrite.terms == saved_terms

        @test !raw_report.su2_moment_quotient
        @test raw_report.su2_moment_raw_count == length(linear.moments)
        @test raw_report.su2_moment_quotient_count == length(linear.moments)
        @test raw_report.su2_moment_quotient_reduction_ratio == 1.0
        @test isempty(raw_report.su2_moment_singlet_channel_support_counts)
        @test isempty(raw_report.su2_moment_eliminated_zero_feature_histogram)

        integrated, integrated_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop,
                ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=true,
                contiguous_rdm_decomposition=:su2,
                contiguous_rdm_support=:extend,
                su2_symmetry=true,
                base_su2_extend_rdm=true,
                su2_moment_quotient=true,
        )
        @test integrated_report.su2_moment_quotient
        @test integrated_report.su2_moment_raw_count <= length(linear.moments)
        @test integrated_report.block_labels == raw_report.block_labels
        @test integrated_report.su2_moment_quotient_count == length(integrated.moments)
        @test integrated_report.su2_moment_quotient_count <
            integrated_report.su2_moment_raw_count
        @test all(
            !block.meta.origin.transform.descriptor.sign_symmetry
            for block in integrated.psd_blocks_lin
        )
        @test all(
            block.meta.origin.transform.descriptor.coordinate_provenance ==
                :generated_translation_relaxation
            for block in integrated.psd_blocks_lin
        )
        @test integrated_report.su2_moment_quotient_reduction_ratio ≈
            length(integrated.moments) / integrated_report.su2_moment_raw_count
        @test integrated_report.su2_moment_support_orbit_count ==
            descriptor.support_orbit_count
        @test integrated_report.su2_moment_eliminated_zero_row_count <=
            quotient_result.eliminated_zero_row_count
        @test [block.size for block in integrated.psd_blocks_lin] ==
            [block.size for block in linear.psd_blocks_lin]
        @test [block.meta.cone for block in integrated.psd_blocks_lin] ==
            [block.meta.cone for block in linear.psd_blocks_lin]
        @test NCTSSoS.assert_moment_linear_data_invariants(integrated) === nothing
        metrics = translation_report_metrics(integrated_report)
        @test metrics.linear_moment_count == length(integrated.moments)
        @test metrics.su2_moment_quotient
        @test metrics.su2_moment_raw_count == integrated_report.su2_moment_raw_count
        @test metrics.estimated_sos_dual_scalar_equalities_upper_bound ==
            length(integrated.moments) - 1
        @test translation_solve_support(integrated_report).supported

        signed, signed_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=true,
            real_moment_matrix=true,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
            su2_moment_quotient=true,
        )
        @test signed_report.sign_symmetry
        @test signed_report.su2_moment_quotient
        @test signed_report.su2_moment_quotient_count == length(signed.moments)
        @test signed_report.su2_moment_quotient_count <
            signed_report.su2_moment_raw_count
        @test all(
            block.meta.origin.transform.descriptor.sign_symmetry
            for block in signed.psd_blocks_lin
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(signed) === nothing

        higher_order, higher_order_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop,
                ops,
                3;
                sign_symmetry=true,
                real_moment_matrix=false,
                contiguous_rdm_k=6,
                contiguous_rdm_decomposition=:su2,
                contiguous_rdm_support=:extend,
                su2_symmetry=true,
                base_su2_extend_rdm=true,
                linear_state_opt_width=5,
                psd_state_opt_width=2,
                su2_moment_quotient=true,
            )
        higher_order_descriptor =
            first(higher_order.psd_blocks_lin).meta.origin.transform.descriptor
        @test higher_order_descriptor.sign_symmetry
        @test all(
            orbit.source_coordinate_mask in
                (:complete, :sign_trivial, :sign_trivial_superset, :generated_projection)
            for orbit in higher_order_descriptor.orbits
        )
        @test any(
            orbit.source_coordinate_mask in
                (:sign_trivial, :sign_trivial_superset, :generated_projection)
            for orbit in higher_order_descriptor.orbits
        )
        @test !any(
            pair -> first(pair).feature == :pauli_su2_base_zero,
            higher_order_report.su2_moment_eliminated_zero_feature_histogram,
        )
        @test !haskey(
            higher_order_report.construction_stage_times_ns,
            :su2_wigner_spin_stream,
        )
        @test all(
            pair -> first(pair).feature != :pauli_su2_base_zero,
            higher_order_report.zero_constraint_feature_histogram,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(higher_order) === nothing

        @test_throws ArgumentError NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            su2_symmetry=true,
            su2_moment_quotient=true,
        )
        @test_throws ArgumentError NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
            base_su2_extend_rdm=true,
            su2_moment_quotient=true,
        )

        registry4, ops4 = create_pauli_variables(1:4)
        hamiltonian4 = heisenberg_chain_hamiltonian(ops4)
        combined_pop4 = polyopt(
            hamiltonian4,
            registry4;
            moment_eq_constraints=[hamiltonian4 * hamiltonian4],
        )
        combined4, combined_report4 =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                combined_pop4,
                ops4,
                2;
                sign_symmetry=false,
                reflection_symmetry=true,
                real_moment_matrix=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:su2,
                contiguous_rdm_support=:extend,
                su2_symmetry=true,
                base_su2_extend_rdm=true,
                singlet_channel_equalities=true,
                linear_state_opt_width=2,
                psd_state_opt_width=1,
                su2_moment_quotient=true,
            )
        @test combined_report4.su2_moment_quotient
        @test combined_report4.su2_moment_quotient_count == length(combined4.moments)
        @test combined_report4.su2_moment_quotient_count <
            combined_report4.su2_moment_raw_count
        @test any(
            label -> label.feature == :moment_matrix &&
                label.decomposition == :translation_su2_reflection,
            combined_report4.block_labels,
        )
        @test any(
            label -> label isa NamedTuple && haskey(label, :feature) &&
                label.feature == :contiguous_rdm,
            combined_report4.block_labels,
        )
        @test any(
            label -> label isa NamedTuple && haskey(label, :feature) &&
                label.feature == :psd_state_opt,
            combined_report4.block_labels,
        )
        source_zero_histogram4 = vcat(
            combined_report4.zero_constraint_feature_histogram,
            combined_report4.su2_moment_eliminated_zero_feature_histogram,
        )
        source_zero_features4 = Set(first(pair).feature for pair in source_zero_histogram4)
        @test :moment_equality in source_zero_features4
        @test :linear_state_opt in source_zero_features4
        @test !(:pauli_su2_translation_orbit_singlet_channel_equality in
            source_zero_features4)
        @test translation_solve_support(combined_report4).supported
        @test NCTSSoS.assert_moment_linear_data_invariants(combined4) === nothing

        noninvariant_objective = one(hamiltonian4) * first(ops4[1])
        noninvariant_pop = polyopt(noninvariant_objective, registry4)
        @test_throws ArgumentError NCTSSoS._pauli_translation_base_linear_relaxation(
            noninvariant_pop,
            ops4,
            1;
            sign_symmetry=false,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
            su2_moment_quotient=true,
        )
    end

    @testset "small solver and certificate equivalence" begin
        registry, ops = create_pauli_variables(1:4)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
        common_options = (
            sign_symmetry=false,
            real_moment_matrix=true,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
        )
        raw, _ = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            1;
            common_options...,
        )
        quotient, quotient_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop,
                ops,
                1;
                common_options...,
                su2_moment_quotient=true,
            )

        function solve_certificate(linear)
            sos = NCTSSoS.sos_dualize(
                linear;
                hermitian_representation=:native,
            )
            set_optimizer(sos.model, Clarabel.Optimizer)
            set_silent(sos.model)
            set_optimizer_attribute(sos.model, "tol_gap_abs", 1e-9)
            set_optimizer_attribute(sos.model, "tol_gap_rel", 1e-9)
            set_optimizer_attribute(sos.model, "tol_feas", 1e-9)
            optimize!(sos.model)
            residual = NCTSSoS.sos_dual_certificate_residual(linear, sos)
            diagnostics = NCTSSoS.sos_dual_block_diagnostics(linear, sos)
            return (
                status=termination_status(sos.model),
                objective=objective_value(sos.model),
                residual=residual.max_abs_residual,
                min_eigenvalue=minimum(
                    diagnostic.native_min_eigenvalue for diagnostic in diagnostics
                ),
            )
        end

        raw_solution = solve_certificate(raw)
        quotient_solution = solve_certificate(quotient)
        @test raw_solution.status in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test quotient_solution.status in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test abs(raw_solution.objective - quotient_solution.objective) <= 1e-8
        @test raw_solution.residual <= 1e-7
        @test quotient_solution.residual <= 1e-7
        @test raw_solution.min_eigenvalue >= -1e-7
        @test quotient_solution.min_eigenvalue >= -1e-7
        @test quotient_report.su2_moment_quotient_count <
            quotient_report.su2_moment_raw_count
    end
end

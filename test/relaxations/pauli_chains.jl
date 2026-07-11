# Tests for sparse Pauli spin-chain bases and charge-compatible sign symmetry.

using Test, NCTSSoS, LinearAlgebra

function _pauli_support(mono)
    return sort!(Int[NCTSSoS._pauli_site(idx) for idx in mono.word])
end

function _linear_forms_isapprox(left, right; atol=1e-12)
    length(left) == length(right) || return false
    for (left_term, right_term) in zip(left, right)
        first(left_term) == first(right_term) || return false
        isapprox(last(left_term), last(right_term); atol, rtol=0) || return false
    end
    return true
end

function _linear_blocks_isapprox(left, right; atol=1e-12)
    size(left) == size(right) || return false
    return all(
        _linear_forms_isapprox(left[idx], right[idx]; atol) for idx in eachindex(left)
    )
end

function _dict_value_by_key_isequal(dict, key)
    for (candidate, value) in dict
        NCTSSoS.key_isequal(candidate, key) && return value
    end
    error("missing key $(repr(key))")
end

function _linear_form_coefficient_by_key_isequal(form, key)
    for (candidate, coef) in form
        NCTSSoS.key_isequal(candidate, key) && return coef
    end
    return zero(last(first(form)))
end

function _linear_form_is_self_adjoint(linear, form)
    isempty(form) && return true
    rtol = sqrt(eps(float(real(one(last(first(form)))))))
    for (key, coef) in form
        adjoint_key = _dict_value_by_key_isequal(linear.adjoint_key, key)
        adjoint_coef = _linear_form_coefficient_by_key_isequal(form, adjoint_key)
        isapprox(adjoint_coef, conj(coef); atol=0, rtol) || return false
    end
    return true
end

function _is_periodic_contiguous(sites::Vector{Int}, n::Int)
    isempty(sites) && return true
    target = Set(sites)
    width = length(sites)
    return any(Set(mod1(start + offset, n) for offset in 0:(width - 1)) == target for start in 1:n)
end

@testset "Pauli contiguous chain basis" begin
    reg4, (σx4, _, _) = create_pauli_variables(1:4)

    @test length(pauli_contiguous_chain_basis(reg4, 0)) == 1
    @test length(pauli_contiguous_chain_basis(reg4, 1)) == 13
    @test length(pauli_contiguous_chain_basis(reg4, 2)) == 49
    @test length(pauli_contiguous_chain_basis(reg4, 3)) == 157
    @test length(pauli_contiguous_chain_basis(reg4, 4)) == 238
    @test length(pauli_contiguous_chain_basis(reg4, 2; periodic=false)) == 40

    basis = pauli_contiguous_chain_basis(reg4, 3)
    @test one(σx4[1]) in basis
    @test all(mono -> degree(mono) <= 3, basis)
    @test all(mono -> _is_periodic_contiguous(_pauli_support(mono), 4), basis)

    reg10, _ = create_pauli_variables(1:10)
    @test length(pauli_contiguous_chain_basis(reg10, 2)) == 1 + 10 * (3 + 9)

    reg100, _ = create_pauli_variables(1:100)
    @test length(pauli_contiguous_chain_basis(reg100, 4)) == 1 + 100 * (3 + 9 + 27 + 81)

    @test_throws ArgumentError pauli_contiguous_chain_basis(reg4, -1)
end

@testset "Pauli chain basis closure under spatial and sign symmetries" begin
    n = 6
    reg, (σx, σy, σz) = create_pauli_variables(1:n)
    basis = pauli_contiguous_chain_basis(reg, 4)
    lookup = Set(basis)

    translation = pauli_site_permutation([2:n; 1])
    reflection = pauli_site_permutation(reverse(1:n))
    sign = pauli_sign_symmetry(n; integer_type=eltype(σx[1].word))

    for g in (translation, reflection, sign), mono in basis
        _, image = NCTSSoS._act_monomial(g, mono)
        @test image in lookup
    end

    @test NCTSSoS._act_monomial(sign, σx[1]) == (-1, σx[1])
    @test NCTSSoS._act_monomial(sign, σy[1]) == (-1, σy[1])
    @test NCTSSoS._act_monomial(sign, σz[1]) == (1, σz[1])
end

@testset "Pauli axis orbit diagnostics quantify axis-orbit reduction" begin
    _, ops4 = create_pauli_variables(1:4)

    @test isdefined(NCTSSoS, :pauli_axis_orbit_diagnostics)
    if isdefined(NCTSSoS, :pauli_axis_orbit_diagnostics)
        order1 = pauli_axis_orbit_diagnostics(ops4, 1)
        @test order1.translation_orbit_count == 4
        @test order1.axis_orbit_count == 2
        @test order1.axis_orbit_size_histogram == [1 => 1, 3 => 1]
        @test order1.max_axis_orbit_size == 3
        @test order1.axis_group_order == 24

        order2 = pauli_axis_orbit_diagnostics(ops4, 2)
        @test order2.translation_orbit_count == 13
        @test order2.axis_orbit_count == 4
        @test order2.axis_orbit_size_histogram == [1 => 1, 3 => 2, 6 => 1]
        @test order2.max_axis_orbit_size == 6
        @test order2.reduction_ratio == 13 / 4

        _, ops12 = create_pauli_variables(1:12)
        order4 = pauli_axis_orbit_diagnostics(ops12, 4)
        @test order4.translation_orbit_count == 121
        @test order4.axis_orbit_count == 23
        @test order4.axis_orbit_size_histogram == [1 => 1, 3 => 4, 6 => 18]
        @test order4.max_axis_orbit_size == 6
        @test order4.reduction_ratio == 121 / 23
    end
end

@testset "Pauli axis rotations act as signed permutations on translation orbit basis" begin
    _, (σx, σy, σz) = create_pauli_variables(1:4)

    actions = pauli_axis_rotation_action_matrices((σx, σy, σz), 1)
    @test actions.n_sites == 4
    @test actions.order == 1
    @test actions.translation_orbit_count == 4
    @test actions.generator_labels == [:global_h, :global_s]
    @test length(actions.generator_action_matrices) == 2

    basis = actions.translation_orbit_representatives
    indices = Dict(mono => idx for (idx, mono) in pairs(basis))
    identity = Matrix{Int}(I, length(basis), length(basis))
    global_h = Matrix(actions.generator_action_matrices[1])
    global_s = Matrix(actions.generator_action_matrices[2])
    group_actions = Matrix.(actions.group_action_matrices)
    group_keys = Set(Tuple(vec(action)) for action in group_actions)

    for action in (global_h, global_s)
        @test all(j -> count(!iszero, action[:, j]) == 1, axes(action, 2))
        @test all(i -> count(!iszero, action[i, :]) == 1, axes(action, 1))
        @test action' * action == identity
    end

    @test global_h^2 == identity
    @test global_s^4 == identity
    @test global_h[indices[σz[1]], indices[σx[1]]] == 1
    @test global_h[indices[σy[1]], indices[σy[1]]] == -1
    @test global_h[indices[σx[1]], indices[σz[1]]] == 1
    @test global_s[indices[σy[1]], indices[σx[1]]] == 1
    @test global_s[indices[σx[1]], indices[σy[1]]] == -1
    @test global_s[indices[σz[1]], indices[σz[1]]] == 1

    @test actions.axis_group_order == 24
    @test length(group_actions) == actions.axis_group_order
    @test group_actions[1] == identity
    @test length(group_keys) == actions.axis_group_order
    @test Tuple(vec(global_h)) in group_keys
    @test Tuple(vec(global_s)) in group_keys
    @test sort(actions.group_element_orders) == [1; fill(2, 9); fill(3, 8); fill(4, 6)]
    for left in group_actions, right in group_actions
        @test Tuple(vec(left * right)) in group_keys
    end

    projectors = pauli_axis_rotation_irrep_projectors((σx, σy, σz), 1)
    @test projectors.irrep_labels == [
        :trivial,
        :orientation_sign,
        :two_dimensional,
        :standard,
        :axis_vector,
    ]
    @test projectors.irrep_dimensions == [1, 1, 2, 3, 3]
    @test projectors.class_size_histogram == [
        :identity => 1,
        :edge_180 => 6,
        :face_180 => 3,
        :body_120 => 8,
        :face_90 => 6,
    ]
    @test projectors.projector_ranks == [1, 0, 0, 0, 3]
    @test projectors.irrep_multiplicities == [1, 0, 0, 0, 1]
    projector_sum = sum(projectors.projector_matrices)
    @test projector_sum ≈ identity atol = 1e-10
    for projector in projectors.projector_matrices
        @test projector * projector ≈ projector atol = 1e-10
        @test projector' ≈ projector atol = 1e-10
    end

    order2_projectors = pauli_axis_rotation_irrep_projectors((σx, σy, σz), 2)
    @test order2_projectors.translation_orbit_count == 13
    @test order2_projectors.projector_ranks == [2, 0, 2, 3, 6]
    @test order2_projectors.irrep_multiplicities == [2, 0, 1, 1, 2]
    @test sum(order2_projectors.projector_ranks) == order2_projectors.translation_orbit_count

    transform = pauli_axis_rotation_isotypic_transform((σx, σy, σz), 2)
    @test transform.isotypic_block_sizes == order2_projectors.projector_ranks
    @test transform.irrep_multiplicities == order2_projectors.irrep_multiplicities
    @test size(transform.transform_matrix) == (13, 13)
    @test transform.transform_matrix' * transform.transform_matrix ≈ Matrix{Float64}(I, 13, 13) atol = 1e-10
    @test transform.transform_matrix * transform.transform_matrix' ≈ Matrix{Float64}(I, 13, 13) atol = 1e-10

    action2 = pauli_axis_rotation_action_matrices((σx, σy, σz), 2)
    offsets = cumsum([0; transform.isotypic_block_sizes])
    for action in action2.group_action_matrices
        reduced_action = transform.transform_matrix' * Matrix{Float64}(action) * transform.transform_matrix
        for row_block in eachindex(transform.isotypic_block_sizes)
            row_range = (offsets[row_block] + 1):offsets[row_block + 1]
            isempty(row_range) && continue
            for col_block in eachindex(transform.isotypic_block_sizes)
                row_block == col_block && continue
                col_range = (offsets[col_block] + 1):offsets[col_block + 1]
                isempty(col_range) && continue
                @test maximum(abs, reduced_action[row_range, col_range]) <= 1e-10
            end
        end
    end
end

@testset "Pauli chain fast-path profile detects supported symmetry specs" begin
    n = 6
    _, ops = create_pauli_variables(1:n)

    @test isdefined(NCTSSoS, :pauli_chain_fast_path_profile)
    if isdefined(NCTSSoS, :pauli_chain_fast_path_profile)
        spec = heisenberg_chain_symmetry_spec(ops; axis_rotations=false, check_invariance=false)
        profile = pauli_chain_fast_path_profile(ops, spec)

        @test profile.n_sites == n
        @test profile.translation_symmetry
        @test profile.reflection_symmetry
        @test !profile.axis_rotation_symmetry
        @test !profile.sign_symmetry
        @test profile.supported_by_translation_fast_path
        @test isempty(profile.unsupported_features)
        @test isempty(profile.missing_required_features)
        @test profile.unrecognized_clifford_generator_count == 0
        @test !profile.check_invariance
        @test profile.offblock_check == :randomized

        sign_spec = try
            heisenberg_chain_symmetry_spec(
                ops;
                sign=true,
                reflection=false,
                axis_rotations=false,
                check_invariance=false,
            )
        catch err
            err
        end
        @test !(sign_spec isa MethodError)
        if sign_spec isa SymmetrySpec
            sign_profile = pauli_chain_fast_path_profile(ops, sign_spec)
            @test sign_profile.translation_symmetry
            @test !sign_profile.reflection_symmetry
            @test sign_profile.sign_symmetry
            @test sign_profile.supported_by_translation_fast_path
            @test isempty(sign_profile.unsupported_features)
            @test isempty(sign_profile.missing_required_features)
        end

        axis_spec = heisenberg_chain_symmetry_spec(ops; axis_rotations=true, check_invariance=false)
        axis_profile = pauli_chain_fast_path_profile(ops, axis_spec)
        @test axis_profile.translation_symmetry
        @test axis_profile.reflection_symmetry
        @test axis_profile.axis_rotation_symmetry
        @test axis_profile.supported_by_translation_fast_path
        @test isempty(axis_profile.unsupported_features)
        @test isempty(axis_profile.missing_required_features)

        axis_only_spec = heisenberg_chain_symmetry_spec(
            ops;
            reflection=false,
            axis_rotations=true,
            check_invariance=false,
        )
        axis_only_profile = pauli_chain_fast_path_profile(ops, axis_only_spec)
        @test axis_only_profile.translation_symmetry
        @test !axis_only_profile.reflection_symmetry
        @test axis_only_profile.axis_rotation_symmetry
        @test axis_only_profile.supported_by_translation_fast_path
        @test isempty(axis_only_profile.unsupported_features)
        @test isempty(axis_only_profile.missing_required_features)

        axis_generator = first(pauli_global_axis_rotation_generators(ops))
        duplicate_axis_spec = SymmetrySpec(
            [pauli_chain_translation(ops), axis_generator, axis_generator];
            check_invariance=false,
        )
        duplicate_axis_profile = pauli_chain_fast_path_profile(ops, duplicate_axis_spec)
        @test !duplicate_axis_profile.axis_rotation_symmetry
        @test duplicate_axis_profile.axis_rotation_generator_count == 2
        @test !duplicate_axis_profile.supported_by_translation_fast_path
        @test duplicate_axis_profile.unsupported_features == [:axis_rotation]

        local_h_spec = SymmetrySpec(
            [CliffordSymmetry(:H, 1; nqubits=n)];
            check_invariance=false,
        )
        local_h_profile = pauli_chain_fast_path_profile(ops, local_h_spec)
        @test !local_h_profile.translation_symmetry
        @test !local_h_profile.supported_by_translation_fast_path
        @test local_h_profile.unsupported_features == [:unrecognized_clifford]
        @test local_h_profile.missing_required_features == [:translation]
        @test local_h_profile.unrecognized_clifford_generator_count == 1
    end
end

@testset "Sparse Pauli charge words follow the supplied chain basis" begin
    n = 6
    reg, _ = create_pauli_variables(1:n)
    basis = pauli_contiguous_chain_basis(reg, 4)

    charge_groups = NCTSSoS._pauli_charge_transform_groups(
        basis,
        PauliChargeSectorSpec(nqubits=n, max_degree=4),
        nothing,
    )
    blocks = collect(Iterators.flatten(charge_groups))

    @test sum(size(block.row_basis, 1) for block in blocks) == length(basis)
    @test Set(block.label.charge for block in blocks) == Set(-4:4)
    @test all(block -> block.provenance == :charge_sector, blocks)

    sign_group = CliffordSymmetryGroup(
        pauli_sign_symmetry(n; integer_type=eltype(basis[1].word));
        nqubits=n,
        integer_type=eltype(basis[1].word),
    )
    signed_groups = NCTSSoS._pauli_charge_transform_groups(
        basis,
        PauliChargeSectorSpec(nqubits=n, max_degree=4),
        sign_group,
    )
    signed_blocks = collect(Iterators.flatten(signed_groups))

    @test sum(size(block.row_basis, 1) for block in signed_blocks) == length(basis)
    @test all(block -> block.label.group_order == 2, signed_blocks)
    @test all(block -> block.provenance == :charge_sector, signed_blocks)
end


# Tests for translation-invariant Pauli chain relaxations.

using Test, NCTSSoS, JuMP, LinearAlgebra
import Clarabel

if !@isdefined(SOLVER)
    # `test/runtests.jl` includes TestUtils first, so CI still uses COSMO.
    # Direct-file runs use a root-environment solver to avoid requiring test extras.
    const SOLVER = optimizer_with_attributes(
        Clarabel.Optimizer,
        "verbose" => false,
    )
end

if !@isdefined(flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

@testset "Translation-invariant Pauli chain relaxation" begin
    quiet(f) = redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            f()
        end
    end

    @testset "contiguous basis and Heisenberg helpers" begin
        n = 8
        registry, ops = create_pauli_variables(1:n)

        basis = pauli_contiguous_chain_basis(ops, 2)
        hamiltonian = heisenberg_chain_hamiltonian(ops)

        @test length(basis) == 1 + n * (3 + 9)
        @test length(terms(hamiltonian)) == 3n
        @test polyopt(hamiltonian, registry).objective == hamiltonian

        orbit_reps = NCTSSoS._pauli_contiguous_chain_orbit_representatives(ops, 3)
        generic_reps = NCTSSoS._translation_orbit_representatives(
            pauli_contiguous_chain_basis(ops, 3),
            n,
        )
        @test orbit_reps == generic_reps
        @test length(orbit_reps) == 1 + sum(3^ℓ for ℓ in 1:3)
    end

    @testset "SolverConfig overload routes recognized chain specs" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
        optimizer = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
        symmetry = heisenberg_chain_symmetry_spec(
            ops;
            axis_rotations=false,
            check_invariance=false,
        )
        config = SolverConfig(; optimizer, order=1, symmetry)

        @test applicable(pauli_translation_invariant_nctssos, pop, config)
        if applicable(pauli_translation_invariant_nctssos, pop, config)
            config_result = quiet() do
                pauli_translation_invariant_nctssos(pop, config; dualize=false)
            end
            direct_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    ops,
                    1,
                    optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                    dualize=false,
                    sign_symmetry=false,
                    reflection_symmetry=true,
                    check_invariance=false,
                )
            end

            @test config_result.objective ≈ direct_result.objective atol = 1e-6
            @test config_result.report.sign_symmetry == false
            @test config_result.report.psd_block_sizes == direct_result.report.psd_block_sizes
            @test any(label -> haskey(label, :reflection), config_result.report.block_labels)
            @test translation_solve_support(config_result) ==
                translation_solve_support(config_result.report)
            @test translation_solve_support(config_result).supported

            native_symmetry = heisenberg_chain_symmetry_spec(
                ops;
                reflection=false,
                axis_rotations=false,
                check_invariance=false,
            )
            native_config = SolverConfig(; optimizer, order=1, symmetry=native_symmetry)
            native_config_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    native_config;
                    dualize=true,
                    real_moment_matrix=false,
                    sos_hermitian_representation=:native,
                )
            end
            @test any(
                pair -> last(pair) == JuMP.MOI.HermitianPositiveSemidefiniteConeTriangle,
                JuMP.list_of_constraint_types(native_config_result.model),
            )

            config_direct_linear_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                )
            end
            @test config_direct_linear_result isa TranslationInvariantResult
            @test config_direct_linear_result.moment_problem isa NCTSSoS.MomentLinearData
            @test config_direct_linear_result.objective ≈ direct_result.objective atol = 1e-6
            @test config_direct_linear_result.report.psd_block_sizes == direct_result.report.psd_block_sizes

            config_qmb_base_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    qmbcertify_base_construct=true,
                )
            end
            auto_qmb_base_result = quiet() do
                cs_nctssos(
                    pop,
                    config;
                    dualize=false,
                    qmbcertify_base_construct=true,
                )
            end
            @test auto_qmb_base_result isa TranslationInvariantResult
            @test auto_qmb_base_result.moment_problem isa NCTSSoS.MomentLinearData
            @test auto_qmb_base_result.objective ≈ config_qmb_base_result.objective atol = 1e-6
            @test auto_qmb_base_result.report.psd_block_sizes ==
                config_qmb_base_result.report.psd_block_sizes
            @test all(
                label -> label.feature == :qmbcertify_base,
                auto_qmb_base_result.report.block_labels,
            )
            auto_qmb_momentum_err = try
                quiet() do
                    cs_nctssos(
                        pop,
                        config;
                        dualize=false,
                        qmbcertify_base_construct=true,
                        momenta=[0],
                    )
                end
                nothing
            catch err
                err
            end
            @test auto_qmb_momentum_err isa ArgumentError
            @test occursin(
                "`qmbcertify_base_construct=true` does not support",
                sprint(showerror, auto_qmb_momentum_err),
            )
            @test occursin("`momenta`", sprint(showerror, auto_qmb_momentum_err))

            config_direct_rdm_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                    contiguous_rdm_k=2,
                )
            end
            @test config_direct_rdm_result isa TranslationInvariantResult
            @test config_direct_rdm_result.moment_problem isa NCTSSoS.MomentLinearData
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :feature) &&
                    label.feature == :contiguous_rdm,
                config_direct_rdm_result.report.block_labels,
            )

            config_axis_equalities_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    axis_rotation_equalities=true,
                )
            end
            @test config_axis_equalities_result isa TranslationInvariantResult
            @test config_axis_equalities_result.objective ≈ direct_result.objective atol = 1e-6
            @test config_axis_equalities_result.report.zero_constraint_count > 0
            @test all(
                zero -> zero.origin.label.feature == :axis_rotation_equality,
                config_axis_equalities_result.moment_problem.linear.zero_constraints,
            )
            @test translation_symmetry_profile(
                config_axis_equalities_result,
            ).axis_rotation_equalities
            @test translation_report_metrics(
                config_axis_equalities_result,
            ).zero_constraint_count == config_axis_equalities_result.report.zero_constraint_count

            config_axis_quotient_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_equalities=true,
                    axis_rotation_quotient=true,
                )
            end
            @test config_axis_quotient_result isa TranslationInvariantResult
            @test config_axis_quotient_result.moment_problem isa NCTSSoS.MomentLinearData
            @test config_axis_quotient_result.objective ≈
                config_axis_equalities_result.objective atol = 1e-6
            @test config_axis_quotient_result.report.zero_constraint_count == 0
            @test config_axis_quotient_result.report.linear_moment_count <
                config_axis_equalities_result.report.linear_moment_count
            config_axis_quotient_metrics =
                translation_report_metrics(config_axis_quotient_result)
            @test config_axis_quotient_metrics.axis_rotation_quotient
            @test !config_axis_quotient_metrics.axis_rotation_equalities
            @test config_axis_quotient_metrics.axis_rotation_quotient_moment_key_count ==
                config_axis_quotient_result.report.linear_moment_count

            translation_only_config = SolverConfig(
                ;
                optimizer,
                order=1,
                symmetry=heisenberg_chain_symmetry_spec(
                    ops;
                    reflection=false,
                    axis_rotations=false,
                    check_invariance=false,
                ),
            )

            auto_direct_su2_extend_result = quiet() do
                cs_nctssos(
                    pop,
                    translation_only_config;
                    dualize=false,
                    direct_linear=true,
                    su2_symmetry=true,
                    base_su2_extend_rdm=true,
                    contiguous_rdm_k=3,
                    contiguous_rdm_decomposition=:su2,
                    contiguous_rdm_support=:extend,
                )
            end
            @test auto_direct_su2_extend_result isa TranslationInvariantResult
            @test auto_direct_su2_extend_result.moment_problem isa NCTSSoS.MomentLinearData
            @test termination_status(auto_direct_su2_extend_result.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            auto_direct_su2_extend_metrics =
                translation_report_metrics(auto_direct_su2_extend_result)
            @test auto_direct_su2_extend_metrics.su2_moment_symmetry
            @test auto_direct_su2_extend_metrics.su2_rdm_symmetry
            @test auto_direct_su2_extend_metrics.contiguous_rdm_block_count ==
                length(pauli_su2_rdm_blocks(3))
            @test auto_direct_su2_extend_metrics.su2_base_zero_row_count > 0
            @test translation_solve_support(auto_direct_su2_extend_result).supported

            auto_su2_extend_symbolic_err = try
                cs_nctssos(
                    pop,
                    translation_only_config;
                    dualize=false,
                    su2_symmetry=true,
                    base_su2_extend_rdm=true,
                    contiguous_rdm_k=3,
                    contiguous_rdm_decomposition=:su2,
                    contiguous_rdm_support=:extend,
                )
                nothing
            catch err
                err
            end
            @test auto_su2_extend_symbolic_err isa ArgumentError
            @test occursin(
                "requires `direct_linear=true`",
                sprint(showerror, auto_su2_extend_symbolic_err),
            )

            config_axis_quotient_missing_axis_err = try
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_quotient=true,
                )
                nothing
            catch err
                err
            end
            @test config_axis_quotient_missing_axis_err isa ArgumentError
            @test occursin(
                "requires `axis_rotation_equalities=true`",
                sprint(showerror, config_axis_quotient_missing_axis_err),
            )

            config_axis_quotient_su2_err = try
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_equalities=true,
                    axis_rotation_quotient=true,
                    su2_symmetry=true,
                )
                nothing
            catch err
                err
            end
            @test config_axis_quotient_su2_err isa ArgumentError
            @test occursin(
                "already subsumes finite axis rotations",
                sprint(showerror, config_axis_quotient_su2_err),
            )

            config_axis_quotient_complex_err = try
                pauli_translation_invariant_nctssos(
                    pop,
                    translation_only_config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_equalities=true,
                    axis_rotation_quotient=true,
                    real_moment_matrix=false,
                )
                nothing
            catch err
                err
            end
            @test config_axis_quotient_complex_err isa ArgumentError
            @test occursin(
                "requires `real_moment_matrix=true`",
                sprint(showerror, config_axis_quotient_complex_err),
            )

            config_axis_quotient_symbolic_err = try
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    axis_rotation_equalities=true,
                    axis_rotation_quotient=true,
                )
                nothing
            catch err
                err
            end
            @test config_axis_quotient_symbolic_err isa ArgumentError
            @test occursin(
                "direct-linear Pauli translation constructor",
                sprint(showerror, config_axis_quotient_symbolic_err),
            )

            direct_linear_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    ops,
                    1,
                    optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                    dualize=false,
                    sign_symmetry=false,
                    reflection_symmetry=true,
                    check_invariance=false,
                    direct_linear=true,
                )
            end
            @test direct_linear_result.moment_problem isa NCTSSoS.MomentLinearData
            @test direct_linear_result.objective ≈ direct_result.objective atol = 1e-6
            @test direct_linear_result.report.psd_block_sizes == direct_result.report.psd_block_sizes

            auto_result = quiet() do
                cs_nctssos(pop, config; dualize=false)
            end
            @test auto_result isa TranslationInvariantResult
            @test auto_result.objective ≈ direct_result.objective atol = 1e-6
            @test auto_result.report.psd_block_sizes == direct_result.report.psd_block_sizes

            auto_direct_linear_result = quiet() do
                cs_nctssos(pop, config; dualize=false, direct_linear=true)
            end
            @test auto_direct_linear_result isa TranslationInvariantResult
            @test auto_direct_linear_result.moment_problem isa NCTSSoS.MomentLinearData
            @test auto_direct_linear_result.objective ≈ direct_result.objective atol = 1e-6

            auto_direct_rdm_result = quiet() do
                cs_nctssos(pop, config; dualize=false, direct_linear=true, contiguous_rdm_k=2)
            end
            @test auto_direct_rdm_result isa TranslationInvariantResult
            @test auto_direct_rdm_result.moment_problem isa NCTSSoS.MomentLinearData
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :feature) &&
                    label.feature == :contiguous_rdm,
                auto_direct_rdm_result.report.block_labels,
            )

            auto_axis_quotient_symbolic_err = try
                cs_nctssos(
                    pop,
                    config;
                    dualize=false,
                    axis_rotation_equalities=true,
                    axis_rotation_quotient=true,
                )
                nothing
            catch err
                err
            end
            @test auto_axis_quotient_symbolic_err isa ArgumentError
            @test occursin(
                "direct-linear Pauli translation constructor",
                sprint(showerror, auto_axis_quotient_symbolic_err),
            )

            auto_qmbcertify_rdm_err = try
                cs_nctssos(
                    pop,
                    config;
                    dualize=false,
                    contiguous_rdm_k=2,
                    contiguous_rdm_decomposition=:qmbcertify,
                )
                nothing
            catch err
                err
            end
            @test auto_qmbcertify_rdm_err isa ArgumentError
            @test occursin(
                "structural-target-only",
                sprint(showerror, auto_qmbcertify_rdm_err),
            )
            @test occursin(
                "shared RDM PSD variables",
                sprint(showerror, auto_qmbcertify_rdm_err),
            )

            auto_direct_qmbcertify_rdm_err = try
                cs_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                    contiguous_rdm_k=2,
                    contiguous_rdm_decomposition=:qmbcertify,
                )
                nothing
            catch err
                err
            end
            @test auto_direct_qmbcertify_rdm_err isa ArgumentError
            @test occursin(
                "structural-target-only",
                sprint(showerror, auto_direct_qmbcertify_rdm_err),
            )
            @test occursin(
                "shared RDM PSD variables",
                sprint(showerror, auto_direct_qmbcertify_rdm_err),
            )

            auto_axis_equalities_result = quiet() do
                cs_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_equalities=true,
                )
            end
            @test auto_axis_equalities_result isa TranslationInvariantResult
            @test auto_axis_equalities_result.moment_problem isa NCTSSoS.MomentLinearData
            @test auto_axis_equalities_result.objective ≈ direct_result.objective atol = 1e-6
            @test auto_axis_equalities_result.report.zero_constraint_count > 0
            @test all(
                zero -> zero.origin.label.feature == :axis_rotation_equality,
                auto_axis_equalities_result.moment_problem.zero_constraints,
            )
            @test translation_zero_origin_histogram(auto_axis_equalities_result.moment_problem) == [
                (feature=:axis_rotation_equality, decomposition=nothing, reason=nothing) =>
                    auto_axis_equalities_result.report.zero_constraint_count,
            ]
            @test translation_zero_origin_histogram(auto_axis_equalities_result) ==
                translation_zero_origin_histogram(auto_axis_equalities_result.moment_problem)

            auto_axis_quotient_result = quiet() do
                cs_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_equalities=true,
                    axis_rotation_quotient=true,
                )
            end
            @test auto_axis_quotient_result isa TranslationInvariantResult
            @test auto_axis_quotient_result.moment_problem isa NCTSSoS.MomentLinearData
            @test auto_axis_quotient_result.objective ≈
                auto_axis_equalities_result.objective atol = 1e-6
            @test auto_axis_quotient_result.report.zero_constraint_count == 0
            @test auto_axis_quotient_result.report.linear_moment_count <
                auto_axis_equalities_result.report.linear_moment_count
            @test translation_report_metrics(auto_axis_quotient_result).axis_rotation_quotient

            qmbcertify_rdm_config = SolverConfig(
                ;
                optimizer,
                order=2,
                symmetry=heisenberg_chain_symmetry_spec(
                    ops;
                    axis_rotations=true,
                    check_invariance=false,
                ),
            )
            qmbcertify_rdm_blocks = pauli_qmbcertify_rdm_blocks(4; ambient_sites=n)
            qmbcertify_batch_blocks = NCTSSoS._translation_qmbcertify_rdm_linear_blocks(
                n,
                4,
                qmbcertify_rdm_blocks,
                Vector{UInt8},
                ComplexF64,
            )
            qmbcertify_single_blocks = [
                NCTSSoS._translation_qmbcertify_rdm_linear_block(
                    n,
                    4,
                    rows,
                    Vector{UInt8},
                    ComplexF64,
                ) for rows in qmbcertify_rdm_blocks
            ]
            @test length(qmbcertify_batch_blocks) == length(qmbcertify_single_blocks)
            for block_idx in eachindex(qmbcertify_batch_blocks)
                @test size(qmbcertify_batch_blocks[block_idx]) ==
                    size(qmbcertify_single_blocks[block_idx])
                for idx in eachindex(qmbcertify_batch_blocks[block_idx])
                    @test getfield(qmbcertify_batch_blocks[block_idx][idx], :terms) ==
                        getfield(qmbcertify_single_blocks[block_idx][idx], :terms)
                end
            end
            qmbcertify_rdm_templates =
                NCTSSoS._qmbcertify_rdm_block_templates(4; ambient_sites=n)
            @test qmbcertify_rdm_templates == qmbcertify_rdm_blocks
            qmbcertify_mutated_blocks = pauli_qmbcertify_rdm_blocks(4; ambient_sites=n)
            push!(first(qmbcertify_mutated_blocks), -1)
            @test pauli_qmbcertify_rdm_blocks(4; ambient_sites=n) ==
                qmbcertify_rdm_templates
            qmbcertify_addons = NCTSSoS._translation_qmbcertify_rdm_addons(
                n,
                [4],
                Vector{UInt8},
                ComplexF64,
            )
            @test qmbcertify_addons.row_blocks_by_k[4] == qmbcertify_rdm_templates
            @test length(qmbcertify_addons.linear_blocks_by_k[4]) ==
                length(qmbcertify_batch_blocks)
            for block_idx in eachindex(qmbcertify_batch_blocks)
                @test size(qmbcertify_addons.linear_blocks_by_k[4][block_idx]) ==
                    size(qmbcertify_batch_blocks[block_idx])
                for idx in eachindex(qmbcertify_batch_blocks[block_idx])
                    @test getfield(
                        qmbcertify_addons.linear_blocks_by_k[4][block_idx][idx],
                        :terms,
                    ) == getfield(qmbcertify_batch_blocks[block_idx][idx], :terms)
                end
            end
            qmbcertify_addon_basis_keys = Set(
                NCTSSoS.symmetric_canon(NCTSSoS.expval(mono)) for
                mono in qmbcertify_addons.moment_basis
            )
            qmbcertify_addon_entry_keys = Set{Vector{UInt8}}()
            for block in qmbcertify_addons.linear_blocks_by_k[4]
                for form in block
                    for (key, _) in form
                        push!(qmbcertify_addon_entry_keys, key)
                        @test key in qmbcertify_addon_basis_keys
                    end
                end
            end
            @test length(qmbcertify_addon_basis_keys) ==
                length(qmbcertify_addon_entry_keys)
            qmbcertify_axis_quotient_basis = pauli_contiguous_chain_basis(
                ops,
                2;
                periodic=true,
            )
            qmbcertify_axis_quotient = NCTSSoS._translation_axis_rotation_quotient_map(
                Vector{UInt8},
                eltype(first(ops)),
                ops,
                n,
                qmbcertify_axis_quotient_basis,
            )
            qmbcertify_axis_representatives =
                NCTSSoS._axis_quotient_representative_moment_map(
                    Vector{UInt8},
                    eltype(first(ops)),
                    qmbcertify_axis_quotient,
                )
            @test length(qmbcertify_axis_representatives) ==
                length(qmbcertify_axis_quotient.representative_pairs)
            for (key, mono) in qmbcertify_axis_quotient.representative_pairs
                @test haskey(qmbcertify_axis_representatives, key)
                @test qmbcertify_axis_representatives[key] == mono
            end
            @test all(
                info -> info.forced_zero || haskey(qmbcertify_axis_representatives, info.key),
                values(qmbcertify_axis_quotient.key_map),
            )
            trusted_axis_builder = NCTSSoS.MomentLinearBuilder(
                Vector{UInt8},
                Float64,
                eltype(qmbcertify_axis_quotient_basis),
            )
            trusted_axis_form = NCTSSoS._linear_moment_form_from_owned_pairs!(
                [copy(trusted_axis_builder.identity) => 1.0],
            )
            NCTSSoS._add_zero_constraint_trusted!(
                trusted_axis_builder,
                trusted_axis_form,
                NCTSSoS.TranslationZeroOrigin((feature=:trusted_axis_rewrite,), 1, 1, :scalar);
                trusted_self_adjoint=true,
            )
            NCTSSoS._apply_translation_axis_rotation_quotient!(
                trusted_axis_builder,
                qmbcertify_axis_quotient,
            )
            @test length(trusted_axis_builder.zero_constraints) == 1
            @test only(trusted_axis_builder.zero_constraints).trusted_self_adjoint

            auto_qmbcertify_axis_rdm_result = quiet() do
                cs_nctssos(
                    pop,
                    qmbcertify_rdm_config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_quotient=true,
                    contiguous_rdm_k=4,
                    contiguous_rdm_decomposition=:qmbcertify,
                )
            end
            qmbcertify_rdm_sizes = length.(qmbcertify_rdm_blocks)
            @test auto_qmbcertify_axis_rdm_result isa TranslationInvariantResult
            @test auto_qmbcertify_axis_rdm_result.moment_problem isa
                NCTSSoS.MomentLinearData
            @test auto_qmbcertify_axis_rdm_result.report.psd_block_sizes[
                end - length(qmbcertify_rdm_sizes) + 1:end
            ] == qmbcertify_rdm_sizes
            @test auto_qmbcertify_axis_rdm_result.report.logical_block_sizes[
                end - length(qmbcertify_rdm_sizes) + 1:end
            ] == qmbcertify_rdm_sizes
            @test all(
                label -> label isa NamedTuple &&
                    label.feature == :contiguous_rdm &&
                    label.decomposition == :qmbcertify &&
                    label.k == 4,
                auto_qmbcertify_axis_rdm_result.report.block_labels[
                    end - length(qmbcertify_rdm_sizes) + 1:end
                ],
            )
            @test translation_solve_support(auto_qmbcertify_axis_rdm_result).supported
            auto_qmbcertify_lso_linear, auto_qmbcertify_lso_report =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    pop,
                    ops,
                    2;
                    reflection_symmetry=true,
                    axis_rotation_symmetry=true,
                    axis_rotation_quotient=true,
                    check_invariance=false,
                    contiguous_rdm_k=4,
                    contiguous_rdm_decomposition=:qmbcertify,
                    contiguous_rdm_support=:extend,
                    linear_state_opt_width=1,
                )
            auto_qmbcertify_lso_metrics =
                translation_report_metrics(auto_qmbcertify_lso_report)
            @test auto_qmbcertify_lso_linear isa NCTSSoS.MomentLinearData
            @test auto_qmbcertify_lso_metrics.contiguous_rdm_block_count ==
                length(qmbcertify_rdm_sizes)
            @test auto_qmbcertify_lso_metrics.linear_state_opt_row_count == 0
            @test auto_qmbcertify_lso_report.zero_constraint_count == 0
            @test auto_qmbcertify_lso_report.axis_rotation_quotient
            @test auto_qmbcertify_lso_report.axis_rotation_forced_zero_moment_class_count > 0
            @test translation_solve_support(auto_qmbcertify_lso_report).supported
            auto_qmbcertify_psd_state_result = quiet() do
                cs_nctssos(
                    pop,
                    qmbcertify_rdm_config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_quotient=true,
                    contiguous_rdm_k=4,
                    contiguous_rdm_decomposition=:qmbcertify,
                    psd_state_opt_width=1,
                )
            end
            auto_qmbcertify_psd_state_metrics =
                translation_report_metrics(auto_qmbcertify_psd_state_result)
            @test auto_qmbcertify_psd_state_result isa TranslationInvariantResult
            @test auto_qmbcertify_psd_state_result.moment_problem isa
                NCTSSoS.MomentLinearData
            @test auto_qmbcertify_psd_state_result.report.psd_block_sizes[
                end - length(qmbcertify_rdm_sizes) + 1:end
            ] == qmbcertify_rdm_sizes
            @test auto_qmbcertify_psd_state_metrics.contiguous_rdm_block_count ==
                length(qmbcertify_rdm_sizes)
            @test auto_qmbcertify_psd_state_metrics.psd_state_opt
            @test auto_qmbcertify_psd_state_metrics.psd_state_opt_block_count > 0
            @test translation_solve_support(auto_qmbcertify_psd_state_result).supported
            auto_qmbcertify_psd_state_psd_block_result = quiet() do
                cs_nctssos(
                    pop,
                    qmbcertify_rdm_config;
                    dualize=false,
                    direct_linear=true,
                    axis_rotation_quotient=true,
                    contiguous_rdm_k=4,
                    contiguous_rdm_decomposition=:qmbcertify,
                    contiguous_rdm_support=:extend,
                    psd_state_opt_width=1,
                    formulation=:psd_blocks,
                    representation=:real,
                    orphan_policy=:free_variables,
                )
            end
            @test auto_qmbcertify_psd_state_psd_block_result isa
                TranslationInvariantResult
            @test auto_qmbcertify_psd_state_psd_block_result.moment_problem isa
                NCTSSoS.MomentLinearData
            @test termination_status(auto_qmbcertify_psd_state_psd_block_result.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test auto_qmbcertify_psd_state_psd_block_result.objective ≈
                auto_qmbcertify_psd_state_result.objective atol = 1e-6
            @test translation_report_metrics(
                auto_qmbcertify_psd_state_psd_block_result,
            ).axis_rotation_quotient
            qmbcertify_source_like_registry, qmbcertify_source_like_ops =
                create_pauli_variables(1:8)
            qmbcertify_source_like_pop = polyopt(
                heisenberg_chain_hamiltonian(qmbcertify_source_like_ops),
                qmbcertify_source_like_registry,
            )
            qmbcertify_source_like_linear, qmbcertify_source_like_report =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    qmbcertify_source_like_pop,
                    qmbcertify_source_like_ops,
                    4;
                    reflection_symmetry=true,
                    axis_rotation_symmetry=true,
                    axis_rotation_quotient=true,
                    contiguous_rdm_k=8,
                    contiguous_rdm_decomposition=:qmbcertify,
                    contiguous_rdm_support=:extend,
                    linear_state_opt_width=7,
                    psd_state_opt_width=3,
                )
            qmbcertify_source_like_metrics =
                translation_report_metrics(qmbcertify_source_like_report)
            @test qmbcertify_source_like_metrics.contiguous_rdm_psd_block_sizes ==
                pauli_qmbcertify_rdm_block_sizes(8; ambient_sites=8)
            @test qmbcertify_source_like_metrics.psd_state_opt
            @test qmbcertify_source_like_metrics.psd_state_opt_block_count ==
                length(qmbcertify_source_like_report.momentum_sectors)
            @test qmbcertify_source_like_metrics.psd_state_opt_psd_block_sizes ==
                fill(18, length(qmbcertify_source_like_report.momentum_sectors))
            @test qmbcertify_source_like_metrics.linear_state_opt
            @test qmbcertify_source_like_metrics.linear_state_opt_row_count == 810
            @test qmbcertify_source_like_metrics.zero_constraint_count == 810
            @test qmbcertify_source_like_report.linear_moment_count == 349
            @test qmbcertify_source_like_report.axis_rotation_quotient
            @test qmbcertify_source_like_report.axis_rotation_quotient_moment_key_count ==
                qmbcertify_source_like_report.linear_moment_count
            @test translation_solve_support(qmbcertify_source_like_report).supported
            @test NCTSSoS.assert_moment_linear_data_invariants(
                qmbcertify_source_like_linear,
            ) === nothing

            auto_momentum_result = quiet() do
                cs_nctssos(pop, config; dualize=false, direct_linear=true, momenta=[0, 1, 2])
            end
            @test auto_momentum_result isa TranslationInvariantResult
            @test auto_momentum_result.moment_problem isa NCTSSoS.MomentLinearData
            @test auto_momentum_result.report.momentum_sectors == [0, 1, 2]

            config_complex_reflection_err = try
                pauli_translation_invariant_nctssos(
                    pop,
                    config;
                    dualize=false,
                    real_moment_matrix=false,
                )
                nothing
            catch err
                err
            end
            @test config_complex_reflection_err isa ArgumentError
            @test occursin(
                "reflection-fixed momentum sectors",
                sprint(showerror, config_complex_reflection_err),
            )

            auto_complex_reflection_err = try
                cs_nctssos(
                    pop,
                    config;
                    dualize=false,
                    direct_linear=true,
                    real_moment_matrix=false,
                )
                nothing
            catch err
                err
            end
            @test auto_complex_reflection_err isa ArgumentError
            @test occursin(
                "reflection-fixed momentum sectors",
                sprint(showerror, auto_complex_reflection_err),
            )

            auto_complex_u1_rdm_result = quiet() do
                cs_nctssos(
                    pop,
                    translation_only_config;
                    dualize=false,
                    direct_linear=true,
                    real_moment_matrix=false,
                    u1_symmetry=true,
                    contiguous_rdm_k=2,
                    contiguous_rdm_decomposition=:u1,
                )
            end
            @test auto_complex_u1_rdm_result isa TranslationInvariantResult
            @test auto_complex_u1_rdm_result.moment_problem isa NCTSSoS.MomentLinearData
            @test termination_status(auto_complex_u1_rdm_result.model) == JuMP.MOI.OPTIMAL
            @test auto_complex_u1_rdm_result.report.real_moment_matrix == false
            @test auto_complex_u1_rdm_result.report.block_labels[end - 2:end] == [
                (feature=:contiguous_rdm, k=2, decomposition=:u1, magnetization=0),
                (feature=:contiguous_rdm, k=2, decomposition=:u1, magnetization=1),
                (feature=:contiguous_rdm, k=2, decomposition=:u1, magnetization=2),
            ]
            @test !isempty(auto_complex_u1_rdm_result.moment_problem.zero_constraints)
            @test all(
                zero -> haskey(zero.origin.label, :component) &&
                    zero.origin.label.component == :complex &&
                    zero.origin.part == :scalar,
                auto_complex_u1_rdm_result.moment_problem.zero_constraints,
            )

            auto_complex_su2_rdm_result = quiet() do
                cs_nctssos(
                    pop,
                    translation_only_config;
                    dualize=false,
                    direct_linear=true,
                    real_moment_matrix=false,
                    su2_symmetry=true,
                    contiguous_rdm_k=2,
                    contiguous_rdm_decomposition=:su2,
                )
            end
            @test auto_complex_su2_rdm_result isa TranslationInvariantResult
            @test auto_complex_su2_rdm_result.moment_problem isa NCTSSoS.MomentLinearData
            @test termination_status(auto_complex_su2_rdm_result.model) == JuMP.MOI.OPTIMAL
            @test auto_complex_su2_rdm_result.report.real_moment_matrix == false
            @test auto_complex_su2_rdm_result.report.block_labels[end - 1:end] == [
                (feature=:contiguous_rdm, k=2, decomposition=:su2, spin2=0),
                (feature=:contiguous_rdm, k=2, decomposition=:su2, spin2=2),
            ]
            @test !isempty(auto_complex_su2_rdm_result.moment_problem.zero_constraints)
            @test all(
                zero -> haskey(zero.origin.label, :component) &&
                    zero.origin.label.component == :complex &&
                    zero.origin.part == :scalar,
                auto_complex_su2_rdm_result.moment_problem.zero_constraints,
            )

            auto_complex_extend_err = try
                quiet() do
                    cs_nctssos(
                        pop,
                        translation_only_config;
                        dualize=false,
                        direct_linear=true,
                        real_moment_matrix=false,
                        contiguous_rdm_k=3,
                        contiguous_rdm_support=:extend,
                        formulation=:psd_blocks,
                        representation=:complex,
                    )
                end
            catch err
                err
            end
            @test auto_complex_extend_err isa ArgumentError
            @test occursin(
                "canonical moment",
                sprint(showerror, auto_complex_extend_err),
            )
            @test occursin(
                "orphan_policy=:free_variables",
                sprint(showerror, auto_complex_extend_err),
            )

            auto_complex_extend_result = quiet() do
                cs_nctssos(
                    pop,
                    translation_only_config;
                    dualize=false,
                    direct_linear=true,
                    real_moment_matrix=false,
                    contiguous_rdm_k=3,
                    contiguous_rdm_support=:extend,
                    formulation=:psd_blocks,
                    representation=:complex,
                    orphan_policy=:free_variables,
                )
            end
            @test auto_complex_extend_result isa TranslationInvariantResult
            @test auto_complex_extend_result.moment_problem isa NCTSSoS.MomentLinearData
            @test termination_status(auto_complex_extend_result.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test !isempty(auto_complex_extend_result.moment_problem.free_keys)
            @test auto_complex_extend_result.report.real_moment_matrix == false
            @test auto_complex_extend_result.report.block_labels[end] ==
                (feature=:contiguous_rdm, k=3)
            @test auto_complex_extend_result.report.psd_block_sizes[end] == 8
            @test translation_solve_support(auto_complex_extend_result).supported
            @test isfinite(auto_complex_extend_result.objective)

            order2_config = SolverConfig(; optimizer, order=2, symmetry)

            config_linear_state_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    order2_config;
                    dualize=false,
                    direct_linear=true,
                    linear_state_opt_width=2,
                )
            end
            @test config_linear_state_result isa TranslationInvariantResult
            @test config_linear_state_result.moment_problem isa NCTSSoS.MomentLinearData
            @test config_linear_state_result.report.zero_constraint_count > 0
            @test all(
                zero -> zero.origin.label.feature == :linear_state_opt,
                config_linear_state_result.moment_problem.zero_constraints,
            )

            config_psd_state_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    order2_config;
                    dualize=false,
                    direct_linear=true,
                    psd_state_opt_width=1,
                )
            end
            @test config_psd_state_result isa TranslationInvariantResult
            @test config_psd_state_result.moment_problem isa NCTSSoS.MomentLinearData
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :feature) &&
                    label.feature == :psd_state_opt,
                config_psd_state_result.report.block_labels,
            )

            auto_linear_state_result = quiet() do
                cs_nctssos(
                    pop,
                    order2_config;
                    dualize=false,
                    direct_linear=true,
                    linear_state_opt_width=2,
                )
            end
            @test auto_linear_state_result isa TranslationInvariantResult
            @test auto_linear_state_result.moment_problem isa NCTSSoS.MomentLinearData
            @test auto_linear_state_result.report.zero_constraint_count > 0
            @test all(
                zero -> zero.origin.label.feature == :linear_state_opt,
                auto_linear_state_result.moment_problem.zero_constraints,
            )

            auto_psd_state_result = quiet() do
                cs_nctssos(
                    pop,
                    order2_config;
                    dualize=false,
                    direct_linear=true,
                    psd_state_opt_width=1,
                )
            end
            @test auto_psd_state_result isa TranslationInvariantResult
            @test auto_psd_state_result.moment_problem isa NCTSSoS.MomentLinearData
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :feature) &&
                    label.feature == :psd_state_opt,
                auto_psd_state_result.report.block_labels,
            )

            axis_config = SolverConfig(
                ;
                optimizer,
                order=1,
                symmetry=heisenberg_chain_symmetry_spec(
                    ops;
                    axis_rotations=true,
                    check_invariance=false,
                ),
            )
            axis_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    axis_config;
                    dualize=false,
                )
            end
            @test axis_result isa TranslationInvariantResult
            @test termination_status(axis_result.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            axis_result_profile = translation_symmetry_profile(axis_result)
            @test axis_result_profile.reflection_symmetry
            @test axis_result_profile.axis_rotation_symmetry
            @test translation_solve_support(axis_result).supported
            auto_axis_result = quiet() do
                cs_nctssos(
                    pop,
                    axis_config;
                    dualize=false,
                    direct_linear=true,
                )
            end
            @test auto_axis_result isa TranslationInvariantResult
            @test auto_axis_result.moment_problem isa NCTSSoS.MomentLinearData
            auto_axis_profile = translation_symmetry_profile(auto_axis_result)
            @test auto_axis_profile.reflection_symmetry
            @test auto_axis_profile.axis_rotation_symmetry
            @test auto_axis_result.objective ≈ axis_result.objective atol = 1e-6
            axis_su2_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    axis_config;
                    dualize=false,
                    su2_symmetry=true,
                )
            end
            @test axis_su2_result isa TranslationInvariantResult
            @test termination_status(axis_su2_result.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test translation_symmetry_profile(axis_su2_result).su2_moment_symmetry
            @test translation_solve_support(axis_su2_result).supported
            auto_axis_su2_result = quiet() do
                cs_nctssos(
                    pop,
                    axis_config;
                    dualize=false,
                    direct_linear=true,
                    su2_symmetry=true,
                )
            end
            @test auto_axis_su2_result isa TranslationInvariantResult
            @test auto_axis_su2_result.moment_problem isa NCTSSoS.MomentLinearData
            @test auto_axis_su2_result.objective ≈ axis_su2_result.objective atol = 1e-6
            @test translation_solve_support(auto_axis_su2_result).supported
            axis_only_config = SolverConfig(
                ;
                optimizer,
                order=2,
                symmetry=heisenberg_chain_symmetry_spec(
                    ops;
                    reflection=false,
                    axis_rotations=true,
                    check_invariance=false,
                ),
            )
            axis_only_result = quiet() do
                pauli_translation_invariant_nctssos(
                    pop,
                    axis_only_config;
                    dualize=false,
                )
            end
            @test axis_only_result isa TranslationInvariantResult
            @test termination_status(axis_only_result.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test translation_symmetry_profile(axis_only_result).axis_rotation_symmetry
            @test translation_solve_support(axis_only_result).supported
            auto_axis_only_result = quiet() do
                cs_nctssos(
                    pop,
                    axis_only_config;
                    dualize=false,
                    direct_linear=true,
                )
            end
            @test auto_axis_only_result isa TranslationInvariantResult
            @test auto_axis_only_result.moment_problem isa NCTSSoS.MomentLinearData
            @test translation_symmetry_profile(auto_axis_only_result).axis_rotation_symmetry
            @test auto_axis_only_result.objective ≈ axis_only_result.objective atol = 1e-6
        end
    end

    @testset "closure check covers contiguous RDM moments" begin
        n = 8
        _, ops = create_pauli_variables(1:n)
        basis = pauli_contiguous_chain_basis(ops, 1)
        basis_order2 = pauli_contiguous_chain_basis(ops, 2)
        basis_order3 = pauli_contiguous_chain_basis(ops, 3)
        hamiltonian = heisenberg_chain_hamiltonian(ops)
        M = eltype(first(ops))

        @test closure_check(:contiguous_rdm, basis; n_sites=n, k=2) === nothing
        @test_throws ArgumentError closure_check(:contiguous_rdm, basis; n_sites=n, k=3)
        @test closure_check(
            :linear_state_opt,
            basis;
            n_sites=n,
            hamiltonian,
            test_width=1,
            sign_symmetry=false,
        ) === nothing
        @test_throws ArgumentError closure_check(
            :linear_state_opt,
            basis;
            n_sites=n,
            hamiltonian,
            test_width=2,
            sign_symmetry=false,
        )
        @test closure_check(
            :psd_state_opt,
            basis_order2;
            n_sites=n,
            hamiltonian,
            test_width=1,
            sign_symmetry=false,
        ) === nothing
        @test_throws ArgumentError closure_check(
            :psd_state_opt,
            basis;
            n_sites=n,
            hamiltonian,
            test_width=1,
            sign_symmetry=false,
        ) === nothing
        @test length(NCTSSoS._psd_state_opt_required_moment_set(
            M,
            n,
            hamiltonian,
            2;
            sign_symmetry=false,
        )) == 1252
        @test_throws ArgumentError closure_check(
            :psd_state_opt,
            basis_order2;
            n_sites=n,
            hamiltonian,
            test_width=3,
            sign_symmetry=false,
        )
        @test closure_check(
            :psd_state_opt,
            basis_order3;
            n_sites=n,
            hamiltonian,
            test_width=2,
            sign_symmetry=false,
        ) === nothing
    end

    @testset "normalized momentum blocks and N>2d block structure" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        mp, report = pauli_translation_invariant_moment_relaxation(pop, ops, 1; sign_symmetry=false)
        identity_cross = mp.constraints[1][2][1, 2]
        @test only(coefficients(identity_cross)) ≈ sqrt(n) + 0im atol = 1e-12
        @test report.logical_block_sizes == [4, 3, 3]
        @test report.psd_block_sizes == [8, 6, 6]
        @test occursin("logical_block_sizes", sprint(show, report))
        @test report.product_cache_hits > 0
        @test report.product_cache_misses > 0
        @test report.product_cache_entries == report.product_cache_misses
        @test report.product_cache_hit_rate > 0.0
        @test occursin("product_cache_hit_rate", sprint(show, report))
        @test occursin("product_cache_entries", sprint(show, report))
        @test report.construction_time_ns > 0
        @test all(
            key -> haskey(report.construction_stage_times_ns, key),
            (:setup, :basis, :precompute, :block_assembly, :constraint_append, :linearization),
        )
        @test all(value >= 0 for value in values(report.construction_stage_times_ns))
        @test occursin("construction_time_ns", sprint(show, report))
        @test report.real_moment_matrix
        support = translation_solve_support(report)
        @test support.supported
        @test support.reason == ""
        @test support.blocker === nothing

        n_large = 12
        registry_large, ops_large = create_pauli_variables(1:n_large)
        pop_large = polyopt(heisenberg_chain_hamiltonian(ops_large), registry_large)
        _, report_large = pauli_translation_invariant_moment_relaxation(pop_large, ops_large, 4)

        @test report_large.basis_size == 1 + n_large * sum(3^ℓ for ℓ in 1:4)
        @test report_large.orbit_basis_size == 1 + sum(3^ℓ for ℓ in 1:4)
        @test report_large.axis_orbit_closed
        @test report_large.axis_orbit_basis_size == 23
        @test report_large.axis_orbit_size_histogram == [1 => 1, 3 => 4, 6 => 18]
        @test report_large.axis_reduction_ratio == 121 / 23
        @test maximum(report_large.logical_block_sizes) == 31
        @test maximum(report_large.psd_block_sizes) == 62
        @test length(report_large.psd_block_sizes) == 4 * (fld(n_large, 2) + 1)
        @test isdefined(NCTSSoS, :translation_block_histogram)
        if isdefined(NCTSSoS, :translation_block_histogram)
            @test translation_block_histogram(report_large; kind=:logical) == [30 => 27, 31 => 1]
            @test translation_block_histogram(report_large; kind=:psd) == [60 => 27, 62 => 1]
            @test_throws ArgumentError translation_block_histogram(report_large; kind=:solver)
        end
        @test isdefined(NCTSSoS, :translation_block_feature_histogram)
        if isdefined(NCTSSoS, :translation_block_feature_histogram)
            @test translation_block_feature_histogram(report_large; kind=:logical) == [
                (feature=:moment_matrix, decomposition=:translation, size=30) => 27,
                (feature=:moment_matrix, decomposition=:translation, size=31) => 1,
            ]
            @test translation_block_feature_histogram(report_large; kind=:psd) == [
                (feature=:moment_matrix, decomposition=:translation, size=60) => 27,
                (feature=:moment_matrix, decomposition=:translation, size=62) => 1,
            ]
            @test_throws ArgumentError translation_block_feature_histogram(report_large; kind=:solver)
        end
        @test isdefined(NCTSSoS, :translation_report_metrics)
        if isdefined(NCTSSoS, :translation_report_metrics)
            metrics = translation_report_metrics(report_large)
            @test metrics.n_sites == n_large
            @test metrics.order == 4
            @test metrics.basis_size == report_large.basis_size
            @test metrics.orbit_basis_size == report_large.orbit_basis_size
            @test metrics.axis_orbit_closed == report_large.axis_orbit_closed
            @test metrics.axis_orbit_basis_size == report_large.axis_orbit_basis_size == 23
            @test metrics.axis_orbit_size_histogram == report_large.axis_orbit_size_histogram
            @test metrics.axis_reduction_ratio == report_large.axis_reduction_ratio
            @test metrics.momentum_sector_count == length(report_large.momentum_sectors)
            @test metrics.n_blocks == 28
            @test metrics.moment_count == report_large.n_unique_moment_matrix_elements
            @test metrics.base_moment_count == report_large.n_unique_moment_matrix_elements
            @test metrics.linear_moment_count == report_large.linear_moment_count
            @test metrics.zero_constraint_count == report_large.zero_constraint_count == 0
            @test metrics.estimated_sos_dual_scalar_equalities_upper_bound ==
                metrics.linear_moment_count - 1
            @test metrics.estimated_sos_dual_dense_schur_bytes ==
                8 * metrics.estimated_sos_dual_scalar_equalities_upper_bound^2
            @test translation_report_metrics(
                report_large;
                scalar_bytes=16,
            ).estimated_sos_dual_dense_schur_bytes ==
                16 * metrics.estimated_sos_dual_scalar_equalities_upper_bound^2
            @test metrics.solve_supported
            @test metrics.solve_blocker === nothing
            @test metrics.solve_blocker_reason == ""
            @test metrics.construction_time_ns == report_large.construction_time_ns
            @test metrics.construction_time_seconds == report_large.construction_time_ns / 1e9
            @test metrics.construction_stage_times_ns == report_large.construction_stage_times_ns
            @test metrics.construction_stage_time_seconds ==
                Dict(key => value / 1e9 for (key, value) in report_large.construction_stage_times_ns)
            @test report_large.linear_moment_count == report_large.n_unique_moment_matrix_elements
            @test metrics.logical_max_block == 31
            @test metrics.psd_max_block == 62
            @test metrics.logical_total_block_side == 841
            @test metrics.psd_total_block_side == 1682
            @test metrics.psd_dense_entries == 101_044
            @test metrics.psd_symmetric_entries == 51_363
            @test metrics.psd_dense_bytes == 808_352
            @test metrics.psd_symmetric_bytes == 410_904
            @test metrics.product_cache_entries == report_large.product_cache_entries
            @test metrics.product_cache_entries == report_large.product_cache_misses
            @test metrics.product_cache_lookups ==
                report_large.product_cache_hits + report_large.product_cache_misses
            @test metrics.product_cache_hit_rate == report_large.product_cache_hit_rate
            @test metrics.product_cache_hit_rate > 0.8
            @test metrics.logical_block_histogram == [30 => 27, 31 => 1]
            @test metrics.psd_block_histogram == [60 => 27, 62 => 1]
            @test metrics.logical_block_feature_histogram == [
                (feature=:moment_matrix, decomposition=:translation, size=30) => 27,
                (feature=:moment_matrix, decomposition=:translation, size=31) => 1,
            ]
            @test metrics.psd_block_feature_histogram == [
                (feature=:moment_matrix, decomposition=:translation, size=60) => 27,
                (feature=:moment_matrix, decomposition=:translation, size=62) => 1,
            ]
            @test_throws ArgumentError translation_report_metrics(report_large; scalar_bytes=0)
        end

        @test isdefined(NCTSSoS, :pauli_translation_structural_targets)
        if isdefined(NCTSSoS, :pauli_translation_structural_targets)
            target_large = pauli_translation_structural_targets(100, 4)
            @test target_large.n_sites == 100
            @test target_large.order == 4
            @test target_large.basis_size == 12_001
            @test target_large.orbit_basis_size == 121
            @test target_large.momentum_sector_count == 51
            @test target_large.n_blocks == 204
            @test target_large.logical_max_block == 31
            @test target_large.psd_max_block == 62
            @test target_large.product_cache_misses == 1_891
            @test target_large.product_cache_entries == target_large.product_cache_misses
            @test target_large.product_cache_hits == 93_000
            @test target_large.product_cache_lookups == 94_891
            @test target_large.product_cache_hit_rate > 0.98
            @test target_large.hamiltonian_width == 2
            @test target_large.max_contiguous_rdm_k == 8
            @test target_large.max_linear_state_opt_width == 7
            @test target_large.max_psd_state_opt_width == 3
            @test target_large.axis_orbit_closed
            @test target_large.axis_orbit_basis_size == 23
            @test target_large.axis_orbit_size_histogram == [1 => 1, 3 => 4, 6 => 18]
            @test target_large.axis_reduction_ratio == 121 / 23
            @test target_large.block_coefficient_domains ==
                fill(:cyclotomic_float64, target_large.n_blocks)
            @test target_large.block_exact_coefficient_domains ==
                fill(:cyclotomic, target_large.n_blocks)
            @test target_large.block_coefficient_domain_histogram ==
                [:cyclotomic_float64 => target_large.n_blocks]
            @test target_large.block_exact_coefficient_domain_histogram ==
                [:cyclotomic => target_large.n_blocks]
            @test !target_large.solve_supported
            @test target_large.solve_blocker == :structural_target_only
            @test occursin(
                "Structural target only",
                target_large.solve_blocker_reason,
            )
            @test isempty(target_large.solve_unsupported_block_features)
            @test isempty(target_large.solve_unsupported_zero_features)
            @test !target_large.requires_construction

            for candidate_n in (10, n_large)
                candidate_report = if candidate_n == n_large
                    report_large
                else
                    candidate_registry, candidate_ops = create_pauli_variables(1:candidate_n)
                    candidate_pop = polyopt(
                        heisenberg_chain_hamiltonian(candidate_ops),
                        candidate_registry,
                    )
                    last(pauli_translation_invariant_moment_relaxation(
                        candidate_pop,
                        candidate_ops,
                        4,
                    ))
                end
                target = pauli_translation_structural_targets(candidate_n, 4)
                metrics = translation_report_metrics(candidate_report)
                @test target.basis_size == metrics.basis_size
                @test target.orbit_basis_size == metrics.orbit_basis_size
                @test target.momentum_sector_count == metrics.momentum_sector_count
                @test target.n_blocks == metrics.n_blocks
                @test target.logical_max_block == metrics.logical_max_block
                @test target.psd_max_block == metrics.psd_max_block
                @test target.logical_total_block_side == metrics.logical_total_block_side
                @test target.psd_total_block_side == metrics.psd_total_block_side
                @test target.logical_block_histogram == metrics.logical_block_histogram
                @test target.psd_block_histogram == metrics.psd_block_histogram
                @test target.psd_symmetric_entries == metrics.psd_symmetric_entries
                @test target.product_cache_hits == metrics.product_cache_hits
                @test target.product_cache_misses == metrics.product_cache_misses
                @test target.product_cache_lookups == metrics.product_cache_lookups
                @test target.product_cache_entries == metrics.product_cache_entries
                @test target.product_cache_hit_rate == metrics.product_cache_hit_rate
                @test target.product_cache_hit_rate > 0.8
                @test target.axis_orbit_basis_size == metrics.axis_orbit_basis_size
                @test target.axis_orbit_size_histogram == metrics.axis_orbit_size_histogram
            end

            addon_registry, addon_ops = create_pauli_variables(1:8)
            addon_pop = polyopt(heisenberg_chain_hamiltonian(addon_ops), addon_registry)
            _, addon_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                contiguous_rdm_k=2,
                psd_state_opt_width=1,
            )
            addon_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                contiguous_rdm_k=2,
                psd_state_opt_width=1,
            )
            addon_metrics = translation_report_metrics(addon_report)
            @test addon_target.logical_block_sizes == addon_report.logical_block_sizes
            @test addon_target.psd_block_sizes == addon_report.psd_block_sizes
            @test addon_target.n_blocks == addon_metrics.n_blocks
            @test addon_target.logical_max_block == addon_metrics.logical_max_block
            @test addon_target.psd_max_block == addon_metrics.psd_max_block
            @test addon_target.psd_symmetric_entries == addon_metrics.psd_symmetric_entries
            @test addon_target.logical_block_feature_histogram ==
                addon_metrics.logical_block_feature_histogram
            @test addon_target.psd_block_feature_histogram ==
                addon_metrics.psd_block_feature_histogram
            @test addon_target.psd_state_opt_candidate_count == 3
            @test addon_target.psd_state_opt_logical_block_sizes == fill(3, 5)
            @test addon_target.psd_state_opt_psd_block_sizes == [3, 6, 6, 6, 6]
            @test addon_target.rdm_logical_block_sizes == [4]
            @test addon_target.rdm_psd_block_sizes == [8]
            complex_full_rdm_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                contiguous_rdm_k=2,
            )
            @test complex_full_rdm_target.rdm_logical_block_sizes == [4]
            @test complex_full_rdm_target.rdm_psd_block_sizes == [4]
            @test complex_full_rdm_target.contiguous_rdm_zero_row_count == 0
            @test complex_full_rdm_target.add_on_zero_row_count == 0
            @test isempty(complex_full_rdm_target.known_zero_constraint_feature_histogram)
            @test complex_full_rdm_target.assumptions.complete_add_on_zero_row_count
            _, complex_full_rdm_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                contiguous_rdm_k=2,
            )
            complex_full_rdm_metrics = translation_report_metrics(
                complex_full_rdm_report,
            )
            @test complex_full_rdm_target.logical_block_sizes ==
                complex_full_rdm_report.logical_block_sizes
            @test complex_full_rdm_target.psd_block_sizes ==
                complex_full_rdm_report.psd_block_sizes
            @test complex_full_rdm_target.psd_symmetric_entries ==
                complex_full_rdm_metrics.psd_symmetric_entries
            @test complex_full_rdm_target.psd_dense_entries ==
                complex_full_rdm_metrics.psd_dense_entries
            @test complex_full_rdm_target.psd_symmetric_entries ==
                complex_full_rdm_target.psd_dense_entries
            @test complex_full_rdm_metrics.contiguous_rdm_zero_row_count == 0
            complex_psd_state_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                psd_state_opt_width=1,
            )
            @test complex_psd_state_target.psd_state_opt_logical_block_sizes == fill(3, 8)
            @test complex_psd_state_target.psd_state_opt_psd_block_sizes == fill(3, 8)
            @test complex_psd_state_target.psd_symmetric_entries ==
                complex_psd_state_target.psd_dense_entries
            _, complex_psd_state_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                psd_state_opt_width=1,
            )
            complex_psd_state_metrics = translation_report_metrics(complex_psd_state_report)
            @test complex_psd_state_target.logical_block_sizes ==
                complex_psd_state_report.logical_block_sizes
            @test complex_psd_state_target.psd_block_sizes ==
                complex_psd_state_report.psd_block_sizes
            @test complex_psd_state_target.psd_symmetric_entries ==
                complex_psd_state_metrics.psd_symmetric_entries
            @test complex_psd_state_target.psd_dense_entries ==
                complex_psd_state_metrics.psd_dense_entries

            qmb_reference_target = pauli_translation_structural_targets(
                100,
                4;
                contiguous_rdm_k=8,
            )
            @test length(qmb_reference_target.qmbcertify_rdm_reference_blocks) == 1
            @test only(qmb_reference_target.qmbcertify_rdm_reference_blocks).k == 8
            @test only(qmb_reference_target.qmbcertify_rdm_reference_blocks).block_sizes ==
                [72, 64, 56]
            @test only(qmb_reference_target.qmbcertify_rdm_reference_blocks).symmetric_entries ==
                6_304
            @test !qmb_reference_target.solve_supported
            @test qmb_reference_target.solve_blocker == :structural_target_only
            @test isempty(qmb_reference_target.solve_unsupported_block_features)
            @test isempty(qmb_reference_target.solve_unsupported_zero_features)
            @test isempty(
                pauli_translation_structural_targets(
                    100,
                    4;
                    contiguous_rdm_k=7,
                ).qmbcertify_rdm_reference_blocks,
            )
            qmb_decomposition_target = pauli_translation_structural_targets(
                100,
                4;
                contiguous_rdm_k=8,
                contiguous_rdm_decomposition=:qmbcertify,
                contiguous_rdm_support=:extend,
                psd_state_opt_width=3,
            )
            @test qmb_decomposition_target.contiguous_rdm_decomposition == :qmbcertify
            @test qmb_decomposition_target.rdm_logical_block_sizes == [72, 64, 56]
            @test qmb_decomposition_target.rdm_psd_block_sizes == [72, 64, 56]
            @test qmb_decomposition_target.qmbcertify_rdm_reference_block_sizes == [72, 64, 56]
            @test qmb_decomposition_target.qmbcertify_rdm_reference_max_block == 72
            @test qmb_decomposition_target.qmbcertify_rdm_reference_dense_entries == 12_416
            @test qmb_decomposition_target.qmbcertify_rdm_reference_symmetric_entries == 6_304
            @test qmb_decomposition_target.psd_max_block == 72
            @test qmb_decomposition_target.n_blocks == 258
            @test qmb_decomposition_target.psd_dense_entries == 763_341
            @test qmb_decomposition_target.psd_symmetric_entries == 388_342
            @test qmb_decomposition_target.block_coefficient_domain_histogram ==
                [:cyclotomic_float64 => 255, nothing => 3]
            @test qmb_decomposition_target.block_exact_coefficient_domain_histogram ==
                [:cyclotomic => 255, nothing => 3]
            @test !qmb_decomposition_target.solve_supported
            @test qmb_decomposition_target.solve_blocker == :structural_target_only
            @test isempty(qmb_decomposition_target.solve_unsupported_block_features)
            @test isempty(qmb_decomposition_target.solve_unsupported_zero_features)
            qmb_decomposition_target_all = pauli_translation_structural_targets(
                100,
                5;
                contiguous_rdm_k=[8, 9, 10],
                contiguous_rdm_decomposition=:qmbcertify,
            )
            @test qmb_decomposition_target_all.rdm_logical_block_sizes ==
                [72, 64, 56, 136, 120, 256, 272, 240]
            @test qmb_decomposition_target_all.qmbcertify_rdm_reference_block_sizes ==
                [72, 64, 56, 136, 120, 256, 272, 240]
            @test qmb_decomposition_target_all.qmbcertify_rdm_reference_max_block == 272
            @test qmb_decomposition_target_all.qmbcertify_rdm_reference_dense_entries == 242_432
            @test qmb_decomposition_target_all.qmbcertify_rdm_reference_symmetric_entries == 121_824
            @test !qmb_decomposition_target_all.qmbcertify_rdm_reference_requires_construction
            @test_throws ArgumentError pauli_translation_structural_targets(
                100,
                4;
                real_moment_matrix=false,
                contiguous_rdm_k=8,
                contiguous_rdm_decomposition=:qmbcertify,
            )
            qmb_base_target = pauli_translation_structural_targets(
                20,
                4;
                qmbcertify_base_extra=9,
            )
            @test qmb_base_target.qmbcertify_base_reference_block_sizes == [
                58,
                fill(114, 9)...,
                57,
                48,
                fill(96, 9)...,
                48,
            ]
            @test qmb_base_target.qmbcertify_base_reference_block_histogram ==
                [48 => 2, 57 => 1, 58 => 1, 96 => 9, 114 => 9]
            @test qmb_base_target.qmbcertify_base_reference_n_blocks == 22
            @test qmb_base_target.qmbcertify_base_reference_max_block == 114
            @test qmb_base_target.qmbcertify_base_reference_dense_entries == 211_129
            @test qmb_base_target.qmbcertify_base_reference_symmetric_entries == 106_615
            @test qmb_base_target.qmbcertify_base_reference_support_unique_count == 3_420
            @test qmb_base_target.qmbcertify_base_reference_support_nonzero_row_count == 37_298
            @test qmb_base_target.qmbcertify_base_reference_support_zero_row_count == 18_232
            @test qmb_base_target.qmbcertify_base_reference_support_word_length_histogram ==
                [0 => 1, 2 => 10, 4 => 575, 6 => 1816, 8 => 1018]
            @test qmb_base_target.qmbcertify_base_reference_active == false
            @test qmb_base_target.qmbcertify_base_reference_requires_construction
            qmb_active_target = pauli_translation_structural_targets(
                20,
                4;
                qmbcertify_base_construct=true,
                qmbcertify_base_extra=9,
            )
            @test qmb_active_target.qmbcertify_base_reference_active
            @test qmb_active_target.basis_size == 2101
            @test qmb_active_target.orbit_basis_size == 106
            @test qmb_active_target.base_psd_block_sizes ==
                qmb_base_target.qmbcertify_base_reference_block_sizes
            @test qmb_active_target.psd_block_sizes ==
                qmb_base_target.qmbcertify_base_reference_block_sizes
            @test qmb_active_target.psd_block_histogram ==
                [48 => 2, 57 => 1, 58 => 1, 96 => 9, 114 => 9]
            @test qmb_active_target.estimated_model_size_gate_status ==
                :blocked_missing_scalar_equality_estimate
            qmb_active_rdm_target = pauli_translation_structural_targets(
                20,
                4;
                qmbcertify_base_construct=true,
                qmbcertify_base_extra=9,
                contiguous_rdm_k=8,
                contiguous_rdm_decomposition=:qmbcertify,
                contiguous_rdm_support=:extend,
            )
            @test qmb_active_rdm_target.qmbcertify_base_reference_active
            @test qmb_active_rdm_target.rdm_psd_block_sizes == [72, 64, 56]
            @test qmb_active_rdm_target.psd_block_sizes ==
                vcat(qmb_base_target.qmbcertify_base_reference_block_sizes, [72, 64, 56])
            qmb_active_pso_target = pauli_translation_structural_targets(
                20,
                4;
                qmbcertify_base_construct=true,
                qmbcertify_base_extra=9,
                psd_state_opt_width=3,
            )
            qmb_pso_expected_blocks20 = [
                36,
                fill(72, 9)...,
                36,
                28,
                fill(56, 9)...,
                28,
            ]
            @test qmb_active_pso_target.qmbcertify_base_reference_active
            @test qmb_active_pso_target.psd_state_opt_psd_block_sizes ==
                qmb_pso_expected_blocks20
            @test qmb_active_pso_target.psd_state_opt_candidate_count == 64
            @test qmb_active_pso_target.psd_block_sizes ==
                vcat(qmb_base_target.qmbcertify_base_reference_block_sizes, qmb_pso_expected_blocks20)
            qmb_active_lso_target = pauli_translation_structural_targets(
                20,
                4;
                qmbcertify_base_construct=true,
                qmbcertify_base_extra=9,
                linear_state_opt_width=7,
            )
            @test qmb_active_lso_target.qmbcertify_base_reference_active
            @test qmb_active_lso_target.linear_state_opt_candidate_count == 1_817
            @test qmb_active_lso_target.linear_state_opt_row_count == 170
            @test qmb_active_lso_target.add_on_zero_row_count == 170
            @test qmb_active_lso_target.known_zero_constraint_feature_histogram == [
                (feature=:linear_state_opt, decomposition=nothing, reason=nothing) => 170,
            ]
            @test qmb_active_lso_target.psd_block_sizes ==
                qmb_base_target.qmbcertify_base_reference_block_sizes
            qmb_active_combined_target = pauli_translation_structural_targets(
                20,
                4;
                qmbcertify_base_construct=true,
                qmbcertify_base_extra=9,
                contiguous_rdm_k=8,
                contiguous_rdm_decomposition=:qmbcertify,
                contiguous_rdm_support=:extend,
                linear_state_opt_width=7,
                psd_state_opt_width=3,
            )
            @test qmb_active_combined_target.qmbcertify_base_reference_active
            @test qmb_active_combined_target.rdm_psd_block_sizes == [72, 64, 56]
            @test qmb_active_combined_target.psd_state_opt_psd_block_sizes ==
                qmb_pso_expected_blocks20
            @test qmb_active_combined_target.linear_state_opt_candidate_count == 1_817
            @test qmb_active_combined_target.linear_state_opt_row_count == 336
            @test qmb_active_combined_target.add_on_zero_row_count == 336
            @test qmb_active_combined_target.known_zero_constraint_feature_histogram == [
                (feature=:linear_state_opt, decomposition=nothing, reason=nothing) => 336,
            ]
            @test qmb_active_combined_target.psd_block_sizes ==
                vcat(
                    qmb_base_target.qmbcertify_base_reference_block_sizes,
                    qmb_pso_expected_blocks20,
                    [72, 64, 56],
                )

            linear_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                linear_state_opt_width=2,
            )
            _, linear_target_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                linear_state_opt_width=2,
            )
            linear_target_metrics = translation_report_metrics(linear_target_report)
            @test linear_target.linear_state_opt_candidate_count == 12
            @test linear_target.linear_state_opt_row_count ==
                linear_target_metrics.linear_state_opt_row_count
            @test linear_target.add_on_zero_row_count ==
                linear_target.linear_state_opt_row_count
            @test linear_target.known_zero_constraint_feature_histogram == [
                (feature=:linear_state_opt, decomposition=nothing, reason=nothing) =>
                    linear_target.linear_state_opt_row_count,
            ]
            @test linear_target.assumptions.complete_add_on_zero_row_count
            @test linear_target.logical_block_sizes == linear_target.base_logical_block_sizes
            @test linear_target.psd_block_sizes == linear_target.base_psd_block_sizes
            complex_linear_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                linear_state_opt_width=2,
            )
            _, complex_linear_target_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                linear_state_opt_width=2,
            )
            complex_linear_target_metrics = translation_report_metrics(
                complex_linear_target_report,
            )
            @test complex_linear_target.linear_state_opt_row_count ==
                complex_linear_target_metrics.linear_state_opt_row_count
            @test complex_linear_target.add_on_zero_row_count ==
                complex_linear_target.linear_state_opt_row_count
            @test complex_linear_target.known_zero_constraint_feature_histogram == [
                (feature=:linear_state_opt, decomposition=nothing, reason=nothing) =>
                    complex_linear_target.linear_state_opt_row_count,
            ]
            @test complex_linear_target.logical_block_sizes ==
                complex_linear_target_report.logical_block_sizes
            @test complex_linear_target.psd_block_sizes ==
                complex_linear_target_report.psd_block_sizes
            @test complex_linear_target.psd_symmetric_entries ==
                complex_linear_target_metrics.psd_symmetric_entries
            @test complex_linear_target.psd_dense_entries ==
                complex_linear_target_metrics.psd_dense_entries
            sign_linear_target = pauli_translation_structural_targets(
                8,
                2;
                linear_state_opt_width=2,
            )
            @test sign_linear_target.sign_symmetry
            @test sign_linear_target.linear_state_opt_candidate_count == 3
            @test sign_linear_target.logical_block_sizes ==
                sign_linear_target.base_logical_block_sizes
            @test sign_linear_target.psd_block_sizes ==
                sign_linear_target.base_psd_block_sizes
            M_large = eltype(ops_large[1])
            raw_state_opt_tests = NCTSSoS._contiguous_state_opt_tests(
                M_large,
                8,
                2;
                sign_symmetry=false,
            )
            sign_state_opt_tests = NCTSSoS._contiguous_state_opt_tests(
                M_large,
                8,
                2;
                sign_symmetry=true,
            )
            @test length(raw_state_opt_tests) == linear_target.linear_state_opt_candidate_count
            @test length(sign_state_opt_tests) ==
                sign_linear_target.linear_state_opt_candidate_count
            @test all(
                mono -> NCTSSoS._pauli_sign_signature(mono) == 0x00,
                sign_state_opt_tests,
            )

            h2_hamiltonian = heisenberg_chain_hamiltonian(addon_ops)
            h2_pop = polyopt(
                h2_hamiltonian,
                addon_registry;
                moment_eq_constraints=[h2_hamiltonian * h2_hamiltonian],
            )
            h2_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                moment_eq_h2=true,
            )
            _, h2_target_report = pauli_translation_invariant_moment_relaxation(
                h2_pop,
                addon_ops,
                2;
                sign_symmetry=false,
            )
            h2_target_metrics = translation_report_metrics(h2_target_report)
            @test h2_target.moment_equality_row_count ==
                h2_target_metrics.moment_equality_row_count
            @test h2_target.moment_equality_row_count > 0
            @test h2_target.add_on_zero_row_count == h2_target.moment_equality_row_count
            @test h2_target.known_zero_constraint_feature_histogram == [
                (feature=:moment_equality, decomposition=nothing, reason=nothing) =>
                    h2_target.moment_equality_row_count,
            ]
            @test h2_target.assumptions.complete_add_on_zero_row_count
            complex_h2_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                moment_eq_h2=true,
            )
            _, complex_h2_target_report = pauli_translation_invariant_moment_relaxation(
                h2_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
            )
            complex_h2_target_metrics = translation_report_metrics(complex_h2_target_report)
            @test complex_h2_target.moment_equality_row_count ==
                complex_h2_target_metrics.moment_equality_row_count
            @test complex_h2_target.moment_equality_row_count > 0
            @test complex_h2_target.add_on_zero_row_count ==
                complex_h2_target.moment_equality_row_count
            @test complex_h2_target.known_zero_constraint_feature_histogram == [
                (feature=:moment_equality, decomposition=nothing, reason=nothing) =>
                    complex_h2_target.moment_equality_row_count,
            ]
            @test complex_h2_target.logical_block_sizes ==
                complex_h2_target_report.logical_block_sizes
            @test complex_h2_target.psd_block_sizes ==
                complex_h2_target_report.psd_block_sizes
            @test complex_h2_target.psd_symmetric_entries ==
                complex_h2_target_metrics.psd_symmetric_entries
            @test complex_h2_target.psd_dense_entries ==
                complex_h2_target_metrics.psd_dense_entries
            h2_order3_orbit_reps =
                NCTSSoS._pauli_contiguous_chain_orbit_representatives(
                    addon_ops,
                    3;
                    periodic=true,
                )
            h2_order3_raw_rows, h2_order3_raw_degrees =
                NCTSSoS._translation_moment_eq_row_bases(
                    h2_order3_orbit_reps;
                    sign_symmetry=false,
                )
            h2_order3_sign_rows, h2_order3_sign_degrees =
                NCTSSoS._translation_moment_eq_row_bases(
                    h2_order3_orbit_reps;
                    sign_symmetry=true,
                )
            @test length(h2_order3_sign_rows) < length(h2_order3_raw_rows)
            @test h2_order3_raw_degrees == degree.(h2_order3_raw_rows)
            @test h2_order3_sign_degrees == degree.(h2_order3_sign_rows)
            @test all(
                mono -> NCTSSoS._pauli_sign_signature(mono) == 0x00,
                h2_order3_sign_rows,
            )

            axis_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                axis_rotation_equalities=true,
            )
            _, axis_target_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                axis_rotation_equalities=true,
            )
            axis_target_metrics = translation_report_metrics(axis_target_report)
            @test axis_target.axis_rotation_equalities
            @test axis_target.axis_rotation_equality_row_count ==
                axis_target_metrics.axis_rotation_equality_row_count
            @test axis_target.axis_rotation_equality_row_count > 0
            @test axis_target.axis_rotation_raw_moment_key_count == 379
            @test axis_target.axis_rotation_moment_class_count == 73
            @test axis_target.axis_rotation_quotient_moment_key_count == 22
            @test axis_target.axis_rotation_forced_zero_moment_class_count == 51
            @test axis_target_metrics.axis_rotation_raw_moment_key_count ==
                axis_target.axis_rotation_raw_moment_key_count
            @test axis_target_metrics.axis_rotation_moment_class_count ==
                axis_target.axis_rotation_moment_class_count
            @test axis_target_metrics.axis_rotation_quotient_moment_key_count ==
                axis_target.axis_rotation_quotient_moment_key_count
            @test axis_target_metrics.axis_rotation_forced_zero_moment_class_count ==
                axis_target.axis_rotation_forced_zero_moment_class_count
            @test axis_target.axis_rotation_moment_quotient_reduction_ratio ≈
                axis_target.axis_rotation_raw_moment_key_count /
                axis_target.axis_rotation_quotient_moment_key_count
            @test axis_target.add_on_zero_row_count ==
                axis_target.axis_rotation_equality_row_count
            @test axis_target.known_zero_constraint_feature_histogram == [
                (feature=:axis_rotation_equality, decomposition=nothing, reason=nothing) =>
                    axis_target.axis_rotation_equality_row_count,
            ]
            complex_axis_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                axis_rotation_equalities=true,
            )
            _, complex_axis_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                axis_rotation_equalities=true,
            )
            complex_axis_metrics = translation_report_metrics(complex_axis_report)
            @test complex_axis_target.axis_rotation_equalities
            @test complex_axis_target.axis_rotation_equality_row_count ==
                complex_axis_metrics.axis_rotation_equality_row_count
            @test complex_axis_target.axis_rotation_equality_row_count > 0
            @test complex_axis_metrics.axis_rotation_raw_moment_key_count ==
                complex_axis_target.axis_rotation_raw_moment_key_count
            @test complex_axis_metrics.axis_rotation_moment_class_count ==
                complex_axis_target.axis_rotation_moment_class_count
            @test complex_axis_metrics.axis_rotation_quotient_moment_key_count ==
                complex_axis_target.axis_rotation_quotient_moment_key_count
            @test complex_axis_metrics.axis_rotation_forced_zero_moment_class_count ==
                complex_axis_target.axis_rotation_forced_zero_moment_class_count
            @test complex_axis_target.add_on_zero_row_count ==
                complex_axis_target.axis_rotation_equality_row_count
            @test complex_axis_target.known_zero_constraint_feature_histogram == [
                (feature=:axis_rotation_equality, decomposition=nothing, reason=nothing) =>
                    complex_axis_target.axis_rotation_equality_row_count,
            ]
            @test complex_axis_target.logical_block_sizes ==
                complex_axis_report.logical_block_sizes
            @test complex_axis_target.psd_block_sizes ==
                complex_axis_report.psd_block_sizes
            @test complex_axis_target.psd_symmetric_entries ==
                complex_axis_metrics.psd_symmetric_entries
            @test complex_axis_target.psd_dense_entries ==
                complex_axis_metrics.psd_dense_entries
            default_axis_target = pauli_translation_structural_targets(
                8,
                2;
                axis_rotation_equalities=true,
            )
            _, default_axis_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                axis_rotation_equalities=true,
            )
            default_axis_metrics = translation_report_metrics(default_axis_report)
            @test default_axis_target.sign_symmetry
            @test default_axis_target.axis_rotation_equality_row_count ==
                default_axis_metrics.axis_rotation_equality_row_count
            @test 0 < default_axis_target.axis_rotation_equality_row_count <
                axis_target.axis_rotation_equality_row_count
            @test default_axis_target.axis_rotation_raw_moment_key_count == 100
            @test default_axis_target.axis_rotation_moment_class_count == 22
            @test default_axis_target.axis_rotation_quotient_moment_key_count == 22
            @test default_axis_target.axis_rotation_forced_zero_moment_class_count == 0
            @test default_axis_metrics.axis_rotation_raw_moment_key_count ==
                default_axis_target.axis_rotation_raw_moment_key_count
            @test default_axis_metrics.axis_rotation_moment_class_count ==
                default_axis_target.axis_rotation_moment_class_count
            @test default_axis_metrics.axis_rotation_quotient_moment_key_count ==
                default_axis_target.axis_rotation_quotient_moment_key_count
            @test default_axis_metrics.axis_rotation_forced_zero_moment_class_count ==
                default_axis_target.axis_rotation_forced_zero_moment_class_count

            axis_split_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                axis_rotation_symmetry=true,
            )
            _, axis_split_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                axis_rotation_symmetry=true,
            )
            axis_split_metrics = translation_report_metrics(axis_split_report)
            @test axis_split_target.axis_rotation_symmetry
            @test axis_split_target.axis_rotation_equalities
            @test axis_split_target.base_logical_block_sizes ==
                axis_split_report.logical_block_sizes
            @test axis_split_target.base_psd_block_sizes ==
                axis_split_report.psd_block_sizes
            @test axis_split_target.axis_rotation_equality_row_count ==
                axis_split_metrics.axis_rotation_equality_row_count
            @test axis_split_metrics.axis_rotation_raw_moment_key_count ==
                axis_split_target.axis_rotation_raw_moment_key_count
            @test axis_split_metrics.axis_rotation_quotient_moment_key_count ==
                axis_split_target.axis_rotation_quotient_moment_key_count
            @test axis_split_metrics.axis_rotation_forced_zero_moment_class_count ==
                axis_split_target.axis_rotation_forced_zero_moment_class_count
            @test axis_split_target.add_on_zero_row_count ==
                axis_split_target.axis_rotation_equality_row_count

            reflected_axis_split_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                reflection_symmetry=true,
                axis_rotation_symmetry=true,
            )
            _, reflected_axis_split_report =
                pauli_translation_invariant_moment_relaxation(
                    addon_pop,
                    addon_ops,
                    2;
                    sign_symmetry=false,
                    reflection_symmetry=true,
                    axis_rotation_symmetry=true,
                )
            reflected_axis_split_metrics =
                translation_report_metrics(reflected_axis_split_report)
            @test reflected_axis_split_target.reflection_symmetry
            @test reflected_axis_split_target.axis_rotation_symmetry
            @test reflected_axis_split_target.axis_rotation_equalities
            @test reflected_axis_split_target.base_logical_block_sizes ==
                reflected_axis_split_report.logical_block_sizes
            @test reflected_axis_split_target.base_psd_block_sizes ==
                reflected_axis_split_report.psd_block_sizes
            @test reflected_axis_split_target.axis_rotation_equality_row_count ==
                reflected_axis_split_metrics.axis_rotation_equality_row_count
            @test all(
                pair -> pair.first.decomposition == :translation_reflection_axis_irrep,
                reflected_axis_split_metrics.logical_block_feature_histogram,
            )
            @test all(
                pair -> pair.first.decomposition == :translation_reflection_axis_irrep,
                reflected_axis_split_metrics.psd_block_feature_histogram,
            )
            @test reflected_axis_split_target.add_on_zero_row_count ==
                reflected_axis_split_target.axis_rotation_equality_row_count
            reflected_axis_order4_target = pauli_translation_structural_targets(
                10,
                4;
                reflection_symmetry=true,
                axis_rotation_symmetry=true,
                linear_state_opt_width=7,
                psd_state_opt_width=2,
            )
            @test reflected_axis_order4_target.reflection_symmetry
            @test reflected_axis_order4_target.axis_rotation_symmetry
            @test reflected_axis_order4_target.axis_rotation_equalities
            @test length(reflected_axis_order4_target.base_psd_block_sizes) == 26
            @test maximum(reflected_axis_order4_target.base_psd_block_sizes) == 120
            @test count(==(120), reflected_axis_order4_target.base_psd_block_sizes) == 8
            @test reflected_axis_order4_target.axis_rotation_equality_row_count == 29_444
            @test reflected_axis_order4_target.axis_rotation_raw_moment_key_count == 23_608
            @test reflected_axis_order4_target.axis_rotation_moment_class_count == 3_986
            @test reflected_axis_order4_target.axis_rotation_quotient_moment_key_count == 1_017
            @test reflected_axis_order4_target.axis_rotation_forced_zero_moment_class_count == 2_969
            @test reflected_axis_order4_target.axis_rotation_moment_quotient_reduction_ratio ≈
                reflected_axis_order4_target.axis_rotation_raw_moment_key_count /
                reflected_axis_order4_target.axis_rotation_quotient_moment_key_count
            @test reflected_axis_order4_target.linear_state_opt_row_count == 819
            @test length(reflected_axis_order4_target.psd_state_opt_psd_block_sizes) == 6
            @test reflected_axis_order4_target.block_coefficient_domain_histogram ==
                [:cyclotomic_float64 => 6, :real_algebraic_float64 => 26]

            u1_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
            )
            @test u1_target.rdm_logical_block_sizes == [1, 3, 3, 1]
            @test u1_target.rdm_psd_block_sizes == [1, 6, 6, 1]
            _, u1_target_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
            )
            u1_target_metrics = translation_report_metrics(u1_target_report)
            @test u1_target.contiguous_rdm_zero_row_count ==
                u1_target_metrics.contiguous_rdm_zero_row_count
            @test u1_target.contiguous_rdm_zero_row_count == 88
            @test u1_target.add_on_zero_row_count == u1_target.contiguous_rdm_zero_row_count
            @test u1_target.known_zero_constraint_feature_histogram == [
                (feature=:contiguous_rdm_zero, decomposition=:u1, reason=:magnetization_offblock) =>
                    u1_target.contiguous_rdm_zero_row_count,
            ]
            @test u1_target.assumptions.complete_add_on_zero_row_count
            complex_u1_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
            )
            @test complex_u1_target.rdm_logical_block_sizes == [1, 3, 3, 1]
            @test complex_u1_target.rdm_psd_block_sizes == [1, 3, 3, 1]
            _, complex_u1_target_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
            )
            complex_u1_target_metrics = translation_report_metrics(
                complex_u1_target_report,
            )
            @test complex_u1_target.contiguous_rdm_zero_row_count ==
                complex_u1_target_metrics.contiguous_rdm_zero_row_count
            @test complex_u1_target.contiguous_rdm_zero_row_count == 88
            @test complex_u1_target.add_on_zero_row_count ==
                complex_u1_target.contiguous_rdm_zero_row_count
            @test complex_u1_target.known_zero_constraint_feature_histogram == [
                (feature=:contiguous_rdm_zero, decomposition=:u1, reason=:magnetization_offblock) =>
                    complex_u1_target.contiguous_rdm_zero_row_count,
            ]
            @test complex_u1_target.logical_block_sizes ==
                complex_u1_target_report.logical_block_sizes
            @test complex_u1_target.psd_block_sizes ==
                complex_u1_target_report.psd_block_sizes
            @test complex_u1_target.psd_symmetric_entries ==
                complex_u1_target_metrics.psd_symmetric_entries
            @test complex_u1_target.psd_dense_entries ==
                complex_u1_target_metrics.psd_dense_entries
            @test complex_u1_target.psd_symmetric_entries ==
                complex_u1_target.psd_dense_entries
            @test complex_u1_target.assumptions.complete_add_on_zero_row_count

            combined_zero_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
                axis_rotation_equalities=true,
            )
            @test combined_zero_target.add_on_zero_row_count ==
                combined_zero_target.axis_rotation_equality_row_count +
                combined_zero_target.contiguous_rdm_zero_row_count
            @test combined_zero_target.known_zero_constraint_feature_histogram == [
                (feature=:axis_rotation_equality, decomposition=nothing, reason=nothing) =>
                    combined_zero_target.axis_rotation_equality_row_count,
                (feature=:contiguous_rdm_zero, decomposition=:u1, reason=:magnetization_offblock) =>
                    combined_zero_target.contiguous_rdm_zero_row_count,
            ]
            @test combined_zero_target.assumptions.complete_add_on_zero_row_count
            combined_all_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
                axis_rotation_equalities=true,
                linear_state_opt_width=2,
                moment_eq_h2=true,
            )
            _, combined_all_report = pauli_translation_invariant_moment_relaxation(
                h2_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
                axis_rotation_equalities=true,
                linear_state_opt_width=2,
            )
            combined_all_metrics = translation_report_metrics(combined_all_report)
            @test combined_all_target.axis_rotation_equality_row_count ==
                combined_all_metrics.axis_rotation_equality_row_count
            @test combined_all_target.contiguous_rdm_zero_row_count ==
                combined_all_metrics.contiguous_rdm_zero_row_count
            @test combined_all_target.linear_state_opt_row_count ==
                combined_all_metrics.linear_state_opt_row_count
            @test combined_all_target.moment_equality_row_count ==
                combined_all_metrics.moment_equality_row_count
            @test combined_all_target.add_on_zero_row_count ==
                combined_all_metrics.zero_constraint_count
            @test combined_all_target.known_zero_constraint_feature_histogram ==
                combined_all_metrics.zero_constraint_feature_histogram
            @test combined_all_target.assumptions.complete_add_on_zero_row_count
            complex_combined_all_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
                axis_rotation_equalities=true,
                linear_state_opt_width=2,
                moment_eq_h2=true,
            )
            _, complex_combined_all_report = pauli_translation_invariant_moment_relaxation(
                h2_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
                axis_rotation_equalities=true,
                linear_state_opt_width=2,
            )
            complex_combined_all_metrics =
                translation_report_metrics(complex_combined_all_report)
            @test complex_combined_all_target.axis_rotation_equality_row_count ==
                complex_combined_all_metrics.axis_rotation_equality_row_count
            @test complex_combined_all_target.contiguous_rdm_zero_row_count ==
                complex_combined_all_metrics.contiguous_rdm_zero_row_count
            @test complex_combined_all_target.linear_state_opt_row_count ==
                complex_combined_all_metrics.linear_state_opt_row_count
            @test complex_combined_all_target.moment_equality_row_count ==
                complex_combined_all_metrics.moment_equality_row_count
            @test complex_combined_all_target.add_on_zero_row_count ==
                complex_combined_all_metrics.zero_constraint_count
            @test complex_combined_all_target.known_zero_constraint_feature_histogram ==
                complex_combined_all_metrics.zero_constraint_feature_histogram
            @test complex_combined_all_target.logical_block_sizes ==
                complex_combined_all_report.logical_block_sizes
            @test complex_combined_all_target.psd_block_sizes ==
                complex_combined_all_report.psd_block_sizes
            @test complex_combined_all_target.psd_symmetric_entries ==
                complex_combined_all_metrics.psd_symmetric_entries
            @test complex_combined_all_target.psd_dense_entries ==
                complex_combined_all_metrics.psd_dense_entries
            @test complex_combined_all_target.assumptions.complete_add_on_zero_row_count

            su2_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                su2_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:su2,
            )
            @test su2_target.rdm_logical_block_sizes == [1, 1]
            @test su2_target.rdm_psd_block_sizes == [1, 1]
            @test su2_target.contiguous_rdm_zero_row_count == 0
            @test !su2_target.assumptions.complete_add_on_zero_row_count
            @test isempty(su2_target.known_zero_constraint_feature_histogram)
            complex_su2_target = pauli_translation_structural_targets(
                8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                su2_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:su2,
            )
            @test complex_su2_target.rdm_logical_block_sizes == [1, 1]
            @test complex_su2_target.rdm_psd_block_sizes == [1, 1]
            @test complex_su2_target.contiguous_rdm_zero_row_count == 0
            @test !complex_su2_target.assumptions.complete_add_on_zero_row_count
            @test isempty(complex_su2_target.known_zero_constraint_feature_histogram)
            _, complex_su2_target_report = pauli_translation_invariant_moment_relaxation(
                addon_pop,
                addon_ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                su2_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:su2,
            )
            complex_su2_target_metrics = translation_report_metrics(
                complex_su2_target_report,
            )
            @test complex_su2_target_metrics.su2_moment_symmetry
            @test complex_su2_target_metrics.su2_rdm_symmetry
            @test complex_su2_target_report.real_moment_matrix == false
            @test complex_su2_target.rdm_logical_block_sizes ==
                complex_su2_target_metrics.contiguous_rdm_logical_block_sizes
            @test complex_su2_target.rdm_psd_block_sizes ==
                complex_su2_target_metrics.contiguous_rdm_psd_block_sizes
            @test complex_su2_target_report.block_labels[end - 1:end] == [
                (feature=:contiguous_rdm, k=2, decomposition=:su2, spin2=0),
                (feature=:contiguous_rdm, k=2, decomposition=:su2, spin2=2),
            ]
            @test complex_su2_target_report.psd_block_sizes[end - 1:end] ==
                complex_su2_target.rdm_psd_block_sizes
            @test complex_su2_target_metrics.contiguous_rdm_zero_row_count > 0

            @test_throws ArgumentError pauli_translation_structural_targets(12, 6)
            @test_throws ArgumentError pauli_translation_structural_targets(100, 4; hamiltonian_width=0)
            @test_throws ArgumentError pauli_translation_structural_targets(8, 2; contiguous_rdm_k=5)
            h2_unclosed_error = try
                pauli_translation_structural_targets(8, 3; moment_eq_h2=true)
                nothing
            catch err
                err
            end
            @test h2_unclosed_error isa ArgumentError
            @test occursin(
                "`moment_eq_h2=true` structural target is not closed",
                sprint(showerror, h2_unclosed_error),
            )
            @test pauli_translation_structural_targets(
                8,
                2;
                contiguous_rdm_k=5,
                contiguous_rdm_support=:extend,
            ).rdm_logical_block_sizes == [32]
            @test_throws DomainError pauli_translation_structural_targets(
                8,
                2;
                contiguous_rdm_k=9,
                contiguous_rdm_support=:extend,
            )
            @test_throws ArgumentError pauli_translation_structural_targets(
                8,
                2;
                linear_state_opt_width=4,
            )
            @test_throws ArgumentError pauli_translation_structural_targets(
                8,
                2;
                psd_state_opt_width=2,
            )
            @test_throws ArgumentError pauli_translation_structural_targets(
                100,
                4;
                contiguous_rdm_k=7,
                contiguous_rdm_decomposition=:qmbcertify,
            )
            qmb_decomposition_extend_target = pauli_translation_structural_targets(
                100,
                4;
                contiguous_rdm_k=8,
                contiguous_rdm_decomposition=:qmbcertify,
                contiguous_rdm_support=:extend,
            )
            @test qmb_decomposition_extend_target.rdm_logical_block_sizes ==
                qmb_decomposition_target.rdm_logical_block_sizes
            @test qmb_decomposition_extend_target.rdm_psd_block_sizes ==
                qmb_decomposition_target.rdm_psd_block_sizes
            @test qmb_decomposition_extend_target.qmbcertify_rdm_reference_block_sizes ==
                qmb_decomposition_target.qmbcertify_rdm_reference_block_sizes
            @test qmb_decomposition_extend_target.contiguous_rdm_support == :extend
            finite_axis_qmb_target = pauli_translation_structural_targets(
                10,
                4;
                reflection_symmetry=true,
                axis_rotation_symmetry=true,
                axis_rotation_quotient=true,
                contiguous_rdm_k=8,
                contiguous_rdm_decomposition=:qmbcertify,
                contiguous_rdm_support=:extend,
                psd_state_opt_width=3,
            )
            @test finite_axis_qmb_target.axis_rotation_symmetry
            @test finite_axis_qmb_target.axis_rotation_quotient
            @test finite_axis_qmb_target.psd_state_opt_support_policy == :extend
            @test finite_axis_qmb_target.psd_state_opt_psd_block_sizes == fill(18, 6)
            @test finite_axis_qmb_target.rdm_psd_block_sizes == [72, 64, 56]
            @test finite_axis_qmb_target.n_blocks == 35
            @test finite_axis_qmb_target.psd_max_block == 120
        end

        @test isdefined(NCTSSoS, :translation_symmetry_profile)
        if isdefined(NCTSSoS, :translation_symmetry_profile)
            base_profile = translation_symmetry_profile(report_large)
            @test base_profile.translation_symmetry
            @test base_profile.sign_symmetry
            @test !base_profile.reflection_symmetry
            @test !base_profile.axis_rotation_symmetry
            @test base_profile.axis_orbit_closed
            @test base_profile.axis_orbit_basis_size == 23
            @test base_profile.axis_reduction_ratio == 121 / 23
            @test base_profile.missing_discrete_symmetries == [:reflection, :axis_rotation]
            @test !base_profile.full_discrete_profile

        end

        mp_reflection, report_reflection = pauli_translation_invariant_moment_relaxation(
            pop_large,
            ops_large,
            4;
            reflection_symmetry=true,
        )
        @test translation_block_histogram(report_reflection; kind=:logical) ==
            [9 => 2, 13 => 6, 17 => 6, 21 => 1, 22 => 1, 30 => 40]
        @test translation_block_histogram(report_reflection; kind=:psd) ==
            [18 => 2, 26 => 6, 30 => 40, 34 => 6, 42 => 1, 44 => 1]
        metrics_reflection = translation_report_metrics(report_reflection)
        @test metrics_reflection.n_blocks == 56
        @test metrics_reflection.logical_max_block == 30
        @test metrics_reflection.psd_max_block <= 44
        @test metrics_reflection.psd_dense_entries == 51_340
        @test metrics_reflection.psd_symmetric_entries == 26_511

        registry_reflected_sign, ops_reflected_sign = create_pauli_variables(1:6)
        pop_reflected_sign = polyopt(
            heisenberg_chain_hamiltonian(ops_reflected_sign),
            registry_reflected_sign,
        )
        mp_reflected_sign, _ = pauli_translation_invariant_moment_relaxation(
            pop_reflected_sign,
            ops_reflected_sign,
            2;
            sign_symmetry=true,
            reflection_symmetry=true,
        )
        model_reflected_sign, _ = build_jump_model(
            mp_reflected_sign.linear;
            formulation=:moment_variables,
            representation=:real,
        )
        @test model_reflected_sign !== nothing

        direct_reflected_sign, _ = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop_reflected_sign,
            ops_reflected_sign,
            2;
            sign_symmetry=true,
            reflection_symmetry=true,
        )
        direct_model_reflected_sign, _ = build_jump_model(
            direct_reflected_sign;
            formulation=:moment_variables,
            representation=:real,
        )
        @test direct_model_reflected_sign !== nothing

        if isdefined(NCTSSoS, :pauli_translation_structural_targets)
            target_reflection = pauli_translation_structural_targets(100, 4; reflection_symmetry=true)
            @test target_reflection.n_blocks == 408
            @test target_reflection.logical_max_block == 30
            @test target_reflection.psd_max_block == 44
            @test target_reflection.logical_block_histogram ==
                [9 => 2, 13 => 6, 17 => 6, 21 => 1, 22 => 1, 30 => 392]
            @test target_reflection.psd_block_histogram ==
                [18 => 2, 26 => 6, 30 => 392, 34 => 6, 42 => 1, 44 => 1]
            @test target_reflection.psd_symmetric_entries == 190_191
            @test target_reflection.block_coefficient_domains ==
                fill(:cyclotomic_float64, target_reflection.n_blocks)
            @test target_reflection.block_exact_coefficient_domains ==
                fill(:cyclotomic_sqrt_rational, target_reflection.n_blocks)
            @test target_reflection.block_coefficient_domain_histogram ==
                [:cyclotomic_float64 => target_reflection.n_blocks]
            @test target_reflection.block_exact_coefficient_domain_histogram ==
                [:cyclotomic_sqrt_rational => target_reflection.n_blocks]
            @test target_reflection.estimated_model_size_gate_status ==
                :blocked_missing_scalar_equality_estimate

            target_reflection12 = pauli_translation_structural_targets(
                n_large,
                4;
                reflection_symmetry=true,
            )
            @test target_reflection12.n_blocks == metrics_reflection.n_blocks
            @test target_reflection12.logical_max_block == metrics_reflection.logical_max_block
            @test target_reflection12.psd_max_block == metrics_reflection.psd_max_block
            @test target_reflection12.logical_block_histogram ==
                metrics_reflection.logical_block_histogram
            @test target_reflection12.psd_block_histogram ==
                metrics_reflection.psd_block_histogram
            @test target_reflection12.psd_symmetric_entries ==
                metrics_reflection.psd_symmetric_entries

            candidate_registry, candidate_ops = create_pauli_variables(1:10)
            candidate_pop = polyopt(
                heisenberg_chain_hamiltonian(candidate_ops),
                candidate_registry,
            )
            _, report_reflection10 = pauli_translation_invariant_moment_relaxation(
                candidate_pop,
                candidate_ops,
                4;
                reflection_symmetry=true,
            )
            target_reflection10 = pauli_translation_structural_targets(
                10,
                4;
                reflection_symmetry=true,
            )
            metrics_reflection10 = translation_report_metrics(report_reflection10)
            @test target_reflection10.n_blocks == metrics_reflection10.n_blocks
            @test target_reflection10.logical_max_block ==
                metrics_reflection10.logical_max_block
            @test target_reflection10.psd_max_block == metrics_reflection10.psd_max_block
            @test target_reflection10.logical_block_histogram ==
                metrics_reflection10.logical_block_histogram
            @test target_reflection10.psd_block_histogram ==
                metrics_reflection10.psd_block_histogram
            @test target_reflection10.psd_symmetric_entries ==
                metrics_reflection10.psd_symmetric_entries

            @test_throws ArgumentError pauli_translation_structural_targets(
                10,
                4;
                reflection_symmetry=true,
                real_moment_matrix=false,
            )
            target_reflection_fixed_complex = pauli_translation_structural_targets(
                10,
                4;
                reflection_symmetry=true,
                real_moment_matrix=false,
                momenta=[0],
            )
            @test target_reflection_fixed_complex.momentum_sectors == [0]
            @test target_reflection_fixed_complex.real_moment_matrix == false
            @test target_reflection_fixed_complex.logical_block_sizes ==
                target_reflection_fixed_complex.psd_block_sizes
            @test target_reflection_fixed_complex.logical_max_block ==
                target_reflection_fixed_complex.psd_max_block
            @test target_reflection_fixed_complex.logical_block_histogram ==
                target_reflection_fixed_complex.psd_block_histogram
            @test target_reflection_fixed_complex.psd_symmetric_entries ==
                target_reflection_fixed_complex.psd_dense_entries
        end
        @test any(
            label -> haskey(label, :reflection) && label.momentum == 0,
            report_reflection.block_labels,
        )
        @test any(
            label -> haskey(label, :reflection) && label.momentum == fld(n_large, 2),
            report_reflection.block_labels,
        )
        @test [
            length(block.meta.origin.logical_row_labels)
            for block in mp_reflection.linear.psd_blocks_lin
        ] == report_reflection.logical_block_sizes
        reflection_origins = [
            block.meta.origin
            for block in mp_reflection.linear.psd_blocks_lin
            if haskey(block.meta.origin.label, :reflection)
        ]
        @test !isempty(reflection_origins)
        @test all(
            origin -> origin.transform.family == :translation_dft_reflection,
            reflection_origins,
        )
        @test all(
            origin -> origin.transform.reflection == origin.label.reflection,
            reflection_origins,
        )
        @test any(origin -> origin.transform.antiunitary, reflection_origins)
        @test any(origin -> !origin.transform.antiunitary, reflection_origins)
        @test all(
            origin -> origin.transform.exact_coefficient_domain == :cyclotomic_sqrt_rational,
            reflection_origins,
        )
        if isdefined(NCTSSoS, :translation_symmetry_profile)
            reflection_profile = translation_symmetry_profile(report_reflection)
            @test reflection_profile.translation_symmetry
            @test reflection_profile.sign_symmetry
            @test reflection_profile.reflection_symmetry
            @test !reflection_profile.axis_rotation_symmetry
            @test reflection_profile.axis_orbit_closed
            @test reflection_profile.axis_orbit_basis_size == 23
            @test reflection_profile.axis_reduction_ratio == 121 / 23
            @test reflection_profile.missing_discrete_symmetries == [:axis_rotation]
            @test !reflection_profile.full_discrete_profile
        end
    end

    @testset "axis diagnostics do not constrain custom bases" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        σx, _, _ = ops
        basis = [one(σx[1]); σx]
        pop = polyopt(sum(σx), registry)

        @test_throws ArgumentError pauli_axis_orbit_diagnostics(ops, 1; basis)
        _, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            basis,
            sign_symmetry=false,
        )
        @test !report.axis_orbit_closed
        @test report.axis_orbit_basis_size == 0
        @test isempty(report.axis_orbit_size_histogram)
        @test report.axis_reduction_ratio == 0.0

        profile = translation_symmetry_profile(report)
        @test !profile.axis_orbit_closed
        @test profile.axis_orbit_basis_size == 0
        @test profile.axis_reduction_ratio == 0.0
    end

    @testset "guardrails reject invalid reductions" begin
        registry, ops = create_pauli_variables(1:4)
        heisenberg_pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(heisenberg_pop, ops, 0)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(heisenberg_pop, ops, 1; momenta=[1, 2], real_moment_matrix=false)
        linear_closure_err = try
            pauli_translation_invariant_moment_relaxation(
                heisenberg_pop,
                ops,
                1;
                sign_symmetry=false,
                linear_state_opt_width=2,
            )
            nothing
        catch err
            err
        end
        @test linear_closure_err isa ArgumentError
        @test occursin("linear_state_opt_width=2", sprint(showerror, linear_closure_err))

        direct_psd_closure_err = try
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                psd_state_opt_width=1,
            )
            nothing
        catch err
            err
        end
        @test direct_psd_closure_err isa ArgumentError
        @test occursin("psd_state_opt_width=1", sprint(showerror, direct_psd_closure_err))

        custom_basis = pauli_contiguous_chain_basis(ops, 1)
        custom_closure_err = try
            pauli_translation_invariant_moment_relaxation(
                heisenberg_pop,
                ops,
                1;
                basis=custom_basis,
                sign_symmetry=false,
                contiguous_rdm_k=3,
            )
            nothing
        catch err
            err
        end
        @test custom_closure_err isa ArgumentError
        @test occursin(
            "closure_check(:contiguous_rdm",
            sprint(showerror, custom_closure_err),
        )

        complex_reflection_err = try
            pauli_translation_invariant_moment_relaxation(
                heisenberg_pop,
                ops,
                1;
                reflection_symmetry=true,
                real_moment_matrix=false,
            )
            nothing
        catch err
            err
        end
        @test complex_reflection_err isa ArgumentError
        @test occursin("reflection-fixed momentum sectors", sprint(showerror, complex_reflection_err))
        _, fixed_complex_reflection = pauli_translation_invariant_moment_relaxation(
            heisenberg_pop,
            ops,
            1;
            reflection_symmetry=true,
            real_moment_matrix=false,
            momenta=[0],
        )
        @test fixed_complex_reflection.real_moment_matrix == false
        @test fixed_complex_reflection.momentum_sectors == [0]
        fixed_complex_metrics = translation_report_metrics(fixed_complex_reflection)
        @test fixed_complex_metrics.psd_symmetric_entries ==
            fixed_complex_metrics.psd_dense_entries
        _, sign_axis_symbolic_report = pauli_translation_invariant_moment_relaxation(
            heisenberg_pop,
            ops,
            2;
            axis_rotation_symmetry=true,
        )
        @test sign_axis_symbolic_report.sign_symmetry
        @test sign_axis_symbolic_report.psd_block_sizes ==
            [2, 4, 6, 12, 1, 2, 6, 12, 1, 4, 6, 12]
        @test translation_symmetry_profile(sign_axis_symbolic_report).axis_rotation_symmetry
        sign_axis_symbolic_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                axis_rotation_symmetry=true,
            )
        end
        @test sign_axis_symbolic_result.moment_problem isa NCTSSoS.MomentProblem
        @test termination_status(sign_axis_symbolic_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test sign_axis_symbolic_result.report.sign_symmetry
        @test sign_axis_symbolic_result.report.psd_block_sizes ==
            sign_axis_symbolic_report.psd_block_sizes
        @test sign_axis_symbolic_result.objective ≈ -2 atol = 1e-6
        axis_reflection_mp, axis_reflection_report =
            pauli_translation_invariant_moment_relaxation(
                heisenberg_pop,
                ops,
                2;
                sign_symmetry=false,
                reflection_symmetry=true,
                axis_rotation_symmetry=true,
            )
        axis_reflection_profile = translation_symmetry_profile(axis_reflection_report)
        @test axis_reflection_profile.reflection_symmetry
        @test axis_reflection_profile.axis_rotation_symmetry
        @test axis_reflection_profile.axis_rotation_equalities
        @test axis_reflection_report.zero_constraint_count == 80
        @test axis_reflection_report.psd_block_sizes ==
            [2, 4, 6, 6, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 4, 6, 6, 6]
        @test NCTSSoS.assert_moment_linear_data_invariants(
            axis_reflection_mp.linear,
            axis_reflection_mp.constraints,
        ) === nothing
        axis_reflection_linear, axis_reflection_direct_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                heisenberg_pop,
                ops,
                2;
                sign_symmetry=false,
                reflection_symmetry=true,
                axis_rotation_symmetry=true,
            )
        @test axis_reflection_direct_report.psd_block_sizes ==
            axis_reflection_report.psd_block_sizes
        @test axis_reflection_direct_report.zero_constraint_count ==
            axis_reflection_report.zero_constraint_count
        @test NCTSSoS.assert_moment_linear_data_invariants(
            axis_reflection_linear,
        ) === nothing
        axis_reflection_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                reflection_symmetry=true,
                axis_rotation_symmetry=true,
            )
        end
        direct_axis_reflection_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                reflection_symmetry=true,
                axis_rotation_symmetry=true,
            )
        end
        @test termination_status(axis_reflection_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(direct_axis_reflection_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test axis_reflection_result.objective ≈ -2 atol = 1e-6
        @test direct_axis_reflection_result.objective ≈
            axis_reflection_result.objective atol = 1e-6
        _, base_symbolic_report = pauli_translation_invariant_moment_relaxation(
            heisenberg_pop,
            ops,
            2;
            sign_symmetry=false,
        )
        axis_mp, axis_symbolic_report = pauli_translation_invariant_moment_relaxation(
            heisenberg_pop,
            ops,
            2;
            sign_symmetry=false,
            axis_rotation_symmetry=true,
        )
        @test axis_mp isa NCTSSoS.MomentProblem
        @test base_symbolic_report.psd_block_sizes == [26, 24, 24]
        @test axis_symbolic_report.psd_block_sizes == [2, 4, 6, 12, 1, 2, 6, 12, 1, 4, 6, 12]
        @test axis_symbolic_report.zero_constraint_count == 80
        @test maximum(axis_symbolic_report.psd_block_sizes) <
            maximum(base_symbolic_report.psd_block_sizes)
        @test translation_symmetry_profile(axis_symbolic_report).axis_rotation_symmetry
        @test translation_symmetry_profile(axis_symbolic_report).axis_rotation_equalities
        symbolic_axis_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                axis_rotation_symmetry=true,
            )
        end
        @test symbolic_axis_result.moment_problem isa NCTSSoS.MomentProblem
        @test termination_status(symbolic_axis_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test symbolic_axis_result.objective ≈ -2 atol = 1e-6
        @test symbolic_axis_result.report.psd_block_sizes ==
            axis_symbolic_report.psd_block_sizes
        @test translation_symmetry_profile(symbolic_axis_result).axis_rotation_symmetry
        @test translation_solve_support(symbolic_axis_result).supported
        _, su2_axis_symmetry_report = pauli_translation_invariant_moment_relaxation(
            heisenberg_pop,
            ops,
            1;
            sign_symmetry=false,
            su2_symmetry=true,
            axis_rotation_symmetry=true,
        )
        @test translation_symmetry_profile(su2_axis_symmetry_report).su2_moment_symmetry
        @test translation_solve_support(su2_axis_symmetry_report).supported
        _, base_direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            heisenberg_pop,
            ops,
            2;
            sign_symmetry=false,
        )
        axis_linear, axis_direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            heisenberg_pop,
            ops,
            2;
            sign_symmetry=false,
            axis_rotation_symmetry=true,
        )
        @test axis_linear isa NCTSSoS.MomentLinearData
        axis_labels = [
            label for label in axis_direct_report.block_labels if
            label isa NamedTuple && haskey(label, :axis_irrep)
        ]
        @test !isempty(axis_labels)
        @test all(label.axis_group_order == 24 for label in axis_labels)
        @test maximum(axis_direct_report.psd_block_sizes) <
            maximum(base_direct_report.psd_block_sizes)
        @test base_direct_report.psd_block_sizes == [26, 24, 24]
        @test axis_direct_report.psd_block_sizes == [2, 4, 6, 12, 1, 2, 6, 12, 1, 4, 6, 12]
        @test axis_direct_report.zero_constraint_count == 80
        @test axis_direct_report.linear_moment_count == 70
        axis_direct_profile = translation_symmetry_profile(axis_direct_report)
        @test axis_direct_profile.axis_rotation_symmetry
        @test axis_direct_profile.axis_rotation_equalities
        axis_direct_metrics = translation_report_metrics(axis_direct_report)
        @test any(axis_direct_metrics.logical_block_feature_histogram) do pair
            key = first(pair)
            key isa NamedTuple &&
                haskey(key, :decomposition) &&
                key.decomposition == :translation_axis_irrep
        end
        axis_linear_mp, axis_linear_report = pauli_translation_invariant_moment_relaxation(
            heisenberg_pop,
            ops,
            2;
            axis_rotation_symmetry=true,
            linear_state_opt_width=2,
        )
        @test axis_linear_report.sign_symmetry
        @test axis_linear_report.psd_block_sizes == axis_symbolic_report.psd_block_sizes
        axis_linear_metrics = translation_report_metrics(axis_linear_report)
        @test axis_linear_metrics.axis_rotation_symmetry
        @test axis_linear_metrics.axis_rotation_equalities
        @test axis_linear_metrics.linear_state_opt
        @test axis_linear_metrics.axis_rotation_equality_row_count == 80
        @test axis_linear_metrics.linear_state_opt_row_count > 0
        @test axis_linear_report.zero_constraint_count ==
            axis_linear_metrics.axis_rotation_equality_row_count +
            axis_linear_metrics.linear_state_opt_row_count
        @test any(
            zc -> zc.origin.label.feature == :axis_rotation_equality,
            axis_linear_mp.linear.zero_constraints,
        )
        @test any(
            zc -> zc.origin.label.feature == :linear_state_opt,
            axis_linear_mp.linear.zero_constraints,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(
            axis_linear_mp.linear,
            axis_linear_mp.constraints,
        ) === nothing
        axis_direct_linear, axis_direct_linear_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                heisenberg_pop,
                ops,
                2;
                axis_rotation_symmetry=true,
                linear_state_opt_width=2,
            )
        @test axis_direct_linear_report.sign_symmetry
        @test axis_direct_linear_report.psd_block_sizes == axis_direct_report.psd_block_sizes
        axis_direct_linear_metrics = translation_report_metrics(axis_direct_linear_report)
        @test axis_direct_linear_metrics.axis_rotation_symmetry
        @test axis_direct_linear_metrics.axis_rotation_equalities
        @test axis_direct_linear_metrics.linear_state_opt
        @test axis_direct_linear_metrics.axis_rotation_equality_row_count == 80
        @test axis_direct_linear_metrics.linear_state_opt_row_count > 0
        @test axis_direct_linear_report.zero_constraint_count ==
            axis_direct_linear_metrics.axis_rotation_equality_row_count +
            axis_direct_linear_metrics.linear_state_opt_row_count
        @test any(
            zc -> zc.origin.label.feature == :axis_rotation_equality,
            axis_direct_linear.zero_constraints,
        )
        @test any(
            zc -> zc.origin.label.feature == :linear_state_opt,
            axis_direct_linear.zero_constraints,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(axis_direct_linear) === nothing
        axis_linear_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                axis_rotation_symmetry=true,
                linear_state_opt_width=2,
            )
        end
        @test termination_status(axis_linear_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test axis_linear_result.objective ≈ -2 atol = 1e-6
        @test axis_linear_result.report.zero_constraint_count ==
            axis_linear_report.zero_constraint_count
        direct_axis_linear_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                axis_rotation_symmetry=true,
                linear_state_opt_width=2,
            )
        end
        @test termination_status(direct_axis_linear_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test direct_axis_linear_result.objective ≈ axis_linear_result.objective atol = 1e-6
        @test direct_axis_linear_result.report.zero_constraint_count ==
            axis_direct_linear_report.zero_constraint_count
        axis_psd_state_mp, axis_psd_state_report =
            pauli_translation_invariant_moment_relaxation(
                heisenberg_pop,
                ops,
                2;
                sign_symmetry=false,
                axis_rotation_symmetry=true,
                psd_state_opt_width=1,
            )
        axis_psd_state_metrics = translation_report_metrics(axis_psd_state_report)
        @test axis_psd_state_metrics.axis_rotation_symmetry
        @test axis_psd_state_metrics.axis_rotation_equalities
        @test axis_psd_state_metrics.psd_state_opt
        @test axis_psd_state_metrics.axis_rotation_equality_row_count == 80
        @test axis_psd_state_metrics.psd_state_opt_block_count == 3
        @test axis_psd_state_report.zero_constraint_count == 80
        @test NCTSSoS.assert_moment_linear_data_invariants(
            axis_psd_state_mp.linear,
            axis_psd_state_mp.constraints,
        ) === nothing
        axis_direct_psd_state, axis_direct_psd_state_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                heisenberg_pop,
                ops,
                2;
                sign_symmetry=false,
                axis_rotation_symmetry=true,
                psd_state_opt_width=1,
            )
        axis_direct_psd_state_metrics =
            translation_report_metrics(axis_direct_psd_state_report)
        @test axis_direct_psd_state_metrics.axis_rotation_symmetry
        @test axis_direct_psd_state_metrics.axis_rotation_equalities
        @test axis_direct_psd_state_metrics.psd_state_opt
        @test axis_direct_psd_state_metrics.axis_rotation_equality_row_count == 80
        @test axis_direct_psd_state_metrics.psd_state_opt_block_count == 3
        @test axis_direct_psd_state_report.zero_constraint_count == 80
        @test NCTSSoS.assert_moment_linear_data_invariants(
            axis_direct_psd_state,
        ) === nothing
        axis_rdm_mp, axis_rdm_report = pauli_translation_invariant_moment_relaxation(
            heisenberg_pop,
            ops,
            2;
            axis_rotation_symmetry=true,
            contiguous_rdm_k=2,
        )
        axis_rdm_metrics = translation_report_metrics(axis_rdm_report)
        @test axis_rdm_metrics.axis_rotation_symmetry
        @test axis_rdm_metrics.axis_rotation_equalities
        @test axis_rdm_metrics.contiguous_rdm
        @test axis_rdm_metrics.axis_rotation_equality_row_count == 80
        @test axis_rdm_metrics.contiguous_rdm_block_count == 1
        @test axis_rdm_report.zero_constraint_count == 80
        @test NCTSSoS.assert_moment_linear_data_invariants(
            axis_rdm_mp.linear,
            axis_rdm_mp.constraints,
        ) === nothing
        axis_direct_rdm, axis_direct_rdm_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                heisenberg_pop,
                ops,
                2;
                axis_rotation_symmetry=true,
                contiguous_rdm_k=2,
            )
        axis_direct_rdm_metrics = translation_report_metrics(axis_direct_rdm_report)
        @test axis_direct_rdm_metrics.axis_rotation_symmetry
        @test axis_direct_rdm_metrics.axis_rotation_equalities
        @test axis_direct_rdm_metrics.contiguous_rdm
        @test axis_direct_rdm_metrics.axis_rotation_equality_row_count == 80
        @test axis_direct_rdm_metrics.contiguous_rdm_block_count == 1
        @test axis_direct_rdm_report.zero_constraint_count == 80
        @test NCTSSoS.assert_moment_linear_data_invariants(axis_direct_rdm) === nothing
        axis_psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                axis_rotation_symmetry=true,
                psd_state_opt_width=1,
            )
        end
        direct_axis_psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                axis_rotation_symmetry=true,
                psd_state_opt_width=1,
            )
        end
        @test termination_status(axis_psd_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(direct_axis_psd_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test direct_axis_psd_state_result.objective ≈
            axis_psd_state_result.objective atol = 1e-6
        axis_rdm_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                axis_rotation_symmetry=true,
                contiguous_rdm_k=2,
            )
        end
        direct_axis_rdm_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                axis_rotation_symmetry=true,
                contiguous_rdm_k=2,
            )
        end
        @test termination_status(axis_rdm_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(direct_axis_rdm_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test direct_axis_rdm_result.objective ≈ axis_rdm_result.objective atol = 1e-6
        direct_axis_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                axis_rotation_symmetry=true,
            )
        end
        @test termination_status(direct_axis_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test direct_axis_result.objective ≈ -2 atol = 1e-6
        @test translation_symmetry_profile(direct_axis_result).axis_rotation_symmetry
        @test translation_solve_support(direct_axis_result).supported
        direct_sign_axis_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                axis_rotation_symmetry=true,
            )
        end
        @test termination_status(direct_sign_axis_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test direct_sign_axis_result.report.sign_symmetry
        @test direct_sign_axis_result.report.psd_block_sizes ==
            direct_axis_result.report.psd_block_sizes
        @test direct_sign_axis_result.objective ≈ direct_axis_result.objective atol = 1e-6
        @test direct_sign_axis_result.objective ≈ -2 atol = 1e-6
        direct_su2_axis_symmetry_result = quiet() do
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                su2_symmetry=true,
                axis_rotation_symmetry=true,
            )
        end
        @test termination_status(direct_su2_axis_symmetry_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test direct_su2_axis_symmetry_result.moment_problem isa NCTSSoS.MomentLinearData
        @test translation_symmetry_profile(
            direct_su2_axis_symmetry_result,
        ).su2_moment_symmetry
        @test translation_solve_support(direct_su2_axis_symmetry_result).supported

        field_pop = polyopt(sum(ops[1]), registry)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(field_pop, ops, 1)
        @test pauli_translation_invariant_moment_relaxation(field_pop, ops, 1; sign_symmetry=false)[2].psd_block_sizes == [8, 6, 6]
        qmbcertify_rdm_err = try
            pauli_translation_invariant_moment_relaxation(
                heisenberg_pop,
                ops,
                1;
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:qmbcertify,
            )
            nothing
        catch err
            err
        end
        @test qmbcertify_rdm_err isa ArgumentError
        @test occursin("structural-target-only", sprint(showerror, qmbcertify_rdm_err))
        @test occursin("shared RDM PSD variables", sprint(showerror, qmbcertify_rdm_err))

        qmbcertify_bad_support_err = try
            pauli_translation_invariant_moment_relaxation(
                heisenberg_pop,
                ops,
                1;
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:qmbcertify,
                contiguous_rdm_support=:bogus,
            )
            nothing
        catch err
            err
        end
        @test qmbcertify_bad_support_err isa ArgumentError
        @test occursin("contiguous_rdm_support", sprint(showerror, qmbcertify_bad_support_err))
        @test occursin("bogus", sprint(showerror, qmbcertify_bad_support_err))

        direct_qmbcertify_rdm_err = try
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:qmbcertify,
            )
            nothing
        catch err
            err
        end
        @test direct_qmbcertify_rdm_err isa ArgumentError
        @test occursin("structural-target-only", sprint(showerror, direct_qmbcertify_rdm_err))
        @test occursin("shared RDM PSD variables", sprint(showerror, direct_qmbcertify_rdm_err))

        direct_qmbcertify_bad_support_err = try
            pauli_translation_invariant_nctssos(
                heisenberg_pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:qmbcertify,
                contiguous_rdm_support=:bogus,
            )
            nothing
        catch err
            err
        end
        @test direct_qmbcertify_bad_support_err isa ArgumentError
        @test occursin(
            "contiguous_rdm_support",
            sprint(showerror, direct_qmbcertify_bad_support_err),
        )
        @test occursin("bogus", sprint(showerror, direct_qmbcertify_bad_support_err))

        σx, σy, σz = ops
        @test_throws ArgumentError pauli_contiguous_chain_basis((σy, σx, σz), 1)

        registry8, ops8 = create_pauli_variables(1:8)
        mismatched_pop = polyopt(ops8[1][6] * ops8[1][7], registry8)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(mismatched_pop, ops, 2; sign_symmetry=false)

        registry5, ops5 = create_pauli_variables(1:5)
        σx5, σy5, _ = ops5
        oriented_pop = polyopt(sum(σx5[i] * σy5[mod1(i + 1, 5)] for i in 1:5), registry5)
        @test pauli_translation_invariant_moment_relaxation(
            oriented_pop,
            ops5,
            2;
            sign_symmetry=false,
            reflection_symmetry=false,
        )[2].basis_size > 0
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            oriented_pop,
            ops5,
            2;
            sign_symmetry=false,
            reflection_symmetry=true,
        )

        registry100, ops100 = create_pauli_variables(1:100)
        type_mismatched_pop = polyopt(ops100[1][1] * ops100[1][2], registry100)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(type_mismatched_pop, ops, 2; sign_symmetry=false)
    end

    @testset "translation path accepts invariant scalar constraints" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        _, _, σz = ops
        objective = heisenberg_chain_hamiltonian(ops)
        eq = sum(σz)
        ineq = one(objective) + objective
        pop = polyopt(objective, registry; eq_constraints=[eq], ineq_constraints=[ineq])

        mp, report = pauli_translation_invariant_moment_relaxation(pop, ops, 1; sign_symmetry=false)

        @test report.psd_block_sizes == [8, 6, 6, 1]
        @test report.logical_block_sizes == [4, 3, 3, 1]
        @test report.block_labels[end] == (feature=:scalar_inequality, index=1)
        @test count(con -> con[1] == :Zero, mp.constraints) == 1
        @test count(con -> con[1] == :PSD, mp.constraints) == length(report.psd_block_sizes)
        @test only(mp.linear.zero_constraints).origin.label == (feature=:scalar_equality, index=1)
        eq_entry = only(mp.constraints[end - 1][2])
        ineq_entry = only(mp.constraints[end][2])
        @test eq_entry == convert(typeof(eq_entry), NCTSSoS._translation_orbit_reduce_polynomial(eq, n))
        @test ineq_entry == convert(typeof(ineq_entry), NCTSSoS._translation_orbit_reduce_polynomial(ineq, n))

        symbolic_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
            )
        end
        direct_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                direct_linear=true,
            )
        end
        @test direct_result.moment_problem isa NCTSSoS.MomentLinearData
        @test direct_result.objective ≈ symbolic_result.objective atol = 1e-6
        @test direct_result.report.block_labels == symbolic_result.report.block_labels

        complex_symbolic_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                real_moment_matrix=false,
            )
        end
        complex_direct_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                real_moment_matrix=false,
                direct_linear=true,
            )
        end
        @test complex_direct_result.moment_problem isa NCTSSoS.MomentLinearData
        @test complex_direct_result.objective ≈ complex_symbolic_result.objective atol = 1e-6
        @test complex_direct_result.report.block_labels == complex_symbolic_result.report.block_labels
    end

    @testset "translation SOS zero duals keep provenance labels" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        _, _, σz = ops
        objective = heisenberg_chain_hamiltonian(ops)
        pop = polyopt(objective, registry; eq_constraints=[sum(σz)])

        mp, _ = pauli_translation_invariant_moment_relaxation(pop, ops, 1; sign_symmetry=false)
        sos = NCTSSoS.sos_dualize(mp)

        @test isdefined(NCTSSoS, :sos_zero_duals)
        @test isdefined(NCTSSoS, :sos_zero_dual_values)
        @test isdefined(NCTSSoS, :translation_zero_origin_histogram)
        if isdefined(NCTSSoS, :sos_zero_duals)
            zero_duals = NCTSSoS.sos_zero_duals(mp, sos)
            @test length(zero_duals) == length(mp.linear.zero_constraints) == 1
            @test only(zero_duals).origin.label == (feature=:scalar_equality, index=1)
            @test only(zero_duals).label == (feature=:scalar_equality, index=1)
            @test only(zero_duals).feature == :scalar_equality
            @test only(zero_duals).decomposition === nothing
            @test only(zero_duals).reason === nothing
            @test only(zero_duals).coefficient_domain === nothing
            @test only(zero_duals).exact_coefficient_domain === nothing
            @test only(zero_duals).term_count == length(only(mp.linear.zero_constraints).form)
            @test only(zero_duals).form == only(mp.linear.zero_constraints).form
        end
        if isdefined(NCTSSoS, :translation_zero_origin_histogram)
            @test translation_zero_origin_histogram(mp) == [
                (feature=:scalar_equality, decomposition=nothing, reason=nothing) => 1,
            ]
        end
    end

    @testset "translation path accepts invariant moment equality constraints" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        _, _, σz = ops
        objective = heisenberg_chain_hamiltonian(ops)
        meq = sum(one(objective) - σz[i] * σz[mod1(i + 1, n)] for i in 1:n)
        pop = polyopt(objective, registry; moment_eq_constraints=[meq])

        mp, report = pauli_translation_invariant_moment_relaxation(pop, ops, 1; sign_symmetry=false)

        @test report.psd_block_sizes == [8, 6, 6]
        @test count(con -> con[1] == :Zero, mp.constraints) > 0
        @test all(size(mat) == (1, 1) for (cone, mat) in mp.constraints if cone == :Zero)
        @test all(
            zc -> zc.origin.label.feature == :moment_equality && zc.origin.label.index == 1,
            mp.linear.zero_constraints,
        )
        meq_metrics = translation_report_metrics(report)
        @test meq_metrics.moment_equality_row_count == report.zero_constraint_count
        @test meq_metrics.axis_rotation_equality_row_count == 0
        @test meq_metrics.linear_state_opt_row_count == 0
        @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
    end

    @testset "translation path can add axis-rotation moment equalities" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        mp, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            axis_rotation_equalities=true,
        )

        @test report.psd_block_sizes == [8, 6, 6]
        @test report.zero_constraint_count == length(mp.linear.zero_constraints)
        @test report.zero_constraint_count > 0
        @test all(
            zc -> zc.origin.label.feature == :axis_rotation_equality,
            mp.linear.zero_constraints,
        )
        @test any(zc -> zc.origin.label.coefficient == -1, mp.linear.zero_constraints)
        @test translation_zero_origin_histogram(mp) == [
            (feature=:axis_rotation_equality, decomposition=nothing, reason=nothing) =>
                report.zero_constraint_count,
        ]
        axis_profile = translation_symmetry_profile(report)
        @test !axis_profile.axis_rotation_symmetry
        @test axis_profile.axis_rotation_equalities
        axis_metrics = translation_report_metrics(report)
        @test axis_metrics.axis_rotation_equality_row_count == report.zero_constraint_count
        @test axis_metrics.moment_equality_row_count == 0
        @test axis_metrics.linear_state_opt_row_count == 0
        @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing

        direct_linear, direct_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            axis_rotation_equalities=true,
        )
        @test direct_report.psd_block_sizes == report.psd_block_sizes
        @test direct_report.zero_constraint_count == report.zero_constraint_count
        @test all(
            zc -> zc.origin.label.feature == :axis_rotation_equality,
            direct_linear.zero_constraints,
        )
        @test translation_zero_origin_histogram(direct_linear) == [
            (feature=:axis_rotation_equality, decomposition=nothing, reason=nothing) =>
                direct_report.zero_constraint_count,
        ]
        @test isdefined(NCTSSoS, :translation_linear_provenance)
        if isdefined(NCTSSoS, :translation_linear_provenance)
            direct_axis_provenance = translation_linear_provenance(direct_linear)
            @test length(direct_axis_provenance.zero_constraints) ==
                direct_report.zero_constraint_count
            @test all(
                row -> row.kind == :zero &&
                    row.feature == :axis_rotation_equality &&
                    row.term_count > 0,
                direct_axis_provenance.zero_constraints,
            )
        end
        direct_axis_metrics = translation_report_metrics(direct_report)
        @test direct_axis_metrics.axis_rotation_equality_row_count ==
            direct_report.zero_constraint_count
        @test direct_axis_metrics.moment_equality_row_count == 0
        @test NCTSSoS.assert_moment_linear_data_invariants(direct_linear) === nothing

        mp_complex, report_complex = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            real_moment_matrix=false,
            axis_rotation_equalities=true,
        )
        direct_complex, direct_complex_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            real_moment_matrix=false,
            axis_rotation_equalities=true,
        )
        @test report_complex.real_moment_matrix == false
        @test direct_complex_report.real_moment_matrix == false
        @test report_complex.zero_constraint_count == direct_complex_report.zero_constraint_count
        @test translation_zero_origin_histogram(mp_complex) ==
            translation_zero_origin_histogram(direct_complex)
        @test all(
            zc -> zc.origin.label.feature == :axis_rotation_equality,
            direct_complex.zero_constraints,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(
            mp_complex.linear,
            mp_complex.constraints,
        ) === nothing
        @test NCTSSoS.assert_moment_linear_data_invariants(direct_complex) === nothing

        σx, _, _ = ops
        field_pop = polyopt(sum(σx), registry)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            field_pop,
            ops,
            1;
            sign_symmetry=false,
            axis_rotation_equalities=true,
        )
    end

    @testset "translation path carries block labels into linear metadata" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        mp, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            contiguous_rdm_k=2,
        )

        origins = [block.meta.origin for block in mp.linear.psd_blocks_lin]
        psd_cones = [cone for (cone, _) in mp.constraints if cone in (:PSD, :HPSD)]

        @test length(origins) == length(report.block_labels)
        @test [block.meta.cone for block in mp.linear.psd_blocks_lin] == psd_cones
        @test all(origin -> hasproperty(origin, :label), origins)
        if all(origin -> hasproperty(origin, :label), origins)
            @test [origin.label for origin in origins] == report.block_labels
        end
        @test report.block_coefficient_domains == [
            origin.transform === nothing ? nothing : origin.transform.coefficient_domain
            for origin in origins
        ]
        @test report.block_exact_coefficient_domains == [
            origin.transform === nothing ? nothing : origin.transform.exact_coefficient_domain
            for origin in origins
        ]
        @test all(origin -> hasproperty(origin, :logical_row_labels), origins)
        if all(origin -> hasproperty(origin, :logical_row_labels), origins)
            @test [length(origin.logical_row_labels) for origin in origins] == report.logical_block_sizes
        end
        @test all(origin -> hasproperty(origin, :transform), origins)
        @test all(origin -> origin.transform === nothing || origin.transform.family == :translation_dft, origins)
        @test all(origin -> origin.transform === nothing || origin.transform.coefficient_domain == :cyclotomic_float64, origins)
        @test all(origin -> origin.transform === nothing || origin.transform.exact_coefficient_domain == :cyclotomic, origins)
        @test all(
            origin -> origin.transform === nothing || origin.transform.n_sites == n,
            origins,
        )
        @test isdefined(NCTSSoS, :translation_linear_provenance)
        if isdefined(NCTSSoS, :translation_linear_provenance)
            provenance = translation_linear_provenance(mp)
            @test provenance == translation_linear_provenance(mp.linear)
            @test length(provenance.psd_blocks) == length(report.block_labels)
            @test [block.index for block in provenance.psd_blocks] ==
                collect(eachindex(report.block_labels))
            @test [block.cone for block in provenance.psd_blocks] == psd_cones
            @test [block.size for block in provenance.psd_blocks] == report.psd_block_sizes
            @test [block.label for block in provenance.psd_blocks] == report.block_labels
            @test [length(block.row_labels) for block in provenance.psd_blocks] ==
                report.psd_block_sizes
            @test [length(block.logical_row_labels) for block in provenance.psd_blocks] ==
                report.logical_block_sizes
            @test [block.coefficient_domain for block in provenance.psd_blocks] ==
                report.block_coefficient_domains
            @test [block.exact_coefficient_domain for block in provenance.psd_blocks] ==
                report.block_exact_coefficient_domains
            @test all(
                block -> block.transform === nothing ||
                    block.transform_family == :translation_dft,
                provenance.psd_blocks,
            )
            @test isempty(provenance.zero_constraints)
        end
        @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
    end

    @testset "translation SOS dual blocks keep provenance labels" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        mp, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            contiguous_rdm_k=2,
        )

        @test isdefined(NCTSSoS, :sos_dual_blocks)
        if isdefined(NCTSSoS, :sos_dual_blocks)
            sos = NCTSSoS.sos_dualize(mp)
            dual_blocks = NCTSSoS.sos_dual_blocks(sos)

            @test length(dual_blocks) == length(report.block_labels)
            @test [block.meta.origin.label for block in dual_blocks] == report.block_labels
            @test all(block.representation == :psd for block in dual_blocks)
            @test [size(block.variable, 1) for block in dual_blocks] == report.psd_block_sizes
            @test [size(block.variable, 2) for block in dual_blocks] == report.psd_block_sizes
        end

        mp_complex, report_complex = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            real_moment_matrix=false,
        )
        if isdefined(NCTSSoS, :sos_dual_blocks)
            sos_complex = NCTSSoS.sos_dualize(mp_complex)
            dual_blocks_complex = NCTSSoS.sos_dual_blocks(sos_complex)

            @test [block.meta.origin.label for block in dual_blocks_complex] == report_complex.block_labels
            @test all(block.representation == :hermitian_lift for block in dual_blocks_complex)
            @test [size(block.variable, 1) for block in dual_blocks_complex] == 2 .* report_complex.psd_block_sizes
            @test [size(block.variable, 2) for block in dual_blocks_complex] == 2 .* report_complex.psd_block_sizes
        end

        sos_complex = NCTSSoS.sos_dualize(mp_complex)
        set_optimizer(sos_complex.model, Clarabel.Optimizer)
        set_silent(sos_complex.model)
        optimize!(sos_complex.model)
        @test termination_status(sos_complex.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        dual_values_complex = sos_dual_block_values(mp_complex, sos_complex)
        @test all(block.representation == :hermitian_lift for block in dual_values_complex)
        @test [size(block.native_value, 1) for block in dual_values_complex] == report_complex.psd_block_sizes
        @test [size(block.native_value, 2) for block in dual_values_complex] == report_complex.psd_block_sizes
        @test all(
            block -> norm(block.native_value - block.native_value', Inf) <= 1e-6,
            dual_values_complex,
        )
    end

    @testset "translation SOS dual certificate residual is small" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        _, _, σz = ops
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry; eq_constraints=[sum(σz)])

        mp, _ = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
        )
        sos = NCTSSoS.sos_dualize(mp)
        set_optimizer(sos.model, Clarabel.Optimizer)
        set_silent(sos.model)
        optimize!(sos.model)

        @test termination_status(sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        dual_block_values = sos_dual_block_values(sos)
        @test sos_dual_block_values(mp, sos) == dual_block_values
        @test length(dual_block_values) == length(mp.linear.psd_blocks_lin)
        @test [block.label for block in dual_block_values] == [
            block.meta.origin.label for block in mp.linear.psd_blocks_lin
        ]
        @test [block.logical_row_labels for block in dual_block_values] == [
            block.meta.origin.logical_row_labels for block in mp.linear.psd_blocks_lin
        ]
        @test [block.transform for block in dual_block_values] == [
            block.meta.origin.transform for block in mp.linear.psd_blocks_lin
        ]
        @test [block.transform_family for block in dual_block_values] == [
            block.transform === nothing ? nothing : block.transform.family
            for block in dual_block_values
        ]
        @test all(block -> all(isfinite, block.value), dual_block_values)
        residual = sos_dual_certificate_residual(mp, sos)
        @test residual.max_abs_residual <= 1e-6
        @test abs(residual.identity_residual) <= 1e-6
        @test haskey(residual.residual_by_moment, mp.linear.identity)
        zero_values = sos_zero_dual_values(mp, sos)
        @test length(zero_values) == length(mp.linear.zero_constraints) == 1
        @test only(zero_values).origin.label == (feature=:scalar_equality, index=1)
        @test only(zero_values).feature == :scalar_equality
        @test only(zero_values).decomposition === nothing
        @test only(zero_values).reason === nothing
        @test only(zero_values).term_count == length(only(mp.linear.zero_constraints).form)
        @test isfinite(only(zero_values).value)
        @test isdefined(NCTSSoS, :sos_dual_certificate_diagnostics)
        certificate = sos_dual_certificate_diagnostics(mp, sos)
        @test certificate.residual == residual
        @test certificate.psd_blocks == sos_dual_block_diagnostics(mp, sos)
        @test certificate.zero_duals == sos_zero_dual_diagnostics(mp, sos)
        @test certificate.moment_count == residual.moment_count
        @test certificate.psd_block_count == length(mp.linear.psd_blocks_lin)
        @test certificate.zero_dual_count == length(mp.linear.zero_constraints)
        @test certificate.max_abs_residual == residual.max_abs_residual
        @test certificate.identity_residual == residual.identity_residual
        @test certificate.max_residual_moment == residual.max_residual_moment
        @test certificate.max_residual_value == residual.max_residual_value
    end

    @testset "translation SOS dual block diagnostics expose native cone health" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        @test isdefined(NCTSSoS, :sos_dual_block_diagnostics)
        if isdefined(NCTSSoS, :sos_dual_block_diagnostics)
            for real_moment_matrix in (true, false)
                mp, report = pauli_translation_invariant_moment_relaxation(
                    pop,
                    ops,
                    1;
                    sign_symmetry=false,
                    real_moment_matrix,
                )
                sos = NCTSSoS.sos_dualize(mp)
                set_optimizer(sos.model, Clarabel.Optimizer)
                set_silent(sos.model)
                optimize!(sos.model)

                @test termination_status(sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
                diagnostics = sos_dual_block_diagnostics(mp, sos)
                @test sos_dual_block_diagnostics(sos) == diagnostics
                @test sos_dual_block_diagnostics(sos_dual_block_values(mp, sos)) == diagnostics
                @test length(diagnostics) == length(report.block_labels)
                @test [diag.label for diag in diagnostics] == report.block_labels
                @test [diag.native_value_size for diag in diagnostics] ==
                    [(size, size) for size in report.psd_block_sizes]
                @test [diag.coefficient_domain for diag in diagnostics] ==
                    report.block_coefficient_domains
                @test [diag.exact_coefficient_domain for diag in diagnostics] ==
                    report.block_exact_coefficient_domains
                @test [diag.transform_family for diag in diagnostics] == [
                    diag.transform === nothing ? nothing : diag.transform.family
                    for diag in diagnostics
                ]
                @test all(diag -> diag.native_hermitian_residual <= 1e-6, diagnostics)
                @test all(diag -> isfinite(diag.native_min_eigenvalue), diagnostics)
                @test all(diag -> isfinite(diag.lifted_min_eigenvalue), diagnostics)
                @test all(diag -> diag.native_min_eigenvalue >= -1e-5, diagnostics)

                if real_moment_matrix
                    @test all(diag -> diag.representation == :psd, diagnostics)
                    @test [diag.value_size for diag in diagnostics] ==
                        [(size, size) for size in report.psd_block_sizes]
                else
                    @test all(diag -> diag.representation == :hermitian_lift, diagnostics)
                    @test [diag.value_size for diag in diagnostics] ==
                        [(2size, 2size) for size in report.psd_block_sizes]
                end
            end

            mp_reflection, report_reflection = pauli_translation_invariant_moment_relaxation(
                pop,
                ops,
                1;
                reflection_symmetry=true,
            )
            sos_reflection = NCTSSoS.sos_dualize(mp_reflection)
            set_optimizer(sos_reflection.model, Clarabel.Optimizer)
            set_silent(sos_reflection.model)
            optimize!(sos_reflection.model)

            @test termination_status(sos_reflection.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            reflection_diagnostics = sos_dual_block_diagnostics(mp_reflection, sos_reflection)
            @test [diag.label for diag in reflection_diagnostics] == report_reflection.block_labels
            @test [diag.native_value_size for diag in reflection_diagnostics] ==
                [(size, size) for size in report_reflection.psd_block_sizes]
            @test [diag.coefficient_domain for diag in reflection_diagnostics] ==
                report_reflection.block_coefficient_domains
            @test [diag.exact_coefficient_domain for diag in reflection_diagnostics] ==
                report_reflection.block_exact_coefficient_domains
            @test any(
                diag -> diag.transform_family == :translation_dft_reflection,
                reflection_diagnostics,
            )
            @test any(
                diag -> diag.transform.family == :translation_dft_reflection,
                reflection_diagnostics,
            )
            @test all(diag -> diag.native_hermitian_residual <= 1e-6, reflection_diagnostics)
            @test all(diag -> isfinite(diag.native_min_eigenvalue), reflection_diagnostics)
        end
    end

    @testset "translation path appends contiguous two-site RDM block" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        mp, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            contiguous_rdm_k=2,
        )

        @test report.block_labels[end] == (feature=:contiguous_rdm, k=2)
        @test report.logical_block_sizes[end] == 4
        @test report.psd_block_sizes[end] == 8
        @test mp.constraints[end][1] == :PSD
        @test size(mp.constraints[end][2]) == (8, 8)
        rdm_metrics = translation_report_metrics(report)
        @test rdm_metrics.contiguous_rdm_block_count == 1
        @test rdm_metrics.contiguous_rdm_logical_block_sizes == [4]
        @test rdm_metrics.contiguous_rdm_psd_block_sizes == [8]
        @test rdm_metrics.contiguous_rdm_zero_row_count == 0
        @test rdm_metrics.linear_state_opt_row_count == 0
        @test rdm_metrics.psd_state_opt_block_count == 0
        @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing

        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            contiguous_rdm_k=3,
        )

        registry8, ops8 = create_pauli_variables(1:8)
        pop8 = polyopt(heisenberg_chain_hamiltonian(ops8), registry8)

        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            pop8,
            ops8,
            2;
            sign_symmetry=false,
            contiguous_rdm_k=5,
        )

        mp_ext, report_ext = pauli_translation_invariant_moment_relaxation(
            pop8,
            ops8,
            2;
            sign_symmetry=false,
            contiguous_rdm_k=5,
            contiguous_rdm_support=:extend,
        )

        @test report_ext.block_labels[end] == (feature=:contiguous_rdm, k=5)
        @test report_ext.logical_block_sizes[end] == 32
        @test report_ext.psd_block_sizes[end] == 64
        @test length(mp_ext.linear.moments) > report_ext.n_unique_moment_matrix_elements
        @test report_ext.linear_moment_count == length(mp_ext.linear.moments)
        @test report_ext.linear_moment_count > report_ext.n_unique_moment_matrix_elements
        metrics_ext = translation_report_metrics(report_ext)
        @test metrics_ext.base_moment_count == report_ext.n_unique_moment_matrix_elements
        @test metrics_ext.linear_moment_count == report_ext.linear_moment_count
        @test metrics_ext.contiguous_rdm_block_count == 1
        @test metrics_ext.contiguous_rdm_logical_block_sizes == [32]
        @test metrics_ext.contiguous_rdm_psd_block_sizes == [64]
        @test metrics_ext.contiguous_rdm_zero_row_count == 0
        @test translation_block_feature_histogram(report; kind=:logical) == [
            (feature=:contiguous_rdm, decomposition=nothing, size=4) => 1,
            (feature=:moment_matrix, decomposition=:translation, size=3) => 2,
            (feature=:moment_matrix, decomposition=:translation, size=4) => 1,
        ]
        @test translation_report_metrics(report).psd_block_feature_histogram == [
            (feature=:contiguous_rdm, decomposition=nothing, size=8) => 1,
            (feature=:moment_matrix, decomposition=:translation, size=6) => 2,
            (feature=:moment_matrix, decomposition=:translation, size=8) => 1,
        ]
        @test NCTSSoS.assert_moment_linear_data_invariants(mp_ext.linear, mp_ext.constraints) === nothing

        model_ext, _ = build_jump_model(mp_ext)
        @test num_variables(model_ext) > 0

        direct_ext, direct_report_ext = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop8,
            ops8,
            2;
            sign_symmetry=false,
            contiguous_rdm_k=5,
            contiguous_rdm_support=:extend,
        )
        @test direct_ext isa NCTSSoS.MomentLinearData
        @test direct_report_ext.block_labels[end] == (feature=:contiguous_rdm, k=5)
        @test direct_report_ext.logical_block_sizes[end] == 32
        @test direct_report_ext.psd_block_sizes[end] == 64
        @test direct_report_ext.n_unique_moment_matrix_elements ==
            report_ext.n_unique_moment_matrix_elements
        @test direct_report_ext.linear_moment_count == length(direct_ext.moments)
        @test direct_report_ext.linear_moment_count > direct_report_ext.n_unique_moment_matrix_elements
        direct_ext_metrics = translation_report_metrics(direct_report_ext)
        @test direct_ext_metrics.contiguous_rdm_block_count == 1
        @test direct_ext_metrics.contiguous_rdm_logical_block_sizes == [32]
        @test direct_ext_metrics.contiguous_rdm_psd_block_sizes == [64]
        @test direct_ext_metrics.base_moment_count ==
            direct_report_ext.n_unique_moment_matrix_elements
        @test direct_ext_metrics.linear_moment_count == direct_report_ext.linear_moment_count
        @test NCTSSoS.assert_moment_linear_data_invariants(direct_ext) === nothing

        direct_su2_ext, direct_su2_ext_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop,
                ops,
                1;
                sign_symmetry=false,
                su2_symmetry=true,
                base_su2_extend_rdm=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:su2,
                contiguous_rdm_support=:extend,
            )
        su2_blocks3 = pauli_su2_rdm_blocks(3)
        expected_su2_labels3 = [
            (feature=:contiguous_rdm, k=3, decomposition=:su2, spin2=block.spin2)
            for block in su2_blocks3
        ]
        direct_su2_ext_metrics = translation_report_metrics(direct_su2_ext_report)
        direct_su2_ext_target = pauli_translation_structural_targets(
            4,
            1;
            sign_symmetry=false,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
        )
        @test direct_su2_ext_report.block_labels[end - length(su2_blocks3) + 1:end] ==
            expected_su2_labels3
        @test direct_su2_ext_metrics.su2_moment_symmetry
        @test direct_su2_ext_metrics.su2_rdm_symmetry
        @test direct_su2_ext_metrics.contiguous_rdm_block_count == length(su2_blocks3)
        @test direct_su2_ext_metrics.contiguous_rdm_logical_block_sizes ==
            [block.multiplicity for block in su2_blocks3]
        @test direct_su2_ext_metrics.contiguous_rdm_zero_row_count == 0
        @test direct_su2_ext_metrics.su2_base_zero_row_count > 0
        @test direct_su2_ext_metrics.zero_constraint_count ==
            direct_su2_ext_metrics.su2_base_zero_row_count
        @test any(
            key -> startswith(String(key), "su2_wigner_spin_pair_"),
            keys(direct_su2_ext_report.construction_stage_times_ns),
        )
        @test any(
            key ->
                startswith(String(key), "su2_wigner_spin_pair_") &&
                    endswith(String(key), "_form_build"),
            keys(direct_su2_ext_report.construction_stage_times_ns),
        )
        @test any(
            key ->
                startswith(String(key), "su2_wigner_spin_pair_") &&
                    endswith(String(key), "_zero_append"),
            keys(direct_su2_ext_report.construction_stage_times_ns),
        )
        @test direct_su2_ext_metrics.linear_moment_count >
            direct_su2_ext_metrics.base_moment_count
        direct_su2_ext_transforms = [
            block.meta.origin.transform
            for block in direct_su2_ext.psd_blocks_lin
        ]
        @test direct_su2_ext_report.block_coefficient_domains ==
            [transform.coefficient_domain for transform in direct_su2_ext_transforms]
        @test direct_su2_ext_report.block_exact_coefficient_domains ==
            [transform.exact_coefficient_domain for transform in direct_su2_ext_transforms]
        base_su2_ext_block_count =
            length(direct_su2_ext_report.block_labels) - length(su2_blocks3)
        @test all(
            transform -> transform isa PauliSU2BasisTransformBlock,
            direct_su2_ext_transforms[1:base_su2_ext_block_count],
        )
        @test all(
            transform -> transform.coefficient_domain == :complex_algebraic_float64,
            direct_su2_ext_transforms[1:base_su2_ext_block_count],
        )
        @test all(
            transform -> transform.exact_coefficient_domain == :complex_sqrt_rational,
            direct_su2_ext_transforms[1:base_su2_ext_block_count],
        )
        direct_su2_ext_rdm_transforms =
            direct_su2_ext_transforms[(base_su2_ext_block_count + 1):end]
        @test [transform.family for transform in direct_su2_ext_rdm_transforms] ==
            fill(:pauli_su2_rdm, length(su2_blocks3))
        @test [
            (
                transform.k,
                transform.spin2,
                transform.multiplicity,
                transform.reference_m2,
                transform.coefficient_domain,
                transform.exact_coefficient_domain,
            )
            for transform in direct_su2_ext_rdm_transforms
        ] == [
            (
                3,
                block.spin2,
                block.multiplicity,
                -block.spin2,
                :algebraic_float64,
                :sqrt_rational,
            )
            for block in su2_blocks3
        ]
        @test direct_su2_ext_target.base_su2_extend_rdm
        @test direct_su2_ext_target.logical_block_sizes ==
            direct_su2_ext_report.logical_block_sizes
        @test direct_su2_ext_target.psd_block_sizes ==
            direct_su2_ext_report.psd_block_sizes
        @test direct_su2_ext_target.logical_block_feature_histogram ==
            direct_su2_ext_metrics.logical_block_feature_histogram
        @test direct_su2_ext_target.psd_block_feature_histogram ==
            direct_su2_ext_metrics.psd_block_feature_histogram
        @test direct_su2_ext_target.rdm_logical_block_sizes ==
            direct_su2_ext_metrics.contiguous_rdm_logical_block_sizes
        @test direct_su2_ext_target.rdm_psd_block_sizes ==
            direct_su2_ext_metrics.contiguous_rdm_psd_block_sizes
        @test direct_su2_ext_metrics.product_cache_hits > 0
        @test direct_su2_ext_metrics.product_cache_misses > 0

        direct_su2_ext_complex, direct_su2_ext_complex_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop,
                ops,
                1;
                sign_symmetry=false,
                real_moment_matrix=false,
                su2_symmetry=true,
                base_su2_extend_rdm=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:su2,
                contiguous_rdm_support=:extend,
            )
        direct_su2_ext_complex_metrics =
            translation_report_metrics(direct_su2_ext_complex_report)
        @test !direct_su2_ext_complex_report.real_moment_matrix
        @test direct_su2_ext_complex_metrics.su2_moment_symmetry
        @test direct_su2_ext_complex_metrics.su2_rdm_symmetry
        @test direct_su2_ext_complex_metrics.contiguous_rdm_block_count ==
            length(su2_blocks3)
        @test all(
            block -> block.meta.cone == :HPSD,
            direct_su2_ext_complex.psd_blocks_lin,
        )
        @test direct_su2_ext_complex_report.psd_block_sizes ==
            direct_su2_ext_complex_report.logical_block_sizes
        @test direct_su2_ext_complex_metrics.contiguous_rdm_psd_block_sizes ==
            direct_su2_ext_complex_metrics.contiguous_rdm_logical_block_sizes
        @test NCTSSoS.assert_moment_linear_data_invariants(
            direct_su2_ext_complex,
        ) === nothing
        @test all(
            zc -> zc.trusted_self_adjoint,
            direct_su2_ext_complex.zero_constraints,
        )
        @test all(
            zc -> _linear_form_is_self_adjoint(direct_su2_ext_complex, zc.form),
            direct_su2_ext_complex.zero_constraints,
        )
        @test direct_su2_ext_target.product_cache_hits ==
            direct_su2_ext_metrics.product_cache_hits
        @test direct_su2_ext_target.product_cache_misses ==
            direct_su2_ext_metrics.product_cache_misses
        @test direct_su2_ext_target.product_cache_lookups ==
            direct_su2_ext_metrics.product_cache_lookups
        @test direct_su2_ext_target.product_cache_entries ==
            direct_su2_ext_metrics.product_cache_entries
        @test direct_su2_ext_target.product_cache_hit_rate ==
            direct_su2_ext_metrics.product_cache_hit_rate
        @test direct_su2_ext_target.su2_base_accounting_record_count ==
            direct_su2_ext_target.momentum_sector_count
        @test direct_su2_ext_target.su2_base_singlet_channel_count == 1
        @test direct_su2_ext_target.su2_base_singlet_channel_support_counts == [0 => 1]
        @test direct_su2_ext_target.su2_base_singlet_channel_equality_row_count == 3
        @test direct_su2_ext_target.su2_base_zero_row_budget >=
            direct_su2_ext_metrics.su2_base_zero_row_count > 0
        @test direct_su2_ext_target.su2_base_offblock_zero_row_budget ==
            2 * direct_su2_ext_target.su2_base_offblock_entry_count
        @test direct_su2_ext_target.su2_base_magnetic_copy_zero_row_budget ==
            direct_su2_ext_target.su2_base_copy_entry_count
        @test direct_su2_ext_target.su2_base_zero_row_budget ==
            direct_su2_ext_target.su2_base_offblock_zero_row_budget +
            direct_su2_ext_target.su2_base_magnetic_copy_zero_row_budget
        @test direct_su2_ext_target.su2_base_offblock_zero_row_budget >=
            direct_su2_ext_metrics.su2_base_spin_offblock_row_count +
            direct_su2_ext_metrics.su2_base_magnetic_offdiag_row_count
        @test direct_su2_ext_target.su2_base_magnetic_copy_zero_row_budget >=
            direct_su2_ext_metrics.su2_base_magnetic_copy_row_count
        @test direct_su2_ext_target.su2_base_offblock_entry_count > 0
        @test direct_su2_ext_target.su2_base_copy_entry_count > 0
        @test !direct_su2_ext_target.solve_supported
        @test direct_su2_ext_target.solve_blocker == :structural_target_only
        direct_su2_ext_axis_target = pauli_translation_structural_targets(
            4,
            1;
            sign_symmetry=false,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
            axis_rotation_equalities=true,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
        )
        @test direct_su2_ext_axis_target.base_su2_extend_rdm
        @test direct_su2_ext_axis_target.axis_rotation_equalities
        @test direct_su2_ext_axis_target.axis_rotation_equality_row_count > 0
        @test direct_su2_ext_axis_target.known_zero_constraint_feature_histogram ==
            [
                (feature=:axis_rotation_equality, decomposition=nothing, reason=nothing) =>
                    direct_su2_ext_axis_target.axis_rotation_equality_row_count,
            ]
        @test direct_su2_ext_axis_target.add_on_zero_row_count ==
            direct_su2_ext_axis_target.axis_rotation_equality_row_count
        registry6_ext, ops6_ext = create_pauli_variables(1:6)
        h2_ext_hamiltonian = heisenberg_chain_hamiltonian(ops6_ext)
        h2_ext_target = pauli_translation_structural_targets(
            6,
            2;
            sign_symmetry=false,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
            moment_eq_h2=true,
        )
        h2_ext_pop = polyopt(
            h2_ext_hamiltonian,
            registry6_ext;
            moment_eq_constraints=[h2_ext_hamiltonian * h2_ext_hamiltonian],
        )
        _, h2_ext_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            h2_ext_pop,
            ops6_ext,
            2;
            sign_symmetry=false,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
        )
        h2_ext_metrics = translation_report_metrics(h2_ext_report)
        @test h2_ext_target.base_su2_extend_rdm
        @test h2_ext_target.moment_eq_h2
        @test h2_ext_target.logical_block_sizes == h2_ext_report.logical_block_sizes
        @test h2_ext_target.psd_block_sizes == h2_ext_report.psd_block_sizes
        @test h2_ext_target.moment_equality_row_count ==
            h2_ext_metrics.moment_equality_row_count
        @test h2_ext_target.moment_equality_row_count > 0
        @test h2_ext_target.add_on_zero_row_count ==
            h2_ext_target.moment_equality_row_count
        @test h2_ext_metrics.zero_constraint_count ==
            h2_ext_metrics.su2_base_zero_row_count +
            h2_ext_metrics.moment_equality_row_count
        @test h2_ext_target.known_zero_constraint_feature_histogram ==
            [
                (feature=:moment_equality, decomposition=nothing, reason=nothing) =>
                    h2_ext_target.moment_equality_row_count,
            ]
        direct_su2_reflected_ext, direct_su2_reflected_ext_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop,
                ops,
                1;
                sign_symmetry=false,
                reflection_symmetry=true,
                su2_symmetry=true,
                base_su2_extend_rdm=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:su2,
                contiguous_rdm_support=:extend,
            )
        direct_su2_reflected_ext_metrics =
            translation_report_metrics(direct_su2_reflected_ext_report)
        direct_su2_reflected_ext_target = pauli_translation_structural_targets(
            4,
            1;
            sign_symmetry=false,
            reflection_symmetry=true,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:su2,
            contiguous_rdm_support=:extend,
        )
        @test direct_su2_reflected_ext isa NCTSSoS.MomentLinearData
        direct_su2_reflected_ext_transforms = [
            block.meta.origin.transform
            for block in direct_su2_reflected_ext.psd_blocks_lin
        ]
        @test direct_su2_reflected_ext_report.block_coefficient_domains ==
            [transform.coefficient_domain for transform in direct_su2_reflected_ext_transforms]
        @test direct_su2_reflected_ext_report.block_exact_coefficient_domains ==
            [
                transform.exact_coefficient_domain
                for transform in direct_su2_reflected_ext_transforms
            ]
        @test [
            transform.family
            for transform in direct_su2_reflected_ext_transforms[
                (end - length(su2_blocks3) + 1):end
            ]
        ] == fill(:pauli_su2_rdm, length(su2_blocks3))
        @test direct_su2_reflected_ext_target.logical_block_sizes ==
            direct_su2_reflected_ext_report.logical_block_sizes
        @test direct_su2_reflected_ext_target.psd_block_sizes ==
            direct_su2_reflected_ext_report.psd_block_sizes
        @test direct_su2_reflected_ext_target.su2_base_zero_row_budget >=
            direct_su2_reflected_ext_metrics.su2_base_zero_row_count > 0
        @test direct_su2_reflected_ext_target.su2_base_offblock_zero_row_budget ==
            2 * direct_su2_reflected_ext_target.su2_base_offblock_entry_count
        @test direct_su2_reflected_ext_target.su2_base_magnetic_copy_zero_row_budget ==
            direct_su2_reflected_ext_target.su2_base_copy_entry_count
        @test direct_su2_reflected_ext_target.su2_base_zero_row_budget ==
            direct_su2_reflected_ext_target.su2_base_offblock_zero_row_budget +
            direct_su2_reflected_ext_target.su2_base_magnetic_copy_zero_row_budget
        @test direct_su2_reflected_ext_target.logical_block_feature_histogram ==
            direct_su2_reflected_ext_metrics.logical_block_feature_histogram
        @test direct_su2_reflected_ext_target.psd_block_feature_histogram ==
            direct_su2_reflected_ext_metrics.psd_block_feature_histogram
        @test direct_su2_reflected_ext_target.su2_base_accounting_record_count >
            direct_su2_reflected_ext_target.momentum_sector_count
        @test !direct_su2_reflected_ext_target.solve_supported
        @test NCTSSoS.assert_moment_linear_data_invariants(
            direct_su2_reflected_ext,
        ) === nothing
        @test_throws ArgumentError pauli_translation_structural_targets(
            4,
            1;
            sign_symmetry=false,
            su2_symmetry=true,
            base_su2_extend_rdm=true,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:su2,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(direct_su2_ext) === nothing

        one3 = one(first(first(ops)))
        M3 = typeof(one3)
        K3 = typeof(symmetric_canon(NCTSSoS.expval(one3)))
        full_rdm3 = NCTSSoS._translation_contiguous_rdm_linear_block(
            4,
            3,
            K3,
            ComplexF64,
        )
        schur3, states3 = NCTSSoS._pauli_su2_schur_matrix(3)
        schur_rows3 = NCTSSoS._pauli_sparse_transform_rows(transpose(schur3); atol=1e-12)
        columns3 = NCTSSoS._pauli_su2_state_columns(states3)
        orbit_monos3, orbit_indices3 =
            NCTSSoS._contiguous_rdm_reduced_orbit_data(4, 3, M3)
        for block in su2_blocks3
            reference_rows = [
                columns3[(block.spin2, -block.spin2, mult)]
                for mult in 1:block.multiplicity
            ]
            reference_schur_rows = NCTSSoS._select_sparse_transform_rows(
                schur_rows3,
                reference_rows,
            )
            transformed = NCTSSoS._transform_linear_block(
                full_rdm3,
                reference_schur_rows,
                reference_schur_rows;
                atol=1e-12,
            )
            direct = NCTSSoS._translation_contiguous_rdm_su2_reduced_linear_block(
                orbit_monos3,
                orbit_indices3,
                3,
                reference_schur_rows,
                K3,
                ComplexF64;
                atol=1e-12,
            )
            @test _linear_blocks_isapprox(transformed, direct)
        end

        complex_ext_target = pauli_translation_structural_targets(
            8,
            2;
            sign_symmetry=false,
            real_moment_matrix=false,
            contiguous_rdm_k=5,
            contiguous_rdm_support=:extend,
        )
        complex_direct_ext, complex_direct_report_ext =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop8,
                ops8,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                contiguous_rdm_k=5,
                contiguous_rdm_support=:extend,
            )
        @test complex_direct_ext isa NCTSSoS.MomentLinearData
        @test complex_direct_report_ext.real_moment_matrix == false
        @test complex_direct_report_ext.block_labels[end] ==
            (feature=:contiguous_rdm, k=5)
        @test complex_direct_report_ext.logical_block_sizes[end] == 32
        @test complex_direct_report_ext.psd_block_sizes[end] == 32
        @test complex_direct_report_ext.linear_moment_count ==
            length(complex_direct_ext.moments)
        @test complex_direct_report_ext.linear_moment_count >
            complex_direct_report_ext.n_unique_moment_matrix_elements
        complex_direct_ext_metrics = translation_report_metrics(
            complex_direct_report_ext,
        )
        @test complex_ext_target.logical_block_sizes ==
            complex_direct_report_ext.logical_block_sizes
        @test complex_ext_target.psd_block_sizes ==
            complex_direct_report_ext.psd_block_sizes
        @test complex_ext_target.rdm_logical_block_sizes == [32]
        @test complex_ext_target.rdm_psd_block_sizes == [32]
        @test complex_ext_target.contiguous_rdm_zero_row_count == 0
        @test complex_ext_target.add_on_zero_row_count == 0
        @test isempty(complex_ext_target.known_zero_constraint_feature_histogram)
        @test complex_direct_ext_metrics.contiguous_rdm_block_count == 1
        @test complex_direct_ext_metrics.contiguous_rdm_logical_block_sizes == [32]
        @test complex_direct_ext_metrics.contiguous_rdm_psd_block_sizes == [32]
        @test complex_direct_ext_metrics.contiguous_rdm_zero_row_count == 0
        @test complex_direct_ext_metrics.base_moment_count ==
            complex_direct_report_ext.n_unique_moment_matrix_elements
        @test complex_direct_ext_metrics.linear_moment_count ==
            complex_direct_report_ext.linear_moment_count
        @test complex_ext_target.psd_symmetric_entries ==
            complex_direct_ext_metrics.psd_symmetric_entries
        @test complex_ext_target.psd_dense_entries ==
            complex_direct_ext_metrics.psd_dense_entries
        @test complex_ext_target.psd_symmetric_entries ==
            complex_ext_target.psd_dense_entries
        @test NCTSSoS.assert_moment_linear_data_invariants(
            complex_direct_ext,
        ) === nothing
        @test !isempty(complex_direct_ext.free_keys)
        @test_throws ArgumentError build_jump_model(
            complex_direct_ext;
            formulation=:psd_blocks,
            representation=:complex,
        )
        complex_model_ext, _ = build_jump_model(
            complex_direct_ext;
            formulation=:psd_blocks,
            representation=:complex,
            orphan_policy=:free_variables,
        )
        @test num_variables(complex_model_ext) > 0
    end

    @testset "translation path decomposes contiguous RDMs under U(1) symmetry" begin
        n = 8
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
        nonneutral_pop = polyopt(sum(first(ops)), registry)

        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:u1,
        )
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            nonneutral_pop,
            ops,
            2;
            sign_symmetry=false,
            u1_symmetry=true,
            contiguous_rdm_k=2,
            contiguous_rdm_decomposition=:u1,
        )

        mp, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            u1_symmetry=true,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:u1,
        )

        @test report.block_labels[end - 3:end] == [
            (feature=:contiguous_rdm, k=3, decomposition=:u1, magnetization=0),
            (feature=:contiguous_rdm, k=3, decomposition=:u1, magnetization=1),
            (feature=:contiguous_rdm, k=3, decomposition=:u1, magnetization=2),
            (feature=:contiguous_rdm, k=3, decomposition=:u1, magnetization=3),
        ]
        @test report.logical_block_sizes[end - 3:end] == [1, 3, 3, 1]
        @test [
            length(block.meta.origin.logical_row_labels)
            for block in mp.linear.psd_blocks_lin[end - 3:end]
        ] == [1, 3, 3, 1]
        @test mp.linear.psd_blocks_lin[end - 3].meta.origin.logical_row_labels ==
            [(feature=:contiguous_rdm_state, k=3, state=0)]
        @test report.psd_block_sizes[end - 3:end] == [1, 6, 6, 1]
        u1_metrics = translation_report_metrics(report)
        @test u1_metrics.contiguous_rdm_block_count == 4
        @test u1_metrics.contiguous_rdm_logical_block_sizes == [1, 3, 3, 1]
        @test u1_metrics.contiguous_rdm_psd_block_sizes == [1, 6, 6, 1]
        @test u1_metrics.contiguous_rdm_zero_row_count == u1_metrics.zero_constraint_count
        @test all(cone == :PSD for (cone, _) in mp.constraints[end - 3:end])
        expected_sizes = [1, 6, 6, 1]
        @test all(
            size(mat) == (sz, sz) for ((_, mat), sz) in zip(mp.constraints[end - 3:end], expected_sizes)
        )
        @test !isempty(mp.linear.zero_constraints)
        @test report.zero_constraint_count == length(mp.linear.zero_constraints)
        @test translation_report_metrics(report).zero_constraint_count == length(mp.linear.zero_constraints)
        u1_zero_histogram = translation_zero_origin_histogram(mp)
        @test sum(last, u1_zero_histogram; init=0) == length(mp.linear.zero_constraints)
        @test all(first(pair).feature == :contiguous_rdm_zero for pair in u1_zero_histogram)
        @test all(first(pair).decomposition == :u1 for pair in u1_zero_histogram)
        @test any(first(pair).reason == :magnetization_offblock for pair in u1_zero_histogram)
        @test all(
            zc -> zc.origin.label.feature == :contiguous_rdm_zero &&
                zc.origin.label.decomposition == :u1,
            mp.linear.zero_constraints,
        )
        @test any(
            zc -> zc.origin.label.reason == :magnetization_offblock,
            mp.linear.zero_constraints,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing

        direct_u1, direct_u1_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            u1_symmetry=true,
            contiguous_rdm_k=3,
            contiguous_rdm_decomposition=:u1,
        )
        direct_u1_metrics = translation_report_metrics(direct_u1_report)
        @test direct_u1 isa NCTSSoS.MomentLinearData
        @test direct_u1_report.block_labels[end - 3:end] == report.block_labels[end - 3:end]
        @test direct_u1_report.logical_block_sizes[end - 3:end] == [1, 3, 3, 1]
        @test direct_u1_report.psd_block_sizes[end - 3:end] == [1, 6, 6, 1]
        @test direct_u1_metrics.contiguous_rdm_block_count == 4
        @test direct_u1_metrics.contiguous_rdm_logical_block_sizes == [1, 3, 3, 1]
        @test direct_u1_metrics.contiguous_rdm_psd_block_sizes == [1, 6, 6, 1]
        @test direct_u1_metrics.contiguous_rdm_zero_row_count ==
            direct_u1_metrics.zero_constraint_count
        @test translation_zero_origin_histogram(direct_u1) == u1_zero_histogram
        @test NCTSSoS.assert_moment_linear_data_invariants(direct_u1) === nothing

        direct_complex_u1, direct_complex_u1_report =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop,
                ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                u1_symmetry=true,
                contiguous_rdm_k=3,
                contiguous_rdm_decomposition=:u1,
            )
        direct_complex_u1_metrics = translation_report_metrics(direct_complex_u1_report)
        @test direct_complex_u1 isa NCTSSoS.MomentLinearData
        @test direct_complex_u1_report.real_moment_matrix == false
        @test direct_complex_u1_report.block_labels[end - 3:end] ==
            report.block_labels[end - 3:end]
        @test direct_complex_u1_report.logical_block_sizes[end - 3:end] == [1, 3, 3, 1]
        @test direct_complex_u1_report.psd_block_sizes[end - 3:end] == [1, 3, 3, 1]
        @test direct_complex_u1_metrics.contiguous_rdm_block_count == 4
        @test direct_complex_u1_metrics.contiguous_rdm_logical_block_sizes == [1, 3, 3, 1]
        @test direct_complex_u1_metrics.contiguous_rdm_psd_block_sizes == [1, 3, 3, 1]
        @test direct_complex_u1_metrics.contiguous_rdm_zero_row_count ==
            direct_complex_u1_metrics.zero_constraint_count
        @test translation_zero_origin_histogram(direct_complex_u1) == u1_zero_histogram
        @test all(
            zero -> zero.origin.label.feature == :contiguous_rdm_zero &&
                zero.origin.label.decomposition == :u1 &&
                haskey(zero.origin.label, :component) &&
                zero.origin.label.component == :complex &&
                zero.origin.part == :scalar,
            direct_complex_u1.zero_constraints,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(direct_complex_u1) === nothing

        mp4, report4 = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            u1_symmetry=true,
            contiguous_rdm_k=4,
            contiguous_rdm_decomposition=:u1,
        )

        @test report4.block_labels[end - 4:end] == [
            (feature=:contiguous_rdm, k=4, decomposition=:u1, magnetization=0),
            (feature=:contiguous_rdm, k=4, decomposition=:u1, magnetization=1),
            (feature=:contiguous_rdm, k=4, decomposition=:u1, magnetization=2),
            (feature=:contiguous_rdm, k=4, decomposition=:u1, magnetization=3),
            (feature=:contiguous_rdm, k=4, decomposition=:u1, magnetization=4),
        ]
        @test report4.logical_block_sizes[end - 4:end] == [1, 4, 6, 4, 1]
        @test report4.psd_block_sizes[end - 4:end] == [1, 8, 12, 8, 1]
        u1_metrics4 = translation_report_metrics(report4)
        @test u1_metrics4.contiguous_rdm_block_count == 5
        @test u1_metrics4.contiguous_rdm_logical_block_sizes == [1, 4, 6, 4, 1]
        @test u1_metrics4.contiguous_rdm_psd_block_sizes == [1, 8, 12, 8, 1]
        @test u1_metrics4.contiguous_rdm_zero_row_count == u1_metrics4.zero_constraint_count
        @test all(cone == :PSD for (cone, _) in mp4.constraints[end - 4:end])
        @test NCTSSoS.assert_moment_linear_data_invariants(mp4.linear, mp4.constraints) === nothing

        direct_u1_4, direct_u1_report4 = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            u1_symmetry=true,
            contiguous_rdm_k=4,
            contiguous_rdm_decomposition=:u1,
        )
        direct_u1_metrics4 = translation_report_metrics(direct_u1_report4)
        @test direct_u1_4 isa NCTSSoS.MomentLinearData
        @test direct_u1_report4.block_labels[end - 4:end] == report4.block_labels[end - 4:end]
        @test direct_u1_report4.logical_block_sizes[end - 4:end] == [1, 4, 6, 4, 1]
        @test direct_u1_report4.psd_block_sizes[end - 4:end] == [1, 8, 12, 8, 1]
        @test direct_u1_metrics4.contiguous_rdm_block_count == 5
        @test direct_u1_metrics4.contiguous_rdm_logical_block_sizes == [1, 4, 6, 4, 1]
        @test direct_u1_metrics4.contiguous_rdm_psd_block_sizes == [1, 8, 12, 8, 1]
        @test direct_u1_metrics4.contiguous_rdm_zero_row_count ==
            direct_u1_metrics4.zero_constraint_count
        @test translation_zero_origin_histogram(direct_u1_4) ==
            translation_zero_origin_histogram(mp4)
        @test NCTSSoS.assert_moment_linear_data_invariants(direct_u1_4) === nothing

        direct_complex_u1_4, direct_complex_u1_report4 =
            NCTSSoS._pauli_translation_base_linear_relaxation(
                pop,
                ops,
                2;
                sign_symmetry=false,
                real_moment_matrix=false,
                u1_symmetry=true,
                contiguous_rdm_k=4,
                contiguous_rdm_decomposition=:u1,
            )
        direct_complex_u1_metrics4 = translation_report_metrics(direct_complex_u1_report4)
        @test direct_complex_u1_4 isa NCTSSoS.MomentLinearData
        @test direct_complex_u1_report4.real_moment_matrix == false
        @test direct_complex_u1_report4.block_labels[end - 4:end] ==
            report4.block_labels[end - 4:end]
        @test direct_complex_u1_report4.logical_block_sizes[end - 4:end] == [1, 4, 6, 4, 1]
        @test direct_complex_u1_report4.psd_block_sizes[end - 4:end] == [1, 4, 6, 4, 1]
        @test direct_complex_u1_metrics4.contiguous_rdm_block_count == 5
        @test direct_complex_u1_metrics4.contiguous_rdm_logical_block_sizes == [1, 4, 6, 4, 1]
        @test direct_complex_u1_metrics4.contiguous_rdm_psd_block_sizes == [1, 4, 6, 4, 1]
        @test direct_complex_u1_metrics4.contiguous_rdm_zero_row_count ==
            direct_complex_u1_metrics4.zero_constraint_count
        @test translation_zero_origin_histogram(direct_complex_u1_4) ==
            translation_zero_origin_histogram(mp4)
        @test all(
            zero -> zero.origin.label.feature == :contiguous_rdm_zero &&
                zero.origin.label.decomposition == :u1 &&
                haskey(zero.origin.label, :component) &&
                zero.origin.label.component == :complex &&
                zero.origin.part == :scalar,
            direct_complex_u1_4.zero_constraints,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(direct_complex_u1_4) === nothing

        registry4, ops4 = create_pauli_variables(1:4)
        pop4 = polyopt(heisenberg_chain_hamiltonian(ops4), registry4)
        full_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                contiguous_rdm_k=2,
            )
        end
        u1_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                u1_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:u1,
            )
        end
        @test termination_status(full_result.model) == JuMP.MOI.OPTIMAL
        @test termination_status(u1_result.model) == JuMP.MOI.OPTIMAL
        @test full_result.solve_time_ns > 0
        @test u1_result.solve_time_ns > 0
        @test u1_result.objective ≈ full_result.objective atol = 1e-6
        @test sum(u1_result.report.psd_block_sizes[end - 2:end]) <
            full_result.report.psd_block_sizes[end]
        su2_full_rdm_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                contiguous_rdm_k=2,
            )
        end
        direct_su2_full_rdm_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                su2_symmetry=true,
                contiguous_rdm_k=2,
            )
        end
        @test termination_status(su2_full_rdm_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(direct_su2_full_rdm_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test translation_symmetry_profile(su2_full_rdm_result).su2_moment_symmetry
        @test translation_symmetry_profile(su2_full_rdm_result).contiguous_rdm
        @test !translation_symmetry_profile(su2_full_rdm_result).su2_rdm_symmetry
        @test translation_solve_support(su2_full_rdm_result).supported
        @test direct_su2_full_rdm_result.moment_problem isa NCTSSoS.MomentLinearData
        @test translation_solve_support(direct_su2_full_rdm_result).supported
        @test su2_full_rdm_result.objective ≈ full_result.objective atol = 1e-6
        @test direct_su2_full_rdm_result.objective ≈
            su2_full_rdm_result.objective atol = 1e-6
        su2_u1_rdm_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                u1_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:u1,
            )
        end
        direct_su2_u1_rdm_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                su2_symmetry=true,
                u1_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:u1,
            )
        end
        @test termination_status(su2_u1_rdm_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(direct_su2_u1_rdm_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test translation_symmetry_profile(su2_u1_rdm_result).su2_moment_symmetry
        @test translation_symmetry_profile(su2_u1_rdm_result).u1_rdm_symmetry
        @test !translation_symmetry_profile(su2_u1_rdm_result).su2_rdm_symmetry
        @test translation_solve_support(su2_u1_rdm_result).supported
        @test direct_su2_u1_rdm_result.moment_problem isa NCTSSoS.MomentLinearData
        @test translation_solve_support(direct_su2_u1_rdm_result).supported
        @test su2_u1_rdm_result.objective ≈ u1_result.objective atol = 1e-6
        @test direct_su2_u1_rdm_result.objective ≈
            su2_u1_rdm_result.objective atol = 1e-6
        for su2_rdm_result in (su2_full_rdm_result, su2_u1_rdm_result)
            su2_rdm_sos = NCTSSoS.sos_dualize(su2_rdm_result.moment_problem)
            set_optimizer(su2_rdm_sos.model, Clarabel.Optimizer)
            set_silent(su2_rdm_sos.model)
            optimize!(su2_rdm_sos.model)
            @test termination_status(su2_rdm_sos.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            su2_rdm_residual = sos_dual_certificate_residual(
                su2_rdm_result.moment_problem,
                su2_rdm_sos,
            )
            @test su2_rdm_residual.max_abs_residual <= 1e-6
            su2_rdm_diagnostics = sos_dual_block_diagnostics(
                su2_rdm_result.moment_problem,
                su2_rdm_sos,
            )
            @test [diag.label for diag in su2_rdm_diagnostics] ==
                su2_rdm_result.report.block_labels
            @test [diag.native_value_size for diag in su2_rdm_diagnostics] ==
                [(size, size) for size in su2_rdm_result.report.psd_block_sizes]
            @test all(
                diag -> diag.native_hermitian_residual <= 1e-6,
                su2_rdm_diagnostics,
            )
        end

        complex_u1_direct_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                real_moment_matrix=false,
                direct_linear=true,
                u1_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:u1,
            )
        end
        @test termination_status(complex_u1_direct_result.model) == JuMP.MOI.OPTIMAL
        @test complex_u1_direct_result.moment_problem isa NCTSSoS.MomentLinearData
        @test complex_u1_direct_result.objective ≈ u1_result.objective atol = 1e-6
        @test complex_u1_direct_result.report.block_labels[end - 2:end] ==
            u1_result.report.block_labels[end - 2:end]
        @test all(
            zero -> haskey(zero.origin.label, :component) &&
                zero.origin.label.component == :complex &&
                zero.origin.part == :scalar,
            complex_u1_direct_result.moment_problem.zero_constraints,
        )

        complex_extend_err = try
            quiet() do
                pauli_translation_invariant_nctssos(
                    pop4,
                    ops4,
                    1,
                    optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                    dualize=false,
                    sign_symmetry=false,
                    real_moment_matrix=false,
                    direct_linear=true,
                    contiguous_rdm_k=3,
                    contiguous_rdm_support=:extend,
                    formulation=:psd_blocks,
                    representation=:complex,
                )
            end
        catch err
            err
        end
        @test complex_extend_err isa ArgumentError
        @test occursin("canonical moment", sprint(showerror, complex_extend_err))
        @test occursin(
            "orphan_policy=:free_variables",
            sprint(showerror, complex_extend_err),
        )

        complex_extend_symbolic_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                real_moment_matrix=false,
                contiguous_rdm_k=3,
                contiguous_rdm_support=:extend,
            )
        end
        complex_extend_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                real_moment_matrix=false,
                direct_linear=true,
                contiguous_rdm_k=3,
                contiguous_rdm_support=:extend,
                formulation=:psd_blocks,
                representation=:complex,
                orphan_policy=:free_variables,
            )
        end
        @test termination_status(complex_extend_symbolic_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(complex_extend_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test complex_extend_result.moment_problem isa NCTSSoS.MomentLinearData
        @test !isempty(complex_extend_result.moment_problem.free_keys)
        @test complex_extend_result.report.real_moment_matrix == false
        @test complex_extend_result.report.block_labels ==
            complex_extend_symbolic_result.report.block_labels
        @test complex_extend_result.report.psd_block_sizes ==
            complex_extend_symbolic_result.report.psd_block_sizes
        @test complex_extend_result.report.block_labels[end] ==
            (feature=:contiguous_rdm, k=3)
        @test complex_extend_result.report.psd_block_sizes[end] == 8
        @test translation_solve_support(complex_extend_result).supported
        @test complex_extend_result.objective ≈
            complex_extend_symbolic_result.objective atol = 1e-6

        u1_sos = NCTSSoS.sos_dualize(u1_result.moment_problem)
        set_optimizer(u1_sos.model, Clarabel.Optimizer)
        set_silent(u1_sos.model)
        optimize!(u1_sos.model)
        @test termination_status(u1_sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        u1_residual = sos_dual_certificate_residual(u1_result.moment_problem, u1_sos)
        @test u1_residual.max_abs_residual <= 1e-6
        u1_zero_values = sos_zero_dual_values(u1_result.moment_problem, u1_sos)
        @test length(u1_zero_values) == u1_result.report.zero_constraint_count
        @test all(zero -> zero.origin.label.feature == :contiguous_rdm_zero, u1_zero_values)
        @test all(zero -> zero.origin.label.decomposition == :u1, u1_zero_values)
        @test all(zero -> isfinite(zero.value), u1_zero_values)
        u1_diagnostics = sos_dual_block_diagnostics(u1_result.moment_problem, u1_sos)
        @test [diag.label for diag in u1_diagnostics] == u1_result.report.block_labels
        @test [diag.native_value_size for diag in u1_diagnostics] ==
            [(size, size) for size in u1_result.report.psd_block_sizes]
        @test all(diag -> diag.native_hermitian_residual <= 1e-6, u1_diagnostics)
        @test all(diag -> isfinite(diag.native_min_eigenvalue), u1_diagnostics)
    end

    @testset "translation path decomposes contiguous RDMs under SU(2) symmetry" begin
        expected_su2_blocks = Dict(
            0 => [(0, 1, 1)],
            1 => [(1, 1, 2)],
            2 => [(0, 1, 1), (2, 1, 3)],
            3 => [(1, 2, 2), (3, 1, 4)],
            4 => [(0, 2, 1), (2, 3, 3), (4, 1, 5)],
            5 => [(1, 5, 2), (3, 4, 4), (5, 1, 6)],
            6 => [(0, 5, 1), (2, 9, 3), (4, 5, 5), (6, 1, 7)],
            7 => [(1, 14, 2), (3, 14, 4), (5, 6, 6), (7, 1, 8)],
            8 => [(0, 14, 1), (2, 28, 3), (4, 20, 5), (6, 7, 7), (8, 1, 9)],
        )
        for k in 0:8
            blocks = pauli_su2_rdm_blocks(k)
            @test [(block.spin2, block.multiplicity, block.irrep_dimension) for block in blocks] ==
                expected_su2_blocks[k]
            @test sum(block.multiplicity * block.irrep_dimension for block in blocks) == 1 << k
        end
        rdm_metrics4 = pauli_su2_rdm_metrics(4)
        @test rdm_metrics4.full_side == 16
        @test rdm_metrics4.reduced_block_sizes == [2, 3, 1]
        @test rdm_metrics4.real_psd_block_sizes == [4, 6, 1]
        @test rdm_metrics4.reduced_max_block == 3
        @test rdm_metrics4.real_psd_max_block == 6
        @test rdm_metrics4.reduced_total_block_side == 6
        @test rdm_metrics4.real_psd_total_block_side == 11
        @test rdm_metrics4.full_dense_entries == 16^2
        @test rdm_metrics4.full_symmetric_entries == 16 * 17 ÷ 2
        @test rdm_metrics4.reduced_dense_entries == 2^2 + 3^2 + 1^2
        @test rdm_metrics4.reduced_symmetric_entries == 2 * 3 ÷ 2 + 3 * 4 ÷ 2 + 1
        @test rdm_metrics4.real_psd_dense_entries == 4^2 + 6^2 + 1^2
        @test rdm_metrics4.real_psd_symmetric_entries ==
            4 * 5 ÷ 2 + 6 * 7 ÷ 2 + 1
        @test rdm_metrics4.block_coefficient_domain_histogram ==
            [:algebraic_float64 => length(rdm_metrics4.blocks)]
        @test rdm_metrics4.block_exact_coefficient_domain_histogram ==
            [:sqrt_rational => length(rdm_metrics4.blocks)]
        @test rdm_metrics4.reduced_dense_bytes == rdm_metrics4.reduced_dense_entries * sizeof(Float64)
        @test rdm_metrics4.real_psd_symmetric_bytes ==
            rdm_metrics4.real_psd_symmetric_entries * sizeof(Float64)
        rdm_metrics8 = pauli_su2_rdm_metrics(8)
        @test rdm_metrics8.full_side == 256
        @test rdm_metrics8.reduced_block_sizes == [14, 28, 20, 7, 1]
        @test rdm_metrics8.real_psd_block_sizes == [28, 56, 40, 14, 1]
        @test rdm_metrics8.reduced_dense_entries == 14^2 + 28^2 + 20^2 + 7^2 + 1
        @test rdm_metrics8.real_psd_max_block == 56
        @test rdm_metrics8.real_psd_total_block_side == 139
        @test_throws ArgumentError pauli_su2_rdm_metrics(4; scalar_bytes=0)
        schur_diag4 = pauli_su2_schur_diagnostics(4)
        @test schur_diag4.dimension == 16
        @test schur_diag4.state_count == 16
        @test schur_diag4.spin_multiplicities == [0 => 2, 2 => 3, 4 => 1]
        @test schur_diag4.unitarity_residual <= 1e-12
        @test schur_diag4.sz_residual <= 1e-12
        @test schur_diag4.casimir_residual <= 1e-10
        @test schur_diag4.max_residual <= 1e-10
        @test schur_diag4.coefficient_domain == :algebraic_float64
        @test schur_diag4.exact_coefficient_domain == :sqrt_rational
        schur_diag8 = pauli_su2_schur_diagnostics(8)
        @test schur_diag8.dimension == 256
        @test schur_diag8.state_count == 256
        @test schur_diag8.spin_multiplicities == [0 => 14, 2 => 28, 4 => 20, 6 => 7, 8 => 1]
        @test schur_diag8.unitarity_residual <= 1e-12
        @test schur_diag8.sz_residual <= 1e-12
        @test schur_diag8.casimir_residual <= 1e-10
        @test schur_diag8.max_residual <= 1e-10
        @test schur_diag8.coefficient_domain == :algebraic_float64
        @test schur_diag8.exact_coefficient_domain == :sqrt_rational
        @test_throws DomainError pauli_su2_schur_diagnostics(-1)
        spin_diag4 = pauli_su2_spin_diagnostics(4)
        @test spin_diag4.dimension == 16
        @test spin_diag4.state_count == 16
        @test spin_diag4.sz_residual <= 1e-12
        @test spin_diag4.casimir_residual <= 1e-10
        spin_diag8 = pauli_su2_spin_diagnostics(8)
        @test spin_diag8.dimension == 256
        @test spin_diag8.state_count == 256
        @test spin_diag8.sz_residual <= 1e-12
        @test spin_diag8.casimir_residual <= 1e-10
        @test_throws DomainError pauli_su2_spin_diagnostics(-1)
        rdm_reduction4 = pauli_su2_rdm_reduction_diagnostics(4)
        @test rdm_reduction4.full_dense_entries == 16^2
        @test rdm_reduction4.reduced_dense_entries == 2^2 + 3^2 + 1
        @test rdm_reduction4.offblock_entry_count == 220
        @test rdm_reduction4.copy_entry_count == 22
        @test rdm_reduction4.accounted_entry_count == 16^2
        rdm_reduction8 = pauli_su2_rdm_reduction_diagnostics(8)
        @test rdm_reduction8.reduced_block_sizes == [14, 28, 20, 7, 1]
        @test rdm_reduction8.offblock_entry_count + rdm_reduction8.copy_entry_count +
            rdm_reduction8.reduced_dense_entries == rdm_reduction8.full_dense_entries
        @test_throws DomainError pauli_su2_rdm_reduction_diagnostics(-1)
        qmb_rdm8 = pauli_qmbcertify_rdm_block_metrics(8)
        qmb_rdm8_blocks = pauli_qmbcertify_rdm_blocks(8)
        @test length.(qmb_rdm8_blocks) == [72, 64, 56]
        @test qmb_rdm8_blocks[1][1:4] == [1, 16, 24, 28]
        @test qmb_rdm8_blocks[2][1:4] == [2, 3, 5, 9]
        @test qmb_rdm8_blocks[3][1:4] == [4, 6, 7, 10]
        @test last(qmb_rdm8_blocks[1]) == 256
        @test length(unique(vcat(qmb_rdm8_blocks...))) == sum(length, qmb_rdm8_blocks)
        @test qmb_rdm8.block_sizes == [72, 64, 56]
        @test qmb_rdm8.n_blocks == 3
        @test qmb_rdm8.max_block == 72
        @test qmb_rdm8.total_block_side == 192
        @test qmb_rdm8.dense_entries == 72^2 + 64^2 + 56^2
        @test qmb_rdm8.symmetric_entries == 6_304
        @test !qmb_rdm8.requires_construction
        @test pauli_qmbcertify_rdm_block_sizes(8) == length.(qmb_rdm8_blocks)
        @test pauli_qmbcertify_rdm_block_sizes(9) == [136, 120]
        @test pauli_qmbcertify_rdm_block_sizes(10) == [256, 272, 240]
        qmb_rdm9 = pauli_qmbcertify_rdm_block_metrics(9)
        @test qmb_rdm9.block_sizes == [136, 120]
        @test qmb_rdm9.n_blocks == 2
        @test qmb_rdm9.max_block == 136
        @test qmb_rdm9.total_block_side == 256
        @test qmb_rdm9.dense_entries == 136^2 + 120^2
        @test qmb_rdm9.symmetric_entries == 16_576
        @test !qmb_rdm9.requires_construction
        qmb_rdm10 = pauli_qmbcertify_rdm_block_metrics(10)
        @test qmb_rdm10.block_sizes == [256, 272, 240]
        @test qmb_rdm10.n_blocks == 3
        @test qmb_rdm10.max_block == 272
        @test qmb_rdm10.total_block_side == 768
        @test qmb_rdm10.dense_entries == 256^2 + 272^2 + 240^2
        @test qmb_rdm10.symmetric_entries == 98_944
        @test !qmb_rdm10.requires_construction
        @test_throws DomainError pauli_qmbcertify_rdm_block_sizes(0)
        @test_throws DomainError pauli_qmbcertify_rdm_blocks(0)
        @test_throws ArgumentError pauli_qmbcertify_rdm_block_sizes(4; ambient_sites=3)
        @test_throws ArgumentError pauli_qmbcertify_rdm_blocks(4; ambient_sites=3)
        @test_throws ArgumentError pauli_qmbcertify_rdm_block_metrics(4; scalar_bytes=0)

        registry20, ops20 = create_pauli_variables(1:20)
        qmb_basis20 = pauli_qmbcertify_chain_basis(ops20, 4; extra=9)
        @test length.(qmb_basis20.basis_by_parity) == [1140, 960]
        @test qmb_basis20.family_count_by_parity == [57, 48]
        @test length(qmb_basis20.unique_basis) == 2051
        @test qmb_basis20.source_row_count == 2100
        @test qmb_basis20.unique_nonidentity_row_count == 2050
        @test qmb_basis20.duplicate_source_row_count == 50
        @test qmb_basis20.unique_row_count_by_parity == [1110, 940]
        @test qmb_basis20.duplicate_row_count_by_parity == [30, 20]
        @test qmb_basis20.short_orbit_family_count == 3
        @test qmb_basis20.short_orbit_overcomplete_row_count == 30
        @test qmb_basis20.base_block_sizes == [
            58,
            fill(114, 9)...,
            57,
            48,
            fill(96, 9)...,
            48,
        ]
        @test qmb_basis20.base_block_histogram == [48 => 2, 57 => 1, 58 => 1, 96 => 9, 114 => 9]
        @test qmb_basis20.base_block_count == 22
        @test qmb_basis20.base_max_block == 114
        @test qmb_basis20.base_total_block_side == 2101
        @test qmb_basis20.base_dense_entries == 211_129
        @test qmb_basis20.base_symmetric_entries == 106_615
        @test length(qmb_basis20.base_block_records) == 22
        @test qmb_basis20.base_block_records[1] == (
            parity=0,
            momentum=0,
            block_size=58,
            family_count=57,
            identity_row=true,
            realified=false,
        )
        @test qmb_basis20.base_block_records[11] == (
            parity=0,
            momentum=10,
            block_size=57,
            family_count=57,
            identity_row=false,
            realified=false,
        )
        @test qmb_basis20.base_block_records[12] == (
            parity=1,
            momentum=0,
            block_size=48,
            family_count=48,
            identity_row=false,
            realified=false,
        )
        @test qmb_basis20.base_block_records[13] == (
            parity=1,
            momentum=1,
            block_size=96,
            family_count=48,
            identity_row=false,
            realified=true,
        )
        @test qmb_basis20.base_support_nonzero_row_count == 37_298
        @test qmb_basis20.base_support_zero_row_count == 18_232
        @test qmb_basis20.base_support_diagonal_nonzero_row_count == 958
        @test qmb_basis20.base_support_offdiagonal_nonzero_row_count == 36_340
        @test qmb_basis20.base_support_unique_count == 3_420
        @test qmb_basis20.base_support_diagonal_unique_count == 251
        @test qmb_basis20.base_support_offdiagonal_unique_count == 3_321
        @test qmb_basis20.base_support_word_length_histogram ==
            [0 => 1, 2 => 10, 4 => 575, 6 => 1816, 8 => 1018]
        qmb_h20 = heisenberg_chain_hamiltonian(ops20)
        qmb_pop20 = polyopt(qmb_h20, registry20)
        qmb_linear20, qmb_report20 =
            NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
                qmb_pop20,
                ops20,
                4;
                extra=9,
            )
        @test qmb_report20.psd_block_sizes == qmb_basis20.base_block_sizes
        @test qmb_report20.logical_block_sizes == qmb_basis20.base_block_sizes
        @test length(qmb_linear20.psd_blocks_lin) == qmb_basis20.base_block_count
        @test length(qmb_linear20.moments) == qmb_basis20.base_support_unique_count
        @test length(qmb_linear20.zero_constraints) == 0
        @test all(
            block -> block.meta.origin.label.feature == :qmbcertify_base,
            qmb_linear20.psd_blocks_lin,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(qmb_linear20) === nothing

        qmb_eq_pop20 = polyopt(
            qmb_h20,
            registry20;
            eq_constraints=[qmb_h20],
        )
        qmb_eq_linear20, qmb_eq_report20 =
            NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
                qmb_eq_pop20,
                ops20,
                4;
                extra=9,
            )
        @test qmb_eq_report20.zero_constraint_count == 1
        @test length(qmb_eq_linear20.zero_constraints) == 1
        @test only(qmb_eq_linear20.zero_constraints).origin.label ==
            (feature=:scalar_equality, index=1)
        @test NCTSSoS.assert_moment_linear_data_invariants(qmb_eq_linear20) === nothing

        qmb_ineq_pop20 = polyopt(
            qmb_h20,
            registry20;
            ineq_constraints=[one(qmb_h20) + qmb_h20],
        )
        qmb_ineq_linear20, qmb_ineq_report20 =
            NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
                qmb_ineq_pop20,
                ops20,
                4;
                extra=9,
            )
        @test qmb_ineq_report20.psd_block_sizes == [qmb_basis20.base_block_sizes; 1]
        @test qmb_ineq_report20.logical_block_sizes == [qmb_basis20.base_block_sizes; 1]
        @test qmb_ineq_report20.block_labels[end] == (feature=:scalar_inequality, index=1)
        @test length(qmb_ineq_linear20.psd_blocks_lin) == qmb_basis20.base_block_count + 1
        @test isempty(qmb_ineq_linear20.zero_constraints)
        @test NCTSSoS.assert_moment_linear_data_invariants(qmb_ineq_linear20) === nothing

        @test qmb_basis20.family_word_length_histogram == [1 => 1, 2 => 50, 3 => 13, 4 => 41]
        @test qmb_basis20.family_translation_orbit_length_histogram == [10 => 3, 20 => 102]
        @test [
            record.word for record in qmb_basis20.family_records
            if record.translation_orbit_length == 10
        ] == [UInt16[1, 31], UInt16[2, 32], UInt16[3, 33]]

        qmb_lso_linear20, qmb_lso_report20 =
            NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
                qmb_pop20,
                ops20,
                4;
                extra=9,
                linear_state_opt_width=7,
        )
        qmb_lso_metrics20 = translation_report_metrics(qmb_lso_report20)
        @test qmb_lso_metrics20.linear_state_opt
        # Sign-canonical LSO row keys deduplicate ± source rows.
        @test 2 * qmb_lso_metrics20.linear_state_opt_row_count == 340
        @test qmb_lso_metrics20.linear_state_opt_row_count ==
            length(qmb_lso_linear20.zero_constraints)
        @test all(
            zero -> zero.origin.label.feature == :linear_state_opt &&
                zero.origin.label.mode == :qmbcertify,
            qmb_lso_linear20.zero_constraints,
        )
        @test all(
            block -> block.meta.origin.label.feature == :qmbcertify_base,
            qmb_lso_linear20.psd_blocks_lin,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(qmb_lso_linear20) === nothing

        qmb_rdm_linear20, qmb_rdm_report20 =
            NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
                qmb_pop20,
                ops20,
                4;
                extra=9,
                contiguous_rdm_k=8,
            )
        qmb_rdm_metrics20 = translation_report_metrics(qmb_rdm_report20)
        @test qmb_rdm_metrics20.contiguous_rdm
        @test qmb_rdm_metrics20.contiguous_rdm_psd_block_sizes ==
            pauli_qmbcertify_rdm_block_sizes(8)
        @test qmb_rdm_report20.psd_block_sizes[1:qmb_basis20.base_block_count] ==
            qmb_basis20.base_block_sizes
        @test qmb_rdm_report20.psd_block_sizes[(qmb_basis20.base_block_count + 1):end] ==
            pauli_qmbcertify_rdm_block_sizes(8)
        @test all(
            block -> block.meta.origin.label.feature == :qmbcertify_base,
            qmb_rdm_linear20.psd_blocks_lin[1:qmb_basis20.base_block_count],
        )
        @test all(
            block -> block.meta.origin.label.feature == :contiguous_rdm &&
                block.meta.origin.label.decomposition == :qmbcertify &&
                block.meta.origin.label.k == 8,
            qmb_rdm_linear20.psd_blocks_lin[(qmb_basis20.base_block_count + 1):end],
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(qmb_rdm_linear20) === nothing

        qmb_pso_linear20, qmb_pso_report20 =
            NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
                qmb_pop20,
                ops20,
                4;
                extra=9,
                psd_state_opt_width=3,
            )
        qmb_pso_metrics20 = translation_report_metrics(qmb_pso_report20)
        qmb_pso_expected_blocks20 = [
            36,
            fill(72, 9)...,
            36,
            28,
            fill(56, 9)...,
            28,
        ]
        @test qmb_pso_metrics20.psd_state_opt
        @test qmb_pso_metrics20.psd_state_opt_psd_block_sizes ==
            qmb_pso_expected_blocks20
        @test qmb_pso_report20.psd_block_sizes[1:qmb_basis20.base_block_count] ==
            qmb_basis20.base_block_sizes
        @test qmb_pso_report20.psd_block_sizes[(qmb_basis20.base_block_count + 1):end] ==
            qmb_pso_expected_blocks20
        @test all(
            block -> block.meta.origin.label.feature == :qmbcertify_base,
            qmb_pso_linear20.psd_blocks_lin[1:qmb_basis20.base_block_count],
        )
        @test all(
            block -> block.meta.origin.label.feature == :psd_state_opt &&
                block.meta.origin.label.mode == :qmbcertify &&
                block.meta.origin.label.width == 3,
            qmb_pso_linear20.psd_blocks_lin[(qmb_basis20.base_block_count + 1):end],
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(qmb_pso_linear20) === nothing

        qmb_pso_rows20 = NCTSSoS._qmbcertify_chain_psd_state_opt_rows(
            qmb_basis20,
            NCTSSoS._qmbcertify_chain_family_representatives(qmb_basis20),
            3,
        )
        qmb_pso_block_basis20 = qmb_pso_rows20.rows_by_parity[1]
        qmb_pso_translated20 = Dict(
            rep => [
                NCTSSoS._translate_pauli_monomial(rep, r, length(ops20[1]))
                for r in 0:(length(ops20[1]) - 1)
            ]
            for rep in qmb_pso_block_basis20
        )
        qmb_pso_uncached20 = NCTSSoS._qmbcertify_chain_psd_state_opt_linear_block(
            qmb_pso_block_basis20,
            3,
            length(ops20[1]),
            qmb_pop20.objective,
            qmb_pso_translated20,
            typeof(symmetric_canon(NCTSSoS.expval(one(first(ops20[1]))))),
            ComplexF64,
        )
        qmb_pso_cached_terms20 =
            NCTSSoS._qmbcertify_chain_psd_state_opt_term_cache(
                qmb_pso_block_basis20,
                length(ops20[1]),
                qmb_pop20.objective,
                qmb_pso_translated20,
            )
        qmb_pso_cached20 = NCTSSoS._qmbcertify_chain_psd_state_opt_linear_block(
            qmb_pso_block_basis20,
            3,
            length(ops20[1]),
            qmb_pso_cached_terms20,
            typeof(symmetric_canon(NCTSSoS.expval(one(first(ops20[1]))))),
            ComplexF64,
        )
        @test _linear_blocks_isapprox(qmb_pso_cached20, qmb_pso_uncached20)

        qmb_combined_linear20, qmb_combined_report20 =
            NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
                qmb_pop20,
                ops20,
                4;
                extra=9,
                contiguous_rdm_k=8,
                linear_state_opt_width=7,
                psd_state_opt_width=3,
            )
        qmb_combined_metrics20 = translation_report_metrics(qmb_combined_report20)
        @test qmb_combined_metrics20.contiguous_rdm
        @test qmb_combined_metrics20.linear_state_opt
        @test qmb_combined_metrics20.psd_state_opt
        @test qmb_combined_metrics20.contiguous_rdm_psd_block_sizes ==
            pauli_qmbcertify_rdm_block_sizes(8)
        @test qmb_combined_metrics20.psd_state_opt_psd_block_sizes ==
            qmb_pso_expected_blocks20
        # Sign-canonical LSO row keys deduplicate ± source rows.
        @test 2 * qmb_combined_metrics20.linear_state_opt_row_count == 672
        @test length(qmb_combined_linear20.zero_constraints) ==
            qmb_combined_metrics20.linear_state_opt_row_count
        @test length(qmb_combined_linear20.psd_blocks_lin) ==
            qmb_basis20.base_block_count +
            length(pauli_qmbcertify_rdm_block_sizes(8)) +
            length(qmb_pso_expected_blocks20)
        qmb_combined_provenance20 =
            translation_linear_provenance(qmb_combined_linear20)
        @test [block.label for block in qmb_combined_provenance20.psd_blocks] ==
            qmb_combined_report20.block_labels
        @test [block.size for block in qmb_combined_provenance20.psd_blocks] ==
            qmb_combined_report20.psd_block_sizes
        @test [
            length(block.logical_row_labels) for block in qmb_combined_provenance20.psd_blocks
        ] == qmb_combined_report20.logical_block_sizes
        qmb_combined_labels20 =
            [block.label for block in qmb_combined_provenance20.psd_blocks]
        qmb_combined_rdm_range20 = (qmb_basis20.base_block_count + 1):(
            qmb_basis20.base_block_count + length(pauli_qmbcertify_rdm_block_sizes(8))
        )
        qmb_combined_pso_range20 = (last(qmb_combined_rdm_range20) + 1):length(
            qmb_combined_labels20,
        )
        @test all(
            label -> label.feature == :qmbcertify_base,
            qmb_combined_labels20[1:qmb_basis20.base_block_count],
        )
        @test all(
            label -> label.feature == :contiguous_rdm &&
                label.decomposition == :qmbcertify &&
                label.k == 8,
            qmb_combined_labels20[qmb_combined_rdm_range20],
        )
        @test all(
            label -> label.feature == :psd_state_opt &&
                label.mode == :qmbcertify &&
                label.width == 3,
            qmb_combined_labels20[qmb_combined_pso_range20],
        )
        @test length(qmb_combined_provenance20.zero_constraints) ==
            qmb_combined_metrics20.linear_state_opt_row_count
        @test all(
            zero -> zero.kind == :zero &&
                zero.feature == :linear_state_opt &&
                zero.label.mode == :qmbcertify &&
                zero.label.width == 7 &&
                zero.term_count > 0,
            qmb_combined_provenance20.zero_constraints,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(qmb_combined_linear20) === nothing
        @test_throws DomainError pauli_qmbcertify_chain_basis(ops20, -1)
        @test_throws ArgumentError pauli_qmbcertify_chain_basis(ops20, 4; extra=-1)

        expected_operator_blocks = Dict(
            0 => [(0, 1, 1)],
            1 => [(0, 1, 1), (2, 1, 3)],
            2 => [(0, 2, 1), (2, 3, 3), (4, 1, 5)],
            3 => [(0, 5, 1), (2, 9, 3), (4, 5, 5), (6, 1, 7)],
            4 => [(0, 14, 1), (2, 28, 3), (4, 20, 5), (6, 7, 7), (8, 1, 9)],
        )
        for k in 0:4
            blocks = pauli_su2_operator_blocks(k)
            @test [(block.spin2, block.multiplicity, block.irrep_dimension) for block in blocks] ==
                expected_operator_blocks[k]
            @test sum(block.multiplicity * block.irrep_dimension for block in blocks) == 4^k
        end
        @test_throws DomainError pauli_su2_operator_blocks(-1)
        @test [(block.spin2, block.multiplicity) for block in pauli_su2_operator_blocks(4)] ==
            [(block.spin2, block.multiplicity) for block in pauli_su2_rdm_blocks(8)]
        operator_metrics4 = pauli_su2_operator_metrics(4)
        @test operator_metrics4.full_side == 4^4
        @test operator_metrics4.reduced_block_sizes == [14, 28, 20, 7, 1]
        @test operator_metrics4.reduced_max_block == 28
        @test operator_metrics4.reduced_total_block_side == 70
        @test operator_metrics4.reduced_dense_entries == 14^2 + 28^2 + 20^2 + 7^2 + 1
        @test operator_metrics4.full_symmetric_entries == 4^4 * (4^4 + 1) ÷ 2
        @test_throws ArgumentError pauli_su2_operator_metrics(4; scalar_bytes=0)
        operator_reduction4 = pauli_su2_operator_reduction_diagnostics(4)
        @test operator_reduction4.full_dense_entries == 4^8
        @test operator_reduction4.reduced_dense_entries == 14^2 + 28^2 + 20^2 + 7^2 + 1
        @test operator_reduction4.offblock_entry_count == 60_636
        @test operator_reduction4.copy_entry_count == 3_470
        @test operator_reduction4.accounted_entry_count == 4^8
        @test_throws DomainError pauli_su2_operator_reduction_diagnostics(-1)
        @test_throws ArgumentError pauli_su2_operator_reduction_diagnostics(4; scalar_bytes=0)
        operator_spin_diag2 = pauli_su2_operator_spin_diagnostics(2)
        @test operator_spin_diag2.dimension == 16
        @test operator_spin_diag2.state_count == 16
        @test operator_spin_diag2.spin_multiplicities == [0 => 2, 2 => 3, 4 => 1]
        @test operator_spin_diag2.unitarity_residual <= 1e-12
        @test operator_spin_diag2.sz_residual <= 1e-12
        @test operator_spin_diag2.casimir_residual <= 1e-10
        @test operator_spin_diag2.offblock_residual <= 1e-10
        operator_spin_diag4 = pauli_su2_operator_spin_diagnostics(4)
        @test operator_spin_diag4.dimension == 256
        @test operator_spin_diag4.state_count == 256
        @test operator_spin_diag4.spin_multiplicities ==
            [0 => 14, 2 => 28, 4 => 20, 6 => 7, 8 => 1]
        @test operator_spin_diag4.unitarity_residual <= 1e-12
        @test operator_spin_diag4.sz_residual <= 1e-12
        @test operator_spin_diag4.casimir_residual <= 1e-10
        @test operator_spin_diag4.offblock_residual <= 1e-10
        @test_throws DomainError pauli_su2_operator_spin_diagnostics(-1)
        @test_throws ArgumentError pauli_su2_operator_spin_diagnostics(4; max_dimension=1)

        expected_word_blocks = Dict(
            0 => [(0, 1, 1)],
            1 => [(2, 1, 3)],
            2 => [(0, 1, 1), (2, 1, 3), (4, 1, 5)],
            3 => [(0, 1, 1), (2, 3, 3), (4, 2, 5), (6, 1, 7)],
            4 => [(0, 3, 1), (2, 6, 3), (4, 6, 5), (6, 3, 7), (8, 1, 9)],
        )
        for support_size in 0:4
            blocks = pauli_su2_word_blocks(support_size)
            @test [(block.spin2, block.multiplicity, block.irrep_dimension) for block in blocks] ==
                expected_word_blocks[support_size]
            @test sum(block.multiplicity * block.irrep_dimension for block in blocks) == 3^support_size
        end
        @test_throws DomainError pauli_su2_word_blocks(-1)
        word_metrics4 = pauli_su2_word_metrics(4)
        @test word_metrics4.full_side == 3^4
        @test word_metrics4.reduced_block_sizes == [3, 6, 6, 3, 1]
        @test word_metrics4.reduced_max_block == 6
        @test word_metrics4.reduced_total_block_side == 19
        @test word_metrics4.reduced_dense_entries == 3^2 + 6^2 + 6^2 + 3^2 + 1
        @test word_metrics4.full_symmetric_entries == 3^4 * (3^4 + 1) ÷ 2
        @test_throws ArgumentError pauli_su2_word_metrics(4; scalar_bytes=0)
        word_reduction4 = pauli_su2_word_reduction_diagnostics(4)
        @test word_reduction4.full_dense_entries == 3^8
        @test word_reduction4.reduced_dense_entries == 3^2 + 6^2 + 6^2 + 3^2 + 1
        @test word_reduction4.offblock_entry_count == 6_192
        @test word_reduction4.copy_entry_count == 278
        @test word_reduction4.accounted_entry_count == 3^8
        @test_throws DomainError pauli_su2_word_reduction_diagnostics(-1)
        @test_throws ArgumentError pauli_su2_word_reduction_diagnostics(4; scalar_bytes=0)
        word_spin_diag4 = pauli_su2_word_spin_diagnostics(4)
        @test word_spin_diag4.dimension == 3^4
        @test word_spin_diag4.state_count == 3^4
        @test word_spin_diag4.spin_multiplicities == [0 => 3, 2 => 6, 4 => 6, 6 => 3, 8 => 1]
        @test word_spin_diag4.unitarity_residual <= 1e-12
        @test word_spin_diag4.sz_residual <= 1e-12
        @test word_spin_diag4.casimir_residual <= 1e-10
        @test word_spin_diag4.offblock_residual <= 1e-10
        @test_throws DomainError pauli_su2_word_spin_diagnostics(-1)
        @test_throws ArgumentError pauli_su2_word_spin_diagnostics(4; max_dimension=1)
        word_transform1 = pauli_su2_word_transform(1)
        @test word_transform1.dimension == 3
        @test word_transform1.spin_multiplicities == [2 => 1]
        @test word_transform1.coefficient_domain == :complex_algebraic_float64
        @test word_transform1.exact_coefficient_domain == :complex_sqrt_rational
        @test word_transform1.transform ≈ ComplexF64[
            inv(sqrt(2)) -im * inv(sqrt(2)) 0
            0 0 1
            -inv(sqrt(2)) -im * inv(sqrt(2)) 0
        ] atol = 1e-12
        @test word_transform1.unitarity_residual <= 1e-12
        word_transform2 = pauli_su2_word_transform(2)
        @test word_transform2.dimension == 9
        @test word_transform2.spin_multiplicities == [0 => 1, 2 => 1, 4 => 1]
        @test word_transform2.unitarity_residual <= 1e-12
        @test word_transform2.transform * word_transform2.transform' ≈ I atol = 1e-12
        word_singlets1 = pauli_su2_word_singlet_channels(1)
        @test word_singlets1.channel_count == 0
        @test size(word_singlets1.transform) == (0, 3)
        @test word_singlets1.singlet_orthonormality_residual == 0.0
        word_singlets2 = pauli_su2_word_singlet_channels(2)
        @test word_singlets2.channel_count == 1
        @test word_singlets2.rows == findall(label -> label.spin2 == 0, word_transform2.labels)
        @test all(label -> label.spin2 == 0 && label.m2 == 0, word_singlets2.row_labels)
        @test word_singlets2.coefficient_domain == :complex_algebraic_float64
        @test word_singlets2.exact_coefficient_domain == :complex_sqrt_rational
        @test word_singlets2.singlet_orthonormality_residual <= 1e-12
        singlet2_abs = abs.(vec(word_singlets2.transform))
        @test singlet2_abs[[1, 5, 9]] ≈ fill(inv(sqrt(3)), 3) atol = 1e-12
        @test maximum(singlet2_abs[setdiff(1:9, [1, 5, 9])]) <= 1e-12
        word_equalities2 = pauli_su2_word_singlet_channel_equalities(2)
        @test word_equalities2.equality_count == 8
        @test word_equalities2.rows == findall(label -> label.spin2 != 0, word_transform2.labels)
        @test all(label -> label.spin2 != 0, word_equalities2.row_labels)
        @test word_equalities2.coefficient_domain == :complex_algebraic_float64
        @test word_equalities2.exact_coefficient_domain == :complex_sqrt_rational
        @test size(word_equalities2.transform) == (8, 9)
        @test word_equalities2.transform * word_equalities2.transform' ≈ I atol = 1e-12
        @test word_equalities2.transform * word_singlets2.transform' ≈ zeros(ComplexF64, 8, 1) atol = 1e-12
        @test word_equalities2.equality_orthonormality_residual <= 1e-12
        @test word_equalities2.singlet_orthogonality_residual <= 1e-12
        @test length(word_equalities2.column_forms) == word_equalities2.equality_count
        @test all(!isempty, word_equalities2.column_forms)
        @test all(
            form -> all(1 <= first(term) <= 9 && abs(last(term)) > 1e-12 for term in form),
            word_equalities2.column_forms,
        )
        word_singlets4 = pauli_su2_word_singlet_channels(4)
        @test word_singlets4.channel_count == 3
        @test all(label -> label.spin2 == 0 && label.m2 == 0, word_singlets4.row_labels)
        @test word_singlets4.transform * word_singlets4.transform' ≈ I atol = 1e-12
        @test word_singlets4.singlet_orthonormality_residual <= 1e-12
        @test_throws DomainError pauli_su2_word_singlet_channels(-1)
        @test_throws ArgumentError pauli_su2_word_singlet_channels(4; max_dimension=1)
        @test_throws DomainError pauli_su2_word_singlet_channel_equalities(-1)
        @test_throws ArgumentError pauli_su2_word_singlet_channel_equalities(4; max_dimension=1)
        @test_throws DomainError pauli_su2_word_singlet_channel_equalities(2; atol=-1)
        word_reflection1 = pauli_su2_word_reflection_diagnostics(1)
        @test word_reflection1.dimension == 3
        @test [
            (diag.spin2, diag.plus_multiplicity, diag.minus_multiplicity)
            for diag in word_reflection1.spin_parities
        ] == [(2, 1, 0)]
        @test word_reflection1.max_reflection_residual <= 1e-12
        word_reflection2 = pauli_su2_word_reflection_diagnostics(2)
        @test word_reflection2.dimension == 9
        @test [
            (diag.spin2, diag.plus_multiplicity, diag.minus_multiplicity, diag.trace)
            for diag in word_reflection2.spin_parities
        ] == [(0, 1, 0, 1), (2, 0, 1, -1), (4, 1, 0, 1)]
        @test word_reflection2.hermitian_residual <= 1e-12
        @test word_reflection2.involution_residual <= 1e-12
        @test word_reflection2.spin_magnetic_offblock_residual <= 1e-12
        @test word_reflection2.magnetic_trace_copy_residual <= 1e-12
        @test word_reflection2.trace_round_residual <= 1e-12
        @test word_reflection2.max_reflection_residual <= 1e-12
        @test_throws DomainError pauli_su2_word_transform(-1)
        @test_throws ArgumentError pauli_su2_word_transform(4; max_dimension=1)
        @test_throws DomainError pauli_su2_word_reflection_diagnostics(-1)
        @test_throws ArgumentError pauli_su2_word_reflection_diagnostics(4; max_dimension=1)
        word_counts = Dict{Int,Int}()
        for support_size in 0:4
            for block in pauli_su2_word_blocks(support_size)
                word_counts[block.spin2] =
                    get(word_counts, block.spin2, 0) + binomial(4, support_size) * block.multiplicity
            end
        end
        @test [
            (block.spin2, block.multiplicity)
            for block in pauli_su2_operator_blocks(4)
        ] == [(spin2, multiplicity) for (spin2, multiplicity) in sort!(collect(word_counts); by=first)]

        _, ops4 = create_pauli_variables(1:4)
        σx4, σy4, σz4 = ops4
        basis_order2 = pauli_contiguous_chain_basis(ops4, 2)
        summary_order2 = pauli_su2_basis_summary(basis_order2)
        @test summary_order2.support_counts == [0 => 1, 1 => 4, 2 => 4]
        @test [(block.spin2, block.multiplicity, block.irrep_dimension) for block in summary_order2.blocks] ==
            [(0, 5, 1), (2, 8, 3), (4, 4, 5)]
        @test sum(block.multiplicity * block.irrep_dimension for block in summary_order2.blocks) ==
            length(basis_order2)
        metrics_order2 = pauli_su2_basis_metrics(basis_order2)
        @test metrics_order2.full_side == 49
        @test metrics_order2.reduced_block_sizes == [5, 8, 4]
        @test metrics_order2.reduced_max_block == 8
        @test metrics_order2.reduced_total_block_side == 17
        @test metrics_order2.full_dense_entries == 49^2
        @test metrics_order2.full_symmetric_entries == 49 * 50 ÷ 2
        @test metrics_order2.reduced_dense_entries == 5^2 + 8^2 + 4^2
        @test metrics_order2.reduced_symmetric_entries == 5 * 6 ÷ 2 + 8 * 9 ÷ 2 + 4 * 5 ÷ 2
        @test metrics_order2.reduced_dense_bytes == metrics_order2.reduced_dense_entries * sizeof(Float64)
        @test_throws ArgumentError pauli_su2_basis_metrics(basis_order2; scalar_bytes=0)
        @test isdefined(NCTSSoS, :pauli_su2_contiguous_chain_structural_targets)
        if isdefined(NCTSSoS, :pauli_su2_contiguous_chain_structural_targets)
            registry6, ops6 = create_pauli_variables(1:6)
            basis6_order2 = pauli_contiguous_chain_basis(ops6, 2)
            actual_target6 = pauli_su2_basis_reduction_diagnostics(basis6_order2)
            analytic_target6 = pauli_su2_contiguous_chain_structural_targets(6, 2)
            @test analytic_target6.basis_size == length(basis6_order2) == 73
            @test analytic_target6.word_count == actual_target6.full_side
            @test analytic_target6.support_counts == actual_target6.support_counts == [0 => 1, 1 => 6, 2 => 6]
            @test analytic_target6.singlet_channel_count == 7
            @test analytic_target6.singlet_channel_support_counts == [0 => 1, 2 => 6]
            @test analytic_target6.singlet_channel_equality_row_count == 66
            @test [(block.spin2, block.multiplicity) for block in analytic_target6.blocks] ==
                [(block.spin2, block.multiplicity) for block in actual_target6.blocks]
            @test analytic_target6.reduced_block_sizes == actual_target6.reduced_block_sizes == [7, 12, 6]
            @test analytic_target6.reduced_dense_entries == actual_target6.reduced_dense_entries
            @test analytic_target6.reduced_symmetric_entries == actual_target6.reduced_symmetric_entries
            @test analytic_target6.offblock_entry_count == actual_target6.offblock_entry_count
            @test analytic_target6.copy_entry_count == actual_target6.copy_entry_count
            @test !analytic_target6.translation_combined
            @test !analytic_target6.solve_supported
            @test analytic_target6.solve_blocker == :structural_target_only
            @test occursin(
                "Structural target only",
                analytic_target6.solve_blocker_reason,
            )
            @test isempty(analytic_target6.solve_unsupported_block_features)
            @test isempty(analytic_target6.solve_unsupported_zero_features)
            @test !analytic_target6.requires_construction

            su2_target100 = pauli_su2_contiguous_chain_structural_targets(100, 4)
            @test su2_target100.basis_size == 12_001
            @test su2_target100.support_counts == [0 => 1, 1 => 100, 2 => 100, 3 => 100, 4 => 100]
            @test [(block.spin2, block.multiplicity) for block in su2_target100.blocks] ==
                [(0, 501), (2, 1_100), (4, 900), (6, 400), (8, 100)]
            @test su2_target100.reduced_block_sizes == [501, 1_100, 900, 400, 100]
            @test su2_target100.reduced_max_block == 1_100
            @test su2_target100.reduced_total_block_side == 3_001
            @test su2_target100.block_coefficient_domain_histogram ==
                [:complex_algebraic_float64 => length(su2_target100.blocks)]
            @test su2_target100.block_exact_coefficient_domain_histogram ==
                [:complex_sqrt_rational => length(su2_target100.blocks)]
            @test su2_target100.full_dense_entries == 144_024_001
            @test su2_target100.reduced_dense_entries == 2_441_001
            @test su2_target100.full_symmetric_entries == 72_018_001
            @test su2_target100.reduced_symmetric_entries == 1_222_001
            @test su2_target100.offblock_entry_count == 134_883_000
            @test su2_target100.copy_entry_count == 6_700_000
            @test su2_target100.accounted_entry_count == su2_target100.full_dense_entries
            @test !su2_target100.translation_combined
            @test !su2_target100.solve_supported
            @test su2_target100.solve_blocker == :structural_target_only
            @test isempty(su2_target100.solve_unsupported_block_features)
            @test isempty(su2_target100.solve_unsupported_zero_features)
            @test isdefined(NCTSSoS, :pauli_su2_translation_orbit_structural_targets)
            translation_su2_target6 = pauli_su2_translation_orbit_structural_targets(6, 2)
            @test translation_su2_target6.orbit_basis_size == 13
            @test translation_su2_target6.singlet_channel_count == 2
            @test translation_su2_target6.singlet_channel_support_counts == [0 => 1, 2 => 1]
            @test translation_su2_target6.singlet_channel_equality_row_count == 11
            @test translation_su2_target6.momentum_sector_count == 4
            @test translation_su2_target6.n_blocks == 12
            @test translation_su2_target6.logical_block_histogram == [1 => 7, 2 => 5]
            @test translation_su2_target6.psd_block_histogram == [2 => 7, 4 => 5]
            @test translation_su2_target6.logical_max_block == 2
            @test translation_su2_target6.psd_max_block == 4
            @test translation_su2_target6.logical_total_block_side == 17
            @test translation_su2_target6.psd_total_block_side == 34
            @test translation_su2_target6.translation_combined
            @test translation_su2_target6.su2_full_dense_entries == 601
            @test translation_su2_target6.su2_active_dense_entries == 75
            @test translation_su2_target6.su2_reduced_dense_entries == 27
            @test translation_su2_target6.offblock_entry_count == 526
            @test translation_su2_target6.copy_entry_count == 48
            @test translation_su2_target6.wigner_offblock_zero_row_budget ==
                2 * translation_su2_target6.offblock_entry_count
            @test translation_su2_target6.wigner_magnetic_copy_zero_row_budget ==
                translation_su2_target6.copy_entry_count
            @test translation_su2_target6.wigner_zero_row_budget ==
                translation_su2_target6.wigner_offblock_zero_row_budget +
                translation_su2_target6.wigner_magnetic_copy_zero_row_budget
            @test translation_su2_target6.accounted_entry_count ==
                translation_su2_target6.su2_full_dense_entries
            target6_accounting_records = translation_su2_target6.su2_accounting_records
            @test translation_su2_target6.su2_accounting_record_count ==
                translation_su2_target6.momentum_sector_count
            @test length(target6_accounting_records) ==
                translation_su2_target6.su2_accounting_record_count
            @test first(target6_accounting_records).label == (
                feature=:moment_matrix,
                decomposition=:translation_su2,
                momentum=0,
            )
            @test first(target6_accounting_records).su2_full_dense_entries == 169
            @test first(target6_accounting_records).su2_active_dense_entries == 21
            @test first(target6_accounting_records).su2_reduced_dense_entries == 9
            @test first(target6_accounting_records).offblock_entry_count == 148
            @test first(target6_accounting_records).copy_entry_count == 12
            @test first(target6_accounting_records).accounted_entry_count == 169
            @test sum(
                record.su2_full_dense_entries for record in target6_accounting_records;
                init=0,
            ) == translation_su2_target6.su2_full_dense_entries
            @test sum(
                record.su2_active_dense_entries for record in target6_accounting_records;
                init=0,
            ) == translation_su2_target6.su2_active_dense_entries
            @test sum(
                record.su2_reduced_dense_entries for record in target6_accounting_records;
                init=0,
            ) == translation_su2_target6.su2_reduced_dense_entries
            @test sum(
                record.offblock_entry_count for record in target6_accounting_records;
                init=0,
            ) == translation_su2_target6.offblock_entry_count
            @test sum(
                record.copy_entry_count for record in target6_accounting_records;
                init=0,
            ) == translation_su2_target6.copy_entry_count
            @test !translation_su2_target6.solve_supported
            @test translation_su2_target6.solve_blocker == :structural_target_only
            @test occursin(
                "Structural target only",
                translation_su2_target6.solve_blocker_reason,
            )
            @test isempty(translation_su2_target6.solve_unsupported_block_features)
            @test isempty(translation_su2_target6.solve_unsupported_zero_features)
            @test !translation_su2_target6.reflection_combined
            @test translation_su2_target6.sign_symmetry_subsumed
            @test length(translation_su2_target6.block_labels) ==
                length(translation_su2_target6.logical_block_sizes)
            @test translation_su2_target6.block_labels[1:3] == [
                (feature=:moment_matrix, decomposition=:translation_su2, momentum=0, spin2=0),
                (feature=:moment_matrix, decomposition=:translation_su2, momentum=0, spin2=2),
                (feature=:moment_matrix, decomposition=:translation_su2, momentum=0, spin2=4),
            ]
            @test translation_su2_target6.block_labels[4:6] == [
                (feature=:moment_matrix, decomposition=:translation_su2, momentum=1, spin2=0),
                (feature=:moment_matrix, decomposition=:translation_su2, momentum=1, spin2=2),
                (feature=:moment_matrix, decomposition=:translation_su2, momentum=1, spin2=4),
            ]
            @test sum(
                pair.second for pair in translation_su2_target6.logical_block_feature_histogram;
                init=0,
            ) == length(translation_su2_target6.logical_block_sizes)
            @test (
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2,
                    momentum=0,
                    spin2=0,
                    size=2,
                ) => 1
            ) in translation_su2_target6.logical_block_feature_histogram
            @test (
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2,
                    momentum=3,
                    spin2=4,
                    size=1,
                ) => 1
            ) in translation_su2_target6.logical_block_feature_histogram
            @test sum(
                pair.second for pair in translation_su2_target6.psd_block_feature_histogram;
                init=0,
            ) == length(translation_su2_target6.psd_block_sizes)
            @test (
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2,
                    momentum=0,
                    spin2=0,
                    size=4,
                ) => 1
            ) in translation_su2_target6.psd_block_feature_histogram
            orbit_reps6 = NCTSSoS._pauli_contiguous_chain_orbit_representatives(
                ops6,
                2;
                periodic=true,
            )
            @test length(orbit_reps6) == translation_su2_target6.orbit_basis_size
            @test_throws ArgumentError pauli_su2_basis_summary(orbit_reps6)
            orbit_summary6 = pauli_su2_translation_orbit_basis_summary(orbit_reps6, 6)
            @test orbit_summary6.word_count == translation_su2_target6.orbit_basis_size
            @test orbit_summary6.support_counts == translation_su2_target6.support_counts
            @test [
                (block.spin2, block.multiplicity)
                for block in orbit_summary6.blocks
            ] == [
                (block.spin2, block.multiplicity)
                for block in translation_su2_target6.zero_momentum_blocks
            ]
            orbit_plan6 = pauli_su2_translation_orbit_basis_reduction_plan(
                orbit_reps6,
                6,
            )
            @test orbit_plan6.support_counts == translation_su2_target6.support_counts
            @test orbit_plan6.full_side == translation_su2_target6.orbit_basis_size
            @test orbit_plan6.transformed_block_sizes == [2, 6, 5]
            @test orbit_plan6.transformed_max_block == 6
            @test orbit_plan6.transformed_total_block_side == orbit_plan6.full_side
            @test orbit_plan6.block_coefficient_domains ==
                fill(:complex_algebraic_float64, length(orbit_plan6.blocks))
            @test orbit_plan6.block_exact_coefficient_domains ==
                fill(:complex_sqrt_rational, length(orbit_plan6.blocks))
            @test all(
                block -> block.coefficient_domain == :complex_algebraic_float64,
                orbit_plan6.blocks,
            )
            @test all(
                block -> block.exact_coefficient_domain == :complex_sqrt_rational,
                orbit_plan6.blocks,
            )
            @test [
                (block.spin2, block.multiplicity)
                for block in orbit_plan6.blocks
            ] == [
                (block.spin2, block.multiplicity)
                for block in translation_su2_target6.zero_momentum_blocks
            ]
            @test first(orbit_plan6.blocks[1].logical_row_labels).feature ==
                :pauli_su2_translation_orbit_multiplicity
            @test first(orbit_plan6.blocks[1].logical_row_labels).support_orbit == ()
            @test first(orbit_plan6.blocks[2].logical_row_labels).support_orbit == (1,)
            @test isdefined(NCTSSoS, :pauli_su2_translation_orbit_basis_transform_blocks)
            orbit_transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
                orbit_reps6,
                6,
            )
            @test [
                (block.spin2, block.multiplicity, block.irrep_dimension, size(block.transform))
                for block in orbit_transform_blocks
            ] == [
                (0, 2, 1, (2, length(orbit_reps6))),
                (2, 2, 3, (6, length(orbit_reps6))),
                (4, 1, 5, (5, length(orbit_reps6))),
            ]
            orbit_transform = reduce(vcat, [block.transform for block in orbit_transform_blocks])
            @test orbit_transform * orbit_transform' ≈ I atol = 1e-12
            @test [block.multiplicity for block in orbit_transform_blocks] ==
                orbit_plan6.reduced_block_sizes
            @test first(orbit_transform_blocks[1].row_labels).feature ==
                :pauli_su2_translation_orbit_basis_state
            @test first(orbit_transform_blocks[1].row_labels).support_orbit == ()
            @test count(
                label -> label.support_orbit == (1, 2) && label.m2 == 0,
                orbit_transform_blocks[2].row_labels,
            ) == 1
            orbit_singlets6 = pauli_su2_translation_orbit_singlet_channels(
                orbit_reps6,
                6,
            )
            @test orbit_singlets6.n_sites == 6
            @test orbit_singlets6.basis_size == length(orbit_reps6)
            @test orbit_singlets6.momentum === nothing
            @test orbit_singlets6.channel_count ==
                translation_su2_target6.singlet_channel_count
            @test orbit_singlets6.row_labels == orbit_transform_blocks[1].row_labels
            @test orbit_singlets6.channel_support_orbits == [(), (1, 2)]
            @test orbit_singlets6.coefficient_domain == :complex_algebraic_float64
            @test orbit_singlets6.exact_coefficient_domain == :complex_sqrt_rational
            @test orbit_singlets6.transform == orbit_transform_blocks[1].transform
            @test orbit_singlets6.singlet_orthonormality_residual <= 1e-12
            orbit_singlets6_k1 = pauli_su2_translation_orbit_singlet_channels(
                orbit_reps6,
                6;
                momentum=1,
            )
            @test orbit_singlets6_k1.momentum == 1
            @test orbit_singlets6_k1.channel_count == orbit_singlets6.channel_count
            @test orbit_singlets6_k1.channel_support_orbits ==
                orbit_singlets6.channel_support_orbits
            @test orbit_singlets6_k1.transform * orbit_singlets6_k1.transform' ≈ I atol = 1e-12
            @test orbit_singlets6_k1.singlet_orthonormality_residual <= 1e-12
            orbit_equalities6 = pauli_su2_translation_orbit_singlet_channel_equalities(
                orbit_reps6,
                6,
            )
            @test orbit_equalities6.n_sites == 6
            @test orbit_equalities6.basis_size == length(orbit_reps6)
            @test orbit_equalities6.momentum === nothing
            @test orbit_equalities6.equality_count ==
                length(orbit_reps6) - orbit_singlets6.channel_count
            @test orbit_equalities6.row_labels == reduce(
                vcat,
                [
                    block.row_labels
                    for block in orbit_transform_blocks
                    if block.spin2 != 0
                ],
            )
            @test count(label -> label.support_orbit == (1, 2), orbit_equalities6.row_labels) == 8
            @test orbit_equalities6.transform * orbit_equalities6.transform' ≈ I atol = 1e-12
            @test orbit_equalities6.transform * orbit_singlets6.transform' ≈
                zeros(ComplexF64, orbit_equalities6.equality_count, orbit_singlets6.channel_count) atol = 1e-12
            @test orbit_equalities6.equality_orthonormality_residual <= 1e-12
            @test orbit_equalities6.singlet_orthogonality_residual <= 1e-12
            @test length(orbit_equalities6.column_forms) == orbit_equalities6.equality_count
            @test length(orbit_equalities6.basis_forms) == orbit_equalities6.equality_count
            @test all(
                form -> all(first(term) in orbit_reps6 && abs(last(term)) > 1e-12 for term in form),
                orbit_equalities6.basis_forms,
            )
            orbit_equalities6_k1 = pauli_su2_translation_orbit_singlet_channel_equalities(
                orbit_reps6,
                6;
                momentum=1,
            )
            @test orbit_equalities6_k1.momentum == 1
            @test orbit_equalities6_k1.equality_count == orbit_equalities6.equality_count
            @test orbit_equalities6_k1.transform * orbit_equalities6_k1.transform' ≈ I atol = 1e-12
            @test orbit_equalities6_k1.singlet_orthogonality_residual <= 1e-12
            @test_throws ArgumentError pauli_su2_translation_orbit_singlet_channels(
                orbit_reps6,
            )
            @test_throws ArgumentError pauli_su2_translation_orbit_singlet_channel_equalities(
                orbit_reps6,
            )
            identity_orbit_block = Matrix{ComplexF64}(I, length(orbit_reps6), length(orbit_reps6))
            orbit_reduced_blocks = pauli_su2_translation_orbit_basis_reduced_moment_blocks(
                identity_orbit_block,
                orbit_reps6,
                6,
            )
            @test [
                (
                    block.spin2,
                    block.multiplicity,
                    block.irrep_dimension,
                    block.reference_m2,
                    size(block.matrix),
                    length(block.row_labels),
                )
                for block in orbit_reduced_blocks
            ] == [
                (0, 2, 1, 0, (2, 2), 2),
                (2, 2, 3, -2, (2, 2), 2),
                (4, 1, 5, -4, (1, 1), 1),
            ]
            @test all(block -> isapprox(block.matrix, I; atol=1e-12), orbit_reduced_blocks)
            @test all(
                block -> block.coefficient_domain == :complex_algebraic_float64,
                orbit_reduced_blocks,
            )
            @test all(
                block -> block.exact_coefficient_domain == :complex_sqrt_rational,
                orbit_reduced_blocks,
            )
            orbit_reduction_residuals =
                pauli_su2_translation_orbit_basis_moment_reduction_diagnostics(
                    identity_orbit_block,
                    orbit_reps6,
                    6,
                )
            @test orbit_reduction_residuals.reduced_block_sizes == [2, 2, 1]
            @test orbit_reduction_residuals.transformed_block_sizes == [2, 6, 5]
            @test orbit_reduction_residuals.reference_m2s == [0, -2, -4]
            @test orbit_reduction_residuals.unitarity_residual <= 1e-12
            @test orbit_reduction_residuals.max_reduction_residual <= 1e-12
            @test isdefined(NCTSSoS, :pauli_su2_translation_orbit_basis_moment_problem)
            if isdefined(NCTSSoS, :pauli_su2_translation_orbit_basis_moment_problem)
                T6 = eltype(first(orbit_reps6).word)
                MP_P6 = Polynomial{PauliAlgebra,T6,ComplexF64}
                orbit_poly_identity = Matrix{MP_P6}(
                    undef,
                    length(orbit_reps6),
                    length(orbit_reps6),
                )
                for row in axes(orbit_poly_identity, 1), col in axes(orbit_poly_identity, 2)
                    orbit_poly_identity[row, col] =
                        row == col ? one(MP_P6) : zero(MP_P6)
                end
                orbit_identity_objective = Polynomial(one(first(ops6[1])))
                @test_throws ArgumentError pauli_su2_translation_orbit_basis_moment_problem(
                    orbit_identity_objective,
                    orbit_poly_identity,
                    orbit_reps6,
                    6,
                )
                orbit_su2_mp = pauli_su2_translation_orbit_basis_moment_problem(
                    orbit_identity_objective,
                    orbit_poly_identity,
                    orbit_reps6,
                    6;
                    assume_su2_invariant=true,
                )
                @test [cone for (cone, _) in orbit_su2_mp.constraints] ==
                    [:HPSD, :HPSD, :HPSD]
                @test [size(mat, 1) for (_, mat) in orbit_su2_mp.constraints] == [2, 2, 1]
                @test [block.size for block in orbit_su2_mp.linear.psd_blocks_lin] == [2, 2, 1]
                @test isempty(orbit_su2_mp.linear.zero_constraints)
                @test [
                    block.meta.origin.label.feature
                    for block in orbit_su2_mp.linear.psd_blocks_lin
                ] == fill(:pauli_su2_translation_orbit_basis_moment, 3)
                @test first(
                    orbit_su2_mp.linear.psd_blocks_lin[1].meta.origin.logical_row_labels,
                ).feature == :pauli_su2_translation_orbit_basis_state
                @test all(
                    block -> block.meta.origin.transform.coefficient_domain ==
                        :complex_algebraic_float64,
                    orbit_su2_mp.linear.psd_blocks_lin,
                )
                orbit_su2_model, _ = build_jump_model(orbit_su2_mp)
                @test !isnothing(orbit_su2_model)
                orbit_singlet_equality_mp =
                    pauli_su2_translation_orbit_basis_moment_problem(
                        orbit_identity_objective,
                        orbit_poly_identity,
                        orbit_reps6,
                        6;
                        assume_su2_invariant=true,
                        singlet_channel_equalities=true,
                    )
                orbit_singlet_zero_constraint_count =
                    count(con -> con[1] == :Zero, orbit_singlet_equality_mp.constraints)
                @test orbit_singlet_zero_constraint_count >= orbit_equalities6.equality_count
                @test orbit_singlet_zero_constraint_count ==
                    length(orbit_singlet_equality_mp.linear.zero_constraints)
                orbit_singlet_zero_labels = [
                    zero.origin.label
                    for zero in orbit_singlet_equality_mp.linear.zero_constraints
                ]
                @test all(
                    label -> label.feature ==
                        :pauli_su2_translation_orbit_singlet_channel_equality,
                    orbit_singlet_zero_labels,
                )
                @test Set(label.row for label in orbit_singlet_zero_labels) ==
                    Set(1:orbit_equalities6.equality_count)
                @test all(label -> label.n_sites == 6, orbit_singlet_zero_labels)
                @test all(label -> label.momentum === nothing, orbit_singlet_zero_labels)
                @test all(label -> label.su2_label.spin2 != 0, orbit_singlet_zero_labels)
                @test all(
                    zero -> zero.origin isa NCTSSoS.PauliSU2SingletChannelEqualityOrigin,
                    orbit_singlet_equality_mp.linear.zero_constraints,
                )
                @test NCTSSoS.assert_moment_linear_data_invariants(
                    orbit_singlet_equality_mp.linear,
                    orbit_singlet_equality_mp.constraints,
                ) === nothing
            end
            @test isdefined(NCTSSoS, :pauli_su2_translation_orbit_zero_momentum_blocks)
            if isdefined(NCTSSoS, :pauli_su2_translation_orbit_zero_momentum_blocks)
                orbit_zero_momentum_blocks =
                    pauli_su2_translation_orbit_zero_momentum_blocks(orbit_reps6, 6)
                @test [
                    (
                        block.spin2,
                        block.multiplicity,
                        block.irrep_dimension,
                        block.reference_m2,
                        size(block.matrix),
                    )
                    for block in orbit_zero_momentum_blocks
                ] == [
                    (0, 2, 1, 0, (2, 2)),
                    (2, 2, 3, -2, (2, 2)),
                    (4, 1, 5, -4, (1, 1)),
                ]
                @test [block.label for block in orbit_zero_momentum_blocks] == [
                    (feature=:moment_matrix, decomposition=:translation_su2, momentum=0, spin2=0),
                    (feature=:moment_matrix, decomposition=:translation_su2, momentum=0, spin2=2),
                    (feature=:moment_matrix, decomposition=:translation_su2, momentum=0, spin2=4),
                ]
                @test all(
                    block -> all(entry -> entry isa Polynomial{PauliAlgebra}, block.matrix),
                    orbit_zero_momentum_blocks,
                )
                @test first(orbit_zero_momentum_blocks[1].row_labels).feature ==
                    :pauli_su2_translation_orbit_basis_state
                @test all(
                    block -> block.coefficient_domain == :complex_algebraic_float64,
                    orbit_zero_momentum_blocks,
                )
            end
            @test isdefined(NCTSSoS, :pauli_su2_translation_orbit_momentum_blocks)
            if isdefined(NCTSSoS, :pauli_su2_translation_orbit_momentum_blocks)
                orbit_k1_blocks = pauli_su2_translation_orbit_momentum_blocks(
                    orbit_reps6,
                    6,
                    1,
                )
                @test [
                    (
                        block.spin2,
                        block.multiplicity,
                        block.irrep_dimension,
                        block.reference_m2,
                        size(block.matrix),
                    )
                    for block in orbit_k1_blocks
                ] == [
                    (0, 1, 1, 0, (1, 1)),
                    (2, 2, 3, -2, (2, 2)),
                    (4, 1, 5, -4, (1, 1)),
                ]
                @test [block.label for block in orbit_k1_blocks] == [
                    (feature=:moment_matrix, decomposition=:translation_su2, momentum=1, spin2=0),
                    (feature=:moment_matrix, decomposition=:translation_su2, momentum=1, spin2=2),
                    (feature=:moment_matrix, decomposition=:translation_su2, momentum=1, spin2=4),
                ]
                @test [
                    (block.spin2, block.multiplicity)
                    for block in orbit_k1_blocks
                ] == [
                    (block.spin2, block.multiplicity)
                    for block in translation_su2_target6.nonzero_momentum_blocks
                ]
                @test all(
                    block -> all(entry -> entry isa Polynomial{PauliAlgebra}, block.matrix),
                    orbit_k1_blocks,
                )
                @test first(orbit_k1_blocks[1].row_labels).feature ==
                    :pauli_su2_translation_orbit_basis_state
                @test all(
                    block -> block.coefficient_domain == :complex_algebraic_float64,
                    orbit_k1_blocks,
                )
                translation_su2_reflection_k0_target6 =
                    pauli_su2_translation_orbit_structural_targets(
                        6,
                        2;
                        reflection_symmetry=true,
                        momenta=[0],
                    )
                orbit_reflection_k0_blocks =
                    pauli_su2_translation_orbit_momentum_blocks(
                        orbit_reps6,
                        6,
                        0;
                        reflection_symmetry=true,
                    )
                @test [block.label for block in orbit_reflection_k0_blocks] ==
                    translation_su2_reflection_k0_target6.block_labels
                @test [
                    block.multiplicity for block in orbit_reflection_k0_blocks
                ] == translation_su2_reflection_k0_target6.logical_block_sizes
                @test all(
                    block -> block.label.decomposition == :translation_su2_reflection,
                    orbit_reflection_k0_blocks,
                )
                @test all(
                    block -> block.label.reflection_fixed,
                    orbit_reflection_k0_blocks,
                )
            end
            @test isdefined(
                NCTSSoS,
                :pauli_su2_translation_orbit_real_reflection_momentum_blocks,
            )
            if isdefined(
                NCTSSoS,
                :pauli_su2_translation_orbit_real_reflection_momentum_blocks,
            )
                translation_su2_real_reflection_k1_target6 =
                    pauli_su2_translation_orbit_structural_targets(
                        6,
                        2;
                        reflection_symmetry=true,
                        momenta=[1],
                    )
                orbit_real_reflection_k1_blocks =
                    pauli_su2_translation_orbit_real_reflection_momentum_blocks(
                        orbit_reps6,
                        6,
                        1,
                    )
                @test [
                    block.label for block in orbit_real_reflection_k1_blocks
                ] == translation_su2_real_reflection_k1_target6.block_labels
                @test [
                    block.multiplicity for block in orbit_real_reflection_k1_blocks
                ] == translation_su2_real_reflection_k1_target6.logical_block_sizes
                @test [
                    size(block.matrix, 1) for block in orbit_real_reflection_k1_blocks
                ] == translation_su2_real_reflection_k1_target6.psd_block_sizes
                @test all(
                    block -> block.label.decomposition == :translation_su2_reflection,
                    orbit_real_reflection_k1_blocks,
                )
                @test all(
                    block -> !block.label.reflection_fixed,
                    orbit_real_reflection_k1_blocks,
                )
                @test all(
                    block -> all(entry -> entry isa Polynomial{PauliAlgebra}, block.matrix),
                    orbit_real_reflection_k1_blocks,
                )
                @test_throws ArgumentError pauli_su2_translation_orbit_real_reflection_momentum_blocks(
                    orbit_reps6,
                    6,
                    0,
                )
                registry10, ops10 = create_pauli_variables(1:10)
                orbit_reps10_order3 =
                    NCTSSoS._pauli_contiguous_chain_orbit_representatives(
                        ops10,
                        3;
                        periodic=true,
                    )
                translation_su2_reflection_khalf_target10 =
                    pauli_su2_translation_orbit_structural_targets(
                        10,
                        3;
                        reflection_symmetry=true,
                        momenta=[5],
                    )
                orbit_reflection_khalf_blocks10 =
                    pauli_su2_translation_orbit_momentum_blocks(
                        orbit_reps10_order3,
                        10,
                        5;
                        reflection_symmetry=true,
                    )
                @test [block.label for block in orbit_reflection_khalf_blocks10] ==
                    translation_su2_reflection_khalf_target10.block_labels
                @test [
                    block.multiplicity for block in orbit_reflection_khalf_blocks10
                ] == translation_su2_reflection_khalf_target10.logical_block_sizes
            end
            @test isdefined(NCTSSoS, :pauli_su2_translation_orbit_momentum_block_bundle)
            if isdefined(NCTSSoS, :pauli_su2_translation_orbit_momentum_block_bundle)
                orbit_bundle6 = pauli_su2_translation_orbit_momentum_block_bundle(
                    orbit_reps6,
                    6,
                )
                @test orbit_bundle6.momentum_sectors ==
                    translation_su2_target6.momentum_sectors
                @test orbit_bundle6.n_blocks == translation_su2_target6.n_blocks
                @test orbit_bundle6.block_labels == translation_su2_target6.block_labels
                @test orbit_bundle6.logical_block_sizes ==
                    translation_su2_target6.logical_block_sizes
                @test orbit_bundle6.logical_block_feature_histogram ==
                    translation_su2_target6.logical_block_feature_histogram
                @test all(
                    block -> all(entry -> entry isa Polynomial{PauliAlgebra}, block.matrix),
                    orbit_bundle6.blocks,
                )
                orbit_reflection_bundle6 = pauli_su2_translation_orbit_momentum_block_bundle(
                    orbit_reps6,
                    6;
                    momenta=[0],
                    reflection_symmetry=true,
                )
                @test orbit_reflection_bundle6.block_labels ==
                    translation_su2_reflection_k0_target6.block_labels
                @test orbit_reflection_bundle6.logical_block_sizes ==
                    translation_su2_reflection_k0_target6.logical_block_sizes
                @test orbit_reflection_bundle6.logical_block_feature_histogram ==
                    translation_su2_reflection_k0_target6.logical_block_feature_histogram
                orbit_real_reflection_bundle6 =
                    pauli_su2_translation_orbit_momentum_block_bundle(
                        orbit_reps6,
                        6;
                        momenta=[1],
                        reflection_symmetry=true,
                    )
                @test orbit_real_reflection_bundle6.block_labels ==
                    translation_su2_real_reflection_k1_target6.block_labels
                @test orbit_real_reflection_bundle6.logical_block_sizes ==
                    translation_su2_real_reflection_k1_target6.logical_block_sizes
                @test [
                    size(block.matrix, 1) for block in orbit_real_reflection_bundle6.blocks
                ] == translation_su2_real_reflection_k1_target6.psd_block_sizes
            end
            @test isdefined(NCTSSoS, :pauli_su2_translation_orbit_moment_problem)
            if isdefined(NCTSSoS, :pauli_su2_translation_orbit_moment_problem)
                orbit_auto_objective = heisenberg_chain_hamiltonian(ops6)
                @test_throws ArgumentError pauli_su2_translation_orbit_moment_problem(
                    orbit_auto_objective,
                    orbit_reps6,
                    6,
                )
                orbit_auto_mp = pauli_su2_translation_orbit_moment_problem(
                    orbit_auto_objective,
                    orbit_reps6,
                    6;
                    assume_su2_invariant=true,
                )
                @test [cone for (cone, _) in orbit_auto_mp.constraints] ==
                    fill(:HPSD, translation_su2_target6.n_blocks)
                @test [size(mat, 1) for (_, mat) in orbit_auto_mp.constraints] ==
                    translation_su2_target6.logical_block_sizes
                @test [block.size for block in orbit_auto_mp.linear.psd_blocks_lin] ==
                    translation_su2_target6.logical_block_sizes
                @test [
                    block.meta.origin.label
                    for block in orbit_auto_mp.linear.psd_blocks_lin
                ] == translation_su2_target6.block_labels
                @test first(
                    orbit_auto_mp.linear.psd_blocks_lin[1].meta.origin.logical_row_labels,
                ).feature == :pauli_su2_translation_orbit_basis_state
                @test all(
                    block -> block.meta.origin.transform.coefficient_domain ==
                        :complex_algebraic_float64,
                    orbit_auto_mp.linear.psd_blocks_lin,
                )
                orbit_auto_model, _ = build_jump_model(orbit_auto_mp)
                @test !isnothing(orbit_auto_model)
                orbit_reflection_mp = pauli_su2_translation_orbit_moment_problem(
                    orbit_auto_objective,
                    orbit_reps6,
                    6;
                    assume_su2_invariant=true,
                    real_moment_matrix=false,
                    momenta=[0],
                    reflection_symmetry=true,
                )
                @test [
                    cone for (cone, _) in
                        orbit_reflection_mp.constraints[
                            1:translation_su2_reflection_k0_target6.n_blocks
                        ]
                ] ==
                    fill(:HPSD, translation_su2_reflection_k0_target6.n_blocks)
                @test [
                    size(mat, 1) for (_, mat) in
                        orbit_reflection_mp.constraints[
                            1:translation_su2_reflection_k0_target6.n_blocks
                        ]
                ] ==
                    translation_su2_reflection_k0_target6.logical_block_sizes
                @test count(con -> con[1] == :Zero, orbit_reflection_mp.constraints) ==
                    length(orbit_reflection_mp.linear.zero_constraints)
                @test [
                    block.meta.origin.label
                    for block in orbit_reflection_mp.linear.psd_blocks_lin
                ] == translation_su2_reflection_k0_target6.block_labels
                @test all(
                    block -> block.meta.origin.transform.family ==
                        :pauli_translation_su2_reflection,
                    orbit_reflection_mp.linear.psd_blocks_lin,
                )
                orbit_reflection_all_target6 =
                    pauli_su2_translation_orbit_structural_targets(
                        6,
                        2;
                        reflection_symmetry=true,
                    )
                orbit_reflection_all_mp = pauli_su2_translation_orbit_moment_problem(
                    orbit_auto_objective,
                    orbit_reps6,
                    6;
                    assume_su2_invariant=true,
                    reflection_symmetry=true,
                )
                @test [
                    cone for (cone, _) in
                        orbit_reflection_all_mp.constraints[
                            1:orbit_reflection_all_target6.n_blocks
                        ]
                ] ==
                    fill(:PSD, orbit_reflection_all_target6.n_blocks)
                @test [
                    size(mat, 1) for (_, mat) in
                        orbit_reflection_all_mp.constraints[
                            1:orbit_reflection_all_target6.n_blocks
                        ]
                ] ==
                    orbit_reflection_all_target6.psd_block_sizes
                @test count(con -> con[1] == :Zero, orbit_reflection_all_mp.constraints) ==
                    length(orbit_reflection_all_mp.linear.zero_constraints)
                @test [
                    block.meta.origin.label
                    for block in orbit_reflection_all_mp.linear.psd_blocks_lin
                ] == orbit_reflection_all_target6.block_labels
            end
            translation_su2_complex_target6 =
                pauli_su2_translation_orbit_structural_targets(
                    6,
                    2;
                    real_moment_matrix=false,
                )
            wigner_target6 = pauli_su2_translation_orbit_structural_targets(
                6,
                2;
                real_moment_matrix=false,
                momenta=[0],
            )
            wigner_mp6 = pauli_su2_translation_orbit_wigner_eckart_moment_problem(
                orbit_auto_objective,
                orbit_reps6,
                6,
                0;
                assume_su2_invariant=true,
            )
            @test [
                block.meta.origin.label for block in wigner_mp6.linear.psd_blocks_lin
            ] == wigner_target6.block_labels
            @test [block.size for block in wigner_mp6.linear.psd_blocks_lin] ==
                wigner_target6.logical_block_sizes
            @test count(con -> con[1] == :HPSD, wigner_mp6.constraints) ==
                wigner_target6.n_blocks
            @test count(con -> con[1] == :Zero, wigner_mp6.constraints) ==
                length(wigner_mp6.linear.zero_constraints)
            wigner_zero_labels = [zc.origin.label for zc in wigner_mp6.linear.zero_constraints]
            @test any(label -> label.reason == :spin_offblock, wigner_zero_labels)
            @test any(label -> label.reason == :magnetic_offdiag, wigner_zero_labels)
            @test any(label -> label.reason == :magnetic_copy, wigner_zero_labels)
            @test all(label -> label.decomposition == :translation_su2, wigner_zero_labels)
            @test all(label -> label.component == :complex, wigner_zero_labels)
            @test NCTSSoS.assert_moment_linear_data_invariants(
                wigner_mp6.linear,
                wigner_mp6.constraints,
            ) === nothing
            wigner_all_mp6 = pauli_su2_translation_orbit_wigner_eckart_moment_problem(
                orbit_auto_objective,
                orbit_reps6,
                6;
                assume_su2_invariant=true,
            )
            @test [
                block.meta.origin.label for block in wigner_all_mp6.linear.psd_blocks_lin
            ] == translation_su2_complex_target6.block_labels
            @test [block.size for block in wigner_all_mp6.linear.psd_blocks_lin] ==
                translation_su2_complex_target6.logical_block_sizes
            @test count(con -> con[1] == :HPSD, wigner_all_mp6.constraints) ==
                translation_su2_complex_target6.n_blocks
            @test count(con -> con[1] == :Zero, wigner_all_mp6.constraints) ==
                length(wigner_all_mp6.linear.zero_constraints)
            wigner_all_zero_labels = [
                zc.origin.label for zc in wigner_all_mp6.linear.zero_constraints
            ]
            @test any(label -> label.momentum == 0, wigner_all_zero_labels)
            @test any(label -> label.momentum == 5, wigner_all_zero_labels)
            @test any(label -> label.reason == :spin_offblock, wigner_all_zero_labels)
            @test any(label -> label.reason == :magnetic_offdiag, wigner_all_zero_labels)
            @test any(label -> label.reason == :magnetic_copy, wigner_all_zero_labels)
            @test NCTSSoS.assert_moment_linear_data_invariants(
                wigner_all_mp6.linear,
                wigner_all_mp6.constraints,
            ) === nothing
            public_su2_pop6 = polyopt(heisenberg_chain_hamiltonian(ops6), registry6)
            public_su2_mp6, public_su2_report6 =
                pauli_translation_invariant_moment_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    real_moment_matrix=false,
                )
            @test public_su2_report6.block_labels ==
                translation_su2_complex_target6.block_labels
            @test public_su2_report6.logical_block_sizes ==
                translation_su2_complex_target6.logical_block_sizes
            @test public_su2_report6.psd_block_sizes ==
                translation_su2_complex_target6.psd_block_sizes
            @test count(con -> con[1] == :HPSD, public_su2_mp6.constraints) ==
                translation_su2_complex_target6.n_blocks
            @test count(con -> con[1] == :Zero, public_su2_mp6.constraints) ==
                public_su2_report6.zero_constraint_count
            @test public_su2_report6.zero_constraint_count ==
                length(public_su2_mp6.linear.zero_constraints)
            @test [
                block.meta.origin.label for block in public_su2_mp6.linear.psd_blocks_lin
            ] == translation_su2_complex_target6.block_labels
            public_su2_zero_labels = [
                zc.origin.label for zc in public_su2_mp6.linear.zero_constraints
            ]
            @test any(label -> label.reason == :spin_offblock, public_su2_zero_labels)
            @test any(label -> label.reason == :magnetic_offdiag, public_su2_zero_labels)
            @test any(label -> label.reason == :magnetic_copy, public_su2_zero_labels)
            @test all(
                label -> label.decomposition == :translation_su2,
                public_su2_zero_labels,
            )
            @test all(label -> label.component == :complex, public_su2_zero_labels)
            public_su2_singlet_mp6, public_su2_singlet_report6 =
                pauli_translation_invariant_moment_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    real_moment_matrix=false,
                    singlet_channel_equalities=true,
                )
            @test public_su2_singlet_report6.block_labels ==
                public_su2_report6.block_labels
            @test public_su2_singlet_report6.zero_constraint_count >
                public_su2_report6.zero_constraint_count
            @test public_su2_singlet_report6.zero_constraint_count ==
                length(public_su2_singlet_mp6.linear.zero_constraints)
            public_su2_singlet_zero_labels = [
                zc.origin.label
                for zc in public_su2_singlet_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label.feature ==
                    :pauli_su2_translation_orbit_singlet_channel_equality,
                public_su2_singlet_zero_labels,
            )
            public_su2_singlet_metrics6 =
                translation_report_metrics(public_su2_singlet_report6)
            @test public_su2_singlet_metrics6.su2_singlet_channel_equality_row_count ==
                count(
                    label -> label.feature ==
                        :pauli_su2_translation_orbit_singlet_channel_equality,
                    public_su2_singlet_zero_labels,
                )
            public_su2_singlet_sos6 = NCTSSoS.sos_dualize(public_su2_singlet_mp6)
            public_su2_singlet_zero_duals6 =
                NCTSSoS.sos_zero_duals(public_su2_singlet_mp6, public_su2_singlet_sos6)
            @test length(public_su2_singlet_zero_duals6) ==
                length(public_su2_singlet_mp6.linear.zero_constraints)
            @test count(
                dual -> dual.origin.label.feature ==
                    :pauli_su2_translation_orbit_singlet_channel_equality,
                public_su2_singlet_zero_duals6,
            ) ==
                public_su2_singlet_metrics6.su2_singlet_channel_equality_row_count
            public_su2_singlet_zero_dual_rows6 = filter(
                dual -> dual.origin isa NCTSSoS.PauliSU2SingletChannelEqualityOrigin,
                public_su2_singlet_zero_duals6,
            )
            @test !isempty(public_su2_singlet_zero_dual_rows6)
            @test all(
                dual -> dual.feature ==
                    :pauli_su2_translation_orbit_singlet_channel_equality,
                public_su2_singlet_zero_dual_rows6,
            )
            @test all(
                dual -> dual.decomposition == :translation_su2,
                public_su2_singlet_zero_dual_rows6,
            )
            @test all(dual -> dual.reason === nothing, public_su2_singlet_zero_dual_rows6)
            @test all(dual -> dual.term_count > 0, public_su2_singlet_zero_dual_rows6)
            @test all(
                dual -> dual.coefficient_domain == :complex_algebraic_float64,
                public_su2_singlet_zero_dual_rows6,
            )
            @test all(
                dual -> dual.exact_coefficient_domain == :complex_sqrt_rational,
                public_su2_singlet_zero_dual_rows6,
            )
            @test translation_solve_support(public_su2_singlet_report6).supported
            @test NCTSSoS.assert_moment_linear_data_invariants(
                public_su2_singlet_mp6.linear,
                public_su2_singlet_mp6.constraints,
            ) === nothing
            public_su2_linear6, public_su2_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    real_moment_matrix=false,
                )
            @test public_su2_linear_report6.block_labels ==
                translation_su2_complex_target6.block_labels
            @test public_su2_linear_report6.psd_block_sizes ==
                translation_su2_complex_target6.psd_block_sizes
            @test [
                block.meta.origin.label for block in public_su2_linear6.psd_blocks_lin
            ] == translation_su2_complex_target6.block_labels
            @test public_su2_linear_report6.zero_constraint_count ==
                public_su2_report6.zero_constraint_count
            @test [
                zc.origin.label for zc in public_su2_linear6.zero_constraints
            ] == public_su2_zero_labels
            public_su2_real_mp6, public_su2_real_report6 =
                pauli_translation_invariant_moment_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                )
            @test public_su2_real_report6.real_moment_matrix
            @test public_su2_real_report6.block_labels == translation_su2_target6.block_labels
            @test public_su2_real_report6.logical_block_sizes ==
                translation_su2_target6.logical_block_sizes
            @test public_su2_real_report6.psd_block_sizes ==
                translation_su2_target6.psd_block_sizes
            @test count(con -> con[1] == :PSD, public_su2_real_mp6.constraints) ==
                translation_su2_target6.n_blocks
            @test count(con -> con[1] == :Zero, public_su2_real_mp6.constraints) ==
                public_su2_real_report6.zero_constraint_count
            @test public_su2_real_report6.zero_constraint_count ==
                length(public_su2_real_mp6.linear.zero_constraints)
            @test [
                block.meta.origin.label
                for block in public_su2_real_mp6.linear.psd_blocks_lin
            ] == translation_su2_target6.block_labels
            @test [
                length(block.meta.origin.logical_row_labels)
                for block in public_su2_real_mp6.linear.psd_blocks_lin
            ] == translation_su2_target6.logical_block_sizes
            public_su2_real_zero_labels = [
                zc.origin.label for zc in public_su2_real_mp6.linear.zero_constraints
            ]
            @test any(label -> label.reason == :spin_offblock, public_su2_real_zero_labels)
            @test any(
                label -> label.reason == :magnetic_offdiag,
                public_su2_real_zero_labels,
            )
            @test any(label -> label.reason == :magnetic_copy, public_su2_real_zero_labels)
            @test all(
                label -> label.decomposition == :translation_su2,
                public_su2_real_zero_labels,
            )
            @test all(
                label -> label.component in (:real, :imag),
                public_su2_real_zero_labels,
            )
            public_su2_real_model6, _ = build_jump_model(public_su2_real_mp6)
            @test !isnothing(public_su2_real_model6)
            public_su2_profile6 = translation_symmetry_profile(public_su2_real_report6)
            @test public_su2_profile6.su2_moment_symmetry
            @test !public_su2_profile6.su2_reflection_moment_symmetry
            public_su2_support6 = translation_solve_support(public_su2_real_report6)
            @test public_su2_support6.supported
            @test public_su2_support6.blocker === nothing
            @test public_su2_support6.reason == ""
            public_su2_metrics6 = translation_report_metrics(public_su2_real_report6)
            @test public_su2_metrics6.solve_supported
            @test public_su2_metrics6.solve_blocker === nothing
            @test public_su2_metrics6.solve_blocker_reason == ""
            public_su2_real_linear6, public_su2_real_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                )
            @test public_su2_real_linear_report6.real_moment_matrix
            @test public_su2_real_linear_report6.block_labels ==
                translation_su2_target6.block_labels
            @test public_su2_real_linear_report6.psd_block_sizes ==
                translation_su2_target6.psd_block_sizes
            @test [
                block.meta.origin.label
                for block in public_su2_real_linear6.psd_blocks_lin
            ] == translation_su2_target6.block_labels
            @test public_su2_real_linear_report6.zero_constraint_count ==
                public_su2_real_report6.zero_constraint_count
            @test [
                zc.origin.label for zc in public_su2_real_linear6.zero_constraints
            ] == public_su2_real_zero_labels
            @test translation_solve_support(public_su2_real_linear_report6).supported
            public_su2_reflection_fixed_target6 =
                pauli_su2_translation_orbit_structural_targets(
                    6,
                    2;
                    reflection_symmetry=true,
                    momenta=[0],
                )
            public_su2_reflection_fixed_mp6, public_su2_reflection_fixed_report6 =
                pauli_translation_invariant_moment_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    reflection_symmetry=true,
                    momenta=[0],
                )
            @test public_su2_reflection_fixed_report6.block_labels ==
                public_su2_reflection_fixed_target6.block_labels
            @test public_su2_reflection_fixed_report6.logical_block_sizes ==
                public_su2_reflection_fixed_target6.logical_block_sizes
            @test public_su2_reflection_fixed_report6.psd_block_sizes ==
                public_su2_reflection_fixed_target6.psd_block_sizes
            @test [
                cone for (cone, _) in
                    public_su2_reflection_fixed_mp6.constraints[
                        1:public_su2_reflection_fixed_target6.n_blocks
                    ]
            ] ==
                fill(:PSD, public_su2_reflection_fixed_target6.n_blocks)
            @test public_su2_reflection_fixed_report6.zero_constraint_count > 0
            @test count(
                con -> con[1] == :Zero,
                public_su2_reflection_fixed_mp6.constraints,
            ) == public_su2_reflection_fixed_report6.zero_constraint_count
            public_su2_reflection_fixed_zero_labels = [
                zc.origin.label
                for zc in public_su2_reflection_fixed_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label.reason == :spin_offblock,
                public_su2_reflection_fixed_zero_labels,
            )
            @test any(
                label -> label.reason == :magnetic_offdiag,
                public_su2_reflection_fixed_zero_labels,
            )
            @test any(
                label -> label.reason == :magnetic_copy,
                public_su2_reflection_fixed_zero_labels,
            )
            public_su2_reflection_fixed_profile6 =
                translation_symmetry_profile(public_su2_reflection_fixed_report6)
            @test public_su2_reflection_fixed_profile6.reflection_symmetry
            @test public_su2_reflection_fixed_profile6.su2_reflection_moment_symmetry
            @test translation_solve_support(
                public_su2_reflection_fixed_report6,
            ).supported
            public_su2_reflection_fixed_singlet_mp6,
            public_su2_reflection_fixed_singlet_report6 =
                pauli_translation_invariant_moment_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    reflection_symmetry=true,
                    momenta=[0],
                    singlet_channel_equalities=true,
                )
            @test public_su2_reflection_fixed_singlet_report6.block_labels ==
                public_su2_reflection_fixed_report6.block_labels
            @test public_su2_reflection_fixed_singlet_report6.zero_constraint_count >
                public_su2_reflection_fixed_report6.zero_constraint_count
            public_su2_reflection_fixed_singlet_zero_labels = [
                zc.origin.label
                for zc in
                    public_su2_reflection_fixed_singlet_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label.feature ==
                    :pauli_su2_translation_orbit_singlet_channel_equality,
                public_su2_reflection_fixed_singlet_zero_labels,
            )
            public_su2_reflection_fixed_singlet_metrics6 =
                translation_report_metrics(public_su2_reflection_fixed_singlet_report6)
            @test public_su2_reflection_fixed_singlet_metrics6.su2_singlet_channel_equality_row_count ==
                count(
                    label -> label.feature ==
                        :pauli_su2_translation_orbit_singlet_channel_equality,
                    public_su2_reflection_fixed_singlet_zero_labels,
                )
            @test translation_solve_support(
                public_su2_reflection_fixed_singlet_report6,
            ).supported
            @test NCTSSoS.assert_moment_linear_data_invariants(
                public_su2_reflection_fixed_singlet_mp6.linear,
                public_su2_reflection_fixed_singlet_mp6.constraints,
            ) === nothing
            public_su2_reflection_fixed_linear6,
            public_su2_reflection_fixed_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    reflection_symmetry=true,
                    momenta=[0],
                )
            @test public_su2_reflection_fixed_linear_report6.block_labels ==
                public_su2_reflection_fixed_target6.block_labels
            @test public_su2_reflection_fixed_linear_report6.psd_block_sizes ==
                public_su2_reflection_fixed_target6.psd_block_sizes
            @test [
                block.meta.origin.label
                for block in public_su2_reflection_fixed_linear6.psd_blocks_lin
            ] == public_su2_reflection_fixed_target6.block_labels
            @test public_su2_reflection_fixed_linear_report6.zero_constraint_count ==
                public_su2_reflection_fixed_report6.zero_constraint_count
            @test [
                zc.origin.label
                for zc in public_su2_reflection_fixed_linear6.zero_constraints
            ] == public_su2_reflection_fixed_zero_labels
            public_su2_reflection_fixed_singlet_linear6,
            public_su2_reflection_fixed_singlet_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    reflection_symmetry=true,
                    momenta=[0],
                    singlet_channel_equalities=true,
                )
            @test public_su2_reflection_fixed_singlet_linear_report6.block_labels ==
                public_su2_reflection_fixed_linear_report6.block_labels
            @test public_su2_reflection_fixed_singlet_linear_report6.zero_constraint_count >
                public_su2_reflection_fixed_linear_report6.zero_constraint_count
            @test any(
                zc -> zc.origin.label.feature ==
                    :pauli_su2_translation_orbit_singlet_channel_equality,
                public_su2_reflection_fixed_singlet_linear6.zero_constraints,
            )
            @test any(
                zc -> zc.origin isa NCTSSoS.PauliSU2SingletChannelEqualityOrigin &&
                    zc.origin.label.feature ==
                        :pauli_su2_translation_orbit_singlet_channel_equality,
                public_su2_reflection_fixed_singlet_linear6.zero_constraints,
            )
            @test all(
                zc -> hasproperty(zc.origin.label, :coefficient_domain) &&
                    hasproperty(zc.origin.label, :exact_coefficient_domain) &&
                    zc.origin.label.coefficient_domain == :complex_algebraic_float64 &&
                    zc.origin.label.exact_coefficient_domain == :complex_sqrt_rational,
                filter(
                    zc -> zc.origin isa NCTSSoS.PauliSU2SingletChannelEqualityOrigin,
                    public_su2_reflection_fixed_singlet_linear6.zero_constraints,
                ),
            )
            @test translation_solve_support(
                public_su2_reflection_fixed_singlet_linear_report6,
            ).supported
            @test NCTSSoS.assert_moment_linear_data_invariants(
                public_su2_reflection_fixed_singlet_linear6,
            ) === nothing
            public_su2_reflection_target6 =
                pauli_su2_translation_orbit_structural_targets(
                    6,
                    2;
                    reflection_symmetry=true,
                )
            public_su2_reflection_mp6, public_su2_reflection_report6 =
                pauli_translation_invariant_moment_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    reflection_symmetry=true,
                )
            @test public_su2_reflection_report6.block_labels ==
                public_su2_reflection_target6.block_labels
            @test public_su2_reflection_report6.logical_block_sizes ==
                public_su2_reflection_target6.logical_block_sizes
            @test public_su2_reflection_report6.psd_block_sizes ==
                public_su2_reflection_target6.psd_block_sizes
            @test [
                cone for (cone, _) in
                    public_su2_reflection_mp6.constraints[
                        1:public_su2_reflection_target6.n_blocks
                    ]
            ] ==
                fill(:PSD, public_su2_reflection_target6.n_blocks)
            @test public_su2_reflection_report6.zero_constraint_count > 0
            @test count(con -> con[1] == :Zero, public_su2_reflection_mp6.constraints) ==
                public_su2_reflection_report6.zero_constraint_count
            public_su2_reflection_zero_labels = [
                zc.origin.label for zc in public_su2_reflection_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label.reason == :spin_offblock,
                public_su2_reflection_zero_labels,
            )
            @test any(
                label -> label.reason == :magnetic_offdiag,
                public_su2_reflection_zero_labels,
            )
            @test any(
                label -> label.reason == :magnetic_copy,
                public_su2_reflection_zero_labels,
            )
            public_su2_reflection_profile6 =
                translation_symmetry_profile(public_su2_reflection_report6)
            @test public_su2_reflection_profile6.reflection_symmetry
            @test public_su2_reflection_profile6.su2_reflection_moment_symmetry
            public_su2_reflection_linear6, public_su2_reflection_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                    reflection_symmetry=true,
                )
            @test public_su2_reflection_linear_report6.block_labels ==
                public_su2_reflection_target6.block_labels
            @test public_su2_reflection_linear_report6.psd_block_sizes ==
                public_su2_reflection_target6.psd_block_sizes
            @test [
                block.meta.origin.label
                for block in public_su2_reflection_linear6.psd_blocks_lin
            ] == public_su2_reflection_target6.block_labels
            @test [
                block.size for block in public_su2_reflection_linear6.psd_blocks_lin
            ] == public_su2_reflection_target6.psd_block_sizes
            @test public_su2_reflection_linear_report6.zero_constraint_count ==
                public_su2_reflection_report6.zero_constraint_count
            @test [
                zc.origin.label for zc in public_su2_reflection_linear6.zero_constraints
            ] == public_su2_reflection_zero_labels
            @test translation_solve_support(public_su2_reflection_linear_report6).supported
            registry4_reflected, ops4_reflected = create_pauli_variables(1:4)
            h4_reflected = heisenberg_chain_hamiltonian(ops4_reflected)
            pop4_reflected = polyopt(h4_reflected, registry4_reflected)
            ordinary_mp4, _ = pauli_translation_invariant_moment_relaxation(
                pop4_reflected,
                ops4_reflected,
                2;
                sign_symmetry=false,
            )
            ordinary_reflection_mp4, _ = pauli_translation_invariant_moment_relaxation(
                pop4_reflected,
                ops4_reflected,
                2;
                sign_symmetry=false,
                reflection_symmetry=true,
            )
            ordinary_result4 = NCTSSoS.solve_moment_problem(
                ordinary_mp4,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false),
            )
            ordinary_reflection_result4 = NCTSSoS.solve_moment_problem(
                ordinary_reflection_mp4,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false),
            )
            su2_wigner_result4 = pauli_translation_invariant_nctssos(
                pop4_reflected,
                ops4_reflected,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
            )
            su2_reflection_result4 = pauli_translation_invariant_nctssos(
                pop4_reflected,
                ops4_reflected,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                reflection_symmetry=true,
            )
            su2_wigner_linear_result4 = pauli_translation_invariant_nctssos(
                pop4_reflected,
                ops4_reflected,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                direct_linear=true,
            )
            su2_reflection_linear_result4 = pauli_translation_invariant_nctssos(
                pop4_reflected,
                ops4_reflected,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                reflection_symmetry=true,
                direct_linear=true,
            )
            ordinary_axis_result4 = pauli_translation_invariant_nctssos(
                pop4_reflected,
                ops4_reflected,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                axis_rotation_equalities=true,
            )
            su2_axis_result4 = pauli_translation_invariant_nctssos(
                pop4_reflected,
                ops4_reflected,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                axis_rotation_equalities=true,
            )
            su2_axis_linear_result4 = pauli_translation_invariant_nctssos(
                pop4_reflected,
                ops4_reflected,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                axis_rotation_equalities=true,
                direct_linear=true,
            )
            su2_singlet_config4 = SolverConfig(
                ;
                optimizer=optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false),
                order=2,
                symmetry=heisenberg_chain_symmetry_spec(
                    ops4_reflected;
                    reflection=false,
                    axis_rotations=false,
                    check_invariance=false,
                ),
            )
            su2_singlet_linear_auto_result4 = quiet() do
                cs_nctssos(
                    pop4_reflected,
                    su2_singlet_config4;
                    dualize=false,
                    direct_linear=true,
                    su2_symmetry=true,
                    singlet_channel_equalities=true,
                )
            end
            @test termination_status(ordinary_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(ordinary_reflection_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(su2_wigner_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(su2_reflection_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(su2_wigner_linear_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(su2_reflection_linear_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(ordinary_axis_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(su2_axis_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(su2_axis_linear_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test termination_status(su2_singlet_linear_auto_result4.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test su2_wigner_result4.report.zero_constraint_count > 0
            @test su2_reflection_result4.report.zero_constraint_count > 0
            @test translation_solve_support(su2_wigner_result4).supported
            @test translation_solve_support(su2_reflection_result4).supported
            @test translation_solve_support(su2_wigner_linear_result4).supported
            @test translation_solve_support(su2_reflection_linear_result4).supported
            @test translation_solve_support(su2_axis_result4).supported
            @test translation_solve_support(su2_axis_linear_result4).supported
            @test translation_solve_support(su2_singlet_linear_auto_result4).supported
            @test su2_singlet_linear_auto_result4.moment_problem isa
                NCTSSoS.MomentLinearData
            @test su2_wigner_result4.objective ≈ ordinary_result4.objective atol = 1e-6
            @test su2_reflection_result4.objective ≈
                ordinary_reflection_result4.objective atol = 1e-6
            @test su2_wigner_linear_result4.objective ≈
                ordinary_result4.objective atol = 1e-6
            @test su2_reflection_linear_result4.objective ≈
                ordinary_reflection_result4.objective atol = 1e-6
            @test su2_axis_result4.objective ≈ ordinary_axis_result4.objective atol = 1e-6
            @test su2_axis_linear_result4.objective ≈
                ordinary_axis_result4.objective atol = 1e-6
            @test su2_singlet_linear_auto_result4.objective ≈
                su2_wigner_linear_result4.objective atol = 1e-6
            su2_singlet_auto_metrics4 =
                translation_report_metrics(su2_singlet_linear_auto_result4)
            @test su2_singlet_auto_metrics4.su2_singlet_channel_equality_row_count > 0
            @test any(
                zc -> zc.origin isa NCTSSoS.PauliSU2SingletChannelEqualityOrigin &&
                    zc.origin.label.feature ==
                        :pauli_su2_translation_orbit_singlet_channel_equality,
                su2_singlet_linear_auto_result4.moment_problem.zero_constraints,
            )
            su2_axis_hist4 = translation_zero_origin_histogram(
                su2_axis_result4.moment_problem,
            )
            @test any(
                pair -> first(pair).feature == :axis_rotation_equality,
                su2_axis_hist4,
            )
            @test any(
                pair -> first(pair).feature == :pauli_su2_base_zero,
                su2_axis_hist4,
            )
            @test translation_zero_origin_histogram(
                su2_axis_linear_result4.moment_problem,
            ) == su2_axis_hist4
            for su2_result4 in (
                su2_wigner_result4,
                su2_wigner_linear_result4,
                su2_axis_result4,
                su2_axis_linear_result4,
            )
                su2_sos4 = NCTSSoS.sos_dualize(su2_result4.moment_problem)
                set_optimizer(su2_sos4.model, Clarabel.Optimizer)
                set_silent(su2_sos4.model)
                optimize!(su2_sos4.model)
                @test termination_status(su2_sos4.model) in
                    (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
                su2_residual4 = sos_dual_certificate_residual(
                    su2_result4.moment_problem,
                    su2_sos4,
                )
                @test su2_residual4.max_abs_residual <= 1e-6
                su2_diagnostics4 =
                    sos_dual_block_diagnostics(su2_result4.moment_problem, su2_sos4)
                @test [diag.label for diag in su2_diagnostics4] ==
                    su2_result4.report.block_labels
                @test [diag.native_value_size for diag in su2_diagnostics4] ==
                    [(size, size) for size in su2_result4.report.psd_block_sizes]
                @test all(
                    diag -> diag.native_hermitian_residual <= 1e-6,
                    su2_diagnostics4,
                )
            end
            scalar_cases4 = (
                (
                    name=:ineq,
                    pop=polyopt(
                        h4_reflected,
                        registry4_reflected;
                        ineq_constraints=[one(h4_reflected) + h4_reflected],
                    ),
                ),
                (
                    name=:eq,
                    pop=polyopt(
                        h4_reflected,
                        registry4_reflected;
                        eq_constraints=[one(h4_reflected) + h4_reflected],
                    ),
                ),
            )
            for scalar_case4 in scalar_cases4
                ordinary_scalar_mp4, _ = pauli_translation_invariant_moment_relaxation(
                    scalar_case4.pop,
                    ops4_reflected,
                    2;
                    sign_symmetry=false,
                )
                ordinary_scalar_result4 = NCTSSoS.solve_moment_problem(
                    ordinary_scalar_mp4,
                    optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false),
                )
                su2_scalar_result4 = pauli_translation_invariant_nctssos(
                    scalar_case4.pop,
                    ops4_reflected,
                    2,
                    optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                    dualize=false,
                    sign_symmetry=false,
                    su2_symmetry=true,
                )
                @test termination_status(ordinary_scalar_result4.model) in
                    (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
                @test termination_status(su2_scalar_result4.model) in
                    (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
                @test translation_solve_support(su2_scalar_result4).supported
                @test su2_scalar_result4.objective ≈
                    ordinary_scalar_result4.objective atol = 1e-6
            end
            public_su2_default_mp6, public_su2_default_report6 =
                pauli_translation_invariant_moment_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                )
            @test public_su2_default_report6.sign_symmetry
            @test public_su2_default_report6.block_labels ==
                translation_su2_target6.block_labels
            @test public_su2_default_report6.psd_block_sizes ==
                translation_su2_target6.psd_block_sizes
            @test [
                block.meta.origin.label
                for block in public_su2_default_mp6.linear.psd_blocks_lin
            ] == translation_su2_target6.block_labels
            public_su2_unsigned_mp6, public_su2_unsigned_report6 =
                pauli_translation_invariant_moment_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    sign_symmetry=false,
                    su2_symmetry=true,
                )
            @test !public_su2_unsigned_report6.sign_symmetry
            @test public_su2_unsigned_report6.block_labels ==
                public_su2_default_report6.block_labels
            @test public_su2_unsigned_report6.psd_block_sizes ==
                public_su2_default_report6.psd_block_sizes
            @test public_su2_unsigned_report6.zero_constraint_count ==
                public_su2_default_report6.zero_constraint_count
            @test [
                block.meta.origin.label
                for block in public_su2_unsigned_mp6.linear.psd_blocks_lin
            ] == translation_su2_target6.block_labels
            public_su2_default_linear6, public_su2_default_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    public_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                )
            @test public_su2_default_linear_report6.sign_symmetry
            @test public_su2_default_linear_report6.block_labels ==
                translation_su2_target6.block_labels
            @test public_su2_default_linear_report6.psd_block_sizes ==
                translation_su2_target6.psd_block_sizes
            @test [
                block.meta.origin.label
                for block in public_su2_default_linear6.psd_blocks_lin
            ] == translation_su2_target6.block_labels
            scalar_su2_objective6 = heisenberg_chain_hamiltonian(ops6)
            scalar_su2_ineq6 = one(scalar_su2_objective6) + scalar_su2_objective6
            scalar_su2_pop6 = polyopt(
                scalar_su2_objective6,
                registry6;
                eq_constraints=[scalar_su2_objective6],
                ineq_constraints=[scalar_su2_ineq6],
            )
            scalar_su2_mp6, scalar_su2_report6 =
                pauli_translation_invariant_moment_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                )
            @test scalar_su2_report6.block_labels[1:(end - 1)] ==
                translation_su2_target6.block_labels
            @test scalar_su2_report6.block_labels[end] ==
                (feature=:scalar_inequality, index=1)
            @test scalar_su2_report6.psd_block_sizes[1:(end - 1)] ==
                translation_su2_target6.psd_block_sizes
            @test scalar_su2_report6.psd_block_sizes[end] == 1
            @test scalar_su2_report6.zero_constraint_count > 1
            @test count(con -> con[1] == :PSD, scalar_su2_mp6.constraints) ==
                translation_su2_target6.n_blocks + 1
            @test count(con -> con[1] == :Zero, scalar_su2_mp6.constraints) ==
                scalar_su2_report6.zero_constraint_count
            scalar_su2_zero_labels = [
                zc.origin.label for zc in scalar_su2_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label == (feature=:scalar_equality, index=1),
                scalar_su2_zero_labels,
            )
            scalar_su2_wigner_labels = filter(
                label -> label isa NamedTuple && haskey(label, :reason),
                scalar_su2_zero_labels,
            )
            @test any(label -> label.reason == :spin_offblock, scalar_su2_wigner_labels)
            @test any(label -> label.reason == :magnetic_offdiag, scalar_su2_wigner_labels)
            @test any(label -> label.reason == :magnetic_copy, scalar_su2_wigner_labels)
            scalar_su2_support6 = translation_solve_support(scalar_su2_report6)
            @test scalar_su2_support6.supported
            @test scalar_su2_support6.blocker === nothing
            @test scalar_su2_support6.reason == ""
            scalar_su2_eq_constraint_idx6 = findfirst(
                con -> con[1] == :Zero &&
                    only(con[2]) ==
                        convert(
                            typeof(only(con[2])),
                            NCTSSoS._translation_orbit_reduce_polynomial(
                                scalar_su2_objective6,
                                6,
                            ),
                        ),
                scalar_su2_mp6.constraints,
            )
            @test scalar_su2_eq_constraint_idx6 !== nothing
            scalar_su2_ineq_constraint_idx6 = findfirst(
                con -> con[1] == :PSD &&
                    size(con[2]) == (1, 1) &&
                    only(con[2]) ==
                        convert(
                            typeof(only(con[2])),
                            NCTSSoS._translation_orbit_reduce_polynomial(
                                scalar_su2_ineq6,
                                6,
                            ),
                        ),
                scalar_su2_mp6.constraints,
            )
            @test scalar_su2_ineq_constraint_idx6 !== nothing
            scalar_su2_eq_entry6 =
                only(scalar_su2_mp6.constraints[scalar_su2_eq_constraint_idx6][2])
            scalar_su2_ineq_entry6 =
                only(scalar_su2_mp6.constraints[scalar_su2_ineq_constraint_idx6][2])
            @test scalar_su2_eq_entry6 ==
                convert(
                    typeof(scalar_su2_eq_entry6),
                    NCTSSoS._translation_orbit_reduce_polynomial(
                        scalar_su2_objective6,
                        6,
                    ),
                )
            @test scalar_su2_ineq_entry6 ==
                convert(
                    typeof(scalar_su2_ineq_entry6),
                    NCTSSoS._translation_orbit_reduce_polynomial(scalar_su2_ineq6, 6),
                )
            @test NCTSSoS.assert_moment_linear_data_invariants(
                scalar_su2_mp6.linear,
                scalar_su2_mp6.constraints,
            ) === nothing
            scalar_su2_linear6, scalar_su2_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                )
            @test scalar_su2_linear_report6.block_labels ==
                scalar_su2_report6.block_labels
            @test scalar_su2_linear_report6.psd_block_sizes ==
                scalar_su2_report6.psd_block_sizes
            @test scalar_su2_linear_report6.zero_constraint_count ==
                scalar_su2_report6.zero_constraint_count
            @test translation_solve_support(scalar_su2_linear_report6).supported
            @test [
                zc.origin.label for zc in scalar_su2_linear6.zero_constraints
            ] == scalar_su2_zero_labels
            scalar_su2_complex_mp6, scalar_su2_complex_report6 =
                pauli_translation_invariant_moment_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                    real_moment_matrix=false,
                )
            @test !scalar_su2_complex_report6.real_moment_matrix
            @test scalar_su2_complex_report6.block_labels[1:(end - 1)] ==
                translation_su2_complex_target6.block_labels
            @test scalar_su2_complex_report6.block_labels[end] ==
                (feature=:scalar_inequality, index=1)
            @test scalar_su2_complex_report6.psd_block_sizes[1:(end - 1)] ==
                translation_su2_complex_target6.psd_block_sizes
            @test scalar_su2_complex_report6.psd_block_sizes[end] == 1
            @test scalar_su2_complex_report6.zero_constraint_count > 1
            @test count(con -> con[1] == :HPSD, scalar_su2_complex_mp6.constraints) ==
                translation_su2_complex_target6.n_blocks + 1
            @test count(con -> con[1] == :Zero, scalar_su2_complex_mp6.constraints) ==
                scalar_su2_complex_report6.zero_constraint_count
            scalar_su2_complex_zero_labels = [
                zc.origin.label for zc in scalar_su2_complex_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label == (feature=:scalar_equality, index=1),
                scalar_su2_complex_zero_labels,
            )
            scalar_su2_complex_wigner_labels = filter(
                label -> label isa NamedTuple && haskey(label, :reason),
                scalar_su2_complex_zero_labels,
            )
            @test any(
                label -> label.reason == :spin_offblock,
                scalar_su2_complex_wigner_labels,
            )
            @test any(
                label -> label.reason == :magnetic_offdiag,
                scalar_su2_complex_wigner_labels,
            )
            @test any(
                label -> label.reason == :magnetic_copy,
                scalar_su2_complex_wigner_labels,
            )
            @test NCTSSoS.assert_moment_linear_data_invariants(
                scalar_su2_complex_mp6.linear,
                scalar_su2_complex_mp6.constraints,
            ) === nothing
            scalar_su2_complex_linear6, scalar_su2_complex_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                    real_moment_matrix=false,
                )
            @test scalar_su2_complex_linear_report6.block_labels ==
                scalar_su2_complex_report6.block_labels
            @test scalar_su2_complex_linear_report6.psd_block_sizes ==
                scalar_su2_complex_report6.psd_block_sizes
            @test scalar_su2_complex_linear_report6.zero_constraint_count ==
                scalar_su2_complex_report6.zero_constraint_count
            @test [
                zc.origin.label for zc in scalar_su2_complex_linear6.zero_constraints
            ] == scalar_su2_complex_zero_labels
            scalar_su2_reflection_fixed_mp6, scalar_su2_reflection_fixed_report6 =
                pauli_translation_invariant_moment_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                    reflection_symmetry=true,
                    momenta=[0],
                )
            @test scalar_su2_reflection_fixed_report6.block_labels[1:(end - 1)] ==
                public_su2_reflection_fixed_target6.block_labels
            @test scalar_su2_reflection_fixed_report6.block_labels[end] ==
                (feature=:scalar_inequality, index=1)
            @test scalar_su2_reflection_fixed_report6.psd_block_sizes[1:(end - 1)] ==
                public_su2_reflection_fixed_target6.psd_block_sizes
            @test scalar_su2_reflection_fixed_report6.psd_block_sizes[end] == 1
            @test scalar_su2_reflection_fixed_report6.zero_constraint_count > 1
            @test [
                cone for (cone, _) in
                    scalar_su2_reflection_fixed_mp6.constraints[
                        1:public_su2_reflection_fixed_target6.n_blocks
                    ]
            ] == fill(:PSD, public_su2_reflection_fixed_target6.n_blocks)
            @test scalar_su2_reflection_fixed_mp6.constraints[end][1] == :PSD
            scalar_su2_reflection_fixed_zero_labels = [
                zc.origin.label
                for zc in scalar_su2_reflection_fixed_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label == (feature=:scalar_equality, index=1),
                scalar_su2_reflection_fixed_zero_labels,
            )
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :reason) &&
                    label.reason == :spin_offblock,
                scalar_su2_reflection_fixed_zero_labels,
            )
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :reason) &&
                    label.reason == :magnetic_copy,
                scalar_su2_reflection_fixed_zero_labels,
            )
            scalar_su2_reflection_fixed_linear6,
            scalar_su2_reflection_fixed_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                    reflection_symmetry=true,
                    momenta=[0],
                )
            @test scalar_su2_reflection_fixed_linear_report6.block_labels ==
                scalar_su2_reflection_fixed_report6.block_labels
            @test scalar_su2_reflection_fixed_linear_report6.psd_block_sizes ==
                scalar_su2_reflection_fixed_report6.psd_block_sizes
            @test scalar_su2_reflection_fixed_linear_report6.zero_constraint_count ==
                scalar_su2_reflection_fixed_report6.zero_constraint_count
            @test [
                zc.origin.label
                for zc in scalar_su2_reflection_fixed_linear6.zero_constraints
            ] == scalar_su2_reflection_fixed_zero_labels
            scalar_su2_reflection_mp6, scalar_su2_reflection_report6 =
                pauli_translation_invariant_moment_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                    reflection_symmetry=true,
                )
            @test scalar_su2_reflection_report6.block_labels[1:(end - 1)] ==
                public_su2_reflection_target6.block_labels
            @test scalar_su2_reflection_report6.block_labels[end] ==
                (feature=:scalar_inequality, index=1)
            @test scalar_su2_reflection_report6.psd_block_sizes[1:(end - 1)] ==
                public_su2_reflection_target6.psd_block_sizes
            @test scalar_su2_reflection_report6.psd_block_sizes[end] == 1
            @test scalar_su2_reflection_report6.zero_constraint_count > 1
            @test [
                cone for (cone, _) in scalar_su2_reflection_mp6.constraints[
                    1:public_su2_reflection_target6.n_blocks
                ]
            ] == fill(:PSD, public_su2_reflection_target6.n_blocks)
            @test scalar_su2_reflection_mp6.constraints[end][1] == :PSD
            scalar_su2_reflection_zero_labels = [
                zc.origin.label for zc in scalar_su2_reflection_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label == (feature=:scalar_equality, index=1),
                scalar_su2_reflection_zero_labels,
            )
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :reason) &&
                    label.reason == :spin_offblock,
                scalar_su2_reflection_zero_labels,
            )
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :reason) &&
                    label.reason == :magnetic_copy,
                scalar_su2_reflection_zero_labels,
            )
            scalar_su2_reflection_linear6, scalar_su2_reflection_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                    reflection_symmetry=true,
                )
            @test scalar_su2_reflection_linear_report6.block_labels ==
                scalar_su2_reflection_report6.block_labels
            @test scalar_su2_reflection_linear_report6.psd_block_sizes ==
                scalar_su2_reflection_report6.psd_block_sizes
            @test scalar_su2_reflection_linear_report6.zero_constraint_count ==
                scalar_su2_reflection_report6.zero_constraint_count
            @test [
                zc.origin.label for zc in scalar_su2_reflection_linear6.zero_constraints
            ] == scalar_su2_reflection_zero_labels
            scalar_su2_reflection_fixed_complex_target6 =
                pauli_su2_translation_orbit_structural_targets(
                    6,
                    2;
                    real_moment_matrix=false,
                    reflection_symmetry=true,
                    momenta=[0],
                )
            scalar_su2_reflection_fixed_complex_mp6,
            scalar_su2_reflection_fixed_complex_report6 =
                pauli_translation_invariant_moment_relaxation(
                    scalar_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                    reflection_symmetry=true,
                    real_moment_matrix=false,
                    momenta=[0],
                )
            @test scalar_su2_reflection_fixed_complex_report6.block_labels[1:(end - 1)] ==
                scalar_su2_reflection_fixed_complex_target6.block_labels
            @test scalar_su2_reflection_fixed_complex_report6.psd_block_sizes[1:(end - 1)] ==
                scalar_su2_reflection_fixed_complex_target6.psd_block_sizes
            @test [
                cone for (cone, _) in
                    scalar_su2_reflection_fixed_complex_mp6.constraints[
                        1:scalar_su2_reflection_fixed_complex_target6.n_blocks
                    ]
            ] == fill(:HPSD, scalar_su2_reflection_fixed_complex_target6.n_blocks)
            @test scalar_su2_reflection_fixed_complex_mp6.constraints[end][1] ==
                :HPSD
            @test scalar_su2_reflection_fixed_complex_report6.zero_constraint_count > 1
            scalar_su2_reflection_fixed_complex_zero_labels = [
                zc.origin.label
                for zc in scalar_su2_reflection_fixed_complex_mp6.linear.zero_constraints
            ]
            @test any(
                label -> label == (feature=:scalar_equality, index=1),
                scalar_su2_reflection_fixed_complex_zero_labels,
            )
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :reason) &&
                    label.reason == :spin_offblock,
                scalar_su2_reflection_fixed_complex_zero_labels,
            )
            @test any(
                label -> label isa NamedTuple &&
                    haskey(label, :reason) &&
                    label.reason == :magnetic_copy,
                scalar_su2_reflection_fixed_complex_zero_labels,
            )
            covered_moment_su2_pop6 = polyopt(
                scalar_su2_objective6,
                registry6;
                moment_eq_constraints=[scalar_su2_objective6 * scalar_su2_objective6],
            )
            covered_moment_su2_mp6, covered_moment_su2_report6 =
                pauli_translation_invariant_moment_relaxation(
                    covered_moment_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                )
            @test covered_moment_su2_report6.block_labels ==
                translation_su2_target6.block_labels
            @test covered_moment_su2_report6.psd_block_sizes ==
                translation_su2_target6.psd_block_sizes
            @test covered_moment_su2_report6.zero_constraint_count > 1
            @test count(con -> con[1] == :Zero, covered_moment_su2_mp6.constraints) ==
                covered_moment_su2_report6.zero_constraint_count
            @test any(
                zc -> zc.origin.label.feature == :moment_equality &&
                    zc.origin.label.index == 1,
                covered_moment_su2_mp6.linear.zero_constraints,
            )
            @test any(
                zc -> zc.origin.label isa NamedTuple &&
                    haskey(zc.origin.label, :reason) &&
                    zc.origin.label.reason == :spin_offblock,
                covered_moment_su2_mp6.linear.zero_constraints,
            )
            @test any(
                zc -> zc.origin.label isa NamedTuple &&
                    haskey(zc.origin.label, :reason) &&
                    zc.origin.label.reason == :magnetic_offdiag,
                covered_moment_su2_mp6.linear.zero_constraints,
            )
            @test any(
                zc -> zc.origin.label isa NamedTuple &&
                    haskey(zc.origin.label, :reason) &&
                    zc.origin.label.reason == :magnetic_copy,
                covered_moment_su2_mp6.linear.zero_constraints,
            )
            covered_moment_su2_support6 =
                translation_solve_support(covered_moment_su2_report6)
            @test covered_moment_su2_support6.supported
            @test isnothing(covered_moment_su2_support6.blocker)
            @test isempty(covered_moment_su2_support6.unsupported_block_features)
            @test isempty(covered_moment_su2_support6.unsupported_zero_features)
            covered_moment_su2_metrics6 =
                translation_report_metrics(covered_moment_su2_report6)
            @test covered_moment_su2_metrics6.solve_supported
            @test isnothing(covered_moment_su2_metrics6.solve_blocker)
            @test covered_moment_su2_metrics6.solve_blocker_reason == ""
            @test isempty(covered_moment_su2_metrics6.solve_unsupported_block_features)
            @test isempty(covered_moment_su2_metrics6.solve_unsupported_zero_features)
            @test covered_moment_su2_metrics6.moment_equality_row_count > 0
            @test covered_moment_su2_metrics6.linear_state_opt_row_count == 0
            @test covered_moment_su2_metrics6.su2_base_zero_row_count > 0
            @test covered_moment_su2_metrics6.su2_base_spin_offblock_row_count > 0
            @test covered_moment_su2_metrics6.su2_base_magnetic_offdiag_row_count > 0
            @test covered_moment_su2_metrics6.su2_base_magnetic_copy_row_count > 0
            @test covered_moment_su2_metrics6.su2_base_zero_row_count +
                covered_moment_su2_metrics6.moment_equality_row_count ==
                covered_moment_su2_metrics6.zero_constraint_count
            @test NCTSSoS.assert_moment_linear_data_invariants(
                covered_moment_su2_mp6.linear,
                covered_moment_su2_mp6.constraints,
            ) === nothing
            covered_moment_su2_objectives = Float64[]
            for direct_linear in (false, true)
                covered_moment_su2_result = quiet() do
                    pauli_translation_invariant_nctssos(
                        covered_moment_su2_pop6,
                        ops6,
                        2,
                        optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                        dualize=false,
                        su2_symmetry=true,
                        direct_linear,
                    )
                end
                @test termination_status(covered_moment_su2_result.model) in
                    (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
                @test covered_moment_su2_result.report.block_labels ==
                    covered_moment_su2_report6.block_labels
                @test covered_moment_su2_result.report.zero_constraint_count ==
                    covered_moment_su2_report6.zero_constraint_count
                if direct_linear
                    @test covered_moment_su2_result.moment_problem isa
                        NCTSSoS.MomentLinearData
                else
                    @test covered_moment_su2_result.moment_problem isa NCTSSoS.MomentProblem
                end
                push!(covered_moment_su2_objectives, covered_moment_su2_result.objective)
            end
            @test covered_moment_su2_objectives[1] ≈
                covered_moment_su2_objectives[2] atol = 1e-6
            covered_moment_su2_config6 = SolverConfig(
                ;
                optimizer=optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false),
                order=2,
                symmetry=heisenberg_chain_symmetry_spec(
                    ops6;
                    reflection=false,
                    axis_rotations=false,
                    check_invariance=false,
                ),
            )
            covered_moment_su2_auto_result = quiet() do
                cs_nctssos(
                    covered_moment_su2_pop6,
                    covered_moment_su2_config6;
                    dualize=false,
                    direct_linear=true,
                    su2_symmetry=true,
                )
            end
            @test covered_moment_su2_auto_result.moment_problem isa
                NCTSSoS.MomentLinearData
            @test termination_status(covered_moment_su2_auto_result.model) in
                (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
            @test translation_solve_support(covered_moment_su2_auto_result).supported
            @test covered_moment_su2_auto_result.objective ≈
                covered_moment_su2_objectives[2] atol = 1e-6
            covered_moment_su2_linear6, covered_moment_su2_linear_report6 =
                NCTSSoS._pauli_translation_base_linear_relaxation(
                    covered_moment_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                )
            @test covered_moment_su2_linear_report6.block_labels ==
                covered_moment_su2_report6.block_labels
            @test covered_moment_su2_linear_report6.zero_constraint_count ==
                covered_moment_su2_report6.zero_constraint_count
            @test any(
                zc -> zc.origin.label.feature == :moment_equality &&
                    zc.origin.label.index == 1,
                covered_moment_su2_linear6.zero_constraints,
            )
            σx6, σy6, σz6 = ops6
            long_range_su2_constraint6 = sum(
                σx6[i] * σx6[mod1(i + 3, 6)] +
                σy6[i] * σy6[mod1(i + 3, 6)] +
                σz6[i] * σz6[mod1(i + 3, 6)]
                for i in 1:6
            )
            unsupported_moment_su2_pop6 = polyopt(
                scalar_su2_objective6,
                registry6;
                moment_eq_constraints=[long_range_su2_constraint6],
            )
            unsupported_moment_su2_err = try
                pauli_translation_invariant_moment_relaxation(
                    unsupported_moment_su2_pop6,
                    ops6,
                    2;
                    su2_symmetry=true,
                )
                nothing
            catch err
                err
            end
            @test unsupported_moment_su2_err isa ArgumentError
            @test occursin(
                "requires real Moment equality constraint 1 coefficients",
                sprint(showerror, unsupported_moment_su2_err),
            )
            orbit_diag6 = pauli_su2_translation_orbit_basis_reduction_diagnostics(
                orbit_reps6,
                6,
            )
            @test orbit_diag6.support_counts == translation_su2_target6.support_counts
            @test orbit_diag6.full_side == translation_su2_target6.orbit_basis_size
            @test orbit_diag6.reduced_block_sizes == [2, 2, 1]
            @test orbit_diag6.full_dense_entries == 13^2
            @test orbit_diag6.reduced_dense_entries == 2^2 + 2^2 + 1
            orbit_active_entries = sum(
                block.irrep_dimension * block.multiplicity^2
                for block in orbit_diag6.blocks;
                init=0,
            )
            @test orbit_diag6.offblock_entry_count ==
                orbit_diag6.full_dense_entries - orbit_active_entries
            @test orbit_diag6.copy_entry_count ==
                orbit_active_entries - orbit_diag6.reduced_dense_entries
            @test orbit_diag6.accounted_entry_count == orbit_diag6.full_dense_entries
            malformed_orbit_reps6 = copy(orbit_reps6)
            two_site_orbit_indices = findall(
                mono -> length(NCTSSoS._pauli_support_tuple(mono)) == 2,
                malformed_orbit_reps6,
            )
            malformed_orbit_reps6[last(two_site_orbit_indices)] =
                malformed_orbit_reps6[first(two_site_orbit_indices)]
            @test_throws ArgumentError pauli_su2_translation_orbit_basis_summary(
                malformed_orbit_reps6,
                6,
            )
            @test_throws ArgumentError pauli_su2_translation_orbit_basis_reduction_plan(
                malformed_orbit_reps6,
                6,
            )
            @test_throws ArgumentError pauli_su2_translation_orbit_basis_transform_blocks(
                malformed_orbit_reps6,
                6,
            )
            @test_throws ArgumentError pauli_su2_translation_orbit_basis_summary(
                basis6_order2,
                6,
            )
            translation_su2_reflection6 = pauli_su2_translation_orbit_structural_targets(
                6,
                2;
                reflection_symmetry=true,
            )
            @test translation_su2_reflection6.orbit_basis_size == 13
            @test translation_su2_reflection6.n_blocks == 19
            @test translation_su2_reflection6.logical_max_block == 2
            @test translation_su2_reflection6.psd_max_block == 4
            @test translation_su2_reflection6.logical_total_block_side == 25
            @test translation_su2_reflection6.psd_total_block_side == 34
            @test translation_su2_reflection6.logical_block_histogram == [1 => 13, 2 => 6]
            @test translation_su2_reflection6.psd_block_histogram == [1 => 8, 2 => 9, 4 => 2]
            @test translation_su2_reflection6.psd_dense_entries == 76
            @test translation_su2_reflection6.psd_symmetric_entries == 55
            @test translation_su2_reflection6.reflection_combined
            @test length(translation_su2_reflection6.block_labels) ==
                length(translation_su2_reflection6.logical_block_sizes)
            @test first(translation_su2_reflection6.block_labels).decomposition ==
                :translation_su2_reflection
            @test translation_su2_reflection6.block_labels[1:4] == [
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=0,
                    reflection_fixed=true,
                    spin2=0,
                    reflection_parity=:plus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=0,
                    reflection_fixed=true,
                    spin2=2,
                    reflection_parity=:plus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=0,
                    reflection_fixed=true,
                    spin2=2,
                    reflection_parity=:minus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=0,
                    reflection_fixed=true,
                    spin2=4,
                    reflection_parity=:plus,
                ),
            ]
            @test translation_su2_reflection6.block_labels[5:10] == [
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=1,
                    reflection_fixed=false,
                    spin2=0,
                    reflection_parity=:plus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=1,
                    reflection_fixed=false,
                    spin2=0,
                    reflection_parity=:minus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=1,
                    reflection_fixed=false,
                    spin2=2,
                    reflection_parity=:plus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=1,
                    reflection_fixed=false,
                    spin2=2,
                    reflection_parity=:minus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=1,
                    reflection_fixed=false,
                    spin2=4,
                    reflection_parity=:plus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=1,
                    reflection_fixed=false,
                    spin2=4,
                    reflection_parity=:minus,
                ),
            ]
            @test translation_su2_reflection6.block_labels[end - 2:end] == [
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=3,
                    reflection_fixed=true,
                    spin2=0,
                    reflection_parity=:plus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=3,
                    reflection_fixed=true,
                    spin2=2,
                    reflection_parity=:minus,
                ),
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=3,
                    reflection_fixed=true,
                    spin2=4,
                    reflection_parity=:plus,
                ),
            ]
            @test sum(
                pair.second for pair in
                    translation_su2_reflection6.logical_block_feature_histogram;
                init=0,
            ) == length(translation_su2_reflection6.logical_block_sizes)
            @test (
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=0,
                    reflection_fixed=true,
                    spin2=2,
                    reflection_parity=:minus,
                    size=1,
                ) => 1
            ) in translation_su2_reflection6.logical_block_feature_histogram
            @test (
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=1,
                    reflection_fixed=false,
                    spin2=2,
                    reflection_parity=:minus,
                    size=2,
                ) => 1
            ) in translation_su2_reflection6.logical_block_feature_histogram
            @test sum(
                pair.second for pair in translation_su2_reflection6.psd_block_feature_histogram;
                init=0,
            ) == length(translation_su2_reflection6.psd_block_sizes)
            @test (
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=0,
                    reflection_fixed=true,
                    spin2=0,
                    reflection_parity=:plus,
                    size=4,
                ) => 1
            ) in translation_su2_reflection6.psd_block_feature_histogram
            @test (
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=1,
                    reflection_fixed=false,
                    spin2=2,
                    reflection_parity=:minus,
                    size=2,
                ) => 1
            ) in translation_su2_reflection6.psd_block_feature_histogram

            translation_su2_target100 = pauli_su2_translation_orbit_structural_targets(100, 4)
            @test translation_su2_target100.orbit_basis_size == 121
            @test translation_su2_target100.momentum_sector_count == 51
            @test [(block.spin2, block.multiplicity) for block in translation_su2_target100.zero_momentum_blocks] ==
                [(0, 6), (2, 11), (4, 9), (6, 4), (8, 1)]
            @test [(block.spin2, block.multiplicity) for block in translation_su2_target100.nonzero_momentum_blocks] ==
                [(0, 5), (2, 11), (4, 9), (6, 4), (8, 1)]
            @test translation_su2_target100.n_blocks == 255
            @test translation_su2_target100.logical_max_block == 11
            @test translation_su2_target100.psd_max_block == 22
            @test translation_su2_target100.logical_total_block_side == 1_531
            @test translation_su2_target100.psd_total_block_side == 3_062
            @test translation_su2_target100.logical_block_histogram ==
                [1 => 51, 4 => 51, 5 => 50, 6 => 1, 9 => 51, 11 => 51]
            @test translation_su2_target100.psd_block_histogram ==
                [2 => 51, 8 => 51, 10 => 50, 12 => 1, 18 => 51, 22 => 51]
            @test translation_su2_target100.psd_dense_entries == 49_820
            @test translation_su2_target100.psd_symmetric_entries == 26_441
            @test translation_su2_target100.translation_combined
            @test translation_su2_target100.su2_full_dense_entries == 734_641
            @test translation_su2_target100.su2_active_dense_entries == 46_625
            @test translation_su2_target100.su2_reduced_dense_entries == 12_455
            @test translation_su2_target100.offblock_entry_count == 688_016
            @test translation_su2_target100.copy_entry_count == 34_170
            @test translation_su2_target100.accounted_entry_count ==
                translation_su2_target100.su2_full_dense_entries
            @test translation_su2_target100.block_coefficient_domains ==
                fill(
                    :complex_algebraic_float64,
                    length(translation_su2_target100.block_labels),
                )
            @test translation_su2_target100.block_exact_coefficient_domains ==
                fill(
                    :complex_sqrt_rational,
                    length(translation_su2_target100.block_labels),
                )
            @test translation_su2_target100.block_coefficient_domain_histogram ==
                [:complex_algebraic_float64 => translation_su2_target100.n_blocks]
            @test translation_su2_target100.block_exact_coefficient_domain_histogram ==
                [:complex_sqrt_rational => translation_su2_target100.n_blocks]
            target100_accounting_records = translation_su2_target100.su2_accounting_records
            @test translation_su2_target100.su2_accounting_record_count ==
                translation_su2_target100.momentum_sector_count
            @test translation_su2_target100.estimated_model_size_gate_status ==
                :blocked_missing_scalar_equality_estimate
            @test length(target100_accounting_records) ==
                translation_su2_target100.su2_accounting_record_count
            @test sum(
                record.su2_full_dense_entries for record in target100_accounting_records;
                init=0,
            ) == translation_su2_target100.su2_full_dense_entries
            @test sum(
                record.su2_active_dense_entries for record in target100_accounting_records;
                init=0,
            ) == translation_su2_target100.su2_active_dense_entries
            @test sum(
                record.su2_reduced_dense_entries for record in target100_accounting_records;
                init=0,
            ) == translation_su2_target100.su2_reduced_dense_entries
            @test sum(
                record.offblock_entry_count for record in target100_accounting_records;
                init=0,
            ) == translation_su2_target100.offblock_entry_count
            @test sum(
                record.copy_entry_count for record in target100_accounting_records;
                init=0,
            ) == translation_su2_target100.copy_entry_count
            @test !translation_su2_target100.solve_supported
            @test translation_su2_target100.solve_blocker == :structural_target_only
            @test isempty(translation_su2_target100.solve_unsupported_block_features)
            @test isempty(translation_su2_target100.solve_unsupported_zero_features)
            translation_su2_reflection100 = pauli_su2_translation_orbit_structural_targets(
                100,
                4;
                reflection_symmetry=true,
            )
            @test translation_su2_reflection100.orbit_basis_size == 121
            @test translation_su2_reflection100.n_blocks == 507
            @test translation_su2_reflection100.logical_max_block == 11
            @test translation_su2_reflection100.psd_max_block == 16
            @test translation_su2_reflection100.logical_total_block_side == 3_001
            @test translation_su2_reflection100.psd_total_block_side == 3_062
            @test translation_su2_reflection100.logical_block_histogram ==
                [1 => 102, 2 => 2, 3 => 4, 4 => 98, 5 => 101, 6 => 3, 8 => 1, 9 => 98, 11 => 98]
            @test translation_su2_reflection100.psd_block_histogram ==
                [1 => 98, 2 => 4, 4 => 100, 5 => 98, 6 => 4, 9 => 98, 10 => 3, 11 => 98, 12 => 3, 16 => 1]
            @test translation_su2_reflection100.psd_dense_entries == 25_092
            @test translation_su2_reflection100.psd_symmetric_entries == 14_077
            @test translation_su2_reflection100.reflection_combined
            @test translation_su2_reflection100.su2_full_dense_entries == 1_426_033
            @test translation_su2_reflection100.su2_active_dense_entries == 90_619
            @test translation_su2_reflection100.su2_reduced_dense_entries == 24_207
            @test translation_su2_reflection100.offblock_entry_count == 1_335_414
            @test translation_su2_reflection100.copy_entry_count == 66_412
            @test translation_su2_reflection100.accounted_entry_count ==
                translation_su2_reflection100.su2_full_dense_entries
            @test translation_su2_reflection100.block_coefficient_domains ==
                fill(
                    :complex_algebraic_float64,
                    length(translation_su2_reflection100.block_labels),
                )
            @test translation_su2_reflection100.block_exact_coefficient_domains ==
                fill(
                    :complex_sqrt_rational,
                    length(translation_su2_reflection100.block_labels),
                )
            @test translation_su2_reflection100.block_coefficient_domain_histogram ==
                [:complex_algebraic_float64 => translation_su2_reflection100.n_blocks]
            @test translation_su2_reflection100.block_exact_coefficient_domain_histogram ==
                [:complex_sqrt_rational => translation_su2_reflection100.n_blocks]
            reflection100_accounting_records =
                translation_su2_reflection100.su2_accounting_records
            @test translation_su2_reflection100.su2_accounting_record_count == 102
            @test length(reflection100_accounting_records) ==
                translation_su2_reflection100.su2_accounting_record_count
            @test sum(
                record.su2_full_dense_entries for record in reflection100_accounting_records;
                init=0,
            ) == translation_su2_reflection100.su2_full_dense_entries
            @test sum(
                record.su2_active_dense_entries for record in reflection100_accounting_records;
                init=0,
            ) == translation_su2_reflection100.su2_active_dense_entries
            @test sum(
                record.su2_reduced_dense_entries for record in reflection100_accounting_records;
                init=0,
            ) == translation_su2_reflection100.su2_reduced_dense_entries
            @test sum(
                record.offblock_entry_count for record in reflection100_accounting_records;
                init=0,
            ) == translation_su2_reflection100.offblock_entry_count
            @test sum(
                record.copy_entry_count for record in reflection100_accounting_records;
                init=0,
            ) == translation_su2_reflection100.copy_entry_count
            @test !translation_su2_reflection100.solve_supported
            @test translation_su2_reflection100.solve_blocker == :structural_target_only
            @test occursin(
                "Structural target only",
                translation_su2_reflection100.solve_blocker_reason,
            )
            @test isempty(translation_su2_reflection100.solve_unsupported_block_features)
            @test isempty(translation_su2_reflection100.solve_unsupported_zero_features)
            @test_throws ArgumentError pauli_su2_translation_orbit_structural_targets(
                6,
                2;
                real_moment_matrix=false,
                momenta=[1],
                reflection_symmetry=true,
            )
            @test_throws ArgumentError pauli_su2_translation_orbit_structural_targets(6, 3)
            @test_throws ArgumentError pauli_su2_translation_orbit_structural_targets(100, 4; scalar_bytes=0)
            @test_throws ArgumentError pauli_su2_contiguous_chain_structural_targets(6, 3)
            @test_throws ArgumentError pauli_su2_contiguous_chain_structural_targets(100, 4; scalar_bytes=0)
        end
        reduction_plan_order2 = pauli_su2_basis_reduction_plan(basis_order2)
        @test reduction_plan_order2.support_counts == [0 => 1, 1 => 4, 2 => 4]
        @test reduction_plan_order2.full_side == length(basis_order2)
        @test reduction_plan_order2.reduced_block_sizes == [5, 8, 4]
        @test reduction_plan_order2.transformed_block_sizes == [5, 24, 20]
        @test reduction_plan_order2.transformed_max_block == 24
        @test reduction_plan_order2.transformed_total_block_side ==
            reduction_plan_order2.full_side
        @test [
            (block.spin2, block.multiplicity, block.irrep_dimension, length(block.logical_row_labels))
            for block in reduction_plan_order2.blocks
        ] == [(0, 5, 1, 5), (2, 8, 3, 8), (4, 4, 5, 4)]
        @test reduction_plan_order2.block_coefficient_domains ==
            fill(:complex_algebraic_float64, length(reduction_plan_order2.blocks))
        @test reduction_plan_order2.block_exact_coefficient_domains ==
            fill(:complex_sqrt_rational, length(reduction_plan_order2.blocks))
        @test all(
            block -> block.coefficient_domain == :complex_algebraic_float64,
            reduction_plan_order2.blocks,
        )
        @test all(
            block -> block.exact_coefficient_domain == :complex_sqrt_rational,
            reduction_plan_order2.blocks,
        )
        spin0_labels = reduction_plan_order2.blocks[1].logical_row_labels
        @test first(spin0_labels) == (
            feature=:pauli_su2_basis_multiplicity,
            support=(),
            support_size=0,
            spin2=0,
            local_multiplicity=1,
        )
        @test count(label -> label.support_size == 2, spin0_labels) == 4
        spin2_labels = reduction_plan_order2.blocks[2].logical_row_labels
        @test count(label -> label.support_size == 1, spin2_labels) == 4
        @test count(label -> label.support_size == 2, spin2_labels) == 4
        basis_transform_blocks = pauli_su2_basis_transform_blocks(basis_order2)
        @test [
            (block.spin2, block.multiplicity, block.irrep_dimension, size(block.transform))
            for block in basis_transform_blocks
        ] == [
            (0, 5, 1, (5, length(basis_order2))),
            (2, 8, 3, (24, length(basis_order2))),
            (4, 4, 5, (20, length(basis_order2))),
        ]
        basis_transform = reduce(vcat, [block.transform for block in basis_transform_blocks])
        @test basis_transform * basis_transform' ≈ I atol = 1e-12
        @test [block.multiplicity for block in basis_transform_blocks] ==
            reduction_plan_order2.reduced_block_sizes
        @test all(
            block -> block.coefficient_domain == :complex_algebraic_float64,
            basis_transform_blocks,
        )
        @test all(
            block -> block.exact_coefficient_domain == :complex_sqrt_rational,
            basis_transform_blocks,
        )
        @test first(basis_transform_blocks[1].row_labels) == (
            feature=:pauli_su2_basis_state,
            support=(),
            support_size=0,
            spin2=0,
            m2=0,
            local_multiplicity=1,
        )
        basis_singlets_order2 = pauli_su2_basis_singlet_channels(basis_order2)
        @test basis_singlets_order2.basis_size == length(basis_order2)
        @test basis_singlets_order2.channel_count == 5
        @test basis_singlets_order2.row_labels == basis_transform_blocks[1].row_labels
        @test count(==(()), basis_singlets_order2.channel_supports) == 1
        @test count(support -> length(support) == 2, basis_singlets_order2.channel_supports) == 4
        @test basis_singlets_order2.coefficient_domain == :complex_algebraic_float64
        @test basis_singlets_order2.exact_coefficient_domain == :complex_sqrt_rational
        @test basis_singlets_order2.transform == basis_transform_blocks[1].transform
        @test basis_singlets_order2.singlet_orthonormality_residual <= 1e-12
        basis_equalities_order2 = pauli_su2_basis_singlet_channel_equalities(basis_order2)
        @test basis_equalities_order2.basis_size == length(basis_order2)
        @test basis_equalities_order2.equality_count ==
            length(basis_order2) - basis_singlets_order2.channel_count
        @test basis_equalities_order2.row_labels == reduce(
            vcat,
            [
                block.row_labels
                for block in basis_transform_blocks
                if block.spin2 != 0
            ],
        )
        @test count(label -> label.support_size == 1, basis_equalities_order2.row_labels) == 12
        @test count(label -> label.support_size == 2, basis_equalities_order2.row_labels) == 32
        @test basis_equalities_order2.coefficient_domain == :complex_algebraic_float64
        @test basis_equalities_order2.exact_coefficient_domain == :complex_sqrt_rational
        @test basis_equalities_order2.transform * basis_equalities_order2.transform' ≈ I atol = 1e-12
        @test basis_equalities_order2.transform * basis_singlets_order2.transform' ≈
            zeros(
                ComplexF64,
                basis_equalities_order2.equality_count,
                basis_singlets_order2.channel_count,
            ) atol = 1e-12
        @test basis_equalities_order2.equality_orthonormality_residual <= 1e-12
        @test basis_equalities_order2.singlet_orthogonality_residual <= 1e-12
        @test length(basis_equalities_order2.column_forms) ==
            basis_equalities_order2.equality_count
        @test length(basis_equalities_order2.basis_forms) ==
            basis_equalities_order2.equality_count
        @test all(
            form -> all(first(term) in basis_order2 && abs(last(term)) > 1e-12 for term in form),
            basis_equalities_order2.basis_forms,
        )
        @test count(
            label -> label.support_size == 2 && label.m2 == 0,
            basis_transform_blocks[2].row_labels,
        ) == 4
        identity_basis_block = Matrix{ComplexF64}(I, length(basis_order2), length(basis_order2))
        reduced_moment_blocks = pauli_su2_basis_reduced_moment_blocks(
            identity_basis_block,
            basis_order2,
        )
        @test [
            (
                block.spin2,
                block.multiplicity,
                block.irrep_dimension,
                block.reference_m2,
                size(block.matrix),
                length(block.row_labels),
            )
            for block in reduced_moment_blocks
        ] == [
            (0, 5, 1, 0, (5, 5), 5),
            (2, 8, 3, -2, (8, 8), 8),
            (4, 4, 5, -4, (4, 4), 4),
        ]
        @test all(block -> isapprox(block.matrix, I; atol=1e-12), reduced_moment_blocks)
        @test all(
            block -> block.coefficient_domain == :complex_algebraic_float64,
            reduced_moment_blocks,
        )
        @test all(
            block -> block.exact_coefficient_domain == :complex_sqrt_rational,
            reduced_moment_blocks,
        )
        @test first(reduced_moment_blocks[2].row_labels).m2 == -2
        poly_one = Polynomial(one(σx4[1]))
        _, polynomial_basis_block = NCTSSoS._build_constraint_matrix(poly_one, basis_order2, :HPSD)
        polynomial_moment_blocks = pauli_su2_basis_reduced_moment_blocks(
            polynomial_basis_block,
            basis_order2,
        )
        @test [(block.spin2, size(block.matrix)) for block in polynomial_moment_blocks] ==
            [(0, (5, 5)), (2, (8, 8)), (4, (4, 4))]
        @test all(
            block -> all(entry -> entry isa typeof(poly_one), block.matrix),
            polynomial_moment_blocks,
        )
        @test_throws ArgumentError pauli_su2_basis_moment_problem(poly_one, basis_order2)
        checked_su2_basis_mp = pauli_su2_basis_moment_problem(
            poly_one,
            basis_order2;
            ops=ops4,
        )
        @test [block.size for block in checked_su2_basis_mp.linear.psd_blocks_lin] ==
            [5, 8, 4]
        @test isempty(checked_su2_basis_mp.linear.zero_constraints)
        singlet_equality_mp = pauli_su2_basis_moment_problem(
            poly_one,
            basis_order2;
            ops=ops4,
            singlet_channel_equalities=true,
        )
        singlet_equality_zero_constraint_count =
            count(con -> con[1] == :Zero, singlet_equality_mp.constraints)
        @test singlet_equality_zero_constraint_count >=
            basis_equalities_order2.equality_count
        @test singlet_equality_zero_constraint_count ==
            length(singlet_equality_mp.linear.zero_constraints)
        @test [block.size for block in singlet_equality_mp.linear.psd_blocks_lin] ==
            [5, 8, 4]
        @test length(singlet_equality_mp.linear.zero_constraints) >=
            basis_equalities_order2.equality_count
        singlet_equality_zero_labels = [
            zero.origin.label
            for zero in singlet_equality_mp.linear.zero_constraints
        ]
        @test all(
            label -> label.feature == :pauli_su2_singlet_channel_equality,
            singlet_equality_zero_labels,
        )
        @test Set(label.row for label in singlet_equality_zero_labels) ==
            Set(1:basis_equalities_order2.equality_count)
        @test all(label -> label.su2_label.spin2 != 0, singlet_equality_zero_labels)
        @test all(
            zero -> zero.origin isa NCTSSoS.PauliSU2SingletChannelEqualityOrigin,
            singlet_equality_mp.linear.zero_constraints,
        )
        singlet_equality_zero_diagnostic = only(
            NCTSSoS.sos_zero_dual_diagnostics([
                (
                    origin=first(singlet_equality_mp.linear.zero_constraints).origin,
                    form=first(singlet_equality_mp.linear.zero_constraints).form,
                    kind=:zero,
                    value=0.0,
                ),
            ]),
        )
        @test singlet_equality_zero_diagnostic.coefficient_domain ==
            :complex_algebraic_float64
        @test singlet_equality_zero_diagnostic.exact_coefficient_domain ==
            :complex_sqrt_rational
        @test singlet_equality_zero_diagnostic.feature ==
            :pauli_su2_singlet_channel_equality
        @test singlet_equality_zero_diagnostic.decomposition === nothing
        @test singlet_equality_zero_diagnostic.reason === nothing
        @test NCTSSoS.assert_moment_linear_data_invariants(
            singlet_equality_mp.linear,
            singlet_equality_mp.constraints,
        ) === nothing
        noninvariant_su2_err = try
            pauli_su2_basis_moment_problem(
                Polynomial(σx4[1]),
                basis_order2;
                ops=ops4,
            )
            nothing
        catch err
            err
        end
        @test noninvariant_su2_err isa ArgumentError
        @test occursin("SU(2) basis objective", sprint(showerror, noninvariant_su2_err))
        basis_order1 = pauli_contiguous_chain_basis(ops4, 1)
        uncovered_objective = Polynomial(σx4[1] * σx4[2] * σx4[3])
        uncovered_err = try
            pauli_su2_basis_moment_problem(
                uncovered_objective,
                basis_order1;
                assume_su2_invariant=true,
            )
            nothing
        catch err
            err
        end
        @test uncovered_err isa ArgumentError
        @test occursin("SU(2) basis objective", sprint(showerror, uncovered_err))

        su2_basis_mp = pauli_su2_basis_moment_problem(
            poly_one,
            basis_order2;
            assume_su2_invariant=true,
        )
        @test [cone for (cone, _) in su2_basis_mp.constraints] == [:HPSD, :HPSD, :HPSD]
        @test [size(mat, 1) for (_, mat) in su2_basis_mp.constraints] == [5, 8, 4]
        @test [block.size for block in su2_basis_mp.linear.psd_blocks_lin] == [5, 8, 4]
        @test [
            block.meta.origin.label.spin2
            for block in su2_basis_mp.linear.psd_blocks_lin
        ] == [0, 2, 4]
        @test first(su2_basis_mp.linear.psd_blocks_lin[1].meta.origin.logical_row_labels).feature ==
            :pauli_su2_basis_state
        @test all(
            block -> block.meta.origin.transform.coefficient_domain == :complex_algebraic_float64,
            su2_basis_mp.linear.psd_blocks_lin,
        )
        su2_basis_sos = NCTSSoS.sos_dualize(su2_basis_mp)
        @test [
            block.meta.origin.label.spin2
            for block in sos_dual_blocks(su2_basis_sos)
        ] == [0, 2, 4]
        @test all(
            block -> block.meta.origin.transform.exact_coefficient_domain ==
                :complex_sqrt_rational,
            sos_dual_blocks(su2_basis_sos),
        )
        su2_basis_model, _ = build_jump_model(su2_basis_mp)
        @test !isnothing(su2_basis_model)

        T = eltype(first(basis_order2).word)
        M = eltype(basis_order2)
        MP_P = Polynomial{PauliAlgebra,T,ComplexF64}
        heisenberg_obj = convert(MP_P, heisenberg_chain_hamiltonian(ops4))
        _, full_basis_block = NCTSSoS._build_constraint_matrix(
            one(MP_P),
            basis_order2,
            :HPSD,
        )
        full_constraints = [(:HPSD, full_basis_block)]
        full_basis = NCTSSoS._collect_pauli_su2_moment_basis(
            heisenberg_obj,
            full_constraints,
        )
        full_mp = NCTSSoS.MomentProblem{PauliAlgebra,T,M,MP_P}(
            heisenberg_obj,
            full_constraints,
            full_basis,
            length(full_basis);
            real_moments=false,
        )
        reduced_heisenberg_mp = pauli_su2_basis_moment_problem(
            heisenberg_obj,
            basis_order2;
            assume_su2_invariant=true,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(
            reduced_heisenberg_mp.linear,
            reduced_heisenberg_mp.constraints,
        ) === nothing
        @test [block.size for block in reduced_heisenberg_mp.linear.psd_blocks_lin] == [5, 8, 4]
        full_psd_entries = sum(
            block.size^2 for block in full_mp.linear.psd_blocks_lin;
            init=0,
        )
        reduced_psd_entries = sum(
            block.size^2 for block in reduced_heisenberg_mp.linear.psd_blocks_lin;
            init=0,
        )
        @test reduced_psd_entries < full_psd_entries

        @test_throws DimensionMismatch pauli_su2_basis_reduced_moment_blocks(
            identity_basis_block[1:(end - 1), :],
            basis_order2,
        )
        reduction_diag = pauli_su2_basis_moment_reduction_diagnostics(
            identity_basis_block,
            basis_order2,
        )
        @test reduction_diag.full_side == length(basis_order2)
        @test reduction_diag.reduced_block_sizes == [5, 8, 4]
        @test reduction_diag.transformed_block_sizes == [5, 24, 20]
        @test reduction_diag.unitarity_residual <= 1e-12
        @test reduction_diag.spin_offblock_residual <= 1e-12
        @test reduction_diag.magnetic_offdiag_residual <= 1e-12
        @test reduction_diag.magnetic_copy_residual <= 1e-12
        @test reduction_diag.max_reduction_residual <= 1e-12
        reduction_order2 = pauli_su2_basis_reduction_diagnostics(basis_order2)
        @test reduction_order2.full_dense_entries == length(basis_order2)^2
        @test reduction_order2.reduced_dense_entries == 5^2 + 8^2 + 4^2
        @test reduction_order2.offblock_entry_count == 2_104
        @test reduction_order2.copy_entry_count == 192
        @test reduction_order2.accounted_entry_count == length(basis_order2)^2
        @test_throws ArgumentError pauli_su2_basis_reduction_diagnostics(basis_order2; scalar_bytes=0)
        basis_spin_diag_order2 = pauli_su2_basis_spin_diagnostics(basis_order2)
        @test basis_spin_diag_order2.support_counts == [0 => 1, 1 => 4, 2 => 4]
        @test basis_spin_diag_order2.dimension == length(basis_order2)
        @test basis_spin_diag_order2.state_count == length(basis_order2)
        @test basis_spin_diag_order2.spin_multiplicities == [0 => 5, 2 => 8, 4 => 4]
        @test basis_spin_diag_order2.unitarity_residual <= 1e-12
        @test basis_spin_diag_order2.sz_residual <= 1e-12
        @test basis_spin_diag_order2.casimir_residual <= 1e-10
        @test basis_spin_diag_order2.offblock_residual <= 1e-10
        @test_throws ArgumentError pauli_su2_basis_spin_diagnostics(basis_order2; max_dimension=1)
        incomplete_basis = [
            one(σx4[1]),
            σx4[1], σy4[1], σz4[1],
            only(monomials(σx4[1] * σx4[2])),
        ]
        @test_throws ArgumentError pauli_su2_basis_summary(incomplete_basis)
        @test_throws ArgumentError pauli_su2_basis_reduction_plan(incomplete_basis)
        @test_throws ArgumentError pauli_su2_basis_reduction_diagnostics(incomplete_basis)
        @test_throws ArgumentError pauli_su2_basis_spin_diagnostics(incomplete_basis)
        @test_throws ArgumentError pauli_su2_basis_transform_blocks(incomplete_basis)
        @test_throws ArgumentError pauli_su2_basis_singlet_channels(incomplete_basis)
        @test_throws ArgumentError pauli_su2_basis_singlet_channel_equalities(incomplete_basis)

        blocks4 = pauli_su2_rdm_blocks(4)
        @test [(block.spin2, block.multiplicity, block.irrep_dimension) for block in blocks4] ==
            [(0, 2, 1), (2, 3, 3), (4, 1, 5)]
        @test sum(block.multiplicity * block.irrep_dimension for block in blocks4) == 1 << 4
        triplet_cg = NCTSSoS._pauli_su2_spin_half_coupling_coefficients(1, 2, 0)
        @test triplet_cg == (
            denominator=4,
            up_numerator=2,
            up_sign=1,
            down_numerator=2,
            down_sign=1,
        )
        singlet_cg = NCTSSoS._pauli_su2_spin_half_coupling_coefficients(1, 0, 0)
        @test singlet_cg == (
            denominator=4,
            up_numerator=2,
            up_sign=-1,
            down_numerator=2,
            down_sign=1,
        )
        @test_throws ArgumentError NCTSSoS._pauli_su2_spin_half_coupling_coefficients(1, 1, 0)
        @test_throws DomainError NCTSSoS._pauli_su2_spin_half_coupling_coefficients(1, 2, 1)
        transform2, states2 = NCTSSoS._pauli_su2_schur_matrix(2)
        singlet_col = only(findall(state -> state.spin2 == 0, states2))
        triplet_cols = findall(state -> state.spin2 == 2, states2)
        @test transform2[:, singlet_col] ≈ [0.0, -inv(sqrt(2)), inv(sqrt(2)), 0.0] atol = 1e-12
        @test transform2[:, triplet_cols] * transform2[:, triplet_cols]' ≈ I -
            transform2[:, singlet_col] * transform2[:, singlet_col]' atol = 1e-12
        transform4, states4 = NCTSSoS._pauli_su2_schur_matrix(4)
        @test transform4' * transform4 ≈ I atol = 1e-12
        @test [
            count(state -> state.spin2 == block.spin2 && state.m2 == -block.spin2, states4)
            for block in blocks4
        ] == [block.multiplicity for block in blocks4]
        @test_throws DomainError pauli_su2_rdm_blocks(-1)

        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
        field_pop = polyopt(sum(first(ops)), registry)

        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            contiguous_rdm_k=2,
            contiguous_rdm_decomposition=:su2,
        )
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            field_pop,
            ops,
            1;
            sign_symmetry=false,
            su2_symmetry=true,
            contiguous_rdm_k=2,
            contiguous_rdm_decomposition=:su2,
        )

        mp, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            1;
            sign_symmetry=false,
            su2_symmetry=true,
            contiguous_rdm_k=2,
            contiguous_rdm_decomposition=:su2,
        )

        @test report.block_labels[end - 1:end] == [
            (feature=:contiguous_rdm, k=2, decomposition=:su2, spin2=0),
            (feature=:contiguous_rdm, k=2, decomposition=:su2, spin2=2),
        ]
        @test report.block_coefficient_domains[end - 1:end] ==
            [:algebraic_float64, :algebraic_float64]
        @test report.block_exact_coefficient_domains[end - 1:end] ==
            [:sqrt_rational, :sqrt_rational]
        @test report.logical_block_sizes[end - 1:end] == [1, 1]
        @test report.psd_block_sizes[end - 1:end] == [1, 1]
        su2_rdm_metrics = translation_report_metrics(report)
        @test su2_rdm_metrics.contiguous_rdm
        @test su2_rdm_metrics.su2_rdm_symmetry
        @test su2_rdm_metrics.su2_moment_symmetry
        @test !su2_rdm_metrics.su2_reflection_moment_symmetry
        @test su2_rdm_metrics.contiguous_rdm_block_count == 2
        @test su2_rdm_metrics.contiguous_rdm_logical_block_sizes == [1, 1]
        @test su2_rdm_metrics.contiguous_rdm_psd_block_sizes == [1, 1]
        @test su2_rdm_metrics.su2_base_zero_row_count > 0
        @test su2_rdm_metrics.contiguous_rdm_zero_row_count > 0
        @test su2_rdm_metrics.su2_base_zero_row_count +
            su2_rdm_metrics.contiguous_rdm_zero_row_count ==
            su2_rdm_metrics.zero_constraint_count
        @test count(con -> con[1] == :Zero, mp.constraints) > 0
        su2_zero_labels = [zc.origin.label for zc in mp.linear.zero_constraints]
        @test all(
            label -> label.feature in (:pauli_su2_base_zero, :contiguous_rdm_zero),
            su2_zero_labels,
        )
        @test any(
            label -> label.feature == :pauli_su2_base_zero &&
                label.reason == :spin_offblock,
            su2_zero_labels,
        )
        @test any(
            label -> label.feature == :contiguous_rdm_zero &&
                label.decomposition == :su2 &&
                label.reason == :spin_magnetic_offblock,
            su2_zero_labels,
        )
        @test any(
            label -> label.feature == :contiguous_rdm_zero &&
                label.decomposition == :su2 &&
                label.reason == :magnetic_copy,
            su2_zero_labels,
        )
        su2_zero_histogram = translation_zero_origin_histogram(mp)
        @test sum(last, su2_zero_histogram; init=0) == length(mp.linear.zero_constraints)
        @test all(
            first(pair).feature in (:pauli_su2_base_zero, :contiguous_rdm_zero)
            for pair in su2_zero_histogram
        )
        @test any(first(pair).feature == :pauli_su2_base_zero for pair in su2_zero_histogram)
        @test any(
            first(pair).feature == :contiguous_rdm_zero &&
                first(pair).decomposition == :su2 &&
                first(pair).reason == :spin_magnetic_offblock
            for pair in su2_zero_histogram
        )
        @test any(
            first(pair).feature == :contiguous_rdm_zero &&
                first(pair).decomposition == :su2 &&
                first(pair).reason == :magnetic_copy
            for pair in su2_zero_histogram
        )
        @test all(cone == :PSD for (cone, _) in mp.constraints[end - 1:end])
        @test [
            length(block.meta.origin.logical_row_labels)
            for block in mp.linear.psd_blocks_lin[end - 1:end]
        ] == [1, 1]
        @test [
            block.meta.origin.transform.family
            for block in mp.linear.psd_blocks_lin[end - 1:end]
        ] == [:pauli_su2_rdm, :pauli_su2_rdm]
        @test [
            block.meta.origin.transform.coefficient_domain
            for block in mp.linear.psd_blocks_lin[end - 1:end]
        ] == [:algebraic_float64, :algebraic_float64]
        @test [
            block.meta.origin.transform.exact_coefficient_domain
            for block in mp.linear.psd_blocks_lin[end - 1:end]
        ] == [:sqrt_rational, :sqrt_rational]
        su2_rdm_profile = translation_symmetry_profile(report)
        @test su2_rdm_profile.su2_moment_symmetry
        @test !su2_rdm_profile.su2_reflection_moment_symmetry
        @test su2_rdm_profile.su2_rdm_symmetry
        su2_rdm_support = translation_solve_support(report)
        @test su2_rdm_support.supported
        @test isempty(su2_rdm_support.unsupported_block_features)
        @test isempty(su2_rdm_support.unsupported_zero_features)
        @test [
            (block.meta.origin.transform.k, block.meta.origin.transform.spin2)
            for block in mp.linear.psd_blocks_lin[end - 1:end]
        ] == [(2, 0), (2, 2)]
        @test [
            (block.meta.origin.transform.reference_m2, size(block.meta.origin.transform.schur_matrix))
            for block in mp.linear.psd_blocks_lin[end - 1:end]
        ] == [(0, (4, 4)), (-2, (4, 4))]
        sos = NCTSSoS.sos_dualize(mp)
        dual_blocks = sos_dual_blocks(sos)
        @test [
            block.meta.origin.transform.family
            for block in dual_blocks[end - 1:end]
        ] == [:pauli_su2_rdm, :pauli_su2_rdm]
        @test [
            block.meta.origin.transform.coefficient_domain
            for block in dual_blocks[end - 1:end]
        ] == [:algebraic_float64, :algebraic_float64]

        full_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                contiguous_rdm_k=2,
            )
        end
        su2_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:su2,
            )
        end
        @test termination_status(full_result.model) == JuMP.MOI.OPTIMAL
        @test termination_status(su2_result.model) == JuMP.MOI.OPTIMAL
        @test su2_result.objective ≈ full_result.objective atol = 1e-6
        @test sum(su2_result.report.psd_block_sizes[end - 1:end]) <
            full_result.report.psd_block_sizes[end]

        complex_su2_direct_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop,
                ops,
                1,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                real_moment_matrix=false,
                direct_linear=true,
                su2_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:su2,
            )
        end
        @test termination_status(complex_su2_direct_result.model) == JuMP.MOI.OPTIMAL
        @test complex_su2_direct_result.moment_problem isa NCTSSoS.MomentLinearData
        @test complex_su2_direct_result.objective ≈ su2_result.objective atol = 1e-6
        @test complex_su2_direct_result.report.block_labels[end - 1:end] ==
            su2_result.report.block_labels[end - 1:end]
        complex_su2_direct_zero_labels = [
            zero.origin.label
            for zero in complex_su2_direct_result.moment_problem.zero_constraints
        ]
        @test any(
            label -> label.feature == :pauli_su2_base_zero &&
                label.reason == :spin_offblock,
            complex_su2_direct_zero_labels,
        )
        @test any(
            label -> label.feature == :contiguous_rdm_zero &&
                label.decomposition == :su2 &&
                haskey(label, :component) &&
                label.component == :complex,
            complex_su2_direct_zero_labels,
        )
        @test all(
            zero -> zero.origin.label.feature == :pauli_su2_base_zero ||
                (
                    zero.origin.label.feature == :contiguous_rdm_zero &&
                    haskey(zero.origin.label, :component) &&
                    zero.origin.label.component == :complex &&
                    zero.origin.part == :scalar
                ),
            complex_su2_direct_result.moment_problem.zero_constraints,
        )

        su2_sos = NCTSSoS.sos_dualize(su2_result.moment_problem)
        set_optimizer(su2_sos.model, Clarabel.Optimizer)
        set_silent(su2_sos.model)
        optimize!(su2_sos.model)
        @test termination_status(su2_sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        su2_residual = sos_dual_certificate_residual(su2_result.moment_problem, su2_sos)
        @test su2_residual.max_abs_residual <= 1e-6
        su2_diagnostics = sos_dual_block_diagnostics(su2_result.moment_problem, su2_sos)
        @test [diag.label for diag in su2_diagnostics] == su2_result.report.block_labels
        @test [diag.native_value_size for diag in su2_diagnostics] ==
            [(size, size) for size in su2_result.report.psd_block_sizes]
        @test [
            diag.transform.family
            for diag in su2_diagnostics[end - 1:end]
        ] == [:pauli_su2_rdm, :pauli_su2_rdm]
        @test all(diag -> diag.native_hermitian_residual <= 1e-6, su2_diagnostics)
        @test all(diag -> isfinite(diag.native_min_eigenvalue), su2_diagnostics)
        @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
    end

    @testset "translation path appends linear state-opt rows" begin
        n = 8
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
        M = eltype(first(ops))

        mp, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            linear_state_opt_width=2,
        )

        @test count(con -> con[1] == :Zero, mp.constraints) > 0
        @test all(size(mat) == (1, 1) for (cone, mat) in mp.constraints if cone == :Zero)
        @test maximum(report.psd_block_sizes) > 0
        @test all(zc -> zc.origin.label.feature == :linear_state_opt, mp.linear.zero_constraints)
        @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing
        linear_metrics = translation_report_metrics(report)
        @test linear_metrics.linear_state_opt
        @test !linear_metrics.su2_moment_symmetry
        @test !linear_metrics.psd_state_opt
        @test linear_metrics.zero_constraint_count == length(mp.linear.zero_constraints)
        @test linear_metrics.linear_state_opt_row_count ==
            linear_metrics.zero_constraint_count
        @test linear_metrics.psd_state_opt_block_count == 0
        @test isempty(linear_metrics.psd_state_opt_logical_block_sizes)
        @test isempty(linear_metrics.psd_state_opt_psd_block_sizes)
        @test linear_metrics.zero_constraint_feature_histogram == [
            (feature=:linear_state_opt, decomposition=nothing, reason=nothing) =>
                linear_metrics.zero_constraint_count,
        ]
        su2_linear_mp, su2_linear_report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            su2_symmetry=true,
            linear_state_opt_width=2,
        )
        su2_linear_profile = translation_symmetry_profile(su2_linear_report)
        @test su2_linear_profile.su2_moment_symmetry
        @test !su2_linear_profile.su2_reflection_moment_symmetry
        @test su2_linear_profile.linear_state_opt
        su2_linear_metrics = translation_report_metrics(su2_linear_report)
        @test su2_linear_metrics.su2_moment_symmetry
        @test !su2_linear_metrics.su2_reflection_moment_symmetry
        @test su2_linear_metrics.linear_state_opt
        @test !su2_linear_metrics.psd_state_opt
        @test su2_linear_metrics.moment_equality_row_count == 0
        @test su2_linear_metrics.linear_state_opt_row_count > 0
        @test su2_linear_metrics.su2_base_zero_row_count > 0
        @test su2_linear_metrics.su2_base_spin_offblock_row_count > 0
        @test su2_linear_metrics.su2_base_magnetic_offdiag_row_count > 0
        @test su2_linear_metrics.su2_base_magnetic_copy_row_count > 0
        @test su2_linear_metrics.su2_base_zero_row_count +
            su2_linear_metrics.linear_state_opt_row_count ==
            su2_linear_metrics.zero_constraint_count
        su2_linear_zero_labels = [
            zc.origin.label for zc in su2_linear_mp.linear.zero_constraints
        ]
        @test all(
            label -> label.feature in (:pauli_su2_base_zero, :linear_state_opt),
            su2_linear_zero_labels,
        )
        @test any(label -> label.feature == :linear_state_opt, su2_linear_zero_labels)
        @test any(
            label -> label.feature == :pauli_su2_base_zero &&
                label.reason == :spin_offblock,
            su2_linear_zero_labels,
        )
        @test any(
            label -> label.feature == :pauli_su2_base_zero &&
                label.reason == :magnetic_copy,
            su2_linear_zero_labels,
        )
        @test translation_solve_support(su2_linear_report).supported
        @test NCTSSoS.assert_moment_linear_data_invariants(
            su2_linear_mp.linear,
            su2_linear_mp.constraints,
        ) === nothing

        qmb_lso_tests = NCTSSoS._qmbcertify_linear_state_opt_tests(
            M,
            20,
            7;
            sign_symmetry=true,
        )
        @test length(qmb_lso_tests) == 1817
        @test all(test -> NCTSSoS._pauli_sign_signature(test) == 0x00, qmb_lso_tests)

        qmb_linear, qmb_report = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            4;
            reflection_symmetry=true,
            sign_symmetry=true,
            real_moment_matrix=true,
            axis_rotation_symmetry=true,
            axis_rotation_equalities=true,
            axis_rotation_quotient=true,
            contiguous_rdm_k=8,
            contiguous_rdm_decomposition=:qmbcertify,
            contiguous_rdm_support=:extend,
            linear_state_opt_width=7,
            linear_state_opt_mode=:qmbcertify,
        )
        qmb_metrics = translation_report_metrics(qmb_report)
        @test qmb_metrics.linear_state_opt_row_count == 152
        @test length(qmb_linear.moments) == 346
        @test all(
            zero -> zero.origin.label.feature == :linear_state_opt &&
                zero.origin.label.mode == :qmbcertify,
            qmb_linear.zero_constraints,
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(qmb_linear) === nothing

        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            linear_state_opt_width=4,
        )

        registry4, ops4 = create_pauli_variables(1:4)
        pop4 = polyopt(heisenberg_chain_hamiltonian(ops4), registry4)
        linear_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                linear_state_opt_width=2,
            )
        end
        @test termination_status(linear_result.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test linear_result.report.zero_constraint_count > 0
        @test translation_zero_origin_histogram(linear_result.moment_problem) == [
            (feature=:linear_state_opt, decomposition=nothing, reason=nothing) =>
                linear_result.report.zero_constraint_count,
        ]
        direct_linear_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                linear_state_opt_width=2,
            )
        end
        @test direct_linear_state_result.moment_problem isa NCTSSoS.MomentLinearData
        @test termination_status(direct_linear_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test direct_linear_state_result.objective ≈ linear_result.objective atol = 1e-6
        @test direct_linear_state_result.report.block_labels == linear_result.report.block_labels
        @test direct_linear_state_result.report.zero_constraint_count ==
            linear_result.report.zero_constraint_count
        @test translation_zero_origin_histogram(direct_linear_state_result.moment_problem) ==
            translation_zero_origin_histogram(linear_result.moment_problem)
        direct_linear_state_sos =
            NCTSSoS.sos_dualize(direct_linear_state_result.moment_problem)
        set_optimizer(direct_linear_state_sos.model, Clarabel.Optimizer)
        set_silent(direct_linear_state_sos.model)
        optimize!(direct_linear_state_sos.model)
        @test termination_status(direct_linear_state_sos.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        direct_linear_state_residual = sos_dual_certificate_residual(
            direct_linear_state_result.moment_problem,
            direct_linear_state_sos,
        )
        @test direct_linear_state_residual.max_abs_residual <= 1e-6
        direct_linear_state_diagnostics = sos_dual_block_diagnostics(
            direct_linear_state_result.moment_problem,
            direct_linear_state_sos,
        )
        @test [diag.label for diag in direct_linear_state_diagnostics] ==
            direct_linear_state_result.report.block_labels
        su2_linear_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                linear_state_opt_width=2,
            )
        end
        reflection_linear_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                reflection_symmetry=true,
                linear_state_opt_width=2,
            )
        end
        su2_reflection_linear_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                reflection_symmetry=true,
                su2_symmetry=true,
                linear_state_opt_width=2,
            )
        end
        @test termination_status(su2_linear_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(reflection_linear_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(su2_reflection_linear_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test translation_solve_support(su2_linear_result).supported
        @test translation_solve_support(su2_reflection_linear_result).supported
        @test su2_linear_result.report.zero_constraint_count > linear_result.report.zero_constraint_count
        @test su2_linear_result.objective ≈ linear_result.objective atol = 1e-6
        @test su2_reflection_linear_result.objective ≈
            reflection_linear_result.objective atol = 1e-6
        linear_sos = NCTSSoS.sos_dualize(linear_result.moment_problem)
        set_optimizer(linear_sos.model, Clarabel.Optimizer)
        set_silent(linear_sos.model)
        optimize!(linear_sos.model)
        @test termination_status(linear_sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        linear_residual = sos_dual_certificate_residual(linear_result.moment_problem, linear_sos)
        @test linear_residual.max_abs_residual <= 1e-6
        linear_zero_values = sos_zero_dual_values(linear_result.moment_problem, linear_sos)
        @test length(linear_zero_values) == linear_result.report.zero_constraint_count
        @test all(zero -> zero.origin.label.feature == :linear_state_opt, linear_zero_values)
        @test all(zero -> isfinite(zero.value), linear_zero_values)
        linear_diagnostics = sos_dual_block_diagnostics(linear_result.moment_problem, linear_sos)
        @test [diag.label for diag in linear_diagnostics] == linear_result.report.block_labels
        @test [diag.native_value_size for diag in linear_diagnostics] ==
            [(size, size) for size in linear_result.report.psd_block_sizes]
        @test all(diag -> diag.native_hermitian_residual <= 1e-6, linear_diagnostics)
        @test all(diag -> isfinite(diag.native_min_eigenvalue), linear_diagnostics)
    end

    @testset "translation path appends PSD state-opt blocks" begin
        n = 8
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
        M = eltype(first(ops))
        K = typeof(symmetric_canon(NCTSSoS.expval(one(M))))

        mp, report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            psd_state_opt_width=1,
        )

        @test report.block_labels[end - 4:end] == [
            (feature=:psd_state_opt, width=1, momentum=0),
            (feature=:psd_state_opt, width=1, momentum=1),
            (feature=:psd_state_opt, width=1, momentum=2),
            (feature=:psd_state_opt, width=1, momentum=3),
            (feature=:psd_state_opt, width=1, momentum=4),
        ]
        @test report.logical_block_sizes[end - 4:end] == fill(3, 5)
        psd_state_metrics = translation_report_metrics(report)
        @test psd_state_metrics.psd_state_opt
        @test !psd_state_metrics.su2_moment_symmetry
        @test !psd_state_metrics.contiguous_rdm
        @test psd_state_metrics.linear_state_opt_row_count == 0
        @test psd_state_metrics.moment_equality_row_count == 0
        @test psd_state_metrics.axis_rotation_equality_row_count == 0
        @test psd_state_metrics.su2_base_zero_row_count == 0
        @test psd_state_metrics.psd_state_opt_block_count == 5
        @test psd_state_metrics.psd_state_opt_logical_block_sizes == fill(3, 5)
        @test psd_state_metrics.psd_state_opt_psd_block_sizes == [3, 6, 6, 6, 6]
        @test [
            length(block.meta.origin.logical_row_labels)
            for block in mp.linear.psd_blocks_lin[end - 4:end]
        ] == fill(3, 5)
        @test all(
            label -> label isa NormalMonomial{PauliAlgebra},
            mp.linear.psd_blocks_lin[end].meta.origin.logical_row_labels,
        )
        @test all(
            block -> block.meta.origin.transform.family == :translation_dft,
            mp.linear.psd_blocks_lin[end - 4:end],
        )
        @test all(
            block -> block.meta.origin.transform.coefficient_domain == :cyclotomic_float64,
            mp.linear.psd_blocks_lin[end - 4:end],
        )
        @test all(
            block -> block.meta.origin.transform.exact_coefficient_domain == :cyclotomic,
            mp.linear.psd_blocks_lin[end - 4:end],
        )
        @test all(cone == :PSD for (cone, _) in mp.constraints[end - 4:end])
        @test all(
            all(mat[i, j] == mat[j, i] for i in axes(mat, 1), j in axes(mat, 2)) for
            (_, mat) in mp.constraints[end - 4:end]
        )
        @test NCTSSoS.assert_moment_linear_data_invariants(mp.linear, mp.constraints) === nothing

        state_opt_basis = NCTSSoS._contiguous_state_opt_tests(M, n, 1; sign_symmetry=false)
        state_opt_translated = Dict(
            rep => [NCTSSoS._translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
            for rep in state_opt_basis
        )
        old_rep_cache = Dict{M,M}()
        uncached_block = NCTSSoS._translation_psd_state_opt_linear_block(
            state_opt_basis,
            2,
            n,
            pop.objective,
            state_opt_translated,
            old_rep_cache,
            K,
            ComplexF64,
        )
        cached_terms = NCTSSoS._translation_psd_state_opt_term_cache(
            state_opt_basis,
            n,
            pop.objective,
            state_opt_translated,
        )
        cached_block = NCTSSoS._translation_psd_state_opt_linear_block(
            state_opt_basis,
            2,
            n,
            cached_terms,
            K,
            ComplexF64,
        )
        @test _linear_blocks_isapprox(cached_block, uncached_block)
        @test Set(NCTSSoS._translation_psd_state_opt_required_moments(cached_terms)) ==
            NCTSSoS._psd_state_opt_required_moment_set(
                M,
                n,
                pop.objective,
                1;
                sign_symmetry=false,
            )

        su2_psd_mp, su2_psd_report = pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            su2_symmetry=true,
            psd_state_opt_width=1,
        )
        su2_psd_profile = translation_symmetry_profile(su2_psd_report)
        @test su2_psd_profile.su2_moment_symmetry
        @test !su2_psd_profile.su2_reflection_moment_symmetry
        @test su2_psd_profile.psd_state_opt
        su2_psd_metrics = translation_report_metrics(su2_psd_report)
        @test su2_psd_metrics.su2_moment_symmetry
        @test !su2_psd_metrics.su2_reflection_moment_symmetry
        @test su2_psd_metrics.psd_state_opt
        @test su2_psd_metrics.su2_base_zero_row_count > 0
        @test translation_solve_support(su2_psd_report).supported
        @test [
            block.meta.origin.label.feature
            for block in su2_psd_mp.linear.psd_blocks_lin[end - 4:end]
        ] == fill(:psd_state_opt, 5)

        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            2;
            sign_symmetry=false,
            psd_state_opt_width=3,
        )

        registry4, ops4 = create_pauli_variables(1:4)
        pop4 = polyopt(heisenberg_chain_hamiltonian(ops4), registry4)
        psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                psd_state_opt_width=1,
            )
        end
        @test termination_status(psd_state_result.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test psd_state_result.report.block_labels[end - 2:end] == [
            (feature=:psd_state_opt, width=1, momentum=0),
            (feature=:psd_state_opt, width=1, momentum=1),
            (feature=:psd_state_opt, width=1, momentum=2),
        ]
        direct_psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                psd_state_opt_width=1,
            )
        end
        @test direct_psd_state_result.moment_problem isa NCTSSoS.MomentLinearData
        @test termination_status(direct_psd_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test direct_psd_state_result.objective ≈ psd_state_result.objective atol = 1e-6
        @test direct_psd_state_result.report.block_labels == psd_state_result.report.block_labels
        @test direct_psd_state_result.report.psd_block_sizes ==
            psd_state_result.report.psd_block_sizes
        su2_psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                psd_state_opt_width=1,
            )
        end
        @test termination_status(su2_psd_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test translation_solve_support(su2_psd_state_result).supported
        @test translation_symmetry_profile(su2_psd_state_result).su2_moment_symmetry
        @test translation_symmetry_profile(su2_psd_state_result).psd_state_opt
        @test su2_psd_state_result.objective ≈ psd_state_result.objective atol = 1e-6
        direct_su2_psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                su2_symmetry=true,
                psd_state_opt_width=1,
            )
        end
        @test direct_su2_psd_state_result.moment_problem isa NCTSSoS.MomentLinearData
        @test termination_status(direct_su2_psd_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test translation_solve_support(direct_su2_psd_state_result).supported
        @test direct_su2_psd_state_result.objective ≈
            su2_psd_state_result.objective atol = 1e-6
        rdm_psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                contiguous_rdm_k=2,
                psd_state_opt_width=1,
            )
        end
        su2_rdm_psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                sign_symmetry=false,
                su2_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:su2,
                psd_state_opt_width=1,
            )
        end
        direct_su2_rdm_psd_state_result = quiet() do
            pauli_translation_invariant_nctssos(
                pop4,
                ops4,
                2,
                optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false);
                dualize=false,
                direct_linear=true,
                sign_symmetry=false,
                su2_symmetry=true,
                contiguous_rdm_k=2,
                contiguous_rdm_decomposition=:su2,
                psd_state_opt_width=1,
            )
        end
        @test termination_status(rdm_psd_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(su2_rdm_psd_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test termination_status(direct_su2_rdm_psd_state_result.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        su2_rdm_psd_profile = translation_symmetry_profile(su2_rdm_psd_state_result)
        @test su2_rdm_psd_profile.su2_moment_symmetry
        @test su2_rdm_psd_profile.su2_rdm_symmetry
        @test su2_rdm_psd_profile.psd_state_opt
        @test translation_solve_support(su2_rdm_psd_state_result).supported
        @test direct_su2_rdm_psd_state_result.moment_problem isa NCTSSoS.MomentLinearData
        @test translation_solve_support(direct_su2_rdm_psd_state_result).supported
        @test su2_rdm_psd_state_result.objective ≈
            rdm_psd_state_result.objective atol = 1e-6
        @test direct_su2_rdm_psd_state_result.objective ≈
            su2_rdm_psd_state_result.objective atol = 1e-6
        su2_rdm_psd_metrics =
            translation_report_metrics(su2_rdm_psd_state_result.report)
        @test su2_rdm_psd_metrics.su2_moment_symmetry
        @test su2_rdm_psd_metrics.su2_rdm_symmetry
        @test su2_rdm_psd_metrics.psd_state_opt
        @test su2_rdm_psd_metrics.su2_base_zero_row_count > 0
        @test su2_rdm_psd_metrics.contiguous_rdm_zero_row_count > 0
        @test su2_rdm_psd_metrics.psd_state_opt_block_count > 0
        @test su2_rdm_psd_metrics.contiguous_rdm_block_count == 2
        direct_su2_rdm_psd_state_sos =
            NCTSSoS.sos_dualize(direct_su2_rdm_psd_state_result.moment_problem)
        set_optimizer(direct_su2_rdm_psd_state_sos.model, Clarabel.Optimizer)
        set_silent(direct_su2_rdm_psd_state_sos.model)
        optimize!(direct_su2_rdm_psd_state_sos.model)
        @test termination_status(direct_su2_rdm_psd_state_sos.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        direct_su2_rdm_psd_state_residual = sos_dual_certificate_residual(
            direct_su2_rdm_psd_state_result.moment_problem,
            direct_su2_rdm_psd_state_sos,
        )
        @test direct_su2_rdm_psd_state_residual.max_abs_residual <= 1e-6
        direct_su2_rdm_psd_state_diagnostics = sos_dual_block_diagnostics(
            direct_su2_rdm_psd_state_result.moment_problem,
            direct_su2_rdm_psd_state_sos,
        )
        @test [diag.label for diag in direct_su2_rdm_psd_state_diagnostics] ==
            direct_su2_rdm_psd_state_result.report.block_labels
        @test [diag.native_value_size for diag in direct_su2_rdm_psd_state_diagnostics] ==
            [(size, size) for size in direct_su2_rdm_psd_state_result.report.psd_block_sizes]
        @test any(
            diag -> diag.label.feature == :contiguous_rdm &&
                diag.label.decomposition == :su2,
            direct_su2_rdm_psd_state_diagnostics,
        )
        @test any(
            diag -> diag.label.feature == :psd_state_opt,
            direct_su2_rdm_psd_state_diagnostics,
        )
        direct_psd_state_sos = NCTSSoS.sos_dualize(direct_psd_state_result.moment_problem)
        set_optimizer(direct_psd_state_sos.model, Clarabel.Optimizer)
        set_silent(direct_psd_state_sos.model)
        optimize!(direct_psd_state_sos.model)
        @test termination_status(direct_psd_state_sos.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        direct_psd_state_residual = sos_dual_certificate_residual(
            direct_psd_state_result.moment_problem,
            direct_psd_state_sos,
        )
        @test direct_psd_state_residual.max_abs_residual <= 1e-6
        direct_psd_state_diagnostics = sos_dual_block_diagnostics(
            direct_psd_state_result.moment_problem,
            direct_psd_state_sos,
        )
        @test [diag.label for diag in direct_psd_state_diagnostics] ==
            direct_psd_state_result.report.block_labels
        @test [diag.native_value_size for diag in direct_psd_state_diagnostics] ==
            [(size, size) for size in direct_psd_state_result.report.psd_block_sizes]
        psd_state_sos = NCTSSoS.sos_dualize(psd_state_result.moment_problem)
        set_optimizer(psd_state_sos.model, Clarabel.Optimizer)
        set_silent(psd_state_sos.model)
        optimize!(psd_state_sos.model)
        @test termination_status(psd_state_sos.model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        @test [
            block.meta.origin.label
            for block in sos_dual_blocks(psd_state_sos)[end - 2:end]
        ] == psd_state_result.report.block_labels[end - 2:end]
        psd_state_residual = sos_dual_certificate_residual(psd_state_result.moment_problem, psd_state_sos)
        @test psd_state_residual.max_abs_residual <= 1e-6
        psd_state_diagnostics = sos_dual_block_diagnostics(
            psd_state_result.moment_problem,
            psd_state_sos,
        )
        @test [
            diag.label
            for diag in psd_state_diagnostics[end - 2:end]
        ] == psd_state_result.report.block_labels[end - 2:end]
        @test [diag.native_value_size for diag in psd_state_diagnostics] ==
            [(size, size) for size in psd_state_result.report.psd_block_sizes]
        @test all(diag -> diag.native_hermitian_residual <= 1e-6, psd_state_diagnostics)
        @test all(diag -> isfinite(diag.native_min_eigenvalue), psd_state_diagnostics)
    end

    @testset "buffered momentum products own cached monomial keys" begin
        n = 4
        _, ops = create_pauli_variables(1:n)
        σx, σy, _ = ops
        basis = pauli_contiguous_chain_basis(ops, 1)
        orbit_reps = NCTSSoS._translation_orbit_representatives(basis, n)
        nontrivial_reps = [m for m in orbit_reps if !isone(m)]
        one_mono = one(first(σx))

        M = eltype(orbit_reps)
        translated = Dict{M,Vector{M}}()
        for rep in nontrivial_reps
            translated[rep] = [NCTSSoS._translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
        end
        translated[one_mono] = fill(one_mono, n)

        T = eltype(one_mono.word)
        rep_cache = Dict{M,M}()
        buf = T[]
        terms = NCTSSoS._translation_momentum_product_terms(
            σx[1],
            σy[1],
            n,
            translated,
            rep_cache,
            ComplexF64,
            buf,
        )
        legacy_terms = NCTSSoS._translation_momentum_product_terms(
            σx[1],
            σy[1],
            n,
            translated,
            Dict{M,M}(),
            ComplexF64,
        )
        direct_terms = Tuple{Int,ComplexF64,M}[]
        for r in 0:(n - 1)
            for (coef, mono) in (σx[1] * translated[σy[1]][r + 1]).terms
                rep = NCTSSoS._translation_orbit_representative(mono, n)
                push!(direct_terms, (r, ComplexF64(coef), rep))
            end
        end
        words_before = [copy(mono.word) for (_, _, mono) in terms]

        resize!(buf, max(length(buf), 1))
        fill!(buf, zero(T))

        @test !isempty(terms)
        @test terms == legacy_terms
        @test terms == direct_terms
        @test [mono.word for (_, _, mono) in terms] == words_before
    end

    @testset "momentum product cache is reused across sectors" begin
        n = 4
        _, ops = create_pauli_variables(1:n)
        σx, σy, _ = ops
        basis = pauli_contiguous_chain_basis(ops, 1)
        orbit_reps = NCTSSoS._translation_orbit_representatives(basis, n)
        nontrivial_reps = [m for m in orbit_reps if !isone(m)]
        one_mono = one(first(σx))

        M = eltype(orbit_reps)
        translated = Dict{M,Vector{M}}()
        for rep in nontrivial_reps
            translated[rep] = [NCTSSoS._translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
        end
        translated[one_mono] = fill(one_mono, n)

        T = eltype(one_mono.word)
        P = Polynomial{PauliAlgebra,T,ComplexF64}
        rep_cache = Dict{M,M}()
        owning_cache = NCTSSoS._TranslationProductCache(
            Dict{Tuple{M,M},Vector{Tuple{Int,ComplexF64,M}}}(),
            0,
            0,
        )
        scratch = owning_cache.buf
        NCTSSoS._translation_momentum_product_terms!(
            owning_cache,
            σx[1],
            σy[1],
            n,
            translated,
            rep_cache,
            ComplexF64,
        )
        NCTSSoS._translation_momentum_product_terms!(
            owning_cache,
            σx[1],
            σy[1],
            n,
            translated,
            rep_cache,
            ComplexF64,
        )
        @test owning_cache.buf === scratch
        @test owning_cache.misses == 1
        @test owning_cache.hits == 1

        product_cache = Dict{Tuple{M,M},Vector{Tuple{Int,ComplexF64,M}}}()

        block_k1 = NCTSSoS._translation_momentum_block(
            nontrivial_reps,
            1,
            n,
            translated,
            rep_cache,
            P,
            product_cache,
        )
        cache_size = length(product_cache)
        block_k2 = NCTSSoS._translation_momentum_block(
            nontrivial_reps,
            2,
            n,
            translated,
            rep_cache,
            P,
            product_cache,
        )

        @test cache_size > 0
        @test length(product_cache) == cache_size
        @test size(block_k1) == (length(nontrivial_reps), length(nontrivial_reps))
        @test size(block_k2) == size(block_k1)
    end

    @testset "small chain dense/reduced and qmb base route" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        hamiltonian = heisenberg_chain_hamiltonian(ops)
        pop = polyopt(hamiltonian, registry)

        dense = quiet() do
            cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=1); dualize=false)
        end
        reduced = quiet() do
            pauli_translation_invariant_nctssos(pop, ops, 1, SOLVER; dualize=false)
        end
        qmb_base = quiet() do
            pauli_translation_invariant_nctssos(
                pop,
                ops,
                1,
                SOLVER;
                dualize=false,
                qmbcertify_base_construct=true,
            )
        end

        @test termination_status(dense.model) == JuMP.MOI.OPTIMAL
        @test termination_status(reduced.model) == JuMP.MOI.OPTIMAL
        @test termination_status(qmb_base.model) == JuMP.MOI.OPTIMAL
        @test reduced.objective ≈ dense.objective atol = 1e-6
        @test qmb_base.moment_problem isa NCTSSoS.MomentLinearData
        @test length(qmb_base.moment_problem.moments) == 6
        @test qmb_base.report.psd_block_sizes == [4, 6, 3, 3, 6, 3]
        @test all(
            label -> label.feature == :qmbcertify_base,
            qmb_base.report.block_labels,
        )
        qmb_base_sos = NCTSSoS.sos_dualize(qmb_base.moment_problem)
        set_optimizer(qmb_base_sos.model, Clarabel.Optimizer)
        set_silent(qmb_base_sos.model)
        optimize!(qmb_base_sos.model)
        @test termination_status(qmb_base_sos.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        qmb_base_residual = sos_dual_certificate_residual(
            qmb_base.moment_problem,
            qmb_base_sos,
        )
        @test qmb_base_residual.max_abs_residual <= 1e-6
        qmb_base_diagnostics = sos_dual_block_diagnostics(
            qmb_base.moment_problem,
            qmb_base_sos,
        )
        @test [diag.label for diag in qmb_base_diagnostics] == qmb_base.report.block_labels
        @test [diag.native_value_size for diag in qmb_base_diagnostics] ==
            [(size, size) for size in qmb_base.report.psd_block_sizes]
        @test all(diag -> diag.label.feature == :qmbcertify_base, qmb_base_diagnostics)
        @test all(
            diag -> diag.transform_family == :translation_dft,
            qmb_base_diagnostics,
        )
        @test all(
            diag -> diag.coefficient_domain == :cyclotomic_float64,
            qmb_base_diagnostics,
        )
        @test all(
            diag -> diag.exact_coefficient_domain == :cyclotomic,
            qmb_base_diagnostics,
        )
        qmb_base_certificate = sos_dual_certificate_diagnostics(
            qmb_base.moment_problem,
            qmb_base_sos,
        )
        @test qmb_base_certificate.residual == qmb_base_residual
        @test qmb_base_certificate.psd_blocks == qmb_base_diagnostics
        @test qmb_base_certificate.psd_block_count == length(qmb_base.report.block_labels)
        qmb_scalar_pop = polyopt(
            hamiltonian,
            registry;
            eq_constraints=[hamiltonian],
            ineq_constraints=[one(hamiltonian) + hamiltonian],
        )
        qmb_scalar_linear, qmb_scalar_report =
            NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
                qmb_scalar_pop,
                ops,
                1,
            )
        qmb_scalar_sos = NCTSSoS.sos_dualize(qmb_scalar_linear)
        set_optimizer(qmb_scalar_sos.model, Clarabel.Optimizer)
        set_silent(qmb_scalar_sos.model)
        optimize!(qmb_scalar_sos.model)
        @test termination_status(qmb_scalar_sos.model) in
            (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL)
        qmb_scalar_residual = sos_dual_certificate_residual(
            qmb_scalar_linear,
            qmb_scalar_sos,
        )
        @test qmb_scalar_residual.max_abs_residual <= 1e-6
        qmb_scalar_block_diagnostics = sos_dual_block_diagnostics(
            qmb_scalar_linear,
            qmb_scalar_sos,
        )
        qmb_scalar_zero_diagnostics = sos_zero_dual_diagnostics(
            qmb_scalar_linear,
            qmb_scalar_sos,
        )
        @test qmb_scalar_report.block_labels[end] == (feature=:scalar_inequality, index=1)
        @test last(qmb_scalar_block_diagnostics).label ==
            (feature=:scalar_inequality, index=1)
        @test only(qmb_scalar_zero_diagnostics).label ==
            (feature=:scalar_equality, index=1)
        @test only(qmb_scalar_zero_diagnostics).feature == :scalar_equality
        qmb_scalar_certificate = sos_dual_certificate_diagnostics(
            qmb_scalar_linear,
            qmb_scalar_sos,
        )
        @test qmb_scalar_certificate.residual == qmb_scalar_residual
        @test qmb_scalar_certificate.psd_blocks == qmb_scalar_block_diagnostics
        @test qmb_scalar_certificate.zero_duals == qmb_scalar_zero_diagnostics
        @test qmb_scalar_certificate.psd_block_count == length(qmb_scalar_report.block_labels)
        @test qmb_scalar_certificate.zero_dual_count == 1
        @test maximum(reduced.report.psd_block_sizes) < only(flatten_sizes(dense.moment_matrix_sizes))
    end
end

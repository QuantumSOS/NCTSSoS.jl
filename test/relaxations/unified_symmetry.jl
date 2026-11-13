# test/relaxations/unified_symmetry.jl
# Tests: the unified `Symmetry` front reuses the existing backends unchanged,
# reports their checks, and rejects unsupported compositions.

using Test, NCTSSoS, LinearAlgebra
using JuMP, COSMO

if !@isdefined(SOLVER)
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

# Rates of change of the rotations of each listed vector of letters.
function unified_vector_rotations(reg, vectors)
    n = length(first(vectors))
    z = zero(Polynomial(first(first(vectors))))
    return [SuppliedDerivation(reg, Dict(v[i].word[1] =>
        (i == a ? Polynomial(v[b]) : i == b ? -Polynomial(v[a]) : z)
        for v in vectors for i in 1:n)) for a in 1:n for b in a+1:n]
end

function unified_direct_finite_report(pop, config, spec)
    direct_config = SolverConfig(optimizer=config.optimizer, order=config.order,
        moment_basis=config.moment_basis, symmetry=spec)
    sparsity = compute_sparsity(pop, direct_config)
    _, report = NCTSSoS.moment_relax_symmetric(pop, sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities, spec)
    return report
end

@testset "Unified symmetry front" begin
    @testset "finite route equals the SymmetrySpec route" begin
        # CHSH: the value -2sqrt(2) is Tsirelson's bound, independent of this code.
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        chsh = polyopt(-(1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]), reg)
        generators = [
            SignedPermutation(x[1].word[1] => x[2].word[1], x[2].word[1] => x[1].word[1],
                y[2].word[1] => (-1, y[2].word[1])),
            SignedPermutation(x[2].word[1] => (-1, x[2].word[1]), y[1].word[1] => y[2].word[1],
                y[2].word[1] => y[1].word[1]),
            SignedPermutation(x[1].word[1] => y[1].word[1], x[2].word[1] => y[2].word[1],
                y[1].word[1] => x[1].word[1], y[2].word[1] => x[2].word[1]),
        ]
        spec = SymmetrySpec(generators)
        sym = Symmetry(finite=generators)
        config = SolverConfig(optimizer=SOLVER, order=1)

        round_trip = SymmetrySpec(sym)
        @test round_trip.generators == spec.generators
        @test round_trip.check_invariance && round_trip.offblock_check === :randomized
        @test SymmetrySpec(Symmetry(spec)).generators == spec.generators

        direct_report = unified_direct_finite_report(chsh, config, spec)
        report = check_symmetry(chsh, config, sym)
        @test report.route === :finite
        @test report.backend === :moment_relax_symmetric
        @test report.backend_report.group_order == direct_report.group_order == 16
        @test report.psd_block_sizes == direct_report.psd_block_sizes
        @test report.invariant_moment_count == direct_report.invariant_moment_count
        @test maximum(report.psd_block_sizes) < direct_report.basis_half_size
        @test (report.finite_generators, report.continuous_generators) == (3, 0)
        @test report.invariance === :exact
        @test report.block_arithmetic === :float64
        @test report.block_verification === :numeric
        @test report.relation_to_input === :equivalent_numeric
        @test report.compactness === :not_applicable
        @test isempty(report.notes)

        direct = cs_nctssos(chsh, SolverConfig(optimizer=SOLVER, order=1, symmetry=spec))
        solved = cs_nctssos(chsh, config, sym)
        @test solved.result isa PolyOptResult
        @test solved.objective ≈ direct.objective atol = 1e-6
        @test solved.objective ≈ -2sqrt(2) atol = 1e-5
        @test solved.report.psd_block_sizes == direct.symmetry.psd_block_sizes

        # Constrained problem: |x1 x2 + x2 x1| <= x1^2 + x2^2 <= 1, so the minimum is -1.
        free, (w,) = create_noncommutative_variables([("w", 1:2)])
        ball = polyopt(1.0 * (w[1] * w[2] + w[2] * w[1]), free;
            ineq_constraints=[1.0 - w[1]^2 - w[2]^2])
        swap = [SignedPermutation(w[1].word[1] => w[2].word[1], w[2].word[1] => w[1].word[1])]
        plain = cs_nctssos(ball, config)
        reduced = cs_nctssos(ball, config, Symmetry(finite=swap))
        @test plain.objective ≈ -1 atol = 1e-5
        @test reduced.objective ≈ plain.objective atol = 1e-5
        @test reduced.report.psd_block_sizes ==
            unified_direct_finite_report(ball, config, SymmetrySpec(swap)).psd_block_sizes
        @test maximum(reduced.report.psd_block_sizes) < 3

        # A Clifford map is a finite generator; it takes the same route.
        pauli, (sx, sy, sz) = create_pauli_variables(1:2)
        heisenberg = polyopt(sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (sx, sy, sz)), pauli)
        basis_config = SolverConfig(optimizer=nothing, moment_basis=[one(sx[1]); sx; sy; sz])
        clifford_report = check_symmetry(heisenberg, basis_config,
            Symmetry(finite=[CliffordSymmetry(:SWAP, 1, 2)]))
        @test clifford_report.route === :finite
        @test clifford_report.psd_block_sizes == [4, 3]
    end

    @testset "disabled checks and asserted options are reported" begin
        n = 4
        reg, (sx, sy, sz) = create_pauli_variables(1:n)
        ring = polyopt(sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, n)]
            for op in (sx, sy, sz) for i in 1:n), reg)
        config = SolverConfig(optimizer=nothing, order=1)
        translation = pauli_site_permutation([2, 3, 4, 1])

        checked = check_symmetry(ring, config, Symmetry(finite=[translation]))
        @test checked.relation_to_input === :equivalent_numeric

        unchecked_blocks = check_symmetry(ring, config, Symmetry(finite=[translation], offblock_check=:off))
        @test unchecked_blocks.block_verification === :unverified
        @test unchecked_blocks.relation_to_input === :unverified
        @test unchecked_blocks.psd_block_sizes == checked.psd_block_sizes
        @test any(occursin("offblock_check=:off", note) for note in unchecked_blocks.notes)

        unchecked_problem = check_symmetry(ring, config, Symmetry(finite=[translation], check_invariance=false))
        @test unchecked_problem.invariance === :unverified
        @test unchecked_problem.relation_to_input === :unverified

        charge = PauliChargeSectorSpec(nqubits=n)
        charge_report = check_symmetry(ring, SolverConfig(optimizer=nothing, order=2),
            Symmetry(pauli_charge=charge))
        @test charge_report.psd_block_sizes == [6, 16, 23, 16, 6]
        @test charge_report.options == [:pauli_charge]
        @test charge_report.invariance === :numeric
        @test charge_report.relation_to_input === :equivalent_numeric

        singlet = Symmetry(finite=[translation], pauli_charge=charge,
            pauli_singlet=PauliSingletConstraintSpec(nqubits=n))
        singlet_report = check_symmetry(ring, SolverConfig(optimizer=nothing, order=2), singlet)
        @test singlet_report.options == [:pauli_charge, :pauli_singlet]
        @test singlet_report.relation_to_input === :unverified
        @test any(occursin("pauli_singlet", note) for note in singlet_report.notes)
        @test occursin("unverified", sprint(show, singlet_report))
    end

    @testset "continuous route equals continuous_nctssos" begin
        reg, (x, y, z) = create_pauli_variables(1:2)
        actions = unified_vector_rotations(reg, [[x[i], y[i], z[i]] for i in 1:2])
        metric = Matrix{Int}(I, 6, 6)
        pop = polyopt(x[1] * x[2] + y[1] * y[2] + z[1] * z[2], reg)
        config = SolverConfig(optimizer=SOLVER, order=1)
        sym = Symmetry(continuous=actions, compact_metric=metric)

        report = check_symmetry(pop, config, sym)
        @test report.route === :continuous
        @test report.backend === :continuous_moment_relaxation
        @test report.psd_block_sizes == [1, 2]
        @test report.invariant_moment_count == 2
        @test report.compactness === :verified
        @test report.block_arithmetic === :complex_rationals
        @test report.block_verification === :exact
        @test report.relation_to_input === :equivalent
        @test report.joint_blocks === :not_applicable

        # The singlet energy of two spins is -3. The solver's dual is checked by its
        # exact residual in the original algebra, without chopping.
        direct = continuous_nctssos(pop, config, actions; compact_metric=metric, dualize=true)
        solved = cs_nctssos(pop, config, sym; dualize=true)
        @test solved.objective ≈ direct.objective atol = 1e-6
        @test solved.objective ≈ -3 atol = 2e-5
        @test Float64(solved.result.certificate.residual_coefficient_bound) < 2e-5
        # Independent identity h + 3 = sum_axis (sigma_1 + sigma_2)^2 / 2, certified
        # exactly through the reduction the front returns.
        @test solved.report.backend_report === solved.result.reduction
        exact = verify_continuous_certificate(solved.report.backend_report, -3,
            [zeros(Rational{BigInt}, 1, 1), fill(big(3) // 2, 2, 2)])
        @test exact.identity_exact && exact.grams_psd_exact
        @test exact.certified_bound == -3
        @test !verify_continuous_certificate(solved.report.backend_report, -3,
            [zeros(Int, 1, 1), ones(Int, 2, 2)]).identity_exact

        # Primal solve: each original matrix is rebuilt from its blocks.
        primal = cs_nctssos(pop, config, sym)
        @test primal.objective ≈ -3 atol = 2e-5
        for (idx, mat) in primal.result.matrices
            rebuilt = map(p -> sum(c * primal.result.moments[m] for (c, m) in terms(p); init=0),
                primal.result.reduction.original.constraints[idx][2])
            @test mat ≈ rebuilt atol = 1e-12
        end

        moments_only = check_symmetry(pop, config, Symmetry(continuous=actions, compact_metric=metric, compress=false))
        @test moments_only.psd_block_sizes == [7]
        @test any(occursin("compress=false", note) for note in moments_only.notes)
    end

    @testset "finite and continuous generators compose on the checked Pauli route" begin
        reg, (x, y, z) = create_pauli_variables(1:3)
        actions = unified_vector_rotations(reg, [[x[i], y[i], z[i]] for i in 1:3])
        metric = Matrix{Int}(I, 9, 9)
        bond(i, j) = x[i] * x[j] + y[i] * y[j] + z[i] * z[j]
        triangle = polyopt(bond(1, 2) + bond(2, 3) + bond(1, 3), reg)
        config = SolverConfig(optimizer=SOLVER, order=1)
        cycle = pauli_site_permutation([2, 3, 1])
        sym = Symmetry(finite=[cycle], continuous=actions, compact_metric=metric)

        # Degree one gives -4.5 while the ground energy is -3: a non-tight control,
        # so agreement is not explained by the bound being exact.
        plain = cs_nctssos(triangle, config; dualize=false)
        direct = continuous_nctssos(triangle, config, actions; compact_metric=metric, clifford_generators=[cycle])
        solved = cs_nctssos(triangle, config, sym)
        @test plain.objective ≈ -4.5 atol = 2e-5
        @test solved.objective ≈ plain.objective atol = 2e-5
        @test solved.objective ≈ direct.objective atol = 1e-6
        @test solved.report.route === :finite_continuous
        @test solved.report.joint_blocks === :continuous_only
        @test solved.report.relation_to_input === :equivalent
        @test solved.report.psd_block_sizes ==
            [b.size for b in direct.reduction.moment_problem.linear.psd_blocks_lin]
        @test any(occursin("coarser", note) for note in solved.report.notes)
        for (m, value) in solved.result.moments
            sign, image = NCTSSoS._act_monomial(cycle, m)
            @test value ≈ sign * solved.result.moments[image] atol = 1e-12
        end

        # The finite generator removes moment coordinates the rotations alone keep.
        continuous_only = check_symmetry(triangle, config, Symmetry(continuous=actions, compact_metric=metric))
        @test solved.report.invariant_moment_count < continuous_only.invariant_moment_count

        # The normalizer requirement is enforced by the backend, not waived here.
        axes = pauli_global_axis_rotation_generators((x, y, z))
        @test check_symmetry(triangle, config,
            Symmetry(finite=axes, continuous=actions, compact_metric=metric)).route === :finite_continuous
        @test_throws ArgumentError check_symmetry(triangle, config,
            Symmetry(finite=axes, continuous=actions[1:1], compact_metric=metric))
    end

    @testset "unsupported compositions are rejected" begin
        reg, (x, y, z) = create_pauli_variables(1:3)
        actions = unified_vector_rotations(reg, [[x[i], y[i], z[i]] for i in 1:3])
        metric = Matrix{Int}(I, 9, 9)
        bond(i, j) = x[i] * x[j] + y[i] * y[j] + z[i] * z[j]
        chain = bond(1, 2) + bond(2, 3)
        config = SolverConfig(optimizer=nothing, order=1)
        reflection = pauli_site_permutation([3, 2, 1])
        free, (w,) = create_noncommutative_variables([("w", 1:2)])
        swap = SignedPermutation(w[1].word[1] => w[2].word[1], w[2].word[1] => w[1].word[1])

        antiunitary = try
            Symmetry(finite=[PauliAntiunitary([0, 0, 0], [0, 0, 0])])
        catch err
            err
        end
        @test antiunitary isa ArgumentError
        @test occursin("antiunitary_moment_relaxation", antiunitary.msg)

        @test_throws ArgumentError Symmetry()
        @test_throws ArgumentError Symmetry(finite=[reflection, swap])
        @test_throws ArgumentError Symmetry(finite=[:reflection])
        @test_throws ArgumentError Symmetry(finite=[swap], continuous=actions, compact_metric=metric)
        @test_throws ArgumentError Symmetry(continuous=actions, compact_metric=metric,
            pauli_charge=PauliChargeSectorSpec(nqubits=3))
        @test_throws ArgumentError Symmetry(continuous=actions, compact_metric=metric, check_invariance=false)
        @test_throws ArgumentError Symmetry(continuous=actions, compact_metric=metric, offblock_check=:off)
        @test_throws ArgumentError Symmetry(finite=[reflection], compact_metric=metric)
        @test_throws ArgumentError Symmetry(finite=[reflection], compress=false)
        @test_throws ArgumentError Symmetry(finite=[reflection], offblock_check=:sometimes)
        # Existing `SymmetrySpec` combination rules apply at construction.
        @test_throws ArgumentError Symmetry(finite=[swap], pauli_charge=PauliChargeSectorSpec(nqubits=3))
        @test_throws ArgumentError SymmetrySpec(Symmetry(continuous=actions, compact_metric=metric))

        pop = polyopt(chain, reg)
        @test_throws ArgumentError check_symmetry(pop, config, Symmetry(continuous=actions))
        @test_throws ArgumentError check_symmetry(polyopt(chain, reg; eq_constraints=[bond(1, 3) - 1]), config,
            Symmetry(finite=[reflection], continuous=actions, compact_metric=metric))
        doubled = SolverConfig(optimizer=nothing, order=1, symmetry=SymmetrySpec(reflection))
        @test_throws ArgumentError check_symmetry(pop, doubled, Symmetry(finite=[reflection]))
        @test_throws ArgumentError cs_nctssos(pop, doubled, Symmetry(finite=[reflection]))
        @test_throws ArgumentError cs_nctssos(pop, config, Symmetry(continuous=actions, compact_metric=metric);
            compress=false)

        # A map that is not a symmetry of the problem is refused on every route.
        asymmetric = polyopt(bond(1, 2) + 2 * bond(2, 3), reg)
        @test_throws ArgumentError check_symmetry(asymmetric, config, Symmetry(finite=[reflection]))
        @test_throws ArgumentError check_symmetry(asymmetric, config,
            Symmetry(finite=[reflection], continuous=actions, compact_metric=metric))
        @test_throws ArgumentError check_symmetry(polyopt(chain + z[1], reg), config,
            Symmetry(continuous=actions, compact_metric=metric))
        # The finite route keeps its dense-relaxation requirement.
        @test_throws ArgumentError check_symmetry(pop,
            SolverConfig(optimizer=nothing, order=1, ts_algo=MMD()), Symmetry(finite=[reflection]))
    end
end

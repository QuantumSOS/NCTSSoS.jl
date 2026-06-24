# Fermionic symmetry-specific regressions
#
# Coverage:
# - Abelian orbital-irrep sector splitting on a finite fermionic symmetry path
# - H₂/STO-3G molecular-orbital labels and N=2 spin-adapted order-2 regression
# - SU(2) spin-adapted block structure on the Be 1s/2s worked example from the note
# - low-cost Be-derived numerical slice using reviewed STO-6G one-electron integrals

using Test, NCTSSoS, LinearAlgebra, TOML

const FERMIONIC_SYMMETRY_EXPECTATIONS_PATH = "expectations/fermionic_symmetry.toml"

const _FERMIONIC_SYMM_Z = ComplexF64[1 0; 0 -1]
const _FERMIONIC_SYMM_I = Matrix{ComplexF64}(I, 2, 2)
const _FERMIONIC_SYMM_SP = ComplexF64[0 1; 0 0]
const _FERMIONIC_SYMM_SM = ComplexF64[0 0; 1 0]

_z2_irrep_multiply(a, b) = a === b ? :A : :B
_h2_irrep_multiply(a, b) = a === b ? :Ag : :B1u

function _expectation_case(relpath::AbstractString, id::AbstractString)
    data = TOML.parsefile(joinpath(pkgdir(NCTSSoS), "test", "data", relpath))
    return only(filter(case -> case["id"] == id, data["cases"]))
end

function _fermionic_symmetry_kron_all(mats::AbstractVector{<:AbstractMatrix})
    isempty(mats) && return Matrix{ComplexF64}(I, 1, 1)
    return reduce(kron, mats)
end

function _fermionic_symmetry_jw_op(mode::Int, kind::Symbol, nmodes::Int)
    mats = Vector{Matrix{ComplexF64}}(undef, nmodes)
    for site in 1:nmodes
        if site < mode
            mats[site] = _FERMIONIC_SYMM_Z
        elseif site == mode
            mats[site] = kind === :annihilate ? _FERMIONIC_SYMM_SM : _FERMIONIC_SYMM_SP
        else
            mats[site] = _FERMIONIC_SYMM_I
        end
    end
    return _fermionic_symmetry_kron_all(mats)
end

function _fermionic_symmetry_word_matrix(word::AbstractVector{<:Integer}, nmodes::Int)
    dim = 2^nmodes
    acc = Matrix{ComplexF64}(I, dim, dim)
    for op in word
        mode = abs(Int(op))
        kind = op > 0 ? :annihilate : :create
        acc = acc * _fermionic_symmetry_jw_op(mode, kind, nmodes)
    end
    return acc
end

function _fermionic_symmetry_poly_matrix(p::Polynomial{FermionicAlgebra,T,C}, nmodes::Int) where {T<:Integer,C<:Number}
    dim = 2^nmodes
    acc = zeros(ComplexF64, dim, dim)
    for (coef, mono) in zip(coefficients(p), monomials(p))
        acc .+= ComplexF64(coef) * _fermionic_symmetry_word_matrix(mono.word, nmodes)
    end
    return acc
end

function _fermionic_symmetry_sector_ground_state(H::AbstractMatrix{<:Complex}, nelectrons::Int)
    keep = [idx for idx in 1:size(H, 1) if count_ones(UInt(idx - 1)) == nelectrons]
    return eigmin(Hermitian(H[keep, keep]))
end

function _h2_sto3g_layout(c_up, c_dn)
    up_modes = [Int(op.word[1]) for op in c_up]
    dn_modes = [Int(op.word[1]) for op in c_dn]
    return FermionicModeLayout(
        Dict(vcat(
            [up_modes[i] => i for i in eachindex(up_modes)],
            [dn_modes[i] => i for i in eachindex(dn_modes)],
        ));
        spin2_of=Dict(vcat(
            [up_modes[i] => 1 for i in eachindex(up_modes)],
            [dn_modes[i] => -1 for i in eachindex(dn_modes)],
        )),
        irrep_of=Dict(1 => :Ag, 2 => :B1u),
        irrep_table=AbelianIrrepTable(:Ag; multiply=_h2_irrep_multiply, dual=identity),
    )
end

function _h2_sto3g_integrals()
    # Reviewed offline molecular-orbital constants for the standard minimal-basis
    # H₂ example with two spatial orbitals (σ_g, σ_u). These are electronic-only
    # integrals; no nuclear-repulsion constant is added to the polynomial.
    h = [
        -1.252477303982  0.0;
         0.0            -0.475934275367;
    ]

    v = zeros(2, 2, 2, 2)
    v[1, 1, 1, 1] = 0.674493166046
    v[2, 2, 2, 2] = 0.69739794957
    v[1, 2, 1, 2] = 0.663472101055
    v[2, 1, 2, 1] = 0.663472101055
    v[1, 2, 2, 1] = 0.181287518306
    v[2, 1, 1, 2] = 0.181287518306
    v[1, 1, 2, 2] = 0.181287518306
    v[2, 2, 1, 1] = 0.181287518306

    return h, v
end

function _h2_sto3g_hamiltonian(c_up, c_up_dag, c_dn, c_dn_dag)
    h, v = _h2_sto3g_integrals()
    ham = zero(1.0 * one(c_up[1]))

    for p in 1:2, q in 1:2
        abs(h[p, q]) <= 1e-12 && continue
        ham += h[p, q] * (c_up_dag[p] * c_up[q] + c_dn_dag[p] * c_dn[q])
    end

    for p in 1:2, q in 1:2, r in 1:2, s in 1:2
        abs(v[p, q, r, s]) <= 1e-12 && continue
        for (cσ, cσ_dag) in ((c_up, c_up_dag), (c_dn, c_dn_dag)),
            (cτ, cτ_dag) in ((c_up, c_up_dag), (c_dn, c_dn_dag))
            ham += 0.5 * v[p, q, r, s] * (cσ_dag[p] * cτ_dag[q] * cτ[s] * cσ[r])
        end
    end

    return ham
end

function _h2_sto3g_problem()
    registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) =
        create_fermionic_variables([("c_up", 1:2), ("c_dn", 1:2)])

    ham = _h2_sto3g_hamiltonian(c_up, c_up_dag, c_dn, c_dn_dag)
    layout = _h2_sto3g_layout(c_up, c_dn)
    sector = FermionicSectorSpec(
        mode_layout=layout,
        split_parity=true,
        split_number=true,
        split_spin=true,
        split_irrep=true,
    )
    spin = FermionicSpinAdaptationSpec(mode_layout=layout)
    n_total = sum(c_up_dag[i] * c_up[i] + c_dn_dag[i] * c_dn[i] for i in 1:2)
    pop = polyopt(ham, registry; moment_eq_constraints=[n_total - 2.0 * one(ham)])
    exact = _fermionic_symmetry_sector_ground_state(_fermionic_symmetry_poly_matrix(ham, 4), 2)

    return pop, sector, spin, exact
end

function _be_1s_2s_layout(c_up, c_dn)
    up_modes = [Int(op.word[1]) for op in c_up]
    dn_modes = [Int(op.word[1]) for op in c_dn]
    return FermionicModeLayout(
        Dict(vcat(
            [up_modes[i] => i for i in eachindex(up_modes)],
            [dn_modes[i] => i for i in eachindex(dn_modes)],
        ));
        spin2_of=Dict(vcat(
            [up_modes[i] => 1 for i in eachindex(up_modes)],
            [dn_modes[i] => -1 for i in eachindex(dn_modes)],
        )),
        irrep_of=Dict(1 => :Ag, 2 => :Ag),
        irrep_table=AbelianIrrepTable(:Ag; multiply=(a, b) -> :Ag, dual=identity),
    )
end

function _be_1s_2s_one_body_problem()
    h = [
        -7.9386904776328775  0.27373827553230023;
         0.27373827553230023 -1.7707398712114135;
    ]

    registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) =
        create_fermionic_variables([("c_up", 1:2), ("c_dn", 1:2)])

    ham = zero(1.0 * one(c_up[1]))
    for p in 1:2, q in 1:2
        ham += h[p, q] * (c_up_dag[p] * c_up[q] + c_dn_dag[p] * c_dn[q])
    end

    basis = [one(c_up[1])]
    for (c, c_dag) in ((c_up, c_up_dag), (c_dn, c_dn_dag))
        for p in 1:2, q in 1:2
            push!(basis, only(monomials(c_dag[p] * c[q])))
        end
    end

    layout = _be_1s_2s_layout(c_up, c_dn)
    sector = FermionicSectorSpec(
        mode_layout=layout,
        split_parity=true,
        split_number=true,
        split_spin=true,
        split_irrep=true,
    )
    spin = FermionicSpinAdaptationSpec(mode_layout=layout)
    n_total = sum(c_up_dag[i] * c_up[i] + c_dn_dag[i] * c_dn[i] for i in 1:2)
    pop = polyopt(ham, registry; moment_eq_constraints=[n_total - 2.0 * one(ham)])
    exact = 2 * minimum(eigvals(Symmetric(h)))

    return pop, basis, sector, spin, exact
end

@testset "Fermionic symmetry-specific regressions" begin
    @testset "Abelian orbital-irrep sectors agree with the ordinary path" begin
        reg, (a, a_dag) = create_fermionic_variables(1:4)
        layout = FermionicModeLayout(
            Dict(i => i for i in 1:4);
            irrep_of=Dict(1 => :B, 2 => :B, 3 => :A, 4 => :A),
            irrep_table=AbelianIrrepTable(:A; multiply=_z2_irrep_multiply, dual=identity),
        )
        symmetry = SymmetrySpec(
            fermionic_generators=[
                FermionicModePermutation(1 => 2, 2 => 1),
                FermionicModePermutation(3 => 4, 4 => 3),
            ],
            sector=FermionicSectorSpec(
                mode_layout=layout,
                split_parity=true,
                split_number=true,
                split_irrep=true,
            ),
        )
        basis = [one(a[1]), a[1], a[2], a[3], a[4], a_dag[1], a_dag[2], a_dag[3], a_dag[4]]
        ham = -(a_dag[1] * a[2] + a_dag[2] * a[1] + a_dag[3] * a[4] + a_dag[4] * a[3])
        pop = polyopt(ham, reg)

        plain = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                moment_basis=basis,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
            ),
        )
        symmetric = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                moment_basis=basis,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
                symmetry=symmetry,
            ),
        )

        @test symmetric.objective ≈ plain.objective atol = 1e-6
        @test !isnothing(symmetric.symmetry)
        @test maximum(symmetric.symmetry.psd_block_sizes) == 1
        @test any(label -> label isa FermionicSectorLabel && label.irrep == :A, symmetric.symmetry.block_labels)
        @test any(label -> label isa FermionicSectorLabel && label.irrep == :B, symmetric.symmetry.block_labels)
        @test all(label isa FermionicSectorLabel for label in symmetric.symmetry.block_labels)
    end

    @testset "H₂/STO-3G irrep labels distinguish densities from transitions" begin
        case = _expectation_case(FERMIONIC_SYMMETRY_EXPECTATIONS_PATH, "h2_sto3g_sector_labels")
        expected = case["expected"]
        reg, ((c_up, c_up_dag), (c_dn, c_dn_dag)) =
            create_fermionic_variables([("c_up", 1:2), ("c_dn", 1:2)])
        sector = FermionicSectorSpec(
            mode_layout=_h2_sto3g_layout(c_up, c_dn),
            split_parity=true,
            split_number=true,
            split_spin=true,
            split_irrep=true,
        )

        density_g = NCTSSoS._fermionic_sector_label(only(monomials(c_up_dag[1] * c_up[1])), sector)
        density_u = NCTSSoS._fermionic_sector_label(only(monomials(c_dn_dag[2] * c_dn[2])), sector)
        transition_gu = NCTSSoS._fermionic_sector_label(only(monomials(c_up_dag[1] * c_up[2])), sector)
        transition_ug = NCTSSoS._fermionic_sector_label(only(monomials(c_dn_dag[2] * c_dn[1])), sector)

        diagonal_irrep = Symbol(expected["diagonal_irrep"])
        transition_irrep = Symbol(expected["transition_irrep"])

        @test density_g == FermionicSectorLabel(0, 0, 0, diagonal_irrep)
        @test density_u == FermionicSectorLabel(0, 0, 0, diagonal_irrep)
        @test transition_gu == FermionicSectorLabel(0, 0, 0, transition_irrep)
        @test transition_ug == FermionicSectorLabel(0, 0, 0, transition_irrep)
    end

    @testset "H₂/STO-3G N=2 spin-adapted relaxation matches exact diagonalization" begin
        oracle = expectations_oracle(FERMIONIC_SYMMETRY_EXPECTATIONS_PATH, "h2_sto3g_n2_spin_adapted")
        pop, sector, spin, exact = _h2_sto3g_problem()

        plain = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
            ),
        )
        symmetric = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
                symmetry=SymmetrySpec(sector=sector, spin_adaptation=spin),
            ),
        )

        @test exact ≈ oracle.opt atol = 1e-10
        @test plain.objective ≈ exact atol = 5e-6
        @test symmetric.objective ≈ exact atol = 5e-6
        @test symmetric.objective ≈ plain.objective atol = 5e-6
        @test !isnothing(symmetric.symmetry)
        @test sort(symmetric.symmetry.psd_block_sizes) == sort(oracle.sides)
        @test maximum(symmetric.symmetry.psd_block_sizes) < only(flatten_sizes(plain.moment_matrix_sizes))
        @test 0 in Set(label.total_spin2 for label in symmetric.symmetry.block_labels if label isa FermionicSpinBlockLabel)
        @test 2 in Set(label.total_spin2 for label in symmetric.symmetry.block_labels if label isa FermionicSpinBlockLabel)
        @test Set(label.sector.irrep for label in symmetric.symmetry.block_labels if label isa FermionicSpinBlockLabel) == Set([:Ag, :B1u])
    end

    @testset "Be 1s/2s worked example reproduces the [2, 2] spin blocks" begin
        oracle = expectations_oracle(FERMIONIC_SYMMETRY_EXPECTATIONS_PATH, "be_1s_2s_spin_blocks")
        reg, ((c_up, c_up_dag), (c_dn, c_dn_dag)) =
            create_fermionic_variables([("c_up", 1:2), ("c_dn", 1:2)])
        layout = _be_1s_2s_layout(c_up, c_dn)
        sector = FermionicSectorSpec(
            mode_layout=layout,
            split_parity=true,
            split_number=true,
            split_spin=true,
            split_irrep=true,
        )
        spin = FermionicSpinAdaptationSpec(mode_layout=layout)
        active_modes = sort!(vcat([Int(op.word[1]) for op in c_up], [Int(op.word[1]) for op in c_dn]))
        density_basis = [
            only(monomials(c_up_dag[i] * c_up[i])) for i in 1:2
        ]
        append!(density_basis, [only(monomials(c_dn_dag[i] * c_dn[i])) for i in 1:2])
        label = NCTSSoS._fermionic_sector_label(first(density_basis), sector)

        blocks = NCTSSoS._fermionic_sector_transform_blocks(
            density_basis,
            label,
            nothing,
            spin,
            active_modes,
        )

        @test [size(block.row_basis, 1) for block in blocks] == oracle.sides
        @test [block.label.total_spin2 for block in blocks] == [0, 2]
        @test all(block.label.sector.irrep == :Ag for block in blocks)
    end

    @testset "Be-derived 1s/2s one-body slice keeps the objective and shrinks PSD blocks" begin
        oracle = expectations_oracle(FERMIONIC_SYMMETRY_EXPECTATIONS_PATH, "be_1s_2s_one_body_n2_order1")
        pop, basis, sector, spin, exact = _be_1s_2s_one_body_problem()

        plain = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                moment_basis=basis,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
            ),
        )
        symmetric = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                moment_basis=basis,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
                symmetry=SymmetrySpec(sector=sector, spin_adaptation=spin),
            ),
        )

        @test plain.objective ≈ exact atol = 5e-6
        @test plain.objective ≈ oracle.opt atol = 5e-6
        @test symmetric.objective ≈ oracle.opt atol = 5e-6
        @test symmetric.objective ≈ plain.objective atol = 5e-6
        @test !isnothing(symmetric.symmetry)
        @test symmetric.symmetry.psd_block_sizes == oracle.sides
        @test maximum(symmetric.symmetry.psd_block_sizes) < only(flatten_sizes(plain.moment_matrix_sizes))
        @test [label.total_spin2 for label in symmetric.symmetry.block_labels if label isa FermionicSpinBlockLabel] == [0, 2]
        @test all(label.sector.irrep == :Ag for label in symmetric.symmetry.block_labels if label isa FermionicSpinBlockLabel)
    end
end

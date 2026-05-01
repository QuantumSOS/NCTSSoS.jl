#!/usr/bin/env julia
# H2/Nk=1 periodic PQG V2RDM demo in NCTSSoS.jl language.
#
# This is the small sibling of h4_periodic_moment_sos.jl:
#   * 2 H atoms/cell, 1 Å spacing, lattice a = 2 Å
#   * Nk = 1
#   * active space [2 electrons, 4 spatial orbitals]
#
# No SDP solver is called here. The useful artifact is the symbolic
# NCTSSoS MomentProblem with explicit D, Q, and G Hermitian PSD blocks.
#
# Usage from the repository root:
#   julia --project=. demos/h2_periodic_moment_sos.jl
#   julia --project=. demos/h2_periodic_moment_sos.jl --include-1d
#
# The default integral dump was generated on HAI with PySCF 2.13.0. See
# test/data/assets/h2_chain_nk1_reference.toml for the reference HF/MP2/CCSD
# values and the constant shift needed to compare the active Hamiltonian with
# total periodic PySCF energies.

using NCTSSoS
using Printf

const SPINS = (:up, :dn)
const DEFAULT_INTEGRALS = normpath(joinpath(
    @__DIR__, "..", "test", "data", "assets", "h2_chain_nk1_active_2e4o_integrals.txt"))

const REFERENCE = (
    hf_total = -1.3447123650300021,
    mp2_total = -1.3692501780283810,
    ccsd_total = -1.3749762789129079,
    ccsd_t_total = -1.3749762789129079,
    active_hf = -25.689602003407181,
    constant_shift = 24.344889638377179,
)

struct Options
    integrals_path::String
    include_one_d::Bool
end

function parse_options(argv)
    integrals_path = DEFAULT_INTEGRALS
    include_one_d = false

    for arg in argv
        if startswith(arg, "--integrals=")
            integrals_path = split(arg, "=", limit = 2)[2]
        elseif arg == "--include-1d"
            include_one_d = true
        elseif arg == "--no-1d"
            include_one_d = false
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. demos/h2_periodic_moment_sos.jl [options]\n")
            println("Options:")
            println("  --integrals=PATH  text integral dump; default test/data/assets/h2_chain_nk1_active_2e4o_integrals.txt")
            println("  --include-1d      add an extra ¹D PSD block")
            println("  --no-1d           explicitly drop the ¹D PSD block (default)")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    return Options(integrals_path, include_one_d)
end

bytes_to_gib(bytes::Integer) = bytes / 1024.0^3
rss_string() = @sprintf("%.3f GiB", bytes_to_gib(Sys.maxrss()))

@inline complex_entry(fields, re_idx::Int, im_idx::Int) =
    complex(parse(Float64, fields[re_idx]), parse(Float64, fields[im_idx]))

function load_integrals_txt(path::AbstractString; norb::Int)
    h1e = Dict{Int,Matrix{ComplexF64}}()
    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()

    for (line_no, raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        (isempty(line) || startswith(line, '#')) && continue

        fields = split(line)
        tag = first(fields)
        if tag == "h1e"
            length(fields) == 6 || error("malformed h1e line $line_no in $(repr(path))")
            k = parse(Int, fields[2])
            p = parse(Int, fields[3]) + 1
            q = parse(Int, fields[4]) + 1
            block = get!(h1e, k) do
                zeros(ComplexF64, norb, norb)
            end
            block[p, q] = complex_entry(fields, 5, 6)
        elseif tag == "eri"
            length(fields) == 11 || error("malformed eri line $line_no in $(repr(path))")
            k1 = parse(Int, fields[2])
            k2 = parse(Int, fields[3])
            k3 = parse(Int, fields[4])
            k4 = parse(Int, fields[5])
            p = parse(Int, fields[6]) + 1
            r = parse(Int, fields[7]) + 1
            q = parse(Int, fields[8]) + 1
            s = parse(Int, fields[9]) + 1
            block = get!(eri, (k1, k2, k3, k4)) do
                zeros(ComplexF64, norb, norb, norb, norb)
            end
            block[p, r, q, s] = complex_entry(fields, 10, 11)
        else
            error("unexpected integral tag $(repr(tag)) at line $line_no in $(repr(path))")
        end
    end

    return h1e, eri
end

function active_hf_energy(h1e, eri; nk::Int, nelec_per_cell::Int)
    nocc = nelec_per_cell ÷ 2

    e1 = 0.0 + 0.0im
    for k in 0:(nk - 1), i in 1:nocc
        e1 += 2.0 * h1e[k][i, i] / nk
    end

    e2 = 0.0 + 0.0im
    for ((k1, k2, k3, k4), Vraw) in eri
        V = Vraw / nk^2
        for i in 1:nocc, j in 1:nocc
            if k1 == k3 && k2 == k4
                e2 += 2.0 * V[i, i, j, j]
            end
            if k1 == k4 && k2 == k3
                e2 -= V[i, j, j, i]
            end
        end
    end

    return real(e1 + e2)
end

function only_monomial(poly)
    mono = monomials(poly)
    length(mono) == 1 || error("expected a single monomial, got $(length(mono)) terms")
    return only(mono)
end

function build_h2_nk1_hamiltonian(h1e, eri; norb::Int)
    registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
        ("c_up", 1:norb),
        ("c_dn", 1:norb),
    ])

    ann = Dict(:up => c_up, :dn => c_dn)
    dag = Dict(:up => c_up_dag, :dn => c_dn_dag)
    ham = (0.0 + 0.0im) * (c_up_dag[1] * c_up[1])

    hk = h1e[0]
    for p in 1:norb, q in 1:norb, spin in SPINS
        coeff = hk[p, q]
        iszero(coeff) && continue
        ham += coeff * dag[spin][p] * ann[spin][q]
    end

    V = eri[(0, 0, 0, 0)]
    spin_channels = ((:up, :up), (:up, :dn), (:dn, :up), (:dn, :dn))
    for p in 1:norb, r in 1:norb, q in 1:norb, s in 1:norb
        coeff = V[p, r, q, s]
        iszero(coeff) && continue
        for (σ, τ) in spin_channels
            ham += 0.5 * coeff *
                dag[σ][p] * dag[τ][q] * ann[τ][s] * ann[σ][r]
        end
    end

    ham = 0.5 * (ham + adjoint(ham))
    vars = (; ann, dag, c_up, c_up_dag, c_dn, c_dn_dag)
    return registry, vars, ham
end

function build_spin_orbitals(vars; norb::Int)
    spin_orbitals = NamedTuple[]
    for spin in SPINS, orb in 1:norb
        push!(spin_orbitals, (; a = vars.ann[spin][orb],
                              adag = vars.dag[spin][orb],
                              spin,
                              orb))
    end
    return spin_orbitals
end

function build_rdm_bases(spin_orbitals)
    NM = typeof(only_monomial(spin_orbitals[1].adag * spin_orbitals[1].a))
    M = length(spin_orbitals)

    oneD_basis = NM[]
    twoD_basis = NM[]
    twoQ_basis = NM[]
    twoG_basis = NM[]

    for p in 1:M
        push!(oneD_basis, only_monomial(spin_orbitals[p].a))
    end
    for p in 1:(M - 1), q in (p + 1):M
        push!(twoD_basis, only_monomial(spin_orbitals[p].a * spin_orbitals[q].a))
        push!(twoQ_basis, only_monomial(spin_orbitals[p].adag * spin_orbitals[q].adag))
    end
    for p in 1:M, q in 1:M
        push!(twoG_basis, only_monomial(spin_orbitals[p].adag * spin_orbitals[q].a))
    end

    return (; oneD_basis, twoD_basis, twoQ_basis, twoG_basis)
end

function total_electron_constraint(vars, ham; norb::Int, total_electrons::Int)
    n = 1.0 * sum(
        vars.c_up_dag[i] * vars.c_up[i] + vars.c_dn_dag[i] * vars.c_dn[i]
        for i in 1:norb
    )
    return n - Float64(total_electrons) * one(ham)
end

function trace_constraints(spin_orbitals, ham; total_electrons::Int)
    n_modes = length(spin_orbitals)
    n_holes = n_modes - total_electrons

    trD = zero(ham)
    trQ = zero(ham)
    trG = zero(ham)

    for p in 1:(n_modes - 1), q in (p + 1):n_modes
        drow = spin_orbitals[p].a * spin_orbitals[q].a
        qrow = spin_orbitals[p].adag * spin_orbitals[q].adag
        trD += adjoint(drow) * drow
        trQ += adjoint(qrow) * qrow
    end
    for p in 1:n_modes, q in 1:n_modes
        grow = spin_orbitals[p].adag * spin_orbitals[q].a
        trG += adjoint(grow) * grow
    end

    return [
        trD - Float64(total_electrons * (total_electrons - 1) ÷ 2) * one(ham),
        trQ - Float64(n_holes * (n_holes - 1) ÷ 2) * one(ham),
        trG - Float64(total_electrons * (n_holes + 1)) * one(ham),
    ]
end

function single_full_system_clique(
    registry,
    moment_support_basis::Vector{NM},
    objective::P;
    n_modes::Int,
) where {NM<:NormalMonomial,P}
    T = eltype(one(NM).word)
    clique_modes = T.(1:n_modes)

    return NCTSSoS.CorrelativeSparsity{FermionicAlgebra,T,P,NM,Nothing}(
        [clique_modes],
        registry,
        P[],
        [Int[]],
        Int[],
        [moment_support_basis],
        [Vector{Vector{NM}}()],
    )
end

function build_h2_pqg_moment_problem(options::Options)
    nk = 1
    norb = 4
    nelec_per_cell = 2
    total_electrons = nk * nelec_per_cell

    h1e, eri = load_integrals_txt(options.integrals_path; norb)
    registry, vars, h2_ham = build_h2_nk1_hamiltonian(h1e, eri; norb)
    objective = h2_ham
    spin_orbitals = build_spin_orbitals(vars; norb)

    meq_constraints = typeof(h2_ham)[
        total_electron_constraint(vars, h2_ham; norb, total_electrons),
    ]
    append!(meq_constraints, trace_constraints(spin_orbitals, h2_ham; total_electrons))

    pop = polyopt(objective, registry; moment_eq_constraints = meq_constraints)
    row_bases = build_rdm_bases(spin_orbitals)
    NM = eltype(row_bases.oneD_basis)
    identity = one(NM)

    oneD_blocks = options.include_one_d ? Vector{NM}[row_bases.oneD_basis] : Vector{NM}[]
    twoD_blocks = Vector{NM}[row_bases.twoD_basis]
    twoQ_blocks = Vector{NM}[row_bases.twoQ_basis]
    twoG_blocks = Vector{NM}[NM[identity; row_bases.twoG_basis]]

    moment_support_basis = options.include_one_d ? NM[
        identity;
        row_bases.oneD_basis;
        row_bases.twoD_basis;
        row_bases.twoQ_basis;
        row_bases.twoG_basis;
    ] : NM[
        identity;
        row_bases.twoD_basis;
        row_bases.twoQ_basis;
        row_bases.twoG_basis;
    ]

    schouten_blocks = Vector{NM}[]
    append!(schouten_blocks, oneD_blocks)
    append!(schouten_blocks, twoD_blocks)
    append!(schouten_blocks, twoQ_blocks)
    append!(schouten_blocks, twoG_blocks)

    full_clique = single_full_system_clique(registry, moment_support_basis, objective;
        n_modes = 2 * nk * norb)
    term_sparsity = NCTSSoS.TermSparsity{NM}(copy(moment_support_basis), schouten_blocks)
    moment_problem = NCTSSoS.moment_relax(pop, full_clique, [[term_sparsity]])
    hf_active = active_hf_energy(h1e, eri; nk, nelec_per_cell)

    return (; nk, norb, nelec_per_cell, total_electrons, h1e, eri,
             registry, vars, h2_ham, objective, pop, meq_constraints,
             spin_orbitals, row_bases, oneD_blocks, twoD_blocks, twoQ_blocks,
             twoG_blocks, moment_support_basis, schouten_blocks, term_sparsity,
             moment_problem, hf_active)
end

real_lift_triangle_rows(n::Integer) = (2n) * (2n + 1) ÷ 2

function constraint_stats(mp)
    hpsd_sizes = [size(mat, 1) for (cone, mat) in mp.constraints if cone == :HPSD]
    zero_sizes = [size(mat) for (cone, mat) in mp.constraints if cone == :Zero]
    total_canonical_moments = length(NCTSSoS._sorted_symmetric_basis(mp.total_basis))
    return (; hpsd_sizes,
             zero_sizes,
             total_canonical_moments,
             direct_real_moment_variables = 2 * total_canonical_moments,
             real_lift_psd_rows = sum(real_lift_triangle_rows(n) for n in hpsd_sizes),
             complex_zero_entries = sum((prod(size) for size in zero_sizes); init = 0),
             direct_real_zero_rows = 2 * sum((prod(size) for size in zero_sizes); init = 0))
end

shift_to_active(total_energy::Real) = total_energy - REFERENCE.constant_shift

function print_summary(data, options::Options, build_seconds::Real)
    stats = constraint_stats(data.moment_problem)
    n_modes = length(data.spin_orbitals)
    n_holes = n_modes - data.total_electrons

    formulation_label = options.include_one_d ?
        "H2/Nk=1 periodic 1D+PQG V2RDM as an NCTSSoS MomentProblem" :
        "H2/Nk=1 periodic PQG V2RDM (no ¹D) as an NCTSSoS MomentProblem"
    println("== ", formulation_label, " ==")
    @printf("%-48s %s\n", "integrals", options.integrals_path)
    @printf("%-48s %s\n", "extra ¹D PSD block", options.include_one_d ? "included" : "disabled")
    @printf("%-48s %s\n", "objective", "H2 active-space Hamiltonian")
    @printf("%-48s N=%d via moment_eq_constraints\n", "particle constraint", data.total_electrons)
    @printf("%-48s TrD=%d, TrQ=%d, TrG=%d\n", "trace constraints",
            data.total_electrons * (data.total_electrons - 1) ÷ 2,
            n_holes * (n_holes - 1) ÷ 2,
            data.total_electrons * (n_holes + 1))
    @printf("%-48s %s\n", "correlative sparsity", "none; single full-system clique")
    println()

    @printf("%-48s %d\n", "k-points", data.nk)
    @printf("%-48s %d\n", "active spatial orbitals per k", data.norb)
    @printf("%-48s %d\n", "spin-orbital modes", n_modes)
    @printf("%-48s %d\n", "active electrons total", data.total_electrons)
    @printf("%-48s %.12f Ha\n", "HF active energy reconstructed", data.hf_active)
    @printf("%-48s %.12f Ha\n", "PySCF HF total per cell", REFERENCE.hf_total)
    @printf("%-48s %.12f Ha\n", "PySCF MP2 total per cell", REFERENCE.mp2_total)
    @printf("%-48s %.12f Ha\n", "PySCF CCSD total per cell", REFERENCE.ccsd_total)
    @printf("%-48s %.12f Ha\n", "PySCF CCSD(T) total per cell", REFERENCE.ccsd_t_total)
    @printf("%-48s %.12f Ha\n", "active-to-total constant shift", REFERENCE.constant_shift)
    @printf("%-48s %.12f Ha\n", "CCSD shifted to active Hamiltonian", shift_to_active(REFERENCE.ccsd_total))
    @printf("%-48s %d\n", "Hamiltonian monomials", length(monomials(data.h2_ham)))
    println()

    println("-- Moment support basis --")
    @printf("%-48s %d\n", "identity", 1)
    @printf("%-48s %s\n", "¹D rows a_p",
            options.include_one_d ? string(length(data.row_bases.oneD_basis)) : "not included")
    @printf("%-48s %d\n", "²D rows a_p a_q, p<q", length(data.row_bases.twoD_basis))
    @printf("%-48s %d\n", "²Q rows a†_p a†_q, p<q", length(data.row_bases.twoQ_basis))
    @printf("%-48s %d\n", "²G rows a†_p a_q", length(data.row_bases.twoG_basis))
    @printf("%-48s %d\n", "total support rows", length(data.moment_support_basis))
    println()

    println("-- Schouten HPSD blocks --")
    @printf("%-48s %s\n", "¹D", options.include_one_d ? string(length.(data.oneD_blocks)) : "not included")
    @printf("%-48s %s\n", "²D", string(length.(data.twoD_blocks)))
    @printf("%-48s %s\n", "²Q", string(length.(data.twoQ_blocks)))
    @printf("%-48s %s\n", "²G", string(length.(data.twoG_blocks)))
    @printf("%-48s %d\n", "total HPSD blocks", length(stats.hpsd_sizes))
    @printf("%-48s %d\n", "largest HPSD block", maximum(stats.hpsd_sizes))
    println()

    println("-- Symbolic moment problem --")
    @printf("%-48s %d\n", "unique moment-matrix monomials", data.moment_problem.n_unique_moment_matrix_elements)
    @printf("%-48s %d\n", "total canonical monomials incl. constraints", stats.total_canonical_moments)
    @printf("%-48s %d\n", "direct real moment variables", stats.direct_real_moment_variables)
    @printf("%-48s %d\n", "real-lift PSD scalar rows", stats.real_lift_psd_rows)
    @printf("%-48s %d\n", "moment_eq polynomials", length(data.meq_constraints))
    @printf("%-48s %d\n", "Zero matrix constraints", length(stats.zero_sizes))
    @printf("%-48s %d\n", "Zero complex scalar entries", stats.complex_zero_entries)
    @printf("%-48s %d\n", "direct real zero/equality rows", stats.direct_real_zero_rows)
    println()

    @printf("%-48s %.3f s\n", "build walltime", build_seconds)
    @printf("%-48s %s\n", "max RSS after build", rss_string())
    println()
    println("No optimizer was attached and optimize! was not called.")
    println("Artifact: symbolic PQG MomentProblem over one shared fermionic CAR moment map.")
    return nothing
end

function main(argv = ARGS)
    options = parse_options(argv)
    @printf("%-48s %s\n", "max RSS at start", rss_string())
    build_seconds = @elapsed data = build_h2_pqg_moment_problem(options)
    print_summary(data, options, build_seconds)
    return data
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

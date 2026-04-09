module H4PeriodicAssets

using TOML

export H4_EXPECTATIONS_PATH,
       H4_INTEGRALS_PATH,
       H4_REFERENCE_PATH,
       build_k2,
       complete_graph_edge_count,
       compound_spatial_edge_count,
       fermionic_degree_leq_count,
       fermionic_order2_basis_size,
       fermionic_order2_nuniq,
       hf_energy,
       load_expectation,
       load_integrals_txt,
       load_nk2_asset,
       load_reference

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const H4_EXPECTATIONS_PATH = joinpath(REPO_ROOT, "test", "data", "expectations", "h4_periodic_v2rdm.toml")
const H4_REFERENCE_PATH = joinpath(REPO_ROOT, "test", "data", "assets", "h4_chain_reference.toml")
const H4_INTEGRALS_PATH = joinpath(REPO_ROOT, "test", "data", "assets", "h4_chain_nk2_integrals.txt")

const NK2_PREFLIGHT_EXPECTATION_ID = "nk2_asset_preflight"
const NK2_ORDER2_BLOCKER_EXPECTATION_ID = "spin_orbital_order2_blocker"

load_reference() = TOML.parsefile(H4_REFERENCE_PATH)
load_expectations() = TOML.parsefile(H4_EXPECTATIONS_PATH)

function expectation_by_id(expectations, id::AbstractString)
    for case in expectations["cases"]
        case["id"] == id && return case["expected"]
    end
    error("H4 expectation case not found: $(repr(id))")
end

load_expectation(id::AbstractString) = expectation_by_id(load_expectations(), id)

function reference_energy_entry(reference, nk::Int)
    for entry in reference["energies"]
        entry["nk"] == nk && return entry
    end
    error("H4 reference energy entry not found for nk = $nk")
end

function load_nk2_asset()
    expectations = load_expectations()
    preflight = expectation_by_id(expectations, NK2_PREFLIGHT_EXPECTATION_ID)
    blocker = expectation_by_id(expectations, NK2_ORDER2_BLOCKER_EXPECTATION_ID)
    reference = load_reference()

    nk = preflight["nk"]
    n_active_orb = preflight["active_space_orbitals_per_k"]
    n_active_elec = preflight["active_space_electrons_per_cell"]
    total_active_electrons = preflight["total_active_electrons"]
    total_spatial_orbitals = preflight["total_spatial_orbitals"]
    total_spin_orbital_modes = preflight["total_spin_orbital_modes"]
    figure = reference_energy_entry(reference, nk)
    h1e, eri = load_integrals_txt(H4_INTEGRALS_PATH; n_active_orb)

    return (; reference,
             preflight,
             blocker,
             figure,
             nk,
             n_active_orb,
             n_active_elec,
             total_active_electrons,
             total_spatial_orbitals,
             total_spin_orbital_modes,
             h1e,
             eri)
end

@inline compound_orbital_index(k::Int, orbital::Int, n_active_orb::Int) = k * n_active_orb + orbital
@inline complex_entry(fields, re_idx::Int, im_idx::Int) =
    complex(parse(Float64, fields[re_idx]), parse(Float64, fields[im_idx]))

function load_integrals_txt(path::AbstractString = H4_INTEGRALS_PATH; n_active_orb::Int)
    h1e = Dict{Int, Matrix{ComplexF64}}()
    eri = Dict{NTuple{4, Int}, Array{ComplexF64, 4}}()

    for (line_no, raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, '#') && continue

        fields = split(line)
        tag = first(fields)

        if tag == "h1e"
            length(fields) == 6 || error("Malformed h1e line $line_no in $(repr(path))")

            ik = parse(Int, fields[2])
            p = parse(Int, fields[3]) + 1
            q = parse(Int, fields[4]) + 1
            value = complex_entry(fields, 5, 6)

            block = get!(h1e, ik) do
                zeros(ComplexF64, n_active_orb, n_active_orb)
            end
            block[p, q] = value
        elseif tag == "eri"
            length(fields) == 11 || error("Malformed eri line $line_no in $(repr(path))")

            ik1 = parse(Int, fields[2])
            ik2 = parse(Int, fields[3])
            ik3 = parse(Int, fields[4])
            ik4 = parse(Int, fields[5])

            # The text dump uses chemist-style orbital order
            # (pₖ₁ rₖ₃ | qₖ₂ sₖ₄), so after the k-block tuple we parse the
            # orbital indices as (p, r, q, s).
            p = parse(Int, fields[6]) + 1
            r = parse(Int, fields[7]) + 1
            q = parse(Int, fields[8]) + 1
            s = parse(Int, fields[9]) + 1
            value = complex_entry(fields, 10, 11)

            block = get!(eri, (ik1, ik2, ik3, ik4)) do
                zeros(ComplexF64, n_active_orb, n_active_orb, n_active_orb, n_active_orb)
            end
            block[p, r, q, s] = value
        else
            error("Unexpected H4 integral tag $(repr(tag)) at line $line_no")
        end
    end

    return h1e, eri
end

function build_k2(h1e, eri; nk::Int, n_active_orb::Int, n_total_electrons::Int)
    n_total_electrons > 1 || error("K₂ construction requires at least two electrons")

    dim = nk * n_active_orb
    k2 = zeros(ComplexF64, dim, dim, dim, dim)
    one_body_prefactor = 2.0 / (n_total_electrons - 1)

    h1e_norm = Dict(ik => block / nk for (ik, block) in h1e)
    eri_norm = Dict(key => block / nk^2 for (key, block) in eri)

    for ik1 in 0:nk-1
        h_block = h1e_norm[ik1]
        for ik2 in 0:nk-1
            for p in 1:n_active_orb, q in 1:n_active_orb, r in 1:n_active_orb, s in 1:n_active_orb
                i1 = compound_orbital_index(ik1, p, n_active_orb)
                i2 = compound_orbital_index(ik2, q, n_active_orb)

                if q == s
                    i3 = compound_orbital_index(ik1, r, n_active_orb)
                    i4 = compound_orbital_index(ik2, s, n_active_orb)
                    k2[i1, i2, i3, i4] += one_body_prefactor * h_block[p, r]
                end

                if q == r
                    i3_exchange = compound_orbital_index(ik2, r, n_active_orb)
                    i4_exchange = compound_orbital_index(ik1, s, n_active_orb)
                    k2[i1, i2, i3_exchange, i4_exchange] -= one_body_prefactor * h_block[p, s]
                end
            end
        end
    end

    for ((ik1, ik2, ik3, ik4), eri_block) in eri_norm
        for p in 1:n_active_orb, q in 1:n_active_orb, r in 1:n_active_orb, s in 1:n_active_orb
            i1 = compound_orbital_index(ik1, p, n_active_orb)
            i2 = compound_orbital_index(ik2, q, n_active_orb)
            i3 = compound_orbital_index(ik3, r, n_active_orb)
            i4 = compound_orbital_index(ik4, s, n_active_orb)
            k2[i1, i2, i3, i4] += eri_block[p, r, q, s]
        end
    end

    return k2
end

function hf_energy(h1e, eri; nk::Int, n_active_elec::Int)
    nocc = n_active_elec ÷ 2

    e1 = 0.0 + 0.0im
    for ik in 0:nk-1
        h_block = h1e[ik]
        for i in 1:nocc
            e1 += 2.0 * h_block[i, i]
        end
    end
    e1 /= nk

    e2 = 0.0 + 0.0im
    for eri_block in values(eri)
        for i in 1:nocc, j in 1:nocc
            e2 += 2.0 * eri_block[i, i, j, j]
            e2 -= eri_block[i, j, j, i]
        end
    end
    e2 /= nk^2

    return real(e1 + e2)
end

function compound_spatial_edge_count(h1e, eri; n_active_orb::Int, atol::Float64 = 1e-10)
    edges = Set{Tuple{Int, Int}}()

    for (ik, block) in h1e
        for p in 1:n_active_orb, q in 1:n_active_orb
            if p != q && abs(block[p, q]) > atol
                a = compound_orbital_index(ik, p, n_active_orb)
                b = compound_orbital_index(ik, q, n_active_orb)
                push!(edges, a < b ? (a, b) : (b, a))
            end
        end
    end

    for ((ik1, ik2, ik3, ik4), block) in eri
        for idx in CartesianIndices(block)
            abs(block[idx]) <= atol && continue

            p, r, q, s = Tuple(idx)
            nodes = unique(sort([
                compound_orbital_index(ik1, p, n_active_orb),
                compound_orbital_index(ik2, q, n_active_orb),
                compound_orbital_index(ik3, r, n_active_orb),
                compound_orbital_index(ik4, s, n_active_orb),
            ]))
            for i in 1:length(nodes)-1, j in i+1:length(nodes)
                push!(edges, (nodes[i], nodes[j]))
            end
        end
    end

    return length(edges)
end

fermionic_degree_leq_count(n_modes::Int, max_degree::Int) =
    sum(binomial(big(2 * n_modes), degree) for degree in 0:max_degree)

fermionic_order2_basis_size(n_modes::Int) = fermionic_degree_leq_count(n_modes, 2)
fermionic_order2_nuniq(n_modes::Int) = fermionic_degree_leq_count(n_modes, 4)
complete_graph_edge_count(n_vertices::Int) = n_vertices * (n_vertices - 1) ÷ 2

end # module H4PeriodicAssets

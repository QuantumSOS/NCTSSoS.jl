# Included after strengthening.jl by the prototype module.

"""
    workflow_local_density_candidates(H, family)

Audit ordinary rows before symmetry reduction. Require distinct, coefficient-one
Pauli words on the original registry, then use the strengthening module's
checked saturated rewrite to verify the complete local Pauli algebra and build its RDM.
"""
function workflow_local_density_candidates(H::NC.Polynomial{W,C,R}, family) where {W<:NC.PauliWord,C,R}
    candidates = NamedTuple{(:family_index, :sites, :words, :matrix),
        Tuple{Int, Vector{Int}, Vector{W}, Matrix{NC.Polynomial{W,_PauliCQ,R}}}}[]
    for (i, rows) in enumerate(family)
        len = length(rows)
        len >= 4 && ispow2(len) && iseven(trailing_zeros(len)) || continue
        all(row -> row isa NC.Polynomial && row.registry === H.registry &&
            length(row.data) == 1 && row.data.coefficients[1] == 1, rows) || continue
        words = sort!(W[only(row.data.words) for row in rows])
        allunique(words) || continue
        sites = sort!(unique([s for w in words for s in NC.sites(w)]))
        length(sites) == trailing_zeros(len) ÷ 2 && length(sites) <= 6 &&
            all(s -> 1 <= s <= H.registry.nsites, sites) || continue
        rewrite = workflow_density_rewrite(H, words, sites)
        push!(candidates, (family_index=i, sites, words, matrix=rewrite.matrix))
    end
    return candidates
end

function _density_moment_map(table, words, nvars)
    coordinates = Dict{eltype(words), Dict{Int32,_PauliCQ}}()
    for w in words
        terms = NC.moment_terms(table, w)
        isnothing(terms) && return nothing
        row = Dict{Int32,_PauliCQ}()
        for (id, c) in terms
            1 <= id <= nvars || return nothing
            NC._sector_add!(row, id, _PauliCQ(c))
        end
        coordinates[w] = row
    end
    return coordinates
end

function _density_matches(program, block, words, coordinates)
    block.label == 0 && block.kind in (NC.PSD, NC.HPSD) && block.side == length(words) || return false
    actual = Dict{Int32,_PauliCQ}()
    expected = empty(actual)
    for j in eachindex(words), i in j:length(words)
        phase_star, lhs = NC.star(words[i])
        phase_prod, word = NC.mul(lhs, words[j])
        phase = phase_star * phase_prod
        empty!(expected)
        for (id, c) in coordinates[word]
            expected[id] = NC.scale(phase, c)
        end
        empty!(actual)
        f = block.first_form + NC.tri_offset(block.side, i, j)
        for t in NC.term_range(program.forms, f)
            NC._sector_add!(actual, program.forms.ids[t], _PauliCQ(program.forms.coefs[t]))
        end
        expected == actual || return false
    end
    return true
end

function _density_replacement(::Type{C}, matrix, coordinates, nvars, label) where {C}
    density = NC.MomentProgram{C}(nvars)
    side = size(matrix, 1)
    row = Dict{Int32,_PauliCQ}()
    for j in 1:side, i in j:side
        empty!(row)
        for (factor, w) in NC.terms(matrix[i,j]), (id, c) in coordinates[w]
            NC._sector_add!(row, id, factor*c)
        end
        for (id, c) in sort!(collect(row); by=first)
            # Do not silently round a reparameterization coefficient or
            # coerce a genuinely complex affine form into a real program.
            C <: Real && !isreal(c) && return nothing
            stored = C(C <: Real ? real(c) : c)
            _PauliCQ(stored) == c || return nothing
            NC.push_term!(density.forms, stored, id)
        end
        NC.close_form!(density.forms)
    end
    NC.push_block!(density, C <: Complex ? NC.HPSD : NC.PSD, side, label)
    # Preserve even coordinates introduced solely by the removed cone.
    old_ids = Set(id for row in values(coordinates) for id in keys(row))
    Set(density.forms.ids) == old_ids || return nothing
    return density
end

"""
    workflow_local_density(built, candidates)

Replace an eligible cone only if its complete affine Gram matches the audited
raw basis, including product phases and original moment IDs. For the full
local matrix algebra, this Gram is unitarily similar to `2^r` copies of the
density matrix scaled by `2^r` (up to transposition). A rotated or reduced cone
without a verified ambient Gram is left alone. The table is read-only.

Return `(program, report)` with `report.eligible`, `report.applied`, and
`report.skipped`. Negative labels identify replacements to downstream stages.
"""
function workflow_local_density(built, candidates)
    program, table = built.program, built.table
    replacements = Dict{Int, typeof(program)}()
    applied = NamedTuple{(:family_index, :sites, :block_index, :before, :after), Tuple{Int, Vector{Int}, Int, Int, Int}}[]
    skipped = Pair{Int, Symbol}[]
    for candidate in candidates
        coordinates = _density_moment_map(table, candidate.words, program.nvars)
        if isnothing(coordinates)
            push!(skipped, candidate.family_index => :missing_coordinates)
            continue
        end
        index = findfirst(eachindex(program.blocks)) do i
            !haskey(replacements, i) && _density_matches(program, program.blocks[i], candidate.words, coordinates)
        end
        if isnothing(index)
            push!(skipped, candidate.family_index => :no_matching_ambient_cone)
            continue
        end
        replacement = _density_replacement(eltype(program.forms.coefs), candidate.matrix, coordinates, program.nvars, -index)
        if isnothing(replacement)
            push!(skipped, candidate.family_index => :coordinate_audit_failed)
            continue
        end
        replacements[index] = replacement
        push!(applied, (family_index=candidate.family_index, sites=candidate.sites, block_index=index,
            before=program.blocks[index].side, after=only(replacement.blocks).side))
    end
    if !isempty(replacements)
        out = typeof(program)(program.nvars)
        append!(out.objective_ids, program.objective_ids)
        append!(out.objective_coefs, program.objective_coefs)
        for (i, block) in enumerate(program.blocks)
            if haskey(replacements, i)
                replacement = replacements[i]
                NC._copy_cone!(out, replacement, only(replacement.blocks))
            else
                NC._copy_cone!(out, program, block)
            end
        end
        NC.validate(out)
    else
        out = program
    end
    return (program=out, report=(eligible=length(candidates), applied, skipped))
end

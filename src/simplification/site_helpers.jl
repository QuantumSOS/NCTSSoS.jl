# =============================================================================
# Shared Helpers for Site-Encoded Algebras (Unsigned indices)
# =============================================================================

@inline function _stable_sort_by_site!(word::Vector{T}) where {T<:Unsigned}
    sort!(word; alg=Base.Sort.InsertionSort, by=decode_site)
    return word
end

function _validate_site_sorted_word!(
    word::Vector{T};
    algebra_name::AbstractString,
    sorted_hint::AbstractString="",
    forbid_adjacent_duplicates::Bool=false,
    duplicate_rule::AbstractString="",
    duplicate_hint::AbstractString=""
) where {T<:Unsigned}
    length(word) <= 1 && return nothing

    prev_site = decode_site(word[1])
    prev_idx = word[1]
    for i in 2:length(word)
        curr_site = decode_site(word[i])
        curr_idx = word[i]

        if curr_site < prev_site
            msg = "$algebra_name word not sorted by site: site $curr_site at index $i comes after site $prev_site"
            isempty(sorted_hint) || (msg *= ". " * sorted_hint)
            throw(ArgumentError(msg))
        elseif forbid_adjacent_duplicates && curr_idx == prev_idx
            msg = "$algebra_name word has consecutive identical operators"
            isempty(duplicate_rule) || (msg *= " ($duplicate_rule)")
            msg *= " at indices $(i - 1) and $i"
            isempty(duplicate_hint) || (msg *= ". " * duplicate_hint)
            throw(ArgumentError(msg))
        end

        prev_site = curr_site
        prev_idx = curr_idx
    end

    return nothing
end


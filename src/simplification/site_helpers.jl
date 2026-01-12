# =============================================================================
# Shared Helpers for Site-Encoded Algebras (Unsigned indices)
# =============================================================================

"""
    _stable_sort_by_site!(word::Vector{T}) where {T<:Unsigned} -> Vector{T}

Sort `word` in-place by decoded site using stable insertion sort.
Preserves relative order of operators at the same site.
"""
@inline _stable_sort_by_site!(word::Vector{T}) where {T<:Unsigned} = sort!(word; alg=Base.Sort.InsertionSort, by=decode_site)

"""
    _stable_sort_by_site!(word::Vector{T}) where {T<:Signed} -> Vector{T}

Sort `word` in-place by mode index using stable insertion sort.
Groups operators by their mode (absolute value), preserving relative order within each mode.
"""
@inline _stable_sort_by_site!(word::Vector{T}) where {T<:Signed} = sort!(word; alg=Base.Sort.InsertionSort, by=_operator_mode)

"""
    _is_site_sorted(word::Vector{T}) where {T<:Unsigned} -> Bool

Check if `word` is sorted by decoded site (ascending order).
"""
@inline _is_site_sorted(word::Vector{T}) where {T<:Unsigned} = issorted(word, by=decode_site)

"""
    _is_site_sorted(word::Vector{T}) where {T<:Signed} -> Bool

Check if `word` has sign-based ordering (all negatives before positives).
Used for Fermionic algebra where negative = annihilation, positive = creation.
"""
@inline _is_site_sorted(word::Vector{T}) where {T<:Signed} = issorted(word, by=_is_creation)

"""
    _is_site_sorted_pauli(word::Vector{T}) where {T<:Integer} -> Bool

Check if `word` is sorted by Pauli site (using `_pauli_site` encoding).
"""
@inline _is_site_sorted_pauli(word::Vector{T}) where {T<:Integer} = issorted(word, by=_pauli_site)

struct SimplifyAlgorithm
    comm_gps::Dict{Variable,Int}
    n_gps::Int
    is_unipotent::Bool
    is_projective::Bool
    function SimplifyAlgorithm(;
        comm_gps::Vector{Vector{Variable}},
        is_unipotent::Bool=false,
        is_projective::Bool=false,
    )
        return new(
            Dict(var => i for (i, vars) in enumerate(comm_gps) for var in vars),
            length(comm_gps),
            is_unipotent,
            is_projective,
        )
    end
end

function nosimp(sa::SimplifyAlgorithm)
    return isone(sa.n_gps) && !sa.is_unipotent && !sa.is_projective
end

function simplify!(m::Monomial, sa::SimplifyAlgorithm)
    (isone(m) || nosimp(sa)) && return m
    _comm!(m, sa.comm_gps)

    sa.is_unipotent && (_simplify_unipotent!(m); return m)
    sa.is_projective && (_simplify_projective!(m); return m)
    _simplify_standard!(m)
    return m
end

function _simplify_standard!(m::Monomial)
    head_idx = 1
    @inbounds for i in 2:length(m.vars)
        if m.vars[i] == m.vars[head_idx]
            m.z[head_idx] += m.z[i]
        else
            head_idx += 1
            m.vars[head_idx] = m.vars[i]
            m.z[head_idx] = m.z[i]
        end
    end
    if head_idx != length(m.vars)
        deleteat!(m.vars, (head_idx + 1):length(m.vars))
        deleteat!(m.z, (head_idx + 1):length(m.z))
    end
    return nothing
end

function _simplify_projective!(m::Monomial)
    head_idx = 1
    m.z[head_idx] = 1
    @inbounds for i in 2:length(m.vars)
        m.vars[i] == m.vars[head_idx] && continue
        head_idx += 1
        m.vars[head_idx] = m.vars[i]
        m.z[head_idx] = 1
    end
    if head_idx != length(m.vars)
        deleteat!(m.vars, (head_idx + 1):length(m.vars))
        deleteat!(m.z, (head_idx + 1):length(m.z))
    end
    return nothing
end

function _simplify_unipotent!(m::Monomial)
    head_idx = 1
    tail_idx = findfirst(isodd, m.z)
    isnothing(tail_idx) && (empty!(m.vars); empty!(m.z); return nothing)

    m.vars[head_idx] = m.vars[tail_idx]
    m.z[head_idx] = mod(m.z[tail_idx], 2)

    tail_idx += 1
    @inbounds for i in tail_idx:length(m.vars)
        i < tail_idx && continue
        iseven(m.z[i]) && continue
        if m.vars[i] == m.vars[head_idx]
            head_idx -= 1
            head_idx >= 1 && continue

            head_idx = findfirst(isodd, view(m.z, (i + 1):length(m.z)))
            isnothing(head_idx) && (empty!(m.vars); empty!(m.z); return m)
            m.vars[1] = m.vars[i + head_idx]
            m.z[1] = 1

            tail_idx = i + head_idx + 1
            head_idx = 1
        else
            head_idx += 1
            m.vars[head_idx] = m.vars[i]
            m.z[head_idx] = 1
        end
    end
    if head_idx != length(m.vars)
        deleteat!(m.vars, (head_idx + 1):length(m.vars))
        deleteat!(m.z, (head_idx + 1):length(m.z))
    end
    return nothing
end

function simplify(m::Monomial, sa::SimplifyAlgorithm)
    return simplify!(copy(m), sa)
end

function simplify(sw::StateWord{ST}, sa::SimplifyAlgorithm) where {ST}
    return StateWord{ST}(canonicalize.(sw.state_monos, Ref(sa)))
end

function simplify!(sw::StateWord{ST}, sa::SimplifyAlgorithm) where {ST}
    return simplify(sw, sa)
end

function simplify(ncsw::NCStateWord, sa::SimplifyAlgorithm)
    return NCStateWord(simplify(ncsw.sw, sa), simplify(ncsw.nc_word, sa))
end

function simplify!(ncsw::NCStateWord, sa::SimplifyAlgorithm)
    return simplify(ncsw, sa)
end

function simplified_star_min!(mono::Monomial, sa::SimplifyAlgorithm)
    # return mono = min(simplify!(star(mono), sa), simplify!(mono, sa))
    # return mono = min(simplify!(mono, sa), simplify!(star(mono), sa))
    # the real star of x1x2y1y2 is not y2y1x2x1 but x2x1y2y1, considering comm_gps
    simplify!(mono, sa)
    @inbounds for i in 1:(sa.n_gps)
        a_init = a_idx = findfirst(x -> sa.comm_gps[x] == i, mono.vars)
        b_init = b_idx = findlast(x -> sa.comm_gps[x] == i, mono.vars)
        isnothing(a_init) && continue
        isnothing(b_init) && continue
        while a_idx <= b_init && b_idx >= a_init
            iszero(mono.z[a_idx]) && (a_idx += 1; continue)
            iszero(mono.z[b_idx]) && (b_idx -= 1; continue)

            var_cmp = cmp(mono.vars[a_idx], mono.vars[b_idx])
            if var_cmp > 0
                (star!(mono); _comm!(mono, sa.comm_gps); return mono)
            elseif var_cmp < 0
                return mono
            end
            if mono.z[a_idx] > mono.z[b_idx]
                (star!(mono); _comm!(mono, sa.comm_gps); return mono)
            elseif mono.z[a_idx] < mono.z[b_idx]
                return mono
            end
            a_idx += 1
            b_idx -= 1
        end
    end
    return mono
end

"""
    cyclic_canonicalize(mono::Monomial)

Canonicalizes a monomial using both cyclic and symmetric operations.
Finds the minimum among all cyclic shifts and their adjoints.
Chclic canonical is both cyclic and symmetric

# Arguments
- `mono::Monomial`: The mono to canonicalize

# Returns
- `Monomial`: Canonicalized monomial (minimum across all cyclic shifts and their stars)
"""
function cyclic_canonicalize(mono::Monomial, sa::SimplifyAlgorithm)
    isone(mono) && return mono
    flatten_vars = mapreduce(
        idx -> fill(mono.vars[idx], mono.z[idx]), vcat, eachindex(mono.z)
    )
    flatten_z = ones(Int, sum(mono.z))
    return mapreduce(
        min, 1:(sum(mono.z) - 1); init=min(simplify!(star(mono), sa), simplify(mono, sa))
    ) do _
        shifted_mono = monomial(circshift!(flatten_vars, 1), circshift!(flatten_z, 1))
        simplified_star_min!(shifted_mono, sa)
    end
end

"""
    symmetric_canonicalize(mono::Monomial)

Canonicalizes a mono by taking the minimum between itself and its adjoint.

# Arguments
- `mono::Monomial`: The mono to canonicalize

# Returns
- `Monomial`: Canonicalized monomial (minimum of original and its star)
"""
function symmetric_canonicalize(mono::Monomial, sa::SimplifyAlgorithm)
    isone(mono) && return mono
    return simplified_star_min!(copy(mono), sa)
end

function canonicalize(sw::StateWord{MaxEntangled}, sa::SimplifyAlgorithm)
    return StateWord{MaxEntangled}(cyclic_canonicalize.(sw.state_monos, Ref(sa)))
end

canonicalize(m::Monomial, sa::SimplifyAlgorithm) = symmetric_canonicalize(m, sa)

function canonicalize(sw::StateWord{Arbitrary}, sa::SimplifyAlgorithm)
    return StateWord{Arbitrary}(symmetric_canonicalize.(sw.state_monos, Ref(sa)))
end

function canonicalize(ncsw::NCStateWord, sa::SimplifyAlgorithm)
    return NCStateWord(canonicalize(ncsw.sw, sa), symmetric_canonicalize(ncsw.nc_word, sa))
end

"""
    canonicalize(poly::Polynomial)

Canonicalizes a polynomial by applying symmetric canonicalization to all monomials.

It only makes sense of symmetric canonicalize a polynomial or non-commuting
varaibles because trace polynomial is taken care at StatePoly

# Arguments
- `poly::Polynomial`: The polynomial to canonicalize

# Returns
- `Polynomial`: Canonicalized polynomial with conjugated coefficients and canonicalized monomials
"""
function canonicalize(poly::Polynomial, sa::SimplifyAlgorithm)
    return Polynomial(conj.(poly.coeffs), symmetric_canonicalize.(poly.monos, Ref(sa)))
end

# ς(w) = ς(w') stated in https://arxiv.org/abs/2301.12513, Section 2.1
# You DO CANONICALIZE HERER!!!
function get_state_basis(
    ::Type{ST}, variables::Vector{Variable}, d::Int, sa::SimplifyAlgorithm
) where {ST}
    canon_algo = ST == MaxEntangled ? cyclic_canonicalize : symmetric_canonicalize
    return unique!(
        map(
            a -> NCStateWord(StateWord{ST}(a[1]), a[2]),
            mapreduce(vcat, 0:d) do nc_deg
                nc_basis = simplify.(monomials(variables, Val(nc_deg)), Ref(sa))
                cw_deg = d - nc_deg
                cw_basis = unique!([
                    begin
                        interm = sort(filter(!isone, collect(c_word)))
                        isempty(interm) ? [one(variables[1])] : interm
                    end for c_word in Iterators.product(
                        ntuple(
                            _ ->
                                unique!(canon_algo.(get_basis(variables, cw_deg), Ref(sa))),
                            cw_deg,
                        )...,
                        [one(variables[1])],
                    ) if sum(degree.(c_word)) <= cw_deg
                ])
                reshape(collect(Iterators.product(cw_basis, nc_basis)), :)
            end,
        ),
    )
end

function get_basis(
    ::Type{NCStatePolynomial{T,ST}},
    variables::Vector{Variable},
    d::Int,
    sa::SimplifyAlgorithm,
) where {T,ST}
    return get_state_basis(ST, variables, d, sa)
end

function canonicalize(sp::StatePolynomial, sa::SimplifyAlgorithm)
    return StatePolynomial((sp.coeffs), canonicalize.(sp.state_words, Ref(sa)))
end

function canonicalize(ncsp::NCStatePolynomial, sa::SimplifyAlgorithm)
    return NCStatePolynomial((ncsp.coeffs), canonicalize.(ncsp.nc_state_words, Ref(sa)))
end

function get_basis(
    ::Type{Polynomial{T}}, variables::Vector{Variable}, d::Int, sa::SimplifyAlgorithm
) where {T}
    return sorted_unique(simplify.(get_basis(variables, d), Ref(sa)))
end

function is_symmetric(p::Polynomial, sa::SimplifyAlgorithm)
    return iszero(
        sum(c * simplify(m, sa) for (c, m) in terms(p)) -
        sum(conj(c) * simplify(star(m), sa) for (c, m) in terms(p)),
    )
end

# Exact validation only: this module does not assemble a reduced relaxation.
# A derivation specifies infinitesimal images and extends by the Leibniz rule.
const _ActionCoefficient = Complex{Rational{BigInt}}

function _action_coefficient(c::Number)
    c isa Union{Integer,Rational,AbstractFloat,Complex} || throw(ArgumentError(
        "Supplied generators require rational or finite floating-point coefficients."))
    isfinite(c) || throw(ArgumentError("Supplied generators require finite coefficients."))
    return _ActionCoefficient(c)
end

function _action_polynomial(p::Polynomial{A,T}) where {A,T}
    return Polynomial(Tuple{_ActionCoefficient,NormalMonomial{A,T}}[
        (_action_coefficient(c), m) for (c, m) in terms(p)
    ])
end
_action_polynomial(m::NormalMonomial) = Polynomial((one(_ActionCoefficient), m))

# Ordinary polynomial multiplication promotes through coeff_type(A), which is
# floating point. Use the existing normal-form kernels with an exact output type.
function _action_product(p::Polynomial{A,T}, q::Polynomial{A,T}) where {A,T}
    result = Tuple{_ActionCoefficient,NormalMonomial{A,T}}[]
    for (c, m) in terms(p), (d, n) in terms(q)
        _push_product_terms!(result, c * d, m, n, _ActionCoefficient)
    end
    return Polynomial(result)
end

# PBW normal forms can have integer coefficients larger than Float64 can
# represent. Read those integers directly, before _simplified_to_terms casts.
function _action_normal_form(::Type{A}, word) where {A<:PBWAlgebra}
    # The bosonic normal-ordering kernel uses machine Int coefficients. Bound
    # the number of contractions by n! to exclude overflow from an exact claim.
    if A === BosonicAlgebra && factorial(big(length(word))) > typemax(Int)
        throw(ArgumentError("Exact bosonic action validation is unsupported at this word length: normal-ordering coefficients may overflow Int."))
    end
    return simplify(A, word)
end

function _action_product(p::Polynomial{A,T}, q::Polynomial{A,T}) where {A<:PBWAlgebra,T<:Integer}
    result = Tuple{_ActionCoefficient,NormalMonomial{A,T}}[]
    for (c, m) in terms(p), (d, n) in terms(q)
        for (scale, word) in _action_normal_form(A, vcat(m.word, n.word))
            push!(result, (c * d * scale, NormalMonomial{A,T}(word)))
        end
    end
    return Polynomial(result)
end

_action_adjoint(p::Polynomial) = adjoint(p)
function _action_adjoint(p::Polynomial{A,T}) where {A<:PBWAlgebra,T<:Integer}
    result = Tuple{_ActionCoefficient,NormalMonomial{A,T}}[]
    for (c, m) in terms(p), (scale, word) in _action_normal_form(A, -reverse(m.word))
        push!(result, (conj(c) * scale, NormalMonomial{A,T}(word)))
    end
    return Polynomial(result)
end

# Public construction requires every registered letter, including adjoints, to
# have an explicit image (use a zero polynomial for fixed letters). No omitted
# image is silently interpreted as a zero action.
struct SuppliedDerivation{A<:AlgebraType,T<:Integer}
    images::Dict{T,Polynomial{A,T,_ActionCoefficient}}
end

function SuppliedDerivation(reg::VariableRegistry{A,T}, images::AbstractDict) where {A,T}
    A in (NonCommutativeAlgebra, ProjectorAlgebra, UnipotentAlgebra,
        PauliAlgebra, FermionicAlgebra, BosonicAlgebra) || throw(ArgumentError(
        "Defining-relation validation is unsupported for algebra $A."))
    domain = Set(indices(reg))
    Set(keys(images)) == domain || throw(ArgumentError(
        "Supply exactly one image for every registered index, including adjoints."))
    converted = Dict{T,Polynomial{A,T,_ActionCoefficient}}()
    for idx in indices(reg)
        image = images[idx]
        image isa Union{Polynomial{A,T},NormalMonomial{A,T}} || throw(ArgumentError(
            "Image of index $idx must be a polynomial or monomial in the registry algebra."))
        p = _action_polynomial(image)
        all(i in domain for (_, m) in terms(p) for i in m.word) || throw(ArgumentError(
            "Image of index $idx contains an unregistered letter."))
        converted[idx] = p
    end
    return SuppliedDerivation{A,T}(converted)
end

# Differentiate an unsimplified word first. In particular δ(x*x) is computed
# before x*x is replaced by 1, x, or a Pauli phase.
function _differentiate_word(action::SuppliedDerivation{A,T}, word) where {A,T}
    P = Polynomial{A,T,_ActionCoefficient}
    value, derivative = one(P), zero(P)
    for idx in word
        haskey(action.images, idx) || throw(ArgumentError("Unregistered action index $idx."))
        letter = _action_polynomial(NormalMonomial{A,T}(T[idx]))
        derivative = _action_product(derivative, letter) + _action_product(value, action.images[idx])
        value = _action_product(value, letter)
    end
    return derivative
end

function differentiate(action::SuppliedDerivation{A,T}, p::Polynomial{A,T}) where {A,T}
    result = zero(Polynomial{A,T,_ActionCoefficient})
    for (c, m) in terms(_action_polynomial(p))
        result += c * _differentiate_word(action, m.word)
    end
    return result
end
differentiate(action::SuppliedDerivation, m::NormalMonomial) = differentiate(action, _action_polynomial(m))

function _action_zero(p, context)
    iszero(p) || throw(ArgumentError("$context: nonzero exact derivative/residual $p."))
    return nothing
end

function _validate_derivation(action::SuppliedDerivation{A,T}) where {A,T}
    domain = sort!(collect(keys(action.images)))
    letters = [_action_polynomial(NormalMonomial{A,T}(T[i])) for i in domain]
    for (i, x) in zip(domain, letters)
        _action_zero(differentiate(action, _action_adjoint(x)) - _action_adjoint(action.images[i]),
            "Adjoint compatibility at index $i")
    end
    # These six built-in algebras have quadratic defining rewrite relations.
    # Compare δ(raw pair) with δ(its normal form), not δ(reduced zero).
    for (i, x) in zip(domain, letters), (j, y) in zip(domain, letters)
        _action_zero(_differentiate_word(action, (i, j)) -
            differentiate(action, _action_product(x, y)), "Defining relation at indices ($i, $j)")
    end
    return nothing
end

# Exact elimination for polynomial span membership. A dependent input is allowed
# for constraint spans, but positivity rows and metric rows must be bases.
function _action_coordinates(rows, targets; independent=false, context="row space")
    # Most row and letter bases are distinct scaled monomials. Their coordinates
    # are coefficient lookups, not a dense elimination of an identity matrix.
    if all(p -> length(terms(p)) == 1, rows)
        monomial_rows = Dict(m => (i,_action_coefficient(c)) for (i,p) in enumerate(rows) for (c,m) in terms(p))
        if length(monomial_rows) == length(rows)
            coordinates = zeros(_ActionCoefficient,length(rows),length(targets))
            for (j,p) in enumerate(targets), (c,m) in terms(p)
                haskey(monomial_rows,m) || throw(ArgumentError(
                    "$context is not closed under the supplied generator; no automatic completion is performed."))
                i,scale = monomial_rows[m]
                coordinates[i,j] += _action_coefficient(c)/scale
            end
            return coordinates
        end
    end
    words = sort!(unique([m for p in vcat(rows, targets) for (_, m) in terms(p)]))
    index = Dict(m => i for (i, m) in enumerate(words))
    n = length(rows)
    a = zeros(_ActionCoefficient, length(words), n + length(targets))
    for (j, p) in enumerate(vcat(rows, targets)), (c, m) in terms(p)
        a[index[m], j] = c
    end
    pivots = Tuple{Int,Int}[]
    r = 1
    for col in 1:n
        pivot = findfirst(i -> !iszero(a[i, col]), r:size(a, 1))
        isnothing(pivot) && continue
        k = r + pivot - 1
        a[r, :], a[k, :] = copy(a[k, :]), copy(a[r, :])
        a[r, :] ./= a[r, col]
        for i in axes(a, 1)
            (i == r || iszero(a[i,col])) && continue
            a[i, :] -= a[i, col] * a[r, :]
        end
        push!(pivots, (r, col))
        r += 1
    end
    independent && length(pivots) != n && throw(ArgumentError("$context has dependent rows."))
    all(iszero, a[r:end, n+1:end]) || throw(ArgumentError(
        "$context is not closed under the supplied generator; no automatic completion is performed."))
    coordinates = zeros(_ActionCoefficient, n, length(targets))
    for (i, j) in pivots
        coordinates[j, :] = a[i, n+1:end]
    end
    return coordinates
end

function _action_metric(metric, n)
    size(metric) == (n, n) || throw(ArgumentError("Compact metric must be $n by $n."))
    h = map(_action_coefficient, metric)
    h == h' || throw(ArgumentError("Compact metric must be Hermitian."))
    # A diagonal Hermitian form is positive definite exactly when every
    # diagonal entry is positive. Preserve that structure in later K'H+HK
    # checks instead of doing cubic rational arithmetic on structural zeros.
    if isdiag(h)
        d = diag(h)
        all(c -> isreal(c) && real(c) > 0, d) || throw(ArgumentError(
            "Compact metric must be positive definite."))
        return Diagonal(d)
    end
    # Exact Hermitian elimination: positive pivots certify positive definiteness.
    a = copy(h)
    for j in 1:n
        isreal(a[j, j]) && real(a[j, j]) > 0 || throw(ArgumentError(
            "Compact metric must be positive definite."))
        for k in j+1:n, i in j+1:n
            a[i, k] -= a[i, j] * a[j, k] / a[j, j]
        end
    end
    return h
end

# Validation is relative to the explicit positivity row spaces supplied here.
# Each localizer must be listed separately; callers must not pass their union.
# Equalities may transform within their span. Inequalities must be individually
# invariant: cone-mixing and equality-ideal membership are not implemented.
# A metric on the registered-letter span certifies a compact linear action.
# Nonlinear images remain valid inputs but cannot use that compactness witness.
# No moment map, cone decomposition, or equivalence claim is returned.
function validate_continuous_symmetry(
    pop::PolyOpt{A,T,P}, actions::AbstractVector{<:SuppliedDerivation};
    positivity_bases, compact_metric=nothing,
) where {A,T,P<:Polynomial{A,T}}
    isempty(actions) && throw(ArgumentError("Supply at least one derivation."))
    all(a isa SuppliedDerivation{A,T} for a in actions) || throw(ArgumentError(
        "Every supplied generator must use the problem's algebra and index type."))
    isempty(positivity_bases) && throw(ArgumentError("Supply every positivity block's row basis."))
    rows = [[_action_polynomial(p) for p in basis] for basis in positivity_bases]
    any(isempty, rows) && throw(ArgumentError("Positivity row bases must be nonempty."))
    domain = indices(pop.registry)
    letters = [_action_polynomial(NormalMonomial{A,T}(T[i])) for i in domain]
    h = isnothing(compact_metric) ? nothing : _action_metric(compact_metric, length(domain))
    objective = _action_polynomial(pop.objective)
    _action_zero(objective - _action_adjoint(objective), "Objective must be Hermitian")
    eqs = _action_polynomial.(pop.eq_constraints)
    meqs = _action_polynomial.(pop.moment_eq_constraints)
    ineqs = _action_polynomial.(pop.ineq_constraints)
    for (i, p) in enumerate(ineqs)
        _action_zero(p - _action_adjoint(p), "Inequality $i must be Hermitian")
    end
    for (a, action) in enumerate(actions)
        Set(keys(action.images)) == Set(domain) || throw(ArgumentError("Generator $a has a different registry domain."))
        _validate_derivation(action)
        _action_zero(differentiate(action, objective), "Generator $a objective invariance")
        for (label, constraints) in (("equality span", eqs), ("moment equality span", meqs))
            _action_coordinates(constraints, [differentiate(action, p) for p in constraints];
                context="Generator $a $label")
        end
        for (i, p) in enumerate(ineqs)
            _action_zero(differentiate(action, p), "Generator $a inequality $i invariance (cone mixing unsupported)")
        end
        for (b, basis) in enumerate(rows)
            _action_coordinates(basis, [differentiate(action, p) for p in basis];
                independent=true, context="Generator $a positivity block $b")
        end
        if h !== nothing
            k = _action_coordinates(letters, [differentiate(action, p) for p in letters];
                independent=true, context="Compact witness registered-letter span")
            all(iszero, k' * h + h * k) || throw(ArgumentError(
                "Generator $a does not preserve the supplied compact metric (column-image convention)."))
        end
    end
    return (validation=:exact, validation_scope=:supplied_positivity_bases,
        compactness=isnothing(h) ? :unverified : :verified,
        reduction=:not_performed, relation_to_input=:not_applicable, physical_target=:unchanged,
        positivity_blocks=length(rows))
end

function validate_continuous_symmetry(pop, actions; kwargs...)
    throw(ArgumentError("Continuous validation supports ordinary PolyOpt polynomials and supplied derivations in one matching built-in algebra; state/mixed-algebra actions are unsupported."))
end

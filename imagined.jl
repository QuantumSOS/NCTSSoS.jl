struct Monomial{A,T} <: AbstractMonomial
    # always guaranteed to be normal
    # 1. For non-commutative meaning site sorted
    # 2. For unipotent, projective algebra, in addition to site sorted, meaning simplified to no-consecutive
    # 3. For Pauli algebra, in addition to site sorted, meaning at most 1 pauli operator for each site
    # 4. For fermionic and bosonic algebra: meaning all creation operators on the left of anhilation operators and creation/anhilation operators are sorted according to physical site
    normal_word::Vector{T}
end

struct ComposedMonomial{AS} <:AbstractMonomial
    normal_monos::Tuple{AS} # not sure if this is correct
end

# Type needs to be more careful, perhaps only 2 MT <: AbstractMonomial , could be either Monomial{A,T} or ComposedMonomial{AS}
struct Polynomial{A,T,F} <: AbstractPolynomial
    # We lose the notion of Term here
    coeffs::Vector{F}
    monos::Vector{Monomial{A,T}} # Vector{ComposedMonomial{AS}} in general
end

struct PauliMonomial{T}
    phase::ComplexF64
    normal_word::Vector{Monomial{PauliAlgebra,T}}
end

# Suggest a better name, more mathematical
struct PhysicsMonomial{A,T} # here A can only be either FermionicAlgebra or BosonicAlgebra
    weight::Vector{Float64}
    normal_words::Vector{Monomial{A,T}}
end

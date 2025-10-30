# Symmetry Support for NCTSSoS.jl
# Based on TSSOS's symmetry implementation using SymbolicWedderburn.jl

"""
    NCVariablePermutation{V} <: SymbolicWedderburn.ByPermutations

Action type for permuting non-commutative variables.

# Fields
- `variables::V`: Vector of variables that can be permuted

# Description
Implements the `SymbolicWedderburn.ByPermutations` interface for non-commutative
polynomial variables. When a permutation acts on a monomial, it permutes the
variables while preserving their order in the monomial.

# Example
```julia
@ncpolyvar x[1:3]
action = NCVariablePermutation(x)
# Permutation (1,2) maps x[1]*x[2] to x[2]*x[1]
```
"""
struct NCVariablePermutation{V} <: SymbolicWedderburn.ByPermutations
    variables::V
end

"""
    SymmetryData{G,A,W}

Container for symmetry-related data in polynomial optimization.

# Fields
- `group::G`: Permutation group acting on variables
- `action::A`: Action of the group on variables
- `wedderburn::W`: Wedderburn decomposition result
- `adapted_bases::Vector{Vector{Vector{P}}}`: Symmetry-adapted polynomial bases for each irreducible representation

# Description
Stores all symmetry information needed to exploit group actions in
polynomial optimization. The adapted bases are organized by:
- First level: Irreducible representations
- Second level: Constraints (objective, inequalities, equalities)
- Third level: Individual basis polynomials
"""
struct SymmetryData{G,A,W,P}
    group::G
    action::A
    wedderburn::W
    adapted_bases::Vector{Vector{Vector{P}}}
end

"""
    SymbolicWedderburn.action(a::NCVariablePermutation, g::AbstractPermutation, mon::Monomial) -> Monomial

Apply a permutation to a non-commutative monomial.

# Arguments
- `a::NCVariablePermutation`: The action defining how variables are permuted
- `g::AbstractPermutation`: The group element to apply
- `mon::Monomial`: The monomial to permute

# Returns
- `Monomial`: The permuted monomial

# Description
For a monomial built from variables vars and permutation σ, each variable
in the monomial is replaced by its image under the permutation.

For example, if vars = [x₁, x₂, x₃] and we have monomial x₁²x₂x₃,
and σ = (1,2) (swaps x₁ and x₂), the result is x₂²x₁x₃.

# Example
```julia
@ncpolyvar x[1:3]
action = NCVariablePermutation(x)
mon = x[1] * x[2]  # Creates Monomial
g = perm"(1,2)"    # Swap first two variables
result = SymbolicWedderburn.action(action, g, mon)  # Returns x[2]*x[1]
```
"""
function SymbolicWedderburn.action(
    a::NCVariablePermutation,
    g::AbstractPermutation,
    mon::Monomial
)
    # If monomial is empty (constant 1), return it
    isempty(mon.vars) && return mon

    # For each variable in the monomial, find its index in the variable list
    # and map it through the permutation
    permuted_vars = map(mon.vars) do var
        # Find which variable this is
        idx = findfirst(==(var), a.variables)
        if idx === nothing
            # Variable not in our list, keep it unchanged
            return var
        end
        # Apply permutation: g(idx) gives the new index
        new_idx = g(idx)
        return a.variables[new_idx]
    end

    # Construct new monomial with permuted variables and original exponents
    return monomial(permuted_vars, mon.z)
end

"""
    normalform(mon::Monomial, group, action, sa::SimplifyAlgorithm) -> Monomial

Compute the canonical (normal) form of a monomial under group action.

# Arguments
- `mon::Monomial`: The monomial to normalize
- `group`: Permutation group acting on variables
- `action`: Action defining how group elements act on monomials
- `sa::SimplifyAlgorithm`: Simplification algorithm (for commutative groups, unipotent, etc.)

# Returns
- `Monomial`: The canonical representative of the monomial's orbit under the group action

# Description
Given a monomial and a group action, computes the lexicographically minimal
monomial in the orbit of the input monomial. This is essential for exploiting
symmetry: all monomials in the same orbit are identified by their normal form.

The normal form is computed by:
1. Applying all group elements to the monomial
2. Simplifying each result using the SimplifyAlgorithm
3. Returning the lexicographically minimal result

# Example
```julia
@ncpolyvar x[1:3]
mon = x[2] * x[3]
G = PermGroup([perm"(1,2,3)"])  # Cyclic permutation
action = NCVariablePermutation(x)
sa = SimplifyAlgorithm(comm_gps=[], is_unipotent=false, is_projective=false)
canonical = normalform(mon, G, action, sa)  # Returns x[1]*x[2]
```
"""
function normalform(
    mon::Monomial,
    group,
    action::NCVariablePermutation,
    sa::SimplifyAlgorithm
)
    # Apply all group elements and collect the simplified results
    orbit_monomials = [
        simplify(SymbolicWedderburn.action(action, g, mon), sa)
        for g in group
    ]

    # Return the lexicographically minimal monomial
    return minimum(orbit_monomials)
end

"""
    compute_symmetry_adapted_bases(
        vars::Vector{Variable},
        order::Int,
        group,
        action::NCVariablePermutation,
        sa::SimplifyAlgorithm;
        semisimple::Bool=false
    ) -> SymmetryData

Compute symmetry-adapted polynomial bases using Wedderburn decomposition.

# Arguments
- `vars::Vector{Variable}`: Variables in the polynomial optimization problem
- `order::Int`: Maximum degree of polynomials in the basis
- `group`: Permutation group acting on the variables
- `action::NCVariablePermutation`: Action defining how the group acts on variables
- `sa::SimplifyAlgorithm`: Algorithm for simplifying polynomials

# Keyword Arguments
- `semisimple::Bool=false`: If true, compute semisimple representation instead of irreducible

# Returns
- `SymmetryData`: Container with group, action, Wedderburn decomposition, and adapted bases

# Description
This function computes symmetry-adapted bases using Wedderburn decomposition:

1. Generate standard monomial basis of degree d and 2d
2. Convert to format compatible with SymbolicWedderburn
3. Compute Wedderburn decomposition to get irreducible representations
4. Extract symmetry-adapted basis for each irreducible representation
5. Convert back to FastPolynomials format

The symmetry-adapted bases have the property that moment matrices block-diagonalize
according to the irreducible representations, dramatically reducing SDP size.

# Example
```julia
@ncpolyvar x[1:3]
G = PermGroup([perm"(1,2)", perm"(1,2,3)"])  # Symmetric group S₃
action = NCVariablePermutation(x)
sa = SimplifyAlgorithm(comm_gps=[x], is_unipotent=false, is_projective=false)
sym_data = compute_symmetry_adapted_bases(x, 2, G, action, sa)
```
"""
function compute_symmetry_adapted_bases(
    vars::Vector{Variable},
    order::Int,
    group,
    action::NCVariablePermutation,
    sa::SimplifyAlgorithm;
    semisimple::Bool=false
)
    # TODO: Implementation in next step
    # This requires:
    # 1. Generate monomial bases for degree d and 2d
    # 2. Convert FastPolynomials to DynamicPolynomials (for SymbolicWedderburn)
    # 3. Call WedderburnDecomposition
    # 4. Extract adapted bases from decomposition
    # 5. Convert back to FastPolynomials
    error("compute_symmetry_adapted_bases not yet implemented")
end

"""
    supp_multi(
        p1::P,
        p2::P,
        group,
        action::NCVariablePermutation,
        sa::SimplifyAlgorithm;
        g::P=one(P)
    ) where {P<:AbstractPolynomial} -> Vector{Monomial}

Compute support of product p1 * g * p2 under group action.

# Arguments
- `p1, p2::P`: Polynomials to multiply
- `group`: Permutation group
- `action::NCVariablePermutation`: Group action
- `sa::SimplifyAlgorithm`: Simplification algorithm
- `g::P=one(P)`: Optional middle polynomial (for localizing matrices)

# Returns
- `Vector{Monomial}`: Support monomials of p1*g*p2, normalized to canonical form

# Description
Helper function for computing term sparsity with symmetry. Computes all
monomials in the support of p1*g*p2 and reduces them to normal form.
This is used when building moment and localizing matrices.
"""
function supp_multi(
    p1::P,
    p2::P,
    group,
    action::NCVariablePermutation,
    sa::SimplifyAlgorithm;
    g::P=one(P)
) where {P<:AbstractPolynomial}
    # Get all products of monomials from p1, g, p2
    support = Monomial[]

    for m1 in monomials(p1), mg in monomials(g), m2 in monomials(p2)
        # Compute product monomial
        product = neat_dot(neat_dot(m1, mg), m2)
        # Simplify and normalize
        simplified = simplify(product, sa)
        normalized = normalform(simplified, group, action, sa)
        push!(support, normalized)
    end

    return sort!(unique!(support))
end

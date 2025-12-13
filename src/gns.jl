using LinearAlgebra
# Registry-based imports from FastPolynomials for GNS reconstruction
using ..FastPolynomials:
    Monomial, Polynomial, VariableRegistry, AlgebraType,
    get_ncbasis, neat_dot, monomials, NonCommutativeAlgebra,
    indices, subregistry

"""
    extract_monomials_from_basis(basis_polys::Vector{Polynomial{A,T,C}}) where {A,T,C}

Extract monomials from a vector of single-term polynomials.

For NonCommutativeAlgebra, Pauli, Projector, and Unipotent algebras, each basis polynomial
has exactly one term. This function extracts the monomial from each polynomial.

For Fermionic/Bosonic algebras, basis polynomials may have multiple terms due to
normal ordering corrections. In these cases, we extract the leading monomial.

# Returns
Vector of monomials in the same order as the input polynomials.
"""
function extract_monomials_from_basis(basis_polys::Vector{Polynomial{A,T,C}}) where {A,T,C}
    result = Monomial{A,T}[]
    for poly in basis_polys
        poly_monomials = monomials(poly)
        if isempty(poly_monomials)
            # This shouldn't happen for a proper basis, but handle gracefully
            push!(result, Monomial{A}(T[]))
        else
            # Take the first monomial (leading term after simplification)
            push!(result, first(poly_monomials))
        end
    end
    return result
end

"""
    reconstruct(H::Matrix, registry::VariableRegistry{A,TI}, H_deg::Int; atol::Float64=1e-3)

Perform GNS (Gelfand-Naimark-Segal) reconstruction to extract finite-dimensional
matrix representations of non-commuting variables from moment data encoded in a Hankel matrix.

The GNS construction recovers matrix representations X₁, X₂, ..., Xₙ for variables
indexed by the registry from the moment matrix (Hankel matrix) of a linear functional.

# Arguments
- `H::Matrix`: Hankel matrix where H[i,j] = ⟨basis[i], basis[j]⟩ for the moment functional,
  indexed by the full monomial basis up to degree `H_deg`
- `registry::VariableRegistry{A,TI}`: Variable registry containing the non-commuting variables
  to reconstruct matrix representations for
- `H_deg::Int`: Maximum degree of monomials used to index the full Hankel matrix `H`
- `atol::Float64=1e-3`: Absolute tolerance for determining which singular values are considered non-zero.
  Singular values greater than `atol` are retained in the reconstruction.

# Returns
- `Dict{TI, Matrix{T}}`: Dictionary mapping variable indices to their matrix representations.
  The size of each matrix is determined by the number of singular values exceeding `atol`.

# Algorithm
The reconstruction follows these steps:
1. Extract the principal `(H_deg-1)` × `(H_deg-1)` block from `H`
2. Perform SVD: H_block = U S Uᵀ and keep singular values > `atol`
3. For each variable index i, construct localizing matrix Kᵢ where Kᵢ[j,k] = ⟨basis[j], xᵢ·basis[k]⟩
4. Compute matrix representation: Xᵢ = S^(-1/2) Uᵀ Kᵢ U S^(-1/2)

# References
- Klep, Šivic, Volčič (2018): "Minimizer extraction in polynomial optimization is robust"

# Example
```julia
# Create variables and a Hankel matrix from moment data
reg, (x,) = create_noncommutative_variables([("x", 1:2)])
H = [1.0  0.5  0.5;
     0.5  1.0  0.0;
     0.5  0.0  1.0]

# Reconstruct matrix representations, keeping singular values > 1e-3
matrices = reconstruct(H, reg, 1)

# Access matrix for variable x₁ (index 1)
X1_mat = matrices[first(indices(reg))]

# Use custom tolerance
matrices = reconstruct(H, reg, 1; atol=1e-6)
```
"""
function reconstruct(
    H::Matrix{T},
    registry::VariableRegistry{A,TI},
    H_deg::Int;
    atol::Float64=1e-3,
) where {T<:Number, A<:AlgebraType, TI<:Integer}
    H_deg < 0 && throw(ArgumentError("H_deg must be non-negative"))
    H_deg == 0 && throw(ArgumentError("H_deg must be positive for reconstruction"))
    atol < 0 && throw(ArgumentError("atol must be non-negative"))

    # Set hankel_deg to H_deg - 1 as required
    hankel_deg = H_deg - 1

    # Generate bases using registry
    H_basis_polys = get_ncbasis(registry, H_deg)
    hankel_basis_polys = get_ncbasis(registry, hankel_deg)

    # Extract monomials from polynomials for indexing
    H_basis = extract_monomials_from_basis(H_basis_polys)
    hankel_basis = extract_monomials_from_basis(hankel_basis_polys)

    len_H = length(H_basis)
    len_hankel = length(hankel_basis)

    size(H, 1) != len_H && throw(
        ArgumentError(
            "Hankel matrix row size $(size(H, 1)) does not match basis length $len_H",
        ),
    )
    size(H, 2) != len_H && throw(
        ArgumentError(
            "Hankel matrix column size $(size(H, 2)) does not match basis length $len_H",
        ),
    )

    hankel_block = @view H[1:len_hankel, 1:len_hankel]
    U, S, _ = svd(Matrix(hankel_block))

    # Determine output dimension based on singular values exceeding atol
    output_dim = count(s -> s > atol, S)

    if output_dim == 0
        throw(ArgumentError(
            "No singular values exceed tolerance $atol. " *
            "Maximum singular value is $(maximum(S)). " *
            "Consider decreasing atol."
        ))
    end

    # Check for flatness: rank(H) should equal rank(hankel_block)
    # This is the flat extension property
    rank_H = count(s -> s > atol, svd(H).S)
    rank_hankel = output_dim

    println("Rank of full Hankel matrix H: $rank_H (using atol=$atol)")
    println("Rank of hankel_block (degree $hankel_deg): $rank_hankel (using atol=$atol)")

    if rank_H != rank_hankel
        @warn """Flatness condition violated: rank(H) = $rank_H ≠ rank(hankel_block) = $rank_hankel
        The moment matrix is not a flat extension.
        This may lead to incorrect or rank-deficient matrix representations.
        Consider using hankel_deg = H_deg to ensure flatness."""
    else
        println("✓ Flatness condition satisfied: rank(H) = rank(hankel_block) = $rank_H")
    end

    U_trunc = U[:, 1:output_dim]
    S_trunc = S[1:output_dim]
    sqrt_S = sqrt.(S_trunc)
    sqrt_S_inv = 1 ./ sqrt_S

    println(
        "GNS reconstruction: keeping $output_dim singular values > $atol, reconstructed matrices will be $(output_dim)×$(output_dim)",
    )

    hankel_dict = hankel_entries_dict(H, H_basis)

    # Return dictionary mapping variable indices to matrices
    matrices = Dict{TI, Matrix{T}}()
    diag_inv = Diagonal(sqrt_S_inv)

    var_indices = indices(registry)
    for var_idx in var_indices
        K = construct_localizing_matrix(hankel_dict, var_idx, hankel_basis)
        X = diag_inv * (U_trunc' * K * U_trunc) * diag_inv
        matrices[var_idx] = X
        # Get variable name for display
        var_name = registry[var_idx]
        println("Variable $(var_name): constructed $(size(X)) matrix representation")
    end

    return matrices
end

"""
    hankel_entries_dict(hankel::Matrix{T}, basis::Vector{<:Monomial}) where {T <: Number} -> Dict{Monomial, T}

Create a lookup dictionary mapping monomial products to Hankel matrix entries.

This function constructs a dictionary where keys are monomials of the form
`neat_dot(u, v)` (representing the product u†·v) and values are the corresponding
Hankel matrix entries `hankel[i, j]` where `basis[i] = u` and `basis[j] = v`.

The dictionary allows efficient lookup of moment values ⟨u†·v⟩ from the Hankel
matrix structure, which is essential for constructing localizing matrices.

# Arguments
- `hankel::Matrix{T}`: Square Hankel matrix where hankel[i,j] = ⟨basis[i]†, basis[j]⟩
- `basis::Vector{<:Monomial}`: Monomial basis used to index the rows and columns of `hankel`

# Returns
- `Dict{Monomial, T}`: Dictionary mapping monomial products to their moment values

# Notes
- If multiple pairs (i,j) yield the same product u†·v, only one entry is stored
- The Hankel matrix must be square and match the basis length

# Example
```julia
reg, (x,) = create_noncommutative_variables([("x", 1:1)])
basis = extract_monomials_from_basis(get_ncbasis(reg, 2))
H = [1.0 0.5 0.3; 0.5 1.0 0.6; 0.3 0.6 1.0]
dict = hankel_entries_dict(H, basis)
# Access ⟨x†·x⟩ = ⟨x²⟩
# value = dict[neat_dot(x[1], x[1])]
```
"""
function hankel_entries_dict(hankel::Matrix{T}, basis::Vector{<:Monomial}) where {T<:Number}
    size(hankel, 1) == size(hankel, 2) ||
        throw(ArgumentError("Hankel matrix must be square, got size $(size(hankel))"))
    length(basis) == size(hankel, 1) || throw(
        ArgumentError(
            "Basis length $(length(basis)) must match Hankel size $(size(hankel, 1))",
        ),
    )

    M = eltype(basis)
    dict = Dict{M,T}()
    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            key = neat_dot(row_mono, col_mono)
            if !haskey(dict, key)
                dict[key] = hankel[i, j]
            end
        end
    end

    return dict
end

"""
    construct_localizing_matrix(hankel_dict::Dict{M,T}, var_idx::TI, basis::Vector{M}) where {A,TM,M<:Monomial{A,TM},T,TI}

Construct the localizing matrix for a given variable index using precomputed Hankel data.

The localizing matrix K for variable with index `var_idx` is defined by
K[i,j] = ⟨basis[i]†, xᵥ·basis[j]⟩, which encodes the action of right-multiplication
by the variable on the space spanned by the basis.

This matrix is a key component in the GNS construction, where it's used to compute
the matrix representation of the variable in the reconstructed Hilbert space.

# Arguments
- `hankel_dict::Dict{M, T}`: Dictionary mapping monomial products to moment values,
  typically created by `hankel_entries_dict`
- `var_idx::TI`: The variable index for which to construct the localizing matrix
- `basis::Vector{M}`: Monomial basis for indexing rows and columns of the localizing matrix

# Returns
- `Matrix{T}`: The n×n localizing matrix where n = length(basis),
  with entries K[i,j] = ⟨basis[i]†, xᵥ·basis[j]⟩

# Algorithm
For each position (i,j):
1. Create the variable monomial from var_idx
2. Compute the product: var_monomial * basis[j]
3. Form the inner product key: neat_dot(basis[i], var_monomial * basis[j])
4. Look up this key in hankel_dict, using 0 if not found

# Example
```julia
reg, (x,) = create_noncommutative_variables([("x", 1:1)])
basis = extract_monomials_from_basis(get_ncbasis(reg, 1))
H = [1.0 0.5; 0.5 1.0]
dict = hankel_entries_dict(H, basis)
var_idx = first(indices(reg))
K = construct_localizing_matrix(dict, var_idx, basis)
```
"""
function construct_localizing_matrix(
    hankel_dict::Dict{M,T},
    var_idx::TI,
    basis::Vector{M},
) where {A<:AlgebraType,TM<:Integer,M<:Monomial{A,TM},T<:Number,TI<:Integer}
    n = length(basis)
    K = zeros(T, n, n)

    # Create monomial from variable index
    var_monomial = Monomial{A}([TM(var_idx)])

    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            # Compute key = neat_dot(row_mono, var_monomial * col_mono)
            # This is adjoint(row_mono) * var_monomial * col_mono
            product = var_monomial * col_mono
            key = neat_dot(row_mono, product)
            K[i, j] = get(hankel_dict, key, zero(T))
        end
    end

    return K
end

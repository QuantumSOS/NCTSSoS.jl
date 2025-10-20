using LinearAlgebra
using ..FastPolynomials:
    Variable,
    Monomial,
    get_basis,
    monomials,
    monomial,
    neat_dot,
    SimplifyAlgorithm,
    simplify

"""
    reconstruct(H::Matrix, vars::Vector{Variable}, H_deg::Int; atol::Float64=1e-3)

Perform GNS (Gelfand-Naimark-Segal) reconstruction to extract finite-dimensional
matrix representations of non-commuting variables from moment data encoded in a Hankel matrix.

The GNS construction recovers matrix representations X₁, X₂, ..., Xₙ for variables
x₁, x₂, ..., xₙ from the moment matrix (Hankel matrix) of a linear functional.

This version uses no simplification algorithm (all variables in one commutative group).

# Arguments
- `H::Matrix`: Hankel matrix where H[i,j] = ⟨basis[i], basis[j]⟩ for the moment functional,
  indexed by the full monomial basis up to degree `H_deg`
- `vars::Vector{Variable}`: Vector of non-commuting variables to reconstruct matrix representations for
- `H_deg::Int`: Maximum degree of monomials used to index the full Hankel matrix `H`
- `atol::Float64=1e-3`: Absolute tolerance for determining which singular values are considered non-zero.
  Singular values greater than `atol` are retained in the reconstruction.

# Returns
- `Vector{Matrix}`: Vector of matrix representations, one for each variable in `vars`.
  The size of each matrix is determined by the number of singular values exceeding `atol`.

# Algorithm
The reconstruction follows these steps:
1. Extract the principal `(H_deg-1)` × `(H_deg-1)` block from `H`
2. Perform SVD: H_block = U S Uᵀ and keep singular values > `atol`
3. For each variable xᵢ, construct localizing matrix Kᵢ where Kᵢ[j,k] = ⟨basis[j], xᵢ·basis[k]⟩
4. Compute matrix representation: Xᵢ = S^(-1/2) Uᵀ Kᵢ U S^(-1/2)

# References
- Klep, Šivic, Volčič (2018): "Minimizer extraction in polynomial optimization is robust"

# Example
```julia
@ncpolyvar x y
# Construct Hankel matrix from moment data
H = [1.0  0.5  0.5;
     0.5  1.0  0.0;
     0.5  0.0  1.0]

# Reconstruct matrix representations, keeping singular values > 1e-3
X_mat, Y_mat = reconstruct(H, [x, y], 1)

# Use custom tolerance
X_mat, Y_mat = reconstruct(H, [x, y], 1; atol=1e-6)
```
"""
function reconstruct(
    H::Matrix{T},
    vars::Vector{Variable},
    H_deg::Int;
    atol::Float64=1e-3,
) where {T<:Number}
    # Create a no-simplification algorithm (all variables in one commutative group)
    sa = SimplifyAlgorithm(comm_gps=[vars], is_unipotent=false, is_projective=false)
    return reconstruct(H, vars, H_deg, sa; atol=atol)
end

"""
    reconstruct(H::Matrix, vars::Vector{Variable}, H_deg::Int, sa::SimplifyAlgorithm; atol::Float64=1e-3)

Perform GNS (Gelfand-Naimark-Segal) reconstruction to extract finite-dimensional
matrix representations of non-commuting variables from moment data encoded in a Hankel matrix,
with custom simplification rules for basis vectors.

The GNS construction recovers matrix representations X₁, X₂, ..., Xₙ for variables
x₁, x₂, ..., xₙ from the moment matrix (Hankel matrix) of a linear functional.

# Arguments
- `H::Matrix`: Hankel matrix where H[i,j] = ⟨basis[i], basis[j]⟩ for the moment functional,
  indexed by the full monomial basis up to degree `H_deg`
- `vars::Vector{Variable}`: Vector of non-commuting variables to reconstruct matrix representations for
- `H_deg::Int`: Maximum degree of monomials used to index the full Hankel matrix `H`
- `sa::SimplifyAlgorithm`: Simplification algorithm for simplifying basis vectors during reconstruction
- `atol::Float64=1e-3`: Absolute tolerance for determining which singular values are considered non-zero.
  Singular values greater than `atol` are retained in the reconstruction.

# Returns
- `Vector{Matrix}`: Vector of matrix representations, one for each variable in `vars`.
  The size of each matrix is determined by the number of singular values exceeding `atol`.

# Algorithm
The reconstruction follows these steps:
1. Extract the principal `(H_deg-1)` × `(H_deg-1)` block from `H`
2. Perform SVD: H_block = U S Uᵀ and keep singular values > `atol`
3. For each variable xᵢ, construct localizing matrix Kᵢ where Kᵢ[j,k] = ⟨basis[j], xᵢ·basis[k]⟩
4. Compute matrix representation: Xᵢ = S^(-1/2) Uᵀ Kᵢ U S^(-1/2)

Basis vectors are simplified according to the rules in `sa` during the construction of localizing matrices.

# References
- Klep, Šivic, Volčič (2018): "Minimizer extraction in polynomial optimization is robust"

# Example
```julia
@ncpolyvar x y
# Construct Hankel matrix from moment data
H = [1.0  0.5  0.5;
     0.5  1.0  0.0;
     0.5  0.0  1.0]

# Use custom simplification algorithm
sa = SimplifyAlgorithm(comm_gps=[[x], [y]], is_unipotent=true)
X_mat, Y_mat = reconstruct(H, [x, y], 1, sa; atol=1e-6)
```
"""
function reconstruct(
    H::Matrix{T},
    vars::Vector{Variable},
    H_deg::Int,
    sa::SimplifyAlgorithm;
    atol::Float64=1e-3,
) where {T<:Number}
    H_deg < 0 && throw(ArgumentError("H_deg must be non-negative"))
    H_deg == 0 && throw(ArgumentError("H_deg must be positive for reconstruction"))
    atol < 0 && throw(ArgumentError("atol must be non-negative"))

    # Set hankel_deg to H_deg - 1 as required
    hankel_deg = H_deg - 1

    H_basis = get_basis(vars, H_deg)
    hankel_basis = get_basis(vars, hankel_deg)

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

    hankel_dict = hankel_entries_dict(H, H_basis, sa)

    matrices = Matrix{T}[]
    diag_inv = Diagonal(sqrt_S_inv)
    for var in vars
        K = construct_localizing_matrix(hankel_dict, var, hankel_basis, sa)
        X = diag_inv * (U_trunc' * K * U_trunc) * diag_inv
        push!(matrices, X)
        println("Variable $(var.name): constructed $(size(X)) matrix representation")
    end

    return matrices
end

"""
    hankel_entries_dict(hankel::Matrix{T}, basis::Vector{Monomial}, sa::SimplifyAlgorithm) where {T <: Number} -> Dict{Monomial, T}

Create a lookup dictionary mapping monomial products to Hankel matrix entries,
with simplification applied to the keys.

This function constructs a dictionary where keys are simplified monomials of the form
`simplify(neat_dot(u, v), sa)` (representing the simplified product u·v) and values are
the corresponding Hankel matrix entries `hankel[i, j]` where `basis[i] = u` and `basis[j] = v`.

The dictionary allows efficient lookup of moment values ⟨u·v⟩ from the Hankel
matrix structure, which is essential for constructing localizing matrices.

# Arguments
- `hankel::Matrix{T}`: Square Hankel matrix where hankel[i,j] = ⟨basis[i], basis[j]⟩
- `basis::Vector{Monomial}`: Monomial basis used to index the rows and columns of `hankel`
- `sa::SimplifyAlgorithm`: Simplification algorithm to apply to monomial products

# Returns
- `Dict{Monomial, T}`: Dictionary mapping simplified monomial products to their moment values

# Notes
- If multiple pairs (i,j) yield the same simplified product, only one entry is stored
- The Hankel matrix must be square and match the basis length

# Example
```julia
@ncpolyvar x
basis = get_basis([x], 2)  # [1, x, x²]
H = [1.0 0.5 0.3; 0.5 1.0 0.6; 0.3 0.6 1.0]
sa = SimplifyAlgorithm(comm_gps=[[x]], is_unipotent=false, is_projective=false)
dict = hankel_entries_dict(H, basis, sa)
# Access ⟨x·x⟩ = ⟨x²⟩
value = dict[simplify(neat_dot(x, x), sa)]  # Returns H[2,2] = 1.0
```
"""
function hankel_entries_dict(
    hankel::Matrix{T}, basis::Vector{Monomial}, sa::SimplifyAlgorithm
) where {T<:Number}
    size(hankel, 1) == size(hankel, 2) ||
        throw(ArgumentError("Hankel matrix must be square, got size $(size(hankel))"))
    length(basis) == size(hankel, 1) || throw(
        ArgumentError(
            "Basis length $(length(basis)) must match Hankel size $(size(hankel, 1))",
        ),
    )

    dict = Dict{Monomial,T}()
    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            key = simplify(neat_dot(row_mono, col_mono), sa)
            if !haskey(dict, key)
                dict[key] = hankel[i, j]
            end
        end
    end

    return dict
end

"""
    construct_localizing_matrix(hankel_dict::Dict{Monomial, T}, var::Variable, basis::Vector{Monomial}, sa::SimplifyAlgorithm) where {T <: Number} -> Matrix{T}

Construct the localizing matrix for a given variable using precomputed Hankel data,
with simplification applied to basis vector products.

The localizing matrix K for variable `a` is defined by K[i,j] = ⟨basis[i], a·basis[j]⟩,
which encodes the action of right-multiplication by `a` on the space spanned by the basis.

This matrix is a key component in the GNS construction, where it's used to compute
the matrix representation of the variable in the reconstructed Hilbert space.

# Arguments
- `hankel_dict::Dict{Monomial, T}`: Dictionary mapping simplified monomial products to moment values,
  typically created by `hankel_entries_dict`
- `var::Variable`: The variable for which to construct the localizing matrix
- `basis::Vector{Monomial}`: Monomial basis for indexing rows and columns of the localizing matrix
- `sa::SimplifyAlgorithm`: Simplification algorithm to apply to monomial products

# Returns
- `Matrix{T}`: The n×n localizing matrix where n = length(basis), with entries K[i,j] = ⟨basis[i], a·basis[j]⟩

# Algorithm
For each position (i,j):
1. Compute the product a·basis[j]
2. Form the inner product key: neat_dot(basis[i], a·basis[j])
3. Simplify the key according to `sa`
4. Look up the simplified key in hankel_dict, using 0 if not found

# Example
```julia
@ncpolyvar x
basis = get_basis([x], 1)  # [1, x]
H = [1.0 0.5; 0.5 1.0]
sa = SimplifyAlgorithm(comm_gps=[[x]], is_unipotent=false, is_projective=false)
dict = hankel_entries_dict(H, basis, sa)
K = construct_localizing_matrix(dict, x, basis, sa)
# K[1,1] = ⟨1, x·1⟩ = ⟨x⟩
# K[1,2] = ⟨1, x·x⟩ = ⟨x²⟩
```
"""
function construct_localizing_matrix(
    hankel_dict::Dict{Monomial,T},
    var::Variable,
    basis::Vector{Monomial},
    sa::SimplifyAlgorithm,
) where {T<:Number}
    n = length(basis)
    K = zeros(T, n, n)

    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            key = simplify(neat_dot(row_mono, neat_dot(monomial(var), col_mono)), sa)
            K[i, j] = get(hankel_dict, key, zero(T))
        end
    end

    return K
end

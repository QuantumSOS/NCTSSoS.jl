using LinearAlgebra
using ..FastPolynomials:
    Variable, Monomial, get_basis, monomials, monomial, neat_dot, degree

"""
    reconstruct(H::Matrix, vars::Vector{Variable}, H_deg::Int, hankel_deg::Int, output_dim::Int)

Perform GNS (Gelfand-Naimark-Segal) reconstruction to extract finite-dimensional
matrix representations of non-commuting variables from moment data encoded in a Hankel matrix.

The GNS construction recovers matrix representations X₁, X₂, ..., Xₙ for variables
x₁, x₂, ..., xₙ from the moment matrix (Hankel matrix) of a linear functional.

# Arguments
- `H::Matrix`: Hankel matrix where H[i,j] = ⟨basis[i], basis[j]⟩ for the moment functional,
  indexed by the full monomial basis up to degree `H_deg`
- `vars::Vector{Variable}`: Vector of non-commuting variables to reconstruct matrix representations for
- `H_deg::Int`: Maximum degree of monomials used to index the full Hankel matrix `H`
- `hankel_deg::Int`: Degree of the principal Hankel block to use for SVD decomposition and
  localizing matrix construction (typically smaller than `H_deg` to ensure localizing matrices fit)
- `output_dim::Int`: Number of leading singular values to keep after SVD decomposition.
  This determines the size of the reconstructed matrices (output_dim × output_dim).

# Returns
- `Vector{Matrix}`: Vector of matrix representations, one for each variable in `vars`.
  Each matrix has size output_dim × output_dim.

# Algorithm
The reconstruction follows these steps:
1. Extract the principal `hankel_deg` × `hankel_deg` block from `H`
2. Perform SVD: H_block = U S Uᵀ and keep the top `output_dim` singular values
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

# Reconstruct 2×2 matrix representations keeping top 2 singular values
X_mat, Y_mat = reconstruct(H, [x, y], 1, 1, 2)
```
"""
function reconstruct(
    H::Matrix{T},
    vars::Vector{Variable},
    H_deg::Int,
    hankel_deg::Int,
    output_dim::Int,
) where {T<:Number}
    H_deg < 0 && throw(ArgumentError("total_deg must be non-negative"))
    hankel_deg < 0 && throw(ArgumentError("hankel_deg must be non-negative"))
    hankel_deg > H_deg && throw(ArgumentError("hankel_deg cannot exceed total_deg"))

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
    display(S)

    output_dim > length(S) && throw(
        ArgumentError(
            "output_dim ($output_dim) cannot exceed number of singular values ($(length(S)))",
        ),
    )
    output_dim <= 0 && throw(ArgumentError("output_dim must be positive"))

    U_trunc = U[:, 1:output_dim]
    S_trunc = S[1:output_dim]
    sqrt_S = sqrt.(S_trunc)
    sqrt_S_inv = 1 ./ sqrt_S

    println(
        "GNS reconstruction: keeping top $output_dim singular values, reconstructed matrices will be $(output_dim)×$(output_dim)",
    )

    hankel_dict = hankel_entries_dict(H, H_basis)

    matrices = Matrix{T}[]
    diag_inv = Diagonal(sqrt_S_inv)
    for var in vars
        K = construct_localizing_matrix(hankel_dict, var, hankel_basis)
        X = diag_inv * (transpose(U_trunc) * K * U_trunc) * diag_inv
        push!(matrices, X)
        println("Variable $(var.name): constructed $(size(X)) matrix representation")
    end

    return matrices
end

"""
    hankel_entries_dict(hankel::Matrix{T}, basis::Vector{Monomial}) where {T <: Number} -> Dict{Monomial, T}

Create a lookup dictionary mapping monomial products to Hankel matrix entries.

This function constructs a dictionary where keys are monomials of the form 
`neat_dot(u, v)` (representing the product u·v) and values are the corresponding 
Hankel matrix entries `hankel[i, j]` where `basis[i] = u` and `basis[j] = v`.

The dictionary allows efficient lookup of moment values ⟨u·v⟩ from the Hankel 
matrix structure, which is essential for constructing localizing matrices.

# Arguments
- `hankel::Matrix{T}`: Square Hankel matrix where hankel[i,j] = ⟨basis[i], basis[j]⟩
- `basis::Vector{Monomial}`: Monomial basis used to index the rows and columns of `hankel`

# Returns
- `Dict{Monomial, T}`: Dictionary mapping monomial products to their moment values

# Notes
- If multiple pairs (i,j) yield the same product u·v, only one entry is stored
- The Hankel matrix must be square and match the basis length

# Example
```julia
@ncpolyvar x
basis = get_basis([x], 2)  # [1, x, x²]
H = [1.0 0.5 0.3; 0.5 1.0 0.6; 0.3 0.6 1.0]
dict = hankel_entries_dict(H, basis)
# Access ⟨x·x⟩ = ⟨x²⟩
value = dict[neat_dot(x, x)]  # Returns H[2,2] = 1.0
```
"""
function hankel_entries_dict(hankel::Matrix{T}, basis::Vector{Monomial}) where {T<:Number}
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
            key = neat_dot(row_mono, col_mono)
            if !haskey(dict, key)
                dict[key] = hankel[i, j]
            end
        end
    end

    return dict
end

"""
    construct_localizing_matrix(hankel_dict::Dict{Monomial, T}, var::Variable, basis::Vector{Monomial}) where {T <: Number} -> Matrix{T}

Construct the localizing matrix for a given variable using precomputed Hankel data.

The localizing matrix K for variable `a` is defined by K[i,j] = ⟨basis[i], a·basis[j]⟩,
which encodes the action of right-multiplication by `a` on the space spanned by the basis.

This matrix is a key component in the GNS construction, where it's used to compute 
the matrix representation of the variable in the reconstructed Hilbert space.

# Arguments
- `hankel_dict::Dict{Monomial, T}`: Dictionary mapping monomial products to moment values,
  typically created by `hankel_entries_dict`
- `var::Variable`: The variable for which to construct the localizing matrix
- `basis::Vector{Monomial}`: Monomial basis for indexing rows and columns of the localizing matrix

# Returns
- `Matrix{T}`: The n×n localizing matrix where n = length(basis), with entries K[i,j] = ⟨basis[i], a·basis[j]⟩

# Algorithm
For each position (i,j):
1. Compute the product a·basis[j]
2. Form the inner product key: neat_dot(basis[i], a·basis[j])
3. Look up this key in hankel_dict, using 0 if not found

# Example
```julia
@ncpolyvar x
basis = get_basis([x], 1)  # [1, x]
H = [1.0 0.5; 0.5 1.0]
dict = hankel_entries_dict(H, basis)
K = construct_localizing_matrix(dict, x, basis)
# K[1,1] = ⟨1, x·1⟩ = ⟨x⟩
# K[1,2] = ⟨1, x·x⟩ = ⟨x²⟩
```
"""
function construct_localizing_matrix(
    hankel_dict::Dict{Monomial,T},
    var::Variable,
    basis::Vector{Monomial},
) where {T<:Number}
    n = length(basis)
    K = zeros(T, n, n)

    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            key = neat_dot(row_mono, neat_dot(monomial(var), col_mono))
            K[i, j] = get(hankel_dict, key, zero(T))
        end
    end

    return K
end

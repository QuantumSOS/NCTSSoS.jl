using LinearAlgebra
import JuMP
import JuMP.MOI
using ..FastPolynomials:
    Variable, Monomial, get_basis, monomials, monomial, neat_dot, canonicalize

"""
    MomentMatrix{T,M}

Structure containing a moment matrix and its associated basis information.

# Fields
- `matrix::Matrix{T}`: The moment matrix values, where entry (i,j) = ⟨basis[i] * basis[j]⟩
- `basis::Vector{M}`: Monomial basis indexing rows and columns of the matrix
- `clique_vars::Vector{Variable}`: Variables appearing in this clique
- `clique_index::Int`: Index of the clique in the original problem
"""
struct MomentMatrix{T,M}
    matrix::Matrix{T}
    basis::Vector{M}
    clique_vars::Vector{Variable}
    clique_index::Int
end

function Base.show(io::IO, mm::MomentMatrix)
    println(io, "MomentMatrix for clique $(mm.clique_index)")
    println(io, "  Variables: $(mm.clique_vars)")
    println(io, "  Basis size: $(length(mm.basis))")
    println(io, "  Matrix size: $(size(mm.matrix))")
end

"""
    reconstruct(H::Matrix, vars::Vector{Variable}, H_deg::Int; atol::Float64=1e-3)

Perform GNS (Gelfand-Naimark-Segal) reconstruction to extract finite-dimensional
matrix representations of non-commuting variables from moment data encoded in a Hankel matrix.

The GNS construction recovers matrix representations X₁, X₂, ..., Xₙ for variables
x₁, x₂, ..., xₙ from the moment matrix (Hankel matrix) of a linear functional.

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

    hankel_dict = hankel_entries_dict(H, H_basis)

    matrices = Matrix{T}[]
    diag_inv = Diagonal(sqrt_S_inv)
    for var in vars
        K = construct_localizing_matrix(hankel_dict, var, hankel_basis)
        X = diag_inv * (U_trunc' * K * U_trunc) * diag_inv
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

"""
    moment_matrix(result::PolyOptResult{T,P,M}; clique::Int=1) where {T,P,M}

Extract the moment matrix from a polynomial optimization result for a specific clique.

The moment matrix is a Hankel matrix indexed by a monomial basis, where entry (i,j)
corresponds to the moment ⟨basis[i] * basis[j]⟩. This matrix is used as input to
GNS reconstruction for extracting matrix representations of variables.

# Arguments
- `result::PolyOptResult`: Result from `cs_nctssos` containing the solved model
- `clique::Int=1`: Index of the clique (default: 1 for single-clique problems)

# Returns
- `MomentMatrix{T,M}`: Structure containing:
  - `matrix::Matrix{T}`: The moment matrix values
  - `basis::Vector{M}`: Monomial basis indexing the matrix
  - `clique_vars::Vector{Variable}`: Variables in this clique

# Example
```julia
@ncpolyvar x[1:2]
f = x[1]^2 + x[2]^2 - 1
pop = polyopt(f)
result = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer))

# Extract moment matrix
mm = moment_matrix(result)

# Use with GNS reconstruction
X_mats = reconstruct(mm.matrix, mm.clique_vars, 2)
```

# Notes
- For problems with term sparsity, the moment matrix may be block-structured
- The basis is derived from `result.cliques_term_sparsities[clique][1]`
- Extraction method depends on whether the problem was dualized (stored in `result.problem_type`)
"""
function moment_matrix(result, clique::Int=1)
    # Validate inputs
    num_cliques = length(result.corr_sparsity.cliques)
    if clique < 1 || clique > num_cliques
        throw(ArgumentError("Invalid clique index $clique. Problem has $num_cliques cliques."))
    end

    # Check that we have the necessary metadata
    if isnothing(result.monomap) || isnothing(result.sa)
        error("Cannot extract moment matrix: PolyOptResult is missing monomap or SimplifyAlgorithm metadata. " *
              "This may be from an older version of NCTSSoS.jl.")
    end

    # Get the term sparsity for this clique
    # First element is the moment matrix, rest are localizing matrices
    mom_term_sparsity = result.cliques_term_sparsities[clique][1]

    # Get the basis from block_bases
    full_basis = construct_full_basis(mom_term_sparsity.block_bases)

    # Extract matrix values based on problem type
    if result.problem_type == MomentPrimal
        matrix = extract_moment_matrix_primal(result.model, full_basis, result.monomap, result.sa)
    else  # SOSDual
        matrix = extract_moment_matrix_dual(result.model, full_basis, result.monomap, result.sa)
    end

    # Get clique variables
    clique_vars = [Variable(Symbol("x$i")) for i in result.corr_sparsity.cliques[clique]]

    return MomentMatrix(Matrix(matrix), full_basis, clique_vars, clique)
end

"""
    moment_matrices(result::PolyOptResult{T,P,M}) where {T,P,M}

Extract moment matrices for all cliques in the optimization result.

# Returns
- `Vector{MomentMatrix{T,M}}`: One moment matrix per clique

# Example
```julia
# For a problem with multiple cliques
mms = moment_matrices(result)
for (i, mm) in enumerate(mms)
    println("Clique \$i: variables ", mm.clique_vars)
end
```
"""
function moment_matrices(result)
    return [moment_matrix(result; clique=i) for i in 1:length(result.corr_sparsity.cliques)]
end

"""
    construct_full_basis(block_bases::Vector{Vector{M}}) where {M}

Construct full basis from block-diagonal structure.

With term sparsity, the moment matrix is block-diagonal where each block
has its own basis. This function creates the full basis that indexes the
complete block-diagonal matrix.
"""
function construct_full_basis(block_bases::Vector{Vector{M}}) where {M}
    # If no term sparsity (single block), return it directly
    length(block_bases) == 1 && return block_bases[1]

    # Otherwise, concatenate all block bases and sort uniquely
    all_basis = reduce(vcat, block_bases)
    sort!(all_basis)
    unique!(all_basis)
    return all_basis
end

"""
Extract moment matrix values from a Moment (primal) problem model.

For Moment problems, we use the monomap to reconstruct which JuMP variables
correspond to which moment matrix entries, then use their values.
"""
function extract_moment_matrix_primal(model, full_basis, monomap, sa)
    T = value_type(typeof(model))
    n = length(full_basis)
    matrix = zeros(T, n, n)

    # For each entry (i,j) in the moment matrix
    for i = 1:n, j = i:n
        # Compute the monomial: basis[i]† * basis[j]
        # This is expval(basis[i] * basis[j]) after simplification
        mono = simplify(expval(_neat_dot3(full_basis[i], one(first(full_basis)), full_basis[j])), sa)

        # Look up the corresponding JuMP variable and get its value
        if haskey(monomap, mono)
            matrix[i, j] = JuMP.value(monomap[mono])
            matrix[j, i] = matrix[i, j]  # Symmetric
        end
    end

    return Symmetric(matrix)
end

"""
Extract moment matrix values from an SOS (dual) problem model.

For SOS problems, the equality constraints encode the polynomial equality.
The dual values of these constraints correspond to moment values for each monomial
in the symmetric basis (canonicalized and unique monomials).

Following the SOS dualization procedure, we must use the symmetric basis
(canonicalized unique monomials) rather than the full unsymmetrized monomap keys.
"""
function extract_moment_matrix_dual(model, full_basis, monomap, sa)
    T = value_type(typeof(model))
    n = length(full_basis)
    matrix = zeros(T, n, n)

    # Get all equality constraints
    eq_cons = JuMP.all_constraints(model, JuMP.GenericAffExpr{T,JuMP.VariableRef}, MOI.EqualTo{T})

    if isempty(eq_cons)
        error("No equality constraints found in the model. Cannot extract moment matrix from SOS dual problem.")
    end

    # Get dual values - these correspond to moment values for each monomial
    # The dual vector is negated as in the NCTSSOS reference
    moment_values = -JuMP.dual.(eq_cons)

    # Reconstruct the symmetric basis exactly as done in sos_dualize
    # This is the canonicalized unique version of monomap keys
    unsymmetrized_basis = sort(collect(keys(monomap)))
    symmetric_basis = sort(unique(canonicalize.(unsymmetrized_basis, Ref(sa))))

    # Verify sizes match
    if length(symmetric_basis) != length(moment_values)
        error("Mismatch between symmetric basis size $(length(symmetric_basis)) and dual values $(length(moment_values))")
    end

    # For each entry (i,j) in the moment matrix
    for i = 1:n, j = i:n
        # Compute the monomial: basis[i]† * basis[j] and canonicalize it
        mono = simplify(expval(_neat_dot3(full_basis[i], one(first(full_basis)), full_basis[j])), sa)
        mono_canonical = canonicalize(mono, sa)

        # Find its index in the symmetric basis
        idx = searchsortedfirst(symmetric_basis, mono_canonical)
        if idx <= length(symmetric_basis) && symmetric_basis[idx] == mono_canonical
            matrix[i, j] = moment_values[idx]
            matrix[j, i] = matrix[i, j]  # Symmetric
        end
    end

    return Symmetric(matrix)
end

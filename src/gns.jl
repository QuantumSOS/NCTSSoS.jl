using LinearAlgebra
using ..FastPolynomials: Variable, Monomial, get_basis, monomials, monomial, neat_dot, degree

"""
    reconstruct(H::Matrix, vars::Vector{Variable}, deg::Int; rtol=1e-12) -> Vector{Matrix}

Perform GNS (Gelfand-Naimark-Segal) construction to reconstruct matrix representations 
of non-commuting variables from a Hankel matrix.

# Arguments
- `H::Matrix`: The Hankel matrix indexed by monomials up to degree `deg`
- `vars::Vector{Variable}`: Vector of non-commuting variables to reconstruct
- `deg::Int`: Maximum degree of monomial basis used to index the Hankel matrix

# Keyword Arguments  
- `rtol::Real=1e-12`: Relative tolerance for determining the rank of the Hankel matrix

# Returns
- `Vector{Matrix}`: Vector of matrix representations, one for each variable in `vars`

# Description
The GNS construction algorithm reconstructs finite-dimensional matrix representations
of non-commuting variables from moment data encoded in a Hankel matrix. The algorithm:

1. Performs SVD decomposition of the Hankel matrix H = U S U^T
2. Determines the rank and constructs the GNS basis from U * sqrt(S)
3. Constructs localizing matrices for each variable using the Hankel matrix entries
4. Computes matrix representations as X_i = sqrt(S)^(-1) * U^T * K_i * U * sqrt(S)^(-1)

# Example
```julia
@ncpolyvar x y
H = [1.0 0.5 0.5; 0.5 0.25 0.25; 0.5 0.25 0.25]  # Example Hankel matrix
matrices = reconstruct(H, [x, y], 1)
```
"""
function reconstruct(H::Matrix, vars::Vector{Variable}, deg::Int; rtol=1e-12)
    # Step 1: SVD decomposition of the Hankel matrix
    U, S, Vt = svd(H)
    
    # Determine the rank based on singular values above tolerance
    rank_H = count(s -> s > rtol * S[1], S)
    
    if rank_H == 0
        throw(ArgumentError("Hankel matrix has numerical rank 0"))
    end
    
    # Truncate to the numerical rank
    U_trunc = U[:, 1:rank_H]
    S_trunc = S[1:rank_H]
    sqrt_S = sqrt.(S_trunc)
    sqrt_S_inv = 1 ./ sqrt_S
    
    println("GNS reconstruction: Hankel matrix rank = $rank_H, reconstructed matrices will be $(rank_H)×$(rank_H)")
    
    # Step 2: Generate the full monomial basis used to index the Hankel matrix
    full_basis = get_basis(vars, deg)
    
    if length(full_basis) != size(H, 1)
        throw(ArgumentError("Hankel matrix size $(size(H)) does not match basis size $(length(full_basis))"))
    end
    
    # Step 3: Generate the reduced basis (degree ≤ deg-1) for the localizing matrices
    reduced_basis = filter(mono -> degree(mono) <= deg-1, full_basis)
    n_reduced = length(reduced_basis)
    
    # Create mapping from reduced basis indices to full basis indices
    reduced_to_full = Dict(i => findfirst(==(mono), full_basis) for (i, mono) in enumerate(reduced_basis))
    
    # Extract the relevant part of U corresponding to the reduced basis
    U_reduced = U_trunc[[reduced_to_full[i] for i in 1:n_reduced], :]
    
    println("Using reduced basis of size $n_reduced (degree ≤ $(deg-1)) for localizing matrices")
    
    # Step 4: Construct matrix representations for each variable
    matrices = Matrix{Float64}[]
    
    for (var_idx, var) in enumerate(vars)
        # Construct the localizing matrix K_i for variable var
        K_i = construct_localizing_matrix(H, var, full_basis, deg)
        
        # Step 5: Compute the matrix representation 
        # X_i = sqrt(S)^(-1) * U_reduced^T * K_i * U_reduced * sqrt(S)^(-1)
        X_i = Diagonal(sqrt_S_inv) * (transpose(U_reduced) * K_i * U_reduced) * Diagonal(sqrt_S_inv)
        
        push!(matrices, X_i)
        
        println("Variable $(var.name): constructed $(size(X_i)) matrix representation")
    end
    
    return matrices
end

"""
    construct_localizing_matrix(hankel::Matrix, var::Variable, basis::Vector{Monomial}, deg::Int) -> Matrix

Construct the localizing matrix K_i for a given variable from the Hankel matrix.

# Arguments
- `hankel::Matrix`: The Hankel matrix indexed by the full basis up to degree `deg`
- `var::Variable`: The variable to construct the localizing matrix for  
- `basis::Vector{Monomial}`: The full monomial basis used to index the Hankel matrix (up to degree `deg`)
- `deg::Int`: Maximum degree constraint

# Returns
- `Matrix`: The localizing matrix K_i where K_i[x,y] = ⟨x, a*y⟩
  The matrix has size (n_reduced × n_reduced) where n_reduced is the number of monomials of degree ≤ deg-1

# Description
The localizing matrix encodes the action of right multiplication by variable `var`.
Each entry K_i[x,y] represents ⟨x, a*y⟩ where x and y are from the reduced basis
(degree ≤ deg-1), and a is the variable.

We extract this from the Hankel matrix H where H[u,v] = ⟨u,v⟩. So K[x,y] = H[x, a*y]
if the monomial a*y exists in the full basis, otherwise K[x,y] = 0.
"""
function construct_localizing_matrix(hankel::Matrix, var::Variable, basis::Vector{Monomial}, deg::Int)
    # Create reduced basis containing only monomials of degree ≤ deg-1
    basis_reduced = filter(mono -> degree(mono) <= deg-1, basis)
    n_reduced = length(basis_reduced)
    
    K = zeros(Float64, n_reduced, n_reduced)
    
    # Create a mapping from monomials to their indices in the full basis
    basis_index = Dict(mono => idx for (idx, mono) in enumerate(basis))
    
    # Create a mapping from inner products neat_dot(u,v) to Hankel matrix positions
    # This represents the structure H[i,j] = ⟨basis[i], basis[j]⟩
    hankel_index = Dict()
    for (i, u) in enumerate(basis)
        for (j, v) in enumerate(basis)
            inner_product = neat_dot(u, v)
            hankel_index[inner_product] = (i, j)
        end
    end
    
    for x in 1:n_reduced
        for y in 1:n_reduced
            # Compute a*y (var * basis_reduced[y])
            a_times_y = neat_dot(monomial(var), basis_reduced[y])
            
            # Compute the inner product ⟨x, a*y⟩ = neat_dot(basis_reduced[x], a_times_y)
            inner_prod = neat_dot(basis_reduced[x], a_times_y)
            
            # Check if this inner product corresponds to some H[u_idx, v_idx]
            if haskey(hankel_index, inner_prod)
                u_idx, v_idx = hankel_index[inner_prod]
                # K[x,y] = ⟨x, a*y⟩ = H[u_idx, v_idx]
                K[x, y] = hankel[u_idx, v_idx]
            end
            # If the inner product is not representable in the Hankel matrix, K[x,y] remains 0
        end
    end
    
    return K
end

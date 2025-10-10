using LinearAlgebra
using ..FastPolynomials: Variable, Monomial, get_basis, monomials, monomial, neat_dot, degree

"""
    reconstruct(
        H::Matrix,
        vars::Vector{Variable},
        total_deg::Int,
        hankel_deg::Int,
        localizing_deg::Int;
        rtol::Real = 1e-12,
    ) -> Vector{Matrix}

Perform GNS (Gelfand-Naimark-Segal) reconstruction using a moment matrix `H` and
its principal sub-block indexed by monomials up to `hankel_deg`.

# Arguments
- `H::Matrix`: Full Hankel matrix indexed by monomials up to degree `total_deg`.
- `vars::Vector{Variable}`: Non-commuting variables to reconstruct.
- `total_deg::Int`: Maximum degree of monomials used to index `H`.
- `hankel_deg::Int`: Degree of the principal Hankel block used for reconstruction.
- `localizing_deg::Int`: Degree of the monomial basis used for localizing matrices.

# Keyword Arguments
- `rtol::Real = 1e-12`: Relative tolerance that determines the numerical rank.

# Returns
- `Vector{Matrix}`: One reconstructed matrix representation per variable in `vars`.
"""
function reconstruct(
    H::Matrix,
    vars::Vector{Variable},
    total_deg::Int,
    hankel_deg::Int,
    localizing_deg::Int;
    rtol::Real = 1e-12,
)
    total_deg < 0 && throw(ArgumentError("total_deg must be non-negative"))
    hankel_deg < 0 && throw(ArgumentError("hankel_deg must be non-negative"))
    localizing_deg < 0 && throw(ArgumentError("localizing_deg must be non-negative"))
    hankel_deg > total_deg &&
        throw(ArgumentError("hankel_deg cannot exceed total_deg"))
    localizing_deg > hankel_deg &&
        throw(ArgumentError("localizing_deg cannot exceed hankel_deg"))

    total_basis = get_basis(vars, total_deg)
    hankel_basis = get_basis(vars, hankel_deg)
    local_basis = get_basis(vars, localizing_deg)

    len_total = length(total_basis)
    len_hankel = length(hankel_basis)
    len_local = length(local_basis)

    size(H, 1) != len_total &&
        throw(ArgumentError("Hankel matrix row size $(size(H, 1)) does not match basis length $len_total"))
    size(H, 2) != len_total &&
        throw(ArgumentError("Hankel matrix column size $(size(H, 2)) does not match basis length $len_total"))

    hankel_block = @view H[1:len_hankel, 1:len_hankel]
    U, S, _ = svd(Matrix(hankel_block))

    rank_H = count(s -> s > rtol * S[1], S)
    rank_H == 0 && throw(ArgumentError("Hankel matrix has numerical rank 0"))

    U_trunc = U[:, 1:rank_H]
    S_trunc = S[1:rank_H]
    sqrt_S = sqrt.(S_trunc)
    sqrt_S_inv = 1 ./ sqrt_S

    println("GNS reconstruction: Hankel matrix rank = $rank_H, reconstructed matrices will be $(rank_H)×$(rank_H)")

    hankel_dict = hankel_entries_dict(Matrix(hankel_block), hankel_basis)

    matrices = Matrix{Float64}[]
    diag_inv = Diagonal(sqrt_S_inv)
    for var in vars
        K = construct_localizing_matrix(hankel_dict, var, local_basis)
        display(K)
        X = diag_inv * (transpose(U_trunc) * K * U_trunc) * diag_inv
        push!(matrices, X)
        println("Variable $(var.name): constructed $(size(X)) matrix representation")
        display(X)
    end

    return matrices
end

"""
    hankel_entries_dict(hankel::Matrix{T}, basis::Vector{Monomial}) where {T <: Number}
        -> Dict{Monomial, T}

Create a dictionary that maps monomials `neat_dot(u, v)` to the corresponding
entries of the Hankel matrix `hankel[i, j]`, where `basis[i] = u` and `basis[j] = v`.
"""
function hankel_entries_dict(
    hankel::Matrix{T},
    basis::Vector{Monomial}
) where {T<:Number}
    size(hankel, 1) == size(hankel, 2) ||
        throw(ArgumentError("Hankel matrix must be square, got size $(size(hankel))"))
    length(basis) == size(hankel, 1) ||
        throw(ArgumentError("Basis length $(length(basis)) must match Hankel size $(size(hankel, 1))"))

    dict = Dict{Monomial, T}()
    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            key = neat_dot(row_mono, col_mono)
            if !haskey(dict, key)
                dict[key] = hankel[i,j] 
            end
        end
    end

    return dict
end

"""
    construct_localizing_matrix(hankel_dict, var::Variable, basis::Vector{Monomial}) -> Matrix

Construct the localizing matrix associated with `var` using Hankel data stored in
`hankel_dict`.
"""
function construct_localizing_matrix(
    hankel_dict::Dict{Monomial, T},
    var::Variable,
    basis::Vector{Monomial},
) where {T<:Number}
    n = length(basis)
    K = zeros(T, n, n)

    for (i, row_mono) in enumerate(basis)
        for (j, col_mono) in enumerate(basis)
            key = neat_dot(row_mono, var * col_mono)
            K[i, j] = get(hankel_dict, key, zero(T))
        end
    end

    return K
end

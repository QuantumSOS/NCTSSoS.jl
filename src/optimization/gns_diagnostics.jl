"""
    RobustnessReport

Numerical conditioning and distance-to-flatness diagnostics for a dense GNS
reconstruction.
"""
struct RobustnessReport{RT<:Real}
    sigma_min::RT
    sigma_max::RT
    condition_number::RT
    dist_to_flat::RT
    operator_error_bound::RT
end

"""
    VerificationReport

Basic verification summary for a reconstructed dense GNS model.
"""
struct VerificationReport
    is_symmetric::Bool
    moment_max_error::Float64
    objective_error::Float64
    constraint_min_eigenvalues::Vector{Float64}
    ball_contained::Bool
end

@inline function _gns_expectation(op::AbstractMatrix, xi::AbstractVector)
    return dot(xi, op * xi)
end

function _gns_eval_monomial(
    matrices::Dict,
    mono::NormalMonomial,
)
    sample_matrix = first(values(matrices))
    T = eltype(sample_matrix)
    dim = size(sample_matrix, 1)
    result = Matrix{T}(I, dim, dim)
    for idx in mono.word
        result *= matrices[idx]
    end
    return result
end

function _gns_eval_polynomial(
    gns::GNSResult,
    poly::Polynomial,
)
    sample_matrix = first(values(gns.matrices))
    T = promote_type(eltype(sample_matrix), coeff_type(poly))
    dim = size(sample_matrix, 1)
    total = zeros(T, dim, dim)
    for (coef, mono) in terms(poly)
        total .+= coef .* _gns_eval_monomial(gns.matrices, mono)
    end
    return total
end

"""
    robustness_report(gns, hankel, full_basis, basis)
    robustness_report(gns, hankel)

Compute conditioning and canonical-flatness diagnostics for a dense Hankel
matrix feeding a GNS reconstruction.
"""
function robustness_report(
    gns::GNSResult{T,RT},
    hankel::AbstractMatrix,
    full_basis::Vector,
    basis::Vector,
) where {T,RT<:Real}
    sigma_max = isempty(gns.singular_values) ? zero(RT) : gns.singular_values[1]
    sigma_min = gns.rank == 0 ? zero(RT) : gns.singular_values[gns.rank]
    kappa = iszero(sigma_min) ? RT(Inf) : RT(sigma_max / sigma_min)

    hankel_flat = flat_extend(hankel, full_basis, basis)
    dist_to_flat = RT(opnorm(Matrix(hankel) - hankel_flat, 2))
    operator_error_bound = iszero(sigma_min) ? RT(Inf) : RT(dist_to_flat / sigma_min)

    return RobustnessReport{RT}(sigma_min, sigma_max, kappa, dist_to_flat, operator_error_bound)
end

function robustness_report(
    gns::GNSResult,
    hankel::AbstractMatrix,
)
    return robustness_report(gns, hankel, gns.full_basis, gns.basis)
end

"""
    verify_gns(gns, monomap, registry; poly=nothing, f_star=nothing, constraints=nothing, ball=false, atol=1e-8)

Verify moment reproduction, symmetry, and optional objective / constraint checks
for a dense GNS reconstruction.
"""
function verify_gns(
    gns::GNSResult,
    monomap::AbstractDict,
    registry::VariableRegistry;
    poly=nothing,
    f_star=nothing,
    constraints=nothing,
    ball::Bool=false,
    atol::Real=1e-8,
)
    is_symmetric = all(norm(A - A', Inf) <= atol for A in values(gns.matrices))

    moment_max_error = 0.0
    for mono in gns.full_basis
        observed = _gns_expectation(_gns_eval_monomial(gns.matrices, mono), gns.xi)
        expected = _gns_moment_value(monomap, mono)
        moment_max_error = max(moment_max_error, Float64(abs(observed - expected)))
    end

    objective_error = NaN
    if !isnothing(poly) && !isnothing(f_star)
        poly_matrix = _gns_eval_polynomial(gns, poly)
        objective_value = _gns_expectation(poly_matrix, gns.xi)
        objective_error = Float64(abs(objective_value - f_star))
    end

    constraint_min_eigenvalues = Float64[]
    if !isnothing(constraints)
        for constraint in constraints
            constraint_matrix = _gns_eval_polynomial(gns, constraint)
            hermitian_matrix = Hermitian((constraint_matrix + constraint_matrix') / 2)
            push!(constraint_min_eigenvalues, Float64(real(eigmin(hermitian_matrix))))
        end
    end

    ball_contained = true
    if ball
        sample_matrix = first(values(gns.matrices))
        T = eltype(sample_matrix)
        dim = size(sample_matrix, 1)
        ball_matrix = Matrix{T}(I, dim, dim)
        for var_idx in indices(registry)
            ball_matrix .-= gns.matrices[var_idx]^2
        end
        ball_contained = eigmin(Hermitian((ball_matrix + ball_matrix') / 2)) >= -atol
    end

    return VerificationReport(
        is_symmetric,
        moment_max_error,
        objective_error,
        constraint_min_eigenvalues,
        ball_contained,
    )
end

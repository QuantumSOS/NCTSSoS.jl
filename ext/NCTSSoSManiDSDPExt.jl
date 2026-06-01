module NCTSSoSManiDSDPExt

using LinearAlgebra
using NCTSSoS
using ManiDSDP

const MOI = NCTSSoS.MOI

const _MANIDSDP_DEP_TOL = 1e-9

struct _DiagonalGramFactor{R<:Real,V<:AbstractVector{R}}
    invdiag::V
end

Base.:\(factor::_DiagonalGramFactor, x::AbstractVector) = factor.invdiag .* x
Base.:\(factor::_DiagonalGramFactor, X::AbstractMatrix) = reshape(factor.invdiag, :, 1) .* X

# The extension is intentionally narrow: it bypasses JuMP only for Pauli
# SOS-dual relaxations whose PSD blocks are plain moment-matrix blocks.  Extra
# equality/localizing machinery changes the dual shape and must not be shoved
# through ManiDSDP by accident.
function NCTSSoS.solve_sdp(
    mp::NCTSSoS.MomentProblem{A},
    optimizer::Union{Type{ManiDSDP.Optimizer},ManiDSDP.Optimizer};
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
) where {A<:NCTSSoS.PauliAlgebra}
    opt = optimizer isa ManiDSDP.Optimizer ? optimizer : optimizer()
    return _solve_pauli_with_manidsdp(
        mp,
        opt;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
    )
end

function NCTSSoS.solve_sdp(
    mp::NCTSSoS.MomentProblem{A},
    optimizer::MOI.OptimizerWithAttributes;
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
) where {A<:NCTSSoS.PauliAlgebra}
    opt = _manidsdp_optimizer_from_attributes(optimizer)
    return _dispatch_or_fallback(
        mp,
        optimizer,
        opt;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
    )
end

function NCTSSoS.solve_sdp(
    mp::NCTSSoS.MomentProblem{A},
    optimizer::F;
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
) where {A<:NCTSSoS.PauliAlgebra,F<:Function}
    opt = _try_construct_manidsdp_optimizer(optimizer)
    return _dispatch_or_fallback(
        mp,
        optimizer,
        opt;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
    )
end

function _dispatch_or_fallback(
    mp,
    optimizer,
    opt;
    dualize::Bool,
    formulation::Symbol,
    representation::Symbol,
    orphan_policy::Symbol,
)
    if opt === nothing
        return invoke(
            NCTSSoS.solve_sdp,
            Tuple{Any,Any},
            mp,
            optimizer;
            dualize=dualize,
            formulation=formulation,
            representation=representation,
            orphan_policy=orphan_policy,
        )
    end
    return _solve_pauli_with_manidsdp(
        mp,
        opt;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
    )
end

function _solve_pauli_with_manidsdp(
    mp,
    opt::ManiDSDP.Optimizer;
    dualize::Bool,
    formulation::Symbol,
    representation::Symbol,
    orphan_policy::Symbol,
)
    dualize || throw(ArgumentError("ManiDSDP native Pauli path supports only dualize=true."))
    if formulation != :moment_variables || representation != :real || orphan_policy != :error
        throw(ArgumentError("Moment lowering options apply only with dualize=false; ManiDSDP native Pauli path solves the SOS dual directly."))
    end

    problem, constant = _pauli_moment_problem_to_manidsdp(mp, opt)
    objective, status = if problem === nothing
        constant, MOI.OPTIMAL
    else
        result = ManiDSDP.solve(problem; settings=opt.settings, opt.solve_kwargs...)
        constant + result.objective, _manidsdp_status(result.status)
    end

    return (
        objective = objective,
        model = NCTSSoS.GenericModel{typeof(objective)}(),
        n_unique_elements = mp.n_unique_moment_matrix_elements,
        status = status,
    )
end

function _manidsdp_optimizer_from_attributes(spec::MOI.OptimizerWithAttributes)
    opt = _try_construct_manidsdp_optimizer(spec.optimizer_constructor)
    opt === nothing && return nothing
    for (attr, value) in spec.params
        _apply_manidsdp_attribute!(opt, attr, value)
    end
    return opt
end

function _try_construct_manidsdp_optimizer(constructor)
    if constructor isa Type && constructor <: ManiDSDP.Optimizer
        return constructor()
    end
    opt = try
        constructor()
    catch
        return nothing
    end
    return opt isa ManiDSDP.Optimizer ? opt : nothing
end

function _apply_manidsdp_attribute!(opt::ManiDSDP.Optimizer, attr, value)
    if attr isa MOI.RawOptimizerAttribute
        MOI.set(opt, attr, value)
    elseif attr isa MOI.Silent
        Bool(value) && MOI.set(opt, MOI.RawOptimizerAttribute("verbose"), false)
    elseif attr isa AbstractString
        MOI.set(opt, MOI.RawOptimizerAttribute(attr), value)
    elseif attr isa Symbol
        MOI.set(opt, MOI.RawOptimizerAttribute(String(attr)), value)
    else
        throw(ArgumentError("unsupported ManiDSDP optimizer attribute $(attr) in NCTSSoS extension"))
    end
    return opt
end

function _pauli_moment_problem_to_manidsdp(mp, opt::ManiDSDP.Optimizer)
    L = mp.linear
    _check_supported_pauli_manidsdp_shape(L)

    blocks = L.psd_blocks_lin
    offsets = _block_offsets(blocks)
    n = isempty(blocks) ? 0 : offsets[end] + blocks[end].size - 1
    n > 0 || throw(ArgumentError("Pauli ManiDSDP path requires at least one HPSD moment block."))

    moment_count = length(L.moments)
    raw = Vector{Union{Nothing,Matrix{ComplexF64}}}(nothing, moment_count)
    for (block_idx, block) in pairs(blocks)
        offset = offsets[block_idx]
        for j in 1:block.size, i in 1:block.size
            row = offset + i - 1
            col = offset + j - 1
            for (key, coef) in block.entries[i, j]
                idx = L.moment_index[key]
                mat = raw[idx]
                if mat === nothing
                    mat = zeros(ComplexF64, n, n)
                    raw[idx] = mat
                end
                mat[row, col] += ComplexF64(coef)
            end
        end
    end

    objective = zeros(ComplexF64, moment_count)
    for (key, coef) in L.objective_lin
        objective[L.moment_index[key]] += ComplexF64(coef)
    end

    identity = L.identity
    identity_idx = L.moment_index[identity]
    C_id = raw[identity_idx] === nothing ? zeros(ComplexF64, n, n) : raw[identity_idx]::Matrix{ComplexF64}
    B0 = _hermitian_part(conj.(C_id))
    constant = real(objective[identity_idx])

    _assert_unit_diagonal_manidsdp_objective(B0)

    matrices = Matrix{ComplexF64}[]
    rhs = Float64[]

    # The imaginary identity equation is normally zero for Pauli Hamiltonians,
    # but keep it if present and harmless: it is an affine equality, not the
    # ManiDSDP objective matrix.
    B_id_im = _hermitian_part(im .* conj.(C_id))
    _append_constraint_if_nonzero!(matrices, rhs, B_id_im, imag(objective[identity_idx]))

    for idx in eachindex(L.moments)
        idx == identity_idx && continue
        C_key = raw[idx]
        c_key = objective[idx]
        if C_key === nothing
            abs(c_key) <= _MANIDSDP_DEP_TOL || throw(ArgumentError(
                "SOS objective references a moment with no moment-matrix coefficient; native ManiDSDP problem is infeasible."
            ))
            continue
        end
        _append_constraint_if_nonzero!(matrices, rhs, _hermitian_part(conj.(C_key)), real(c_key))
        _append_constraint_if_nonzero!(matrices, rhs, _hermitian_part(im .* conj.(C_key)), imag(c_key))
    end

    isempty(matrices) && return nothing, constant
    _assert_zero_constraint_diagonals(matrices)

    rank = Int(get(opt.solve_kwargs, :rank, min(n, 2)))
    problem = _orthogonal_pauli_dualsdp(-B0, [-A for A in matrices], -rhs; rank=rank)
    return problem, constant
end

function _check_supported_pauli_manidsdp_shape(L)
    isempty(L.zero_constraints) || throw(ArgumentError(
        "ManiDSDP native Pauli path currently supports pure moment-matrix SOS relaxations only; " *
        "found $(length(L.zero_constraints)) scalar zero/equality constraints."
    ))
    isempty(L.psd_blocks_lin) && throw(ArgumentError("ManiDSDP native Pauli path found no PSD blocks."))
    for (idx, block) in pairs(L.psd_blocks_lin)
        block.meta.cone == :HPSD || throw(ArgumentError(
            "ManiDSDP native Pauli path expected HPSD block $idx, got $(block.meta.cone)."
        ))
        block.meta.origin isa NCTSSoS.MomentMatrixOrigin || throw(ArgumentError(
            "ManiDSDP native Pauli path supports moment-matrix blocks only; " *
            "block $idx has origin $(typeof(block.meta.origin))."
        ))
    end
    return nothing
end

function _block_offsets(blocks)
    offsets = Vector{Int}(undef, length(blocks))
    next = 1
    for (idx, block) in pairs(blocks)
        offsets[idx] = next
        next += block.size
    end
    return offsets
end

_hermitian_part(A::AbstractMatrix) = Matrix(Hermitian((A + A') / 2))

function _append_constraint_if_nonzero!(matrices, rhs, A::AbstractMatrix, b::Real)
    norm(A) <= _MANIDSDP_DEP_TOL && abs(b) <= _MANIDSDP_DEP_TOL && return false
    norm(A) <= _MANIDSDP_DEP_TOL && throw(ArgumentError(
        "SOS coefficient equation has zero matrix but nonzero rhs $(b); the ManiDSDP native problem would be infeasible."
    ))
    push!(matrices, Matrix{ComplexF64}(A))
    push!(rhs, Float64(b))
    return true
end

function _orthogonal_pauli_dualsdp(C::AbstractMatrix, A::Vector{Matrix{ComplexF64}}, b::Vector{Float64}; rank::Integer)
    n = size(C, 1)
    1 <= rank <= n || throw(ArgumentError("rank must be in 1:size(C, 1)"))
    length(A) == length(b) || throw(ArgumentError("length(A) must match length(b)"))

    gram_diag = Vector{Float64}(undef, length(A))
    for (idx, Ai) in pairs(A)
        gram_diag[idx] = real(dot(Ai, Ai))
        gram_diag[idx] > _MANIDSDP_DEP_TOL || throw(ArgumentError(
            "ManiDSDP Pauli SOS coefficient matrix $idx is numerically zero."
        ))
    end

    gram = Diagonal(gram_diag)
    factor = _DiagonalGramFactor(inv.(gram_diag))
    D = ManiDSDP.apply_At(A, factor \ b)
    return ManiDSDP.DualSDP{ComplexF64,Float64,typeof(gram),typeof(factor)}(
        Matrix{ComplexF64}(C),
        A,
        b,
        Int(rank),
        gram,
        factor,
        _hermitian_part(D),
    )
end

function _assert_unit_diagonal_manidsdp_objective(B0::AbstractMatrix)
    for i in axes(B0, 1)
        abs(real(B0[i, i]) - 1.0) <= 1e-8 && abs(imag(B0[i, i])) <= 1e-8 ||
            throw(ArgumentError("ManiDSDP requires the SOS identity matrix to have unit real diagonal; entry $i is $(B0[i, i])."))
    end
    return nothing
end

function _assert_zero_constraint_diagonals(matrices)
    for (idx, A) in pairs(matrices)
        for i in axes(A, 1)
            abs(A[i, i]) <= 1e-8 || throw(ArgumentError(
                "ManiDSDP requires non-identity SOS coefficient matrices to have zero diagonal; " *
                "constraint $idx has diagonal entry $i = $(A[i, i])."
            ))
        end
    end
    return nothing
end

function _manidsdp_status(status::Symbol)
    status === :optimal && return MOI.OPTIMAL
    status === :max_iter && return MOI.ITERATION_LIMIT
    return MOI.OTHER_ERROR
end

end # module

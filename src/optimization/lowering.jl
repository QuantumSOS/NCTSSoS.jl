# =============================================================================
# MomentProblem -> JuMP lowering
# =============================================================================

"""
    build_jump_model(mp::MomentProblem;
        formulation=:moment_variables,
        representation=:real,
        orphan_policy=:error,
    ) -> (model, extract_monomap)

Lower a symbolic `MomentProblem` into a JuMP model without attaching an
optimizer. The returned `extract_monomap` closure reads the solved model values
and returns the canonical-moment map used by the existing direct moment path.

The default `formulation=:moment_variables, representation=:real` preserves the
historical lowering exactly: real algebras use one free real moment variable per
canonical key, and complex algebras use split real/imaginary moment variables
with Hermitian blocks embedded into real PSD cones.
"""
function build_jump_model(
    mp::MomentProblem{A,T,M,P};
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P<:Polynomial{A,T}}
    formulation in (:moment_variables, :psd_blocks) ||
        throw(ArgumentError("Unsupported formulation $(repr(formulation)); expected :moment_variables or :psd_blocks"))
    representation in (:real, :complex) ||
        throw(ArgumentError("Unsupported representation $(repr(representation)); expected :real or :complex"))
    orphan_policy in (:error, :aux_psd_free) ||
        throw(ArgumentError("Unsupported orphan_policy $(repr(orphan_policy)); expected :error or :aux_psd_free"))

    if formulation != :moment_variables
        throw(ArgumentError("formulation=:psd_blocks is not implemented yet"))
    end

    if _is_complex_problem(A)
        representation == :real ||
            throw(ArgumentError("formulation=:moment_variables currently supports representation=:real for complex algebras"))
        return _build_complex_moment_variable_model(mp)
    else
        # Real algebras ignore representation; their moment coordinates are real.
        return _build_real_moment_variable_model(mp)
    end
end

function _moment_variable_basis(mp::MomentProblem{A,T,M,P}) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    basis = [symmetric_canon(expval(m)) for m in mp.total_basis]
    sorted_unique!(basis)
    return basis
end

"""
    _build_real_moment_variable_model(mp) -> (model, extract_monomap)

Historical real-algebra direct moment lowering, factored out of
`solve_moment_problem` without semantic changes.
"""
function _build_real_moment_variable_model(
    mp::MomentProblem{A,T,M,P},
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    C = eltype(coefficients(mp.objective))
    model = GenericModel{C}()

    basis = _moment_variable_basis(mp)

    @variable(model, y[1:length(basis)], set_string_name=false)

    one_sym = symmetric_canon(expval(one(M)))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("Expected identity moment to be present in basis")
    @constraint(model, y[idx_one] == 1)  # Normalization

    monovars = Dict(zip(basis, y))

    # Add constraints.
    for (cone, mat) in mp.constraints
        dim = size(mat, 1)
        jump_mat = [
            _substitute_poly(mat[i, j], monovars)
            for i in 1:dim, j in 1:dim
        ]

        if cone == :Zero
            @constraint(model, jump_mat in Zeros())
        elseif cone == :PSD
            @constraint(model, _checked_symmetric(jump_mat; context="real PSD constraint") in PSDCone())
        else
            error("Unexpected cone type $cone for real problem")
        end
    end

    obj_expr = _substitute_poly(mp.objective, monovars)
    @objective(model, Min, obj_expr)

    extract_monomap = function ()
        return Dict(m => value(monovars[m]) for m in basis)
    end

    return model, extract_monomap
end

function _solve_real_moment_problem(
    mp::MomentProblem{A,T,M,P},
    optimizer,
    silent::Bool,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    model, extract_monomap = _build_real_moment_variable_model(mp)
    set_optimizer(model, optimizer)
    silent && set_silent(model)
    optimize!(model)
    return (
        objective=objective_value(model),
        model=model,
        monomap=extract_monomap(),
        n_unique_elements=mp.n_unique_moment_matrix_elements,
    )
end

"""
    _build_complex_moment_variable_model(mp) -> (model, extract_monomap)

Historical complex/Hermitian direct moment lowering, factored out of
`solve_moment_problem` without semantic changes.
"""
function _build_complex_moment_variable_model(
    mp::MomentProblem{A,T,M,P},
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    C = real(eltype(coefficients(mp.objective)))
    model = GenericModel{C}()

    basis = _moment_variable_basis(mp)
    n_basis = length(basis)

    @variable(model, y_re[1:n_basis], set_string_name=false)
    @variable(model, y_im[1:n_basis], set_string_name=false)

    one_sym = symmetric_canon(expval(one(M)))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("Expected identity moment to be present in basis")

    # Normalization: y[1] = 1 (real), y_im[1] = 0.
    @constraint(model, y_re[idx_one] == 1)
    @constraint(model, y_im[idx_one] == 0)

    basis_to_idx = Dict(m => i for (i, m) in enumerate(basis))

    # Add constraints with Hermitian embedding.
    for (cone, mat) in mp.constraints
        dim = size(mat, 1)

        # Keep these typed: JuMP's triangular PSD path dispatches on the matrix
        # element type, and `Matrix{Any}` would miss the cheaper Symmetric path.
        Aff = typeof(zero(C) * y_re[1])
        mat_re = Matrix{Aff}(undef, dim, dim)
        mat_im = Matrix{Aff}(undef, dim, dim)

        for i in 1:dim, j in 1:dim
            re_expr, im_expr = _substitute_complex_poly(mat[i, j], basis_to_idx, y_re, y_im)
            mat_re[i, j] = re_expr
            mat_im[i, j] = im_expr
        end

        if cone == :Zero
            # Zero cone: both real and imaginary parts must be zero.
            @constraint(model, [mat_re[i, j] for i in 1:dim, j in 1:dim] .== 0)
            @constraint(model, [mat_im[i, j] for i in 1:dim, j in 1:dim] .== 0)
        elseif cone == :HPSD
            # Hermitian PSD: embed as real 2n x 2n PSD.
            # [Re(H), -Im(H); Im(H), Re(H)] in PSD.
            embedded = [
                [mat_re[i, j] for i in 1:dim, j in 1:dim]   [-mat_im[i, j] for i in 1:dim, j in 1:dim]
                [mat_im[i, j] for i in 1:dim, j in 1:dim]   [mat_re[i, j] for i in 1:dim, j in 1:dim]
            ]
            @constraint(model, embedded in PSDCone())
        else
            error("Unexpected cone type $cone for complex problem")
        end
    end

    # `polyopt` validates that complex-algebra objectives are Hermitian, so their
    # expectation values are real. Optimize the real part only.
    obj_re, _ = _substitute_complex_poly(mp.objective, basis_to_idx, y_re, y_im)
    @objective(model, Min, obj_re)

    extract_monomap = function ()
        return Dict(m => Complex(value(y_re[i]), value(y_im[i])) for (m, i) in basis_to_idx)
    end

    return model, extract_monomap
end

function _solve_complex_moment_problem(
    mp::MomentProblem{A,T,M,P},
    optimizer,
    silent::Bool,
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P}
    model, extract_monomap = _build_complex_moment_variable_model(mp)
    set_optimizer(model, optimizer)
    silent && set_silent(model)
    optimize!(model)
    return (
        objective=objective_value(model),
        model=model,
        monomap=extract_monomap(),
        n_unique_elements=mp.n_unique_moment_matrix_elements,
    )
end

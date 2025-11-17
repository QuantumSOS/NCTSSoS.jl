# T: type of the coefficients
# monomap: map from monomials in DynamicPolynomials to variables in JuMP
struct MomentProblem{T,M,CR<:ConstraintRef,JS<:AbstractJuMPScalar}
    model::GenericModel{T}
    constraints::Vector{CR}
    monomap::Dict{M,JS}  # TODO: maybe refactor.
    sa::SimplifyAlgorithm
end

function substitute_variables(poly::P, monomap::Dict{M,JS}) where {T1,P<:AbstractPolynomial{T1},M,JS<:AbstractJuMPScalar}
    iszero(poly) ? zero(T1) * monomap[one(M)] : sum(coef * monomap[mono] for (coef, mono) in zip(coefficients(poly), monomials(poly)))
end

# cliques_cons: groups constraints according to cliques,
# global_cons: constraints that are not in any single clique
# cliques_term_sparsities: each clique, each obj/constraint, each ts_clique, each basis needed to index moment matrix
"""
    moment_relax(pop::PolyOpt{P}, corr_sparsity::CorrelativeSparsity, cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}) where {T,P<:AbstractPolynomial{T},M}

Construct a moment relaxation of a polynomial optimization problem using correlative sparsity.

# Arguments
- `pop::PolyOpt{P}`: The polynomial optimization problem to relax
- `corr_sparsity::CorrelativeSparsity`: The correlative sparsity structure defining cliques and global constraints
- `cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}`: Term sparsity information for each clique, containing block bases for moment matrix indexing

# Returns
- `MomentProblem`: A moment relaxation problem containing the JuMP model, constraint references, monomial mapping, and simplification algorithm

# Description
This function creates a semidefinite programming relaxation of the input polynomial optimization problem by:
1. Computing the total basis from all clique term sparsities
2. Creating JuMP variables for each monomial in the basis
3. Constructing moment matrix constraints for each clique and global constraint
4. Setting up the objective function using variable substitution

The relaxation exploits correlative sparsity to reduce the size of the semidefinite program by partitioning constraints into cliques and handling global constraints separately.
"""
function moment_relax(pop::PolyOpt{P}, corr_sparsity::CorrelativeSparsity, cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}) where {T,P<:AbstractPolynomial{T},M}
    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    !(T <: Real) && error("Moment relaxation is not supported for PolyOpt, use CPolyOpt")
    model = GenericModel{T}()

    sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_unipotent=pop.is_unipotent, is_projective=pop.is_projective)
    # the union of clique_total_basis
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        reduce(vcat, [
            map(monomials(poly)) do m
                simplify(expval(_neat_dot3(rol_idx, m, col_idx)), sa)
            end
            for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities) for basis in term_sparsity.block_bases for rol_idx in basis for col_idx in basis
        ])
    end...)

    # map the monomials to JuMP variables, the first variable must be 1
    # TODO: make set_string_name = false to further improve performance
    @variable(model, y[1:length(total_basis)], set_string_name = false)
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))

    # Create constraints with proper naming for later retrieval
    constraint_matrices =
        [mapreduce(vcat, zip(enumerate(cliques_term_sparsities), corr_sparsity.clq_cons)) do ((clq_idx, term_sparsities), cons_idx)
                mapreduce(vcat, enumerate(zip(term_sparsities, [one(pop.objective), corr_sparsity.cons[cons_idx]...]))) do (poly_idx, (term_sparsity, poly))
                    map(enumerate(term_sparsity.block_bases)) do (blk_idx, ts_sub_basis)
                        # poly_idx == 1 means objective (moment matrix)
                        # poly_idx > 1 means constraint poly_idx-1 (localizing matrix)
                        constraint_name = if poly_idx == 1
                            Symbol("mom_mtx_clique_$(clq_idx)_block_$(blk_idx)")
                        else
                            Symbol("loc_mtx_clique_$(clq_idx)_cons_$(poly_idx-1)_block_$(blk_idx)")
                        end

                        add_matrix_constraint!(
                            model,
                            poly,
                            ts_sub_basis,
                            monomap,
                            poly in pop.eq_constraints ? Zeros() : PSDCone(),
                            sa,
                            constraint_name
                        )
                    end
                end
            end
            map(enumerate(corr_sparsity.global_cons)) do (g_idx, global_con)
                add_matrix_constraint!(
                    model,
                    corr_sparsity.cons[global_con],
                    [one(pop.objective)],
                    monomap,
                    corr_sparsity.cons[global_con] in pop.eq_constraints ? Zeros() : PSDCone(),
                    sa,
                    Symbol("global_cons_$(g_idx)")
                )
            end]

    @objective(model, Min, mapreduce(p -> p[1] * monomap[canonicalize(expval(p[2]), sa)], +, terms(pop.objective)))

    return MomentProblem(model, constraint_matrices, monomap, sa)
end

"""
    add_matrix_constraint!(
        model::GenericModel{T1},
        poly::P,
        local_basis::Vector{M1},
        monomap::Dict{M2,JS},
        cone,
        sa::SimplifyAlgorithm,
        name::Symbol
    )

Add a named matrix constraint (moment matrix or localizing matrix) to the JuMP model.

Creates a matrix where entry (i,j) corresponds to the moment ⟨bᵢ† poly bⱼ⟩, where
bᵢ and bⱼ are basis elements. The matrix is constrained to be in the specified cone
(either PSDCone for positive semidefinite constraints, or Zeros for equality constraints).

The constraint is stored with the given name, allowing later retrieval via model[name].

# Arguments
- `model`: The JuMP model to add the constraint to
- `poly`: The polynomial multiplier (1 for moment matrices, constraint polynomial for localizing matrices)
- `local_basis`: The monomial basis for this block
- `monomap`: Dictionary mapping monomials to JuMP variables
- `cone`: The cone constraint (PSDCone() or Zeros())
- `sa`: Simplification algorithm for reducing monomials
- `name`: Name for the constraint (required, enables retrieval via model[name])

# Returns
The constraint reference

# Example
```julia
add_matrix_constraint!(model, one(f), basis, monomap, PSDCone(), sa, :mom_mtx_clique_1_block_1)
# Later retrieve: model[:mom_mtx_clique_1_block_1]
```
"""
function add_matrix_constraint!(
    model::GenericModel{T1},
    poly::P,
    local_basis::Vector{M1}, # M2 should be expval(M1)
    monomap::Dict{M2,JS},
    cone, # FIXME: which type should I use?
    sa::SimplifyAlgorithm,
    name::Symbol
) where {T,T1,P<:AbstractPolynomial{T},M1,M2,JS<:AbstractJuMPScalar}
    T_prom = promote_type(T, T1)
    moment_mtx = [
        sum([T_prom(coef) * monomap[simplify!(expval(_neat_dot3(row_idx, mono, col_idx)), sa)] for (coef, mono) in zip(coefficients(poly), monomials(poly))]) for
        row_idx in local_basis, col_idx in local_basis
    ]

    return model[name] = @constraint(model, moment_mtx in cone)
end

# T: type of the coefficients
# monomap: map from monomials in DynamicPolynomials to variables in JuMP
struct MomentProblem{T,M,CR<:ConstraintRef,JS<:AbstractJuMPScalar}
    model::GenericModel{T}
    constraints::Vector{CR}
    monomap::Dict{M,JS}  # TODO: maybe refactor.
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
function moment_relax(
    pop::PolyOpt{A,P},
    corr_sparsity::CorrelativeSparsity{A,TI,P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {A<:AlgebraType, TI<:Integer, C<:Number, P<:Polynomial{A,TI,C}, M<:Monomial{A,TI}}
    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    !(C <: Real) && error("Moment relaxation is not supported for PolyOpt with complex coefficients, use complex_moment_relax")
    model = GenericModel{C}()

    # the union of clique_total_basis
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        reduce(vcat, [
            map(monomials(poly)) do m
                expval(_neat_dot3(rol_idx, m, col_idx))
            end
            for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities) for basis in term_sparsity.block_bases for rol_idx in basis for col_idx in basis
        ])
    end...)

    # map the monomials to JuMP variables, the first variable must be 1
    # TODO: make set_string_name = false to further improve performance
    @variable(model, y[1:length(total_basis)], set_string_name = false)
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))

    constraint_matrices =
        [mapreduce(vcat, zip(cliques_term_sparsities, corr_sparsity.clq_cons)) do (term_sparsities, cons_idx)
                mapreduce(vcat, zip(term_sparsities, [one(pop.objective), corr_sparsity.cons[cons_idx]...])) do (term_sparsity, poly)
                    map(term_sparsity.block_bases) do ts_sub_basis
                        constrain_moment_matrix!(
                            model,
                            poly,
                            ts_sub_basis,
                            monomap,
                            poly in pop.eq_constraints ? Zeros() : PSDCone())
                    end
                end
            end
            map(corr_sparsity.global_cons) do global_con
                constrain_moment_matrix!(
                    model,
                    corr_sparsity.cons[global_con],
                    [one(pop.objective)],
                    monomap,
                    corr_sparsity.cons[global_con] in pop.eq_constraints ? Zeros() : PSDCone()
                )
            end]

    @objective(model, Min, mapreduce(p -> p[1] * monomap[symmetric_canon(expval(p[2]))], +, terms(pop.objective)))

    return MomentProblem(model, constraint_matrices, monomap)
end

function constrain_moment_matrix!(
    model::GenericModel{T1},
    poly::P,
    local_basis::Vector{M1}, # M2 should be expval(M1)
    monomap::Dict{M2,JS},
    cone # FIXME: which type should I use?
) where {T,T1,P<:AbstractPolynomial{T},M1,M2,JS<:AbstractJuMPScalar}
    T_prom = promote_type(T, T1)
    moment_mtx = [
        sum([T_prom(coef) * monomap[expval(_neat_dot3(row_idx, mono, col_idx))] for (coef, mono) in zip(coefficients(poly), monomials(poly))]) for
        row_idx in local_basis, col_idx in local_basis
    ]
    return @constraint(model, moment_mtx in cone)
end

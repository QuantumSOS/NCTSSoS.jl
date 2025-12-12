struct ComplexMomentProblem{T,M,P<:AbstractPolynomial{T}}
    objective::P
    # constraint matrix + type in Symbol
    constraints::Vector{Tuple{Symbol,Matrix{P}}}
    total_basis::Vector{M}
end

"""
    complex_moment_relax(cpop, corr_sparsity, cliques_term_sparsities)

Construct a complex moment relaxation (Hermitian positive semidefinite constraints).

This is the complex/Hermitian variant of moment_relax. In the new design with unified
PolyOpt{A,P}, dispatch between real and complex is based on traits, not types.

See also: [`moment_relax`](@ref)
"""
function complex_moment_relax(cpop::PolyOpt{A,P}, corr_sparsity::CorrelativeSparsity, cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}) where {A<:AlgebraType, T, P<:Polynomial{A,T}, M}
    # the union of clique_total_basis
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        reduce(vcat, [
            map(monomials(poly)) do m
                expval(_neat_dot3(rol_idx, m, col_idx))
            end
            for (poly, term_sparsity) in zip([one(cpop.objective); corr_sparsity.cons[cons_idx]], term_sparsities) for basis in term_sparsity.block_bases for rol_idx in basis for col_idx in basis
        ])
    end...)

    constraints =
        [mapreduce(vcat, zip(cliques_term_sparsities, corr_sparsity.clq_cons)) do (term_sparsities, cons_idx)
                mapreduce(vcat, zip(term_sparsities, [one(cpop.objective), corr_sparsity.cons[cons_idx]...])) do (term_sparsity, poly)
                    map(term_sparsity.block_bases) do ts_sub_basis
                        constrain_moment_matrix(
                            poly,
                            ts_sub_basis,
                            poly in cpop.eq_constraints ? :Zero : :HPSD)
                    end
                end
            end
            map(corr_sparsity.global_cons) do global_con
                constrain_moment_matrix(
                    corr_sparsity.cons[global_con],
                    [one(eltype(total_basis))],
                    corr_sparsity.cons[global_con] in cpop.eq_constraints ? :Zero : :HPSD
                )
            end]

    return ComplexMomentProblem(cpop.objective, constraints, total_basis)
end

function constrain_moment_matrix(
    poly::P,
    local_basis::Vector{M},
    cone::Symbol
) where {T,P<:AbstractPolynomial{T},M}
    moment_mtx = [
        sum(coef * _neat_dot3(row_idx, mono, col_idx) for (coef, mono) in terms(poly)) for
        row_idx in local_basis, col_idx in local_basis
    ]
    return (cone, moment_mtx)
end

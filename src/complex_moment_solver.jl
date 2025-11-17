struct ComplexMomentProblem{T,M,P<:AbstractPolynomial{T}}
    objective::P
    # constraint matrix + type in Symbol
    constraints::Vector{Tuple{Symbol,Matrix{P}}}
    total_basis::Vector{M}
    sa::SimplifyAlgorithm
end

function moment_relax(cpop::ComplexPolyOpt{P}, corr_sparsity::CorrelativeSparsity, cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}) where {T,P<:AbstractPolynomial{T},M}
    sa = SimplifyAlgorithm(comm_gps=cpop.comm_gps, is_unipotent=cpop.is_unipotent, is_projective=cpop.is_projective)
    # the union of clique_total_basis
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        reduce(vcat, [
            map(monomials(poly)) do m
                simplify!(expval(_neat_dot3(rol_idx, m, col_idx)), sa)
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
                            poly in cpop.eq_constraints ? :Zero : :HPSD, sa)
                    end
                end
            end
            map(corr_sparsity.global_cons) do global_con
                constrain_moment_matrix(
                    corr_sparsity.cons[global_con],
                    [one(eltype(total_basis))],
                    corr_sparsity.cons[global_con] in cpop.eq_constraints ? :Zero : :HPSD,
                    sa
                )
            end]

    simplified_objective = sum(c * simplify(m, sa) for (c, m) in terms(cpop.objective))
    return ComplexMomentProblem(simplified_objective, constraints, total_basis, sa)
end

function constrain_moment_matrix(
    poly::P,
    local_basis::Vector{M},
    cone::Symbol,
    sa::SimplifyAlgorithm
) where {T,P<:AbstractPolynomial{T},M}
    moment_mtx = [
        sum(coef * simplify!(_neat_dot3(row_idx, mono, col_idx), sa) for (coef, mono) in terms(poly)) for
        row_idx in local_basis, col_idx in local_basis
    ]
    return (cone, moment_mtx)
end

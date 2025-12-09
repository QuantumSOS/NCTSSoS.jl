"""
    CorrelativeSparsity

Structure representing the correlative sparsity pattern of a polynomial optimization problem.

# Fields
- `cliques::Vector{Vector{Variable}}`: Groups of variables that form cliques in the sparsity graph
- `cons::Vector{P}`: All constraints in the problem
- `clq_cons::Vector{Vector{Int}}`: Constraint indices assigned to each clique, regardless of equality or inequality
- `global_cons::Vector{Int}`: Constraint indices not captured by any single clique
- `clq_mtx_basis::Vector{Vector{Monomial}}`: Monomial bases for moment/localizing matrices within each clique
"""
struct CorrelativeSparsity{P,M}
    cliques::Vector{Vector{Variable}}
    cons::Vector{P} # making sure context of `Int` in following variables are clear
    clq_cons::Vector{Vector{Int}}
    global_cons::Vector{Int}
    clq_mom_mtx_bases::Vector{Vector{M}}
    clq_localizing_mtx_bases::Vector{Vector{Vector{M}}}
end

function Base.show(io::IO, cs::CorrelativeSparsity)
    max_size = maximum(length.(cs.cliques))
    println(io, "Correlative Sparsity: \n")
    println(io, "   maximum size: $max_size")
    for clique_i in 1:length(cs.cliques)
        println(io, "   Clique $clique_i: ")
        println(io, "       Variables: ", cs.cliques[clique_i])
        println(io, "       Bases length: ", length(cs.clq_mom_mtx_bases[clique_i]))
        println(io, "       Constraints: ")
        for cons_j in eachindex(cs.clq_cons[clique_i])
            println(io, "           ", cs.cons[cons_j], " :  with $(length(cs.clq_localizing_mtx_bases[clique_i][cons_j])) basis monomials")
        end
    end
    println(io, "   Global Constraints: ")
    for geq_cons in cs.global_cons
        println(io, "     $(geq_cons)")
    end
end

"""
    get_correlative_graph(ordered_vars::Vector{Variable}, obj::Polynomial{T}, cons::Vector{Polynomial{T}}, order::Int) where {T}

Constructs a correlative sparsity graph from polynomial optimization problem components.

# Arguments
- `ordered_vars::Vector{Variable}`: Variables in the order to appear in the graph
- `obj::Polynomial{T}`: Objective polynomial
- `cons::Vector{Polynomial{T}}`: Constraint polynomials
- `order::Int`: Order of the moment relaxation

# Returns
- `SimpleGraph`: Graph representing variable correlations
"""
function get_correlative_graph(ordered_vars::Vector{Variable}, obj::P, cons::Vector{P}) where {T,P<:AbstractPolynomial{T}}
    # NOTE: code will be buggy is ordered_vars is not the same as the one reference in other functions
    @assert issorted(ordered_vars) "Variables must be sorted"

    nvars = length(ordered_vars)
    G = SimpleGraph(nvars)

    findvar(v) = searchsortedfirst(ordered_vars, v)

    map(mono -> add_clique!(G, findvar.(variables(mono))), monomials(obj))
    map(con -> add_clique!(G, findvar.(variables(con))), cons)
    return G
end

"""
    assign_constraint(cliques::Vector{Vector{Variable}}, cons::Vector{Polynomial{T}}) where {T}

Assigns constraints to cliques based on variable support.

# Arguments
- `cliques::Vector{Vector{Variable}}`: Variable cliques
- `cons::Vector{P}`: Constraint polynomials

# Returns
- `Tuple{Vector{Vector{Int}}, Vector{Int}}`: Tuple containing:
  - Constraint indices for each clique
  - Global constraint indices not captured by any single clique
"""
function assign_constraint(cliques::Vector{Vector{Variable}}, cons::Vector{P}) where {T,P<:AbstractPolynomial{T}}
    # assign each constraint to a clique
    # there might be constraints that are not captured by any single clique,
    # NOTE: we ignore this constraint. This should only occur at lower order of relaxation.

    # clique_cons: vector of vector of constraints index, each belong to a clique
    clique_cons = map(cliques) do clique
        findall(g -> issubset(unique!(variables(g)), clique), cons)
    end
    return clique_cons, setdiff(1:length(cons), union(clique_cons...))
end

"""
    correlative_sparsity(pop::PolyOpt{T}, order::Int, elim_algo::EliminationAlgorithm) where {T}

Decomposes a polynomial optimization problem into a correlative sparsity pattern by identifying
variable cliques and assigning constraints to cliques, enabling block-structured semidefinite relaxations.

# Arguments
- `pop::PolyOpt{T}`: Polynomial optimization problem containing objective, constraints, and variables
- `order::Int`: Order of the moment relaxation
- `elim_algo::EliminationAlgorithm`: Algorithm for clique tree elimination

# Returns
- `CorrelativeSparsity`: Structure containing:
  - `cliques`: Groups of variables that form cliques in the sparsity graph
  - `cliques_cons`: Constraint indices assigned to each clique
  - `global_cons`: Constraint indices not captured by any single clique
  - `cliques_idcs_bases`: Monomial bases for indexing moment/localizing matrices within each clique
"""
function correlative_sparsity(pop::OP, order::Int, elim_algo::EliminationAlgorithm) where {T,P<:AbstractPolynomial{T},OP<:OptimizationProblem{P}}
    all_cons = vcat(pop.eq_constraints, pop.ineq_constraints)
    cliques = map(x -> sort(pop.variables[x]), clique_decomp(get_correlative_graph(pop.variables, pop.objective, all_cons), elim_algo))

    cliques_cons, global_cons = assign_constraint(cliques, all_cons)

    cliques_moment_matrix_bases = map(cliques) do clique
        sorted_unique(get_basis(clique, order))
    end

    cliques_moment_matrix_bases_dg = map(bs -> NCTSSoS.FastPolynomials.degree.(bs), cliques_moment_matrix_bases)

    cliques_idx_bases = map(zip(eachindex(cliques), cliques_cons)) do (clique_idx, clique_cons)
        if length(clique_cons) == 0
            return eltype(cliques_moment_matrix_bases)[]  # Return empty array of the same type
        end
        # get the basis of the moment matrix in a clique, then sort it
        cur_orders = order .- cld.(maxdegree.(all_cons[clique_cons]), 2)
        cur_lengths = map(o -> searchsortedfirst(cliques_moment_matrix_bases_dg[clique_idx], o + 1) - 1, cur_orders)
        map(cur_lengths) do len
            iszero(len) ? eltype(cliques_moment_matrix_bases)[] : cliques_moment_matrix_bases[clique_idx][1:len]
        end
    end

    return CorrelativeSparsity(cliques, all_cons, cliques_cons, global_cons, cliques_moment_matrix_bases, cliques_idx_bases)
end

"""
    TermSparsity

Structure representing term sparsity information for polynomial optimization.

# Fields
- `term_sparse_graph_supp::Vector{Monomial}`: Support of the term sparsity graph
- `block_bases::Vector{Vector{Monomial}}`: Bases of moment/localizing matrices in each block
"""
struct TermSparsity{M}
    term_sparse_graph_supp::Vector{M}
    block_bases::Vector{Vector{M}}
end

function debug(ts::TermSparsity)
    println("Term Sparse Graph Support", ts.term_sparse_graph_supp)
    println("Block bases", ts.block_bases)
end

function Base.show(io::IO, sparsity::TermSparsity)
    println(io, "Number of Activated supp:   ", length(sparsity.term_sparse_graph_supp))
    println(io, "Number of Bases Activated in each sub-block", length.(sparsity.block_bases))
end

function init_activated_supp(partial_obj::P, cons::Vector{P}, mom_mtx_bases::Vector{M}) where {T,P<:AbstractPolynomial{T},M}
    return sorted_union(symmetric_canon.(monomials(partial_obj)), mapreduce(monomials, vcat, cons; init=M[]), [neat_dot(b, b) for b in mom_mtx_bases])
end

"""
    term_sparsities(initial_activated_supp::Vector{M}, cons::Vector{P}, mom_mtx_bases::Vector{M}, localizing_mtx_bases::Vector{Vector{M}}, ts_algo::EliminationAlgorithm) where {T,P<:AbstractPolynomial{T},M}

Computes term sparsity structures for the moment matrix and all localizing matrices.

# Arguments
- `initial_activated_supp::Vector{M}`: Initial set of activated support monomials
- `cons::Vector{P}`: Vector of constraint polynomials
- `mom_mtx_bases::Vector{M}`: Basis monomials for the moment matrix
- `localizing_mtx_bases::Vector{Vector{M}}`: Basis monomials for each localizing matrix corresponding to constraints
- `ts_algo::EliminationAlgorithm`: Algorithm for clique tree elimination in term sparsity graphs

# Returns
- `Vector{TermSparsity}`: Vector containing term sparsity structures, with the first element corresponding to the moment matrix and subsequent elements corresponding to localizing matrices for each constraint
"""
function term_sparsities(initial_activated_supp::Vector{M}, cons::Vector{P}, mom_mtx_bases::Vector{M}, localizing_mtx_bases::Vector{Vector{M}}, ts_algo::EliminationAlgorithm) where {T,P<:AbstractPolynomial{T},M}
    [
        [iterate_term_sparse_supp(initial_activated_supp, one(P), mom_mtx_bases, ts_algo)];
        map(zip(cons, localizing_mtx_bases)) do (poly, basis)
            iterate_term_sparse_supp(initial_activated_supp, poly, basis, ts_algo)
        end
    ]
end

"""
    get_term_sparsity_graph(cons_support::Vector{Monomial}, activated_supp::Vector{Monomial}, basis::Vector{Monomial})

Constructs a term sparsity graph for polynomial constraints.

# Arguments
- `cons_support::Vector{Monomial}`: Support monomials of constraints
- `activated_supp::Vector{Monomial}`: Support from previous iterations
- `bases::Vector{Monomial}`: Basis used to index the moment matrix

# Returns
- `SimpleGraph`: Term sparsity graph
"""
function get_term_sparsity_graph(cons_support::Vector{M}, activated_supp::Vector{M}, bases::Vector{M}) where {M}
    nterms = length(bases)
    G = SimpleGraph(nterms)
    sorted_activated_supp = sort(activated_supp)
    for i in 1:nterms, j in i+1:nterms
        for supp in cons_support
            connected_mono_lr = _neat_dot3(bases[i], supp, bases[j])
            connected_mono_rl = _neat_dot3(bases[j], supp, bases[i])
            expval_cm_lr = expval(connected_mono_lr) * one(Monomial)
            expval_cm_rl = expval(connected_mono_rl) * one(Monomial)
            if expval_cm_lr in sorted_activated_supp || expval_cm_rl in sorted_activated_supp
                add_edge!(G, i, j)
                continue
            end
        end
    end
    return G
end

"""
    iterate_term_sparse_supp(activated_supp::Vector{Monomial}, poly::Polynomial, basis::Vector{Monomial}, elim_algo::EliminationAlgorithm)

Iteratively computes term sparsity support for a polynomial.

# Arguments
- `activated_supp::Vector{Monomial}`: Currently activated support monomials
- `poly::Polynomial`: Input polynomial
- `basis::Vector{Monomial}`: Basis monomials
- `elim_algo::EliminationAlgorithm`: Elimination algorithm for clique decomposition

# Returns
- `TermSparsity`: Term sparsity structure containing graph support and block bases
"""
function iterate_term_sparse_supp(activated_supp::Vector{M}, poly::P, basis::Vector{M}, elim_algo::EliminationAlgorithm) where {M,P}
    F = get_term_sparsity_graph(monomials(poly), activated_supp, basis)
    blocks = clique_decomp(F, elim_algo)
    map(block -> add_clique!(F, block), blocks)
    return TermSparsity(term_sparsity_graph_supp(F, basis, poly), map(x -> basis[x], blocks))
end

"""
    term_sparsity_graph_supp(G::SimpleGraph, basis::Vector{Monomial}, g::Polynomial)

Computes the support of a term sparsity graph for a given polynomial.

# Arguments
- `G::SimpleGraph`: Term sparsity graph
- `basis::Vector{Monomial}`: Basis monomials
- `g::Polynomial`: Input polynomial

# Returns
- `Vector{Monomial}`: Support monomials for the term sparsity graph
"""
function term_sparsity_graph_supp(G::SimpleGraph, basis::Vector{M}, g::P) where {M,P}
    # following (10.4) in Sparse Polynomial Optimization: Theory and Practise
    # NOTE: Do I need to symmetric canonicalize it?
    gsupp(a, b) = map(g_supp -> _neat_dot3(a, g_supp, b), monomials(g))
    return union([gsupp(basis[v], basis[v]) for v in vertices(G)]..., [gsupp(basis[e.src], basis[e.dst]) for e in edges(G)]...)
end

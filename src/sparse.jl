"""
    CorrelativeSparsity{A<:AlgebraType, T<:Integer, P<:Polynomial{A,T}, M<:Monomial{A,T}}

Structure representing the correlative sparsity pattern of a polynomial optimization problem.

# Type Parameters
- `A<:AlgebraType`: The algebra type (PauliAlgebra, NonCommutativeAlgebra, etc.)
- `T<:Integer`: The index type used in the registry
- `P<:Polynomial{A,T}`: The polynomial type
- `M<:Monomial{A,T}`: The monomial type

# Fields
- `cliques::Vector{Vector{T}}`: Groups of variable indices that form cliques in the sparsity graph
- `registry::VariableRegistry{A,T}`: Variable registry for symbol lookups and sub-registry creation
- `cons::Vector{P}`: All constraints in the problem
- `clq_cons::Vector{Vector{Int}}`: Constraint indices assigned to each clique, regardless of equality or inequality
- `global_cons::Vector{Int}`: Constraint indices not captured by any single clique
- `clq_mom_mtx_bases::Vector{Vector{M}}`: Monomial bases for moment matrices within each clique
- `clq_localizing_mtx_bases::Vector{Vector{Vector{M}}}`: Monomial bases for localizing matrices within each clique

# Notes
- Cliques store variable indices (type `T`), not symbols
- Use `subregistry(cs.registry, cs.cliques[i])` to get a clique-local registry
- Use `cs.registry[idx]` to look up symbol names for indices
"""
struct CorrelativeSparsity{A<:AlgebraType, T<:Integer, P<:Polynomial{A,T}, M<:Monomial{A,T}}
    cliques::Vector{Vector{T}}
    registry::VariableRegistry{A,T}
    cons::Vector{P}
    clq_cons::Vector{Vector{Int}}
    global_cons::Vector{Int}
    clq_mom_mtx_bases::Vector{Vector{M}}
    clq_localizing_mtx_bases::Vector{Vector{Vector{M}}}
end

function Base.show(io::IO, cs::CorrelativeSparsity{A,T,P,M}) where {A,T,P,M}
    max_size = isempty(cs.cliques) ? 0 : maximum(length.(cs.cliques))
    println(io, "Correlative Sparsity ($(nameof(A))): \n")
    println(io, "   maximum clique size: $max_size")
    println(io, "   number of cliques: $(length(cs.cliques))")
    for clique_i in 1:length(cs.cliques)
        println(io, "   Clique $clique_i: ")
        # Display variable symbols from registry
        clique_syms = [cs.registry[idx] for idx in cs.cliques[clique_i]]
        println(io, "       Variables: ", clique_syms)
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
    get_correlative_graph(registry::VariableRegistry{A,T}, obj::P, cons::Vector{P}) where {A<:AlgebraType, T<:Integer, P<:Polynomial{A,T}}

Constructs a correlative sparsity graph from polynomial optimization problem components.

The graph nodes are numbered 1 to N where N is the number of unique variable indices.
A position-based mapping is used to convert variable indices to graph node positions.

# Arguments
- `registry::VariableRegistry{A,T}`: Variable registry containing all problem variables
- `obj::P`: Objective polynomial
- `cons::Vector{P}`: Constraint polynomials

# Returns
- `Tuple{SimpleGraph, Vector{T}, Dict{T,Int}}`: Tuple containing:
  - `SimpleGraph`: Graph representing variable correlations
  - `sorted_indices::Vector{T}`: All unique indices in sorted order
  - `idx_to_node::Dict{T,Int}`: Mapping from variable index to graph node position

# Notes
- Uses `variable_indices()` to extract indices from polynomials
- Graph nodes are 1:N positions, not raw variable indices
- For signed index types (Fermionic/Bosonic), uses `abs(idx)` for consistency
"""
function get_correlative_graph(
    registry::VariableRegistry{A,T},
    obj::P,
    cons::Vector{P}
) where {A<:AlgebraType, T<:Integer, C, P<:Polynomial{A,T,C}}
    # Collect all unique variable indices from objective and constraints
    all_indices = Set{T}()
    union!(all_indices, variable_indices(obj))
    for c in cons
        union!(all_indices, variable_indices(c))
    end

    # Sort indices for consistent ordering
    sorted_indices = sort(collect(all_indices))
    nvars = length(sorted_indices)

    # Create mapping: variable index -> graph node position (1:N)
    idx_to_node = Dict{T,Int}(idx => pos for (pos, idx) in enumerate(sorted_indices))

    # Build correlative graph
    G = SimpleGraph(nvars)

    # Helper to get graph positions from a polynomial
    function get_positions(p::Polynomial)
        idxs = variable_indices(p)
        return [idx_to_node[idx] for idx in idxs]
    end

    # Helper to get graph positions from a monomial
    function get_positions(m::Monomial{A,T}) where {A,T}
        positions = Int[]
        seen = Set{T}()
        for idx in m.word
            abs_idx = T(abs(idx))
            if abs_idx in seen
                continue
            end
            push!(seen, abs_idx)
            if haskey(idx_to_node, abs_idx)
                push!(positions, idx_to_node[abs_idx])
            end
        end
        return positions
    end

    # Add cliques from objective monomials
    for mono in monomials(obj)
        positions = get_positions(mono)
        add_clique!(G, positions)
    end

    # Add cliques from constraint monomials
    for con in cons
        positions = get_positions(con)
        add_clique!(G, positions)
    end

    return G, sorted_indices, idx_to_node
end

"""
    assign_constraint(cliques::Vector{Vector{T}}, cons::Vector{P}, registry::VariableRegistry{A,T}) where {A<:AlgebraType, T<:Integer, C, P<:Polynomial{A,T,C}}

Assigns constraints to cliques based on variable support (index-based).

# Arguments
- `cliques::Vector{Vector{T}}`: Variable index cliques
- `cons::Vector{P}`: Constraint polynomials
- `registry::VariableRegistry{A,T}`: Variable registry (for future extensions)

# Returns
- `Tuple{Vector{Vector{Int}}, Vector{Int}}`: Tuple containing:
  - Constraint indices for each clique
  - Global constraint indices not captured by any single clique

# Notes
- Uses `variable_indices()` to get indices from constraints
- A constraint is assigned to a clique if all its variable indices are in the clique
"""
function assign_constraint(
    cliques::Vector{Vector{T}},
    cons::Vector{P},
    registry::VariableRegistry{A,T}
) where {A<:AlgebraType, T<:Integer, C, P<:Polynomial{A,T,C}}
    # For each clique, find constraints whose variables are all in the clique
    clique_sets = [Set{T}(clq) for clq in cliques]

    clique_cons = map(clique_sets) do clique_set
        findall(cons) do con
            con_indices = variable_indices(con)
            issubset(con_indices, clique_set)
        end
    end

    # Global constraints are those not assigned to any clique
    assigned = isempty(clique_cons) ? Int[] : union(clique_cons...)
    global_cons = setdiff(1:length(cons), assigned)

    return clique_cons, collect(global_cons)
end

"""
    extract_monomials_from_basis(basis_polys::Vector{Polynomial{A,T,C}}) where {A,T,C}

Extract monomials from basis polynomials for use in moment matrix construction.

For "simple" algebras (NonCommutative, Pauli, Projector, Unipotent), each basis polynomial
is a single-term polynomial, so we extract the monomial directly.

For "multi-term" algebras (Fermionic, Bosonic), basis polynomials may have multiple terms
due to normal ordering. In this case, we use the leading monomial.

# Arguments
- `basis_polys::Vector{Polynomial{A,T,C}}`: Basis polynomials from `get_ncbasis`

# Returns
- `Vector{Monomial{A,T}}`: Extracted monomials for moment matrix indexing
"""
function extract_monomials_from_basis(basis_polys::Vector{Polynomial{A,T,C}}) where {A<:AlgebraType, T<:Integer, C<:Number}
    result = Monomial{A,T}[]
    for p in basis_polys
        if isempty(p.terms)
            # Zero polynomial - shouldn't happen but handle gracefully
            continue
        end
        # Use the first (leading) monomial from the polynomial
        # For simple algebras, this is the only monomial
        push!(result, p.terms[1].monomial)
    end
    return result
end

"""
    correlative_sparsity(pop::OP, order::Int, elim_algo::EliminationAlgorithm) where {A<:AlgebraType, T, P<:Polynomial{A,T}, OP<:OptimizationProblem{A,P}}

Decomposes a polynomial optimization problem into a correlative sparsity pattern by identifying
variable cliques and assigning constraints to cliques, enabling block-structured semidefinite relaxations.

# Arguments
- `pop::OptimizationProblem{A,P}`: Polynomial optimization problem containing objective, constraints, and registry
- `order::Int`: Order of the moment relaxation
- `elim_algo::EliminationAlgorithm`: Algorithm for clique tree elimination

# Returns
- `CorrelativeSparsity{A,T,P,M}`: Structure containing:
  - `cliques`: Groups of variable indices forming cliques
  - `registry`: Variable registry for symbol lookups
  - `cons`: All constraints
  - `clq_cons`: Constraint indices assigned to each clique
  - `global_cons`: Constraint indices not captured by any single clique
  - `clq_mom_mtx_bases`: Monomial bases for moment matrices
  - `clq_localizing_mtx_bases`: Monomial bases for localizing matrices

# Notes
- Uses `subregistry()` and `get_ncbasis()` for basis generation
- Extracts monomials from basis polynomials for matrix indexing
"""
function correlative_sparsity(
    pop::OP,
    order::Int,
    elim_algo::EliminationAlgorithm
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, OP<:OptimizationProblem{A,P}}
    registry = pop.registry
    all_cons = vcat(pop.eq_constraints, pop.ineq_constraints)

    # Build correlative graph and get index mappings
    G, sorted_indices, idx_to_node = get_correlative_graph(registry, pop.objective, all_cons)

    # Decompose graph into cliques
    clique_node_sets = clique_decomp(G, elim_algo)

    # Convert clique node positions back to variable indices
    node_to_idx = Dict{Int,T}(pos => idx for (idx, pos) in idx_to_node)
    cliques = map(clique_node_sets) do node_set
        sort([node_to_idx[node] for node in node_set])
    end

    # Assign constraints to cliques
    cliques_cons, global_cons = assign_constraint(cliques, all_cons, registry)

    # Generate moment matrix bases for each clique using subregistry + get_ncbasis
    M = Monomial{A,T}
    cliques_moment_matrix_bases = Vector{Vector{M}}()

    for clique_indices in cliques
        # Create sub-registry for this clique
        sub_reg = subregistry(registry, clique_indices)
        # Generate basis polynomials
        basis_polys = get_ncbasis(sub_reg, order)
        # Extract monomials from polynomials
        basis_monos = extract_monomials_from_basis(basis_polys)
        # Sort and deduplicate
        push!(cliques_moment_matrix_bases, sorted_unique(basis_monos))
    end

    # Compute degrees for localizing matrix basis truncation
    cliques_moment_matrix_bases_dg = [NCTSSoS.FastPolynomials.degree.(bs) for bs in cliques_moment_matrix_bases]

    # Generate localizing matrix bases (truncated based on constraint degree)
    cliques_localizing_bases = map(zip(eachindex(cliques), cliques_cons)) do (clique_idx, clique_cons_indices)
        if isempty(clique_cons_indices)
            return Vector{M}[]
        end

        # For each constraint, compute the reduced order basis
        cur_orders = order .- cld.(maxdegree.(all_cons[clique_cons_indices]), 2)
        cur_lengths = map(cur_orders) do o
            searchsortedfirst(cliques_moment_matrix_bases_dg[clique_idx], o + 1) - 1
        end

        map(cur_lengths) do len
            if iszero(len)
                M[]
            else
                cliques_moment_matrix_bases[clique_idx][1:len]
            end
        end
    end

    return CorrelativeSparsity{A,T,P,M}(
        cliques,
        registry,
        all_cons,
        cliques_cons,
        global_cons,
        cliques_moment_matrix_bases,
        cliques_localizing_bases
    )
end

"""
    TermSparsity{M}

Structure representing term sparsity information for polynomial optimization.

# Type Parameters
- `M`: Monomial type

# Fields
- `term_sparse_graph_supp::Vector{M}`: Support of the term sparsity graph
- `block_bases::Vector{Vector{M}}`: Bases of moment/localizing matrices in each block
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

"""
    init_activated_supp(partial_obj::P, cons::Vector{P}, mom_mtx_bases::Vector{M}) where {A,T,C,P,M}

Initialize the activated support for term sparsity iteration.

# Arguments
- `partial_obj::P`: Partial objective polynomial for this clique
- `cons::Vector{P}`: Constraint polynomials assigned to this clique
- `mom_mtx_bases::Vector{M}`: Moment matrix basis monomials

# Returns
- `Vector{M}`: Sorted union of symmetric-canonicalized objective monomials,
  constraint monomials, and diagonal moment entries
"""
function init_activated_supp(
    partial_obj::P, cons::Vector{P}, mom_mtx_bases::Vector{M}
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:Monomial{A,T}}
    return sorted_union(
        symmetric_canon.(monomials(partial_obj)),
        mapreduce(monomials, vcat, cons; init=M[]),
        [neat_dot(b, b) for b in mom_mtx_bases]
    )
end

"""
    term_sparsities(initial_activated_supp::Vector{M}, cons::Vector{P}, mom_mtx_bases::Vector{M}, localizing_mtx_bases::Vector{Vector{M}}, ts_algo::EliminationAlgorithm) where {A,T,C,P,M}

Computes term sparsity structures for the moment matrix and all localizing matrices.

# Arguments
- `initial_activated_supp::Vector{M}`: Initial set of activated support monomials
- `cons::Vector{P}`: Vector of constraint polynomials
- `mom_mtx_bases::Vector{M}`: Basis monomials for the moment matrix
- `localizing_mtx_bases::Vector{Vector{M}}`: Basis monomials for each localizing matrix corresponding to constraints
- `ts_algo::EliminationAlgorithm`: Algorithm for clique tree elimination in term sparsity graphs

# Returns
- `Vector{TermSparsity{M}}`: Vector containing term sparsity structures, with the first element
  corresponding to the moment matrix and subsequent elements corresponding to localizing matrices
"""
function term_sparsities(
    initial_activated_supp::Vector{M},
    cons::Vector{P},
    mom_mtx_bases::Vector{M},
    localizing_mtx_bases::Vector{Vector{M}},
    ts_algo::EliminationAlgorithm
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:Monomial{A,T}}
    [
        [iterate_term_sparse_supp(initial_activated_supp, one(P), mom_mtx_bases, ts_algo)];
        map(zip(cons, localizing_mtx_bases)) do (poly, basis)
            iterate_term_sparse_supp(initial_activated_supp, poly, basis, ts_algo)
        end
    ]
end

"""
    get_term_sparsity_graph(cons_support::Vector{M}, activated_supp::Vector{M}, bases::Vector{M}) where {A,T,M}

Constructs a term sparsity graph for polynomial constraints.

# Arguments
- `cons_support::Vector{M}`: Support monomials of constraints
- `activated_supp::Vector{M}`: Support from previous iterations
- `bases::Vector{M}`: Basis used to index the moment matrix

# Returns
- `SimpleGraph`: Term sparsity graph where edges connect basis elements
  whose products have support in the activated support set
"""
function get_term_sparsity_graph(
    cons_support::Vector{M}, activated_supp::Vector{M}, bases::Vector{M}
) where {A<:AlgebraType, T<:Integer, M<:Monomial{A,T}}
    nterms = length(bases)
    G = SimpleGraph(nterms)
    sorted_activated_supp = sort(activated_supp)

    for i in 1:nterms, j in i+1:nterms
        for supp in cons_support
            # _neat_dot3 returns Polynomial (with simplified terms)
            connected_poly_lr = _neat_dot3(bases[i], supp, bases[j])
            connected_poly_rl = _neat_dot3(bases[j], supp, bases[i])
            # Check if any monomial from the simplified result is in activated support
            monos_lr = monomials(connected_poly_lr)
            monos_rl = monomials(connected_poly_rl)
            if any(m in sorted_activated_supp for m in monos_lr) ||
               any(m in sorted_activated_supp for m in monos_rl)
                add_edge!(G, i, j)
                break  # Edge already added, no need to continue checking supports
            end
        end
    end
    return G
end

"""
    iterate_term_sparse_supp(activated_supp::Vector{M}, poly::P, basis::Vector{M}, elim_algo::EliminationAlgorithm) where {A,T,C,P,M}

Iteratively computes term sparsity support for a polynomial.

# Arguments
- `activated_supp::Vector{M}`: Currently activated support monomials
- `poly::P`: Input polynomial
- `basis::Vector{M}`: Basis monomials
- `elim_algo::EliminationAlgorithm`: Elimination algorithm for clique decomposition

# Returns
- `TermSparsity{M}`: Term sparsity structure containing graph support and block bases
"""
function iterate_term_sparse_supp(
    activated_supp::Vector{M}, poly::P, basis::Vector{M}, elim_algo::EliminationAlgorithm
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:Monomial{A,T}}
    F = get_term_sparsity_graph(monomials(poly), activated_supp, basis)
    blocks = clique_decomp(F, elim_algo)
    map(block -> add_clique!(F, block), blocks)
    return TermSparsity(term_sparsity_graph_supp(F, basis, poly), map(x -> basis[x], blocks))
end

"""
    term_sparsity_graph_supp(G::SimpleGraph, basis::Vector{M}, g::P) where {A,T,C,P,M}

Computes the support of a term sparsity graph for a given polynomial.

Implements equation (10.4) from "Sparse Polynomial Optimization: Theory and Practice".

# Arguments
- `G::SimpleGraph`: Term sparsity graph
- `basis::Vector{M}`: Basis monomials
- `g::P`: Input polynomial

# Returns
- `Vector{M}`: Support monomials for the term sparsity graph
"""
function term_sparsity_graph_supp(
    G::SimpleGraph, basis::Vector{M}, g::P
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:Monomial{A,T}}
    # Compute products basis[a]† * g_support * basis[b] for all graph edges
    # _neat_dot3 returns Polynomial, extract monomials from each result
    gsupp(a, b) = reduce(vcat, [monomials(_neat_dot3(a, g_supp, b)) for g_supp in monomials(g)])
    return union(
        [gsupp(basis[v], basis[v]) for v in vertices(G)]...,
        [gsupp(basis[e.src], basis[e.dst]) for e in edges(G)]...
    )
end

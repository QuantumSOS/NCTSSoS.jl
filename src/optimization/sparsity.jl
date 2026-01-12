"""
    CorrelativeSparsity{A<:AlgebraType, T<:Integer, P, M, ST}

Structure representing the correlative sparsity pattern of a polynomial optimization problem.

# Type Parameters
- `A<:AlgebraType`: The algebra type (PauliAlgebra, NonCommutativeAlgebra, etc.)
- `T<:Integer`: The index type used in the registry
- `P`: The polynomial type (Polynomial{A,T} or NCStatePolynomial)
- `M`: The monomial/basis type (NormalMonomial{A,T} or NCStateWord)
- `ST`: State type - `Nothing` for regular polynomials, `StateType` subtype for state polynomials

# Fields
- `cliques::Vector{Vector{T}}`: Groups of variable indices that form cliques in the sparsity graph
- `registry::VariableRegistry{A,T}`: Variable registry for symbol lookups and sub-registry creation
- `cons::Vector{P}`: All constraints in the problem
- `clq_cons::Vector{Vector{Int}}`: Constraint indices assigned to each clique, regardless of equality or inequality
- `global_cons::Vector{Int}`: Constraint indices not captured by any single clique
- `clq_mom_mtx_bases::Vector{Vector{M}}`: NormalMonomial bases for moment matrices within each clique
- `clq_localizing_mtx_bases::Vector{Vector{Vector{M}}}`: NormalMonomial bases for localizing matrices within each clique

# Notes
- Cliques store variable indices (type `T`), not symbols
- Use `subregistry(cs.registry, cs.cliques[i])` to get a clique-local registry
- Use `cs.registry[idx]` to look up symbol names for indices
"""
struct CorrelativeSparsity{A<:AlgebraType, T<:Integer, P, M, ST}
    cliques::Vector{Vector{T}}
    registry::VariableRegistry{A,T}
    cons::Vector{P}
    clq_cons::Vector{Vector{Int}}
    global_cons::Vector{Int}
    clq_mom_mtx_bases::Vector{Vector{M}}
    clq_localizing_mtx_bases::Vector{Vector{Vector{M}}}
end

# Show for regular polynomials (ST = Nothing)
function Base.show(io::IO, cs::CorrelativeSparsity{A,T,P,M,Nothing}) where {A,T,P,M}
    _show_correlative_sparsity(io, cs, "Correlative Sparsity ($(nameof(A)))", "monomials")
end

# Show for state polynomials (ST <: StateType)
function Base.show(io::IO, cs::CorrelativeSparsity{A,T,P,M,ST}) where {A,T,P,M,ST<:StateType}
    _show_correlative_sparsity(io, cs, "Correlative Sparsity ($(nameof(A)), $(nameof(ST)))", "NCStateWords")
end

function _show_correlative_sparsity(io::IO, cs::CorrelativeSparsity, header::String, basis_name::String)
    max_size = isempty(cs.cliques) ? 0 : maximum(length.(cs.cliques))
    println(io, "$header: \n")
    println(io, "   maximum clique size: $max_size")
    println(io, "   number of cliques: $(length(cs.cliques))")
    for clique_i in 1:length(cs.cliques)
        println(io, "   Clique $clique_i: ")
        clique_syms = [cs.registry[idx] for idx in cs.cliques[clique_i]]
        println(io, "       Variables: ", clique_syms)
        println(io, "       Bases length: ", length(cs.clq_mom_mtx_bases[clique_i]))
        println(io, "       Constraints: ")
        for (local_j, cons_idx) in enumerate(cs.clq_cons[clique_i])
            println(io, "           ", cs.cons[cons_idx], " :  with $(length(cs.clq_localizing_mtx_bases[clique_i][local_j])) basis $basis_name")
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

    # Helper to get graph positions from a monomial element
    function get_positions(m::NormalMonomial{A,T}) where {A,T}
        idxs = variable_indices(m)
        return sort([idx_to_node[idx] for idx in idxs])
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
  - `clq_mom_mtx_bases`: NormalMonomial bases for moment matrices
  - `clq_localizing_mtx_bases`: NormalMonomial bases for localizing matrices

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
    # get_ncbasis returns Vector{NormalMonomial{A,T}} directly
    cliques_moment_matrix_bases = map(cliques) do clique_indices
        sub_reg = subregistry(registry, clique_indices)
        basis = get_ncbasis(sub_reg, order)
        sorted_unique!(basis)
    end

    # Infer M type from first non-empty clique basis
    M = eltype(first(filter(!isempty, cliques_moment_matrix_bases)))

    # Compute degrees for localizing matrix basis truncation
    cliques_moment_matrix_bases_dg = [degree.(bs) for bs in cliques_moment_matrix_bases]

    # Generate localizing matrix bases (truncated based on constraint degree)
    cliques_localizing_bases = map(zip(eachindex(cliques), cliques_cons)) do (clique_idx, clique_cons_indices)
        if isempty(clique_cons_indices)
            return Vector{M}[]
        end

        # For each constraint, compute the reduced order basis
        cur_orders = order .- cld.(degree.(all_cons[clique_cons_indices]), 2)
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

    return CorrelativeSparsity{A,T,P,M,Nothing}(
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
- `M`: NormalMonomial type

# Fields
- `term_sparse_graph_supp::Vector{M}`: Support of the term sparsity graph
- `block_bases::Vector{Vector{M}}`: Bases of moment/localizing matrices in each block
"""
struct TermSparsity{M}
    term_sparse_graph_supp::Vector{M}
    block_bases::Vector{Vector{M}}
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
- `mom_mtx_bases::Vector{M}`: Moment matrix basis NormalMonomials

# Returns
- `Vector{NormalMonomial{A,T}}`: Sorted union of symmetric-canonicalized objective monomials,
  constraint monomials, and diagonal moment entries
"""
function init_activated_supp(
    partial_obj::P, cons::Vector{P}, mom_mtx_bases::Vector{M}
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:NormalMonomial{A,T}}
    # Compute diagonal entries b†b with bilinear expansion for NormalMonomial bases
    # For each basis NormalMonomial b, expand: Σ_{k,l} conj(c_k) * c_l * (word_k)† * word_l
    NM = NormalMonomial{A,T}
    diagonal_entries = NM[]
    for b in mom_mtx_bases
        for (c1, word1) in b
            for (c2, word2) in b
                diag_word = neat_dot(word1.word, word2.word)
                simplified_poly = Polynomial(simplify(A, diag_word))
                append!(diagonal_entries, last.(simplified_poly.terms))
            end
        end
    end
    return sorted_union(
        symmetric_canon.(last.(partial_obj.terms)),
        mapreduce(p -> last.(p.terms), vcat, cons; init=NM[]),
        diagonal_entries
    )
end

"""
    term_sparsities(initial_activated_supp::Vector{NormalMonomial{A,T}}, cons::Vector{P}, mom_mtx_bases::Vector{M}, localizing_mtx_bases::Vector{Vector{M}}, ts_algo::EliminationAlgorithm)

Computes term sparsity structures for the moment matrix and all localizing matrices.

# Arguments
- `initial_activated_supp::Vector{NormalMonomial}`: Initial set of activated support NormalMonomials
- `cons::Vector{P}`: Vector of constraint polynomials
- `mom_mtx_bases::Vector{M}`: Basis NormalMonomials for the moment matrix
- `localizing_mtx_bases::Vector{Vector{M}}`: Basis Monomials for each localizing matrix corresponding to constraints
- `ts_algo::EliminationAlgorithm`: Algorithm for clique tree elimination in term sparsity graphs

# Returns
- `Vector{TermSparsity{M}}`: Vector containing term sparsity structures, with the first element
  corresponding to the moment matrix and subsequent elements corresponding to localizing matrices
"""
function term_sparsities(
    initial_activated_supp::Vector{NormalMonomial{A,T}},
    cons::Vector{P},
    mom_mtx_bases::Vector{M},
    localizing_mtx_bases::Vector{Vector{M}},
    ts_algo::EliminationAlgorithm
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:NormalMonomial{A,T}}
    [
        [iterate_term_sparse_supp(initial_activated_supp, one(P), mom_mtx_bases, ts_algo)];
        map(zip(cons, localizing_mtx_bases)) do (poly, basis)
            iterate_term_sparse_supp(initial_activated_supp, poly, basis, ts_algo)
        end
    ]
end

"""
    get_term_sparsity_graph(cons_support::Vector{NormalMonomial{A,T}}, activated_supp::Vector{NormalMonomial{A,T}}, bases::Vector{M}) where {A,T,M}

Constructs a term sparsity graph for polynomial constraints.

# Arguments
- `cons_support::Vector{NormalMonomial}`: Support NormalMonomials of constraints
- `activated_supp::Vector{NormalMonomial}`: Support from previous iterations
- `bases::Vector{M}`: Basis Monomials used to index the moment matrix

# Returns
- `SimpleGraph`: Term sparsity graph where edges connect basis elements
  whose products have support in the activated support set
"""
function get_term_sparsity_graph(
    cons_support::Vector{NormalMonomial{A,T}}, activated_supp::Vector{NormalMonomial{A,T}}, bases::Vector{M}
) where {A<:AlgebraType, T<:Integer, M<:NormalMonomial{A,T}}
    nterms = length(bases)
    G = SimpleGraph(nterms)
    sorted_activated_supp = sort(activated_supp)

    for i in 1:nterms, j in i+1:nterms
        edge_found = false
        # Expand bilinear form over all term pairs from NormalMonomial bases
        for (_, word_i) in bases[i]
            edge_found && break
            for (_, word_j) in bases[j]
                edge_found && break
                for supp in cons_support
                    # _neat_dot3 on NormalMonomials returns Vector{T}
                    word_lr = _neat_dot3(word_i, supp, word_j)
                    word_rl = _neat_dot3(word_j, supp, word_i)
                    connected_poly_lr = Polynomial(simplify(A, word_lr))
                    connected_poly_rl = Polynomial(simplify(A, word_rl))
                    # Check if any result monomial is in activated support
                    monos_lr = last.(connected_poly_lr.terms)
                    monos_rl = last.(connected_poly_rl.terms)
                    if any(m in sorted_activated_supp for m in monos_lr) ||
                       any(m in sorted_activated_supp for m in monos_rl)
                        add_edge!(G, i, j)
                        edge_found = true
                        break
                    end
                end
            end
        end
    end
    return G
end

"""
    iterate_term_sparse_supp(activated_supp::Vector{NormalMonomial{A,T}}, poly::P, basis::Vector{M}, elim_algo::EliminationAlgorithm) where {A,T,C,P,M}

Iteratively computes term sparsity support for a polynomial.

# Arguments
- `activated_supp::Vector{NormalMonomial}`: Currently activated support NormalMonomials
- `poly::P`: Input polynomial
- `basis::Vector{M}`: Basis Monomials
- `elim_algo::EliminationAlgorithm`: Elimination algorithm for clique decomposition

# Returns
- `TermSparsity{M}`: Term sparsity structure containing graph support and block bases
"""
function iterate_term_sparse_supp(
    activated_supp::Vector{NormalMonomial{A,T}}, poly::P, basis::Vector{M}, elim_algo::EliminationAlgorithm
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:NormalMonomial{A,T}}
    # poly.terms contains (coef, NormalMonomial) tuples
    F = get_term_sparsity_graph(last.(poly.terms), activated_supp, basis)
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
- `basis::Vector{M}`: Basis Monomials
- `g::P`: Input polynomial

# Returns
- `Vector{NormalMonomial{A,T}}`: Support NormalMonomials for the term sparsity graph
"""
function term_sparsity_graph_supp(
    G::SimpleGraph, basis::Vector{M}, g::P
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:NormalMonomial{A,T}}
    NM = NormalMonomial{A,T}
    # Compute products with bilinear expansion over Monomial terms
    function gsupp(a::M, b::M)
        out = NM[]
        for (_, word_a) in a
            for (_, word_b) in b
                for g_supp in last.(g.terms)
                    word = _neat_dot3(word_a, g_supp, word_b)
                    p = Polynomial(simplify(A, word))
                    append!(out, last.(p.terms))
                end
            end
        end
        return out
    end
    return union(
        [gsupp(basis[v], basis[v]) for v in vertices(G)]...,
        [gsupp(basis[e.src], basis[e.dst]) for e in edges(G)]...
    )
end

# =============================================================================
# State Correlative Sparsity Implementation
# =============================================================================

"""
    get_state_correlative_graph(registry, obj, cons) -> (SimpleGraph, Vector{T}, Dict{T,Int})

Constructs a correlative sparsity graph for state polynomial optimization.

Similar to `get_correlative_graph` but uses `variable_indices` on NCStatePolynomials.

# Arguments
- `registry`: Variable registry
- `obj`: Objective NCStatePolynomial
- `cons`: Vector of constraint NCStatePolynomials

# Returns
- Tuple containing:
  - `SimpleGraph`: Graph representing variable correlations
  - `sorted_indices::Vector{T}`: All unique indices in sorted order
  - `idx_to_node::Dict{T,Int}`: Mapping from variable index to graph node position
"""
function get_state_correlative_graph(
    registry::VariableRegistry{A,T},
    obj::NCStatePolynomial{C,ST,A,T},
    cons::Vector{NCStatePolynomial{C,ST,A,T}}
) where {A<:AlgebraType, ST<:StateType, T<:Integer, C<:Number}
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

    # Helper to get graph positions from an NCStateWord
    function get_positions(ncsw::NCStateWord{ST,A,T})
        positions = Int[]
        for idx in variables(ncsw)
            if haskey(idx_to_node, idx)
                push!(positions, idx_to_node[idx])
            end
        end
        return positions
    end

    # Add cliques from objective NC state words
    for ncsw in monomials(obj)
        positions = get_positions(ncsw)
        add_clique!(G, positions)
    end

    # Add cliques from constraints (whole constraint, not per-monomial)
    # All variables in a constraint must be connected to ensure the constraint
    # can be assigned to a single clique in assign_state_constraint()
    for con in cons
        con_indices = variable_indices(con)
        positions = [idx_to_node[idx] for idx in con_indices if haskey(idx_to_node, idx)]
        add_clique!(G, positions)
    end

    return G, sorted_indices, idx_to_node
end

"""
    assign_state_constraint(cliques, cons, registry) -> (Vector{Vector{Int}}, Vector{Int})

Assigns state polynomial constraints to cliques based on variable support.

# Arguments
- `cliques::Vector{Vector{T}}`: Variable index cliques
- `cons::Vector{NCStatePolynomial}`: Constraint state polynomials
- `registry`: Variable registry

# Returns
- Tuple containing:
  - Constraint indices for each clique
  - Global constraint indices not captured by any single clique
"""
function assign_state_constraint(
    cliques::Vector{Vector{T}},
    cons::Vector{NCStatePolynomial{C,ST,A,T}},
    registry::VariableRegistry{A,T}
) where {A<:AlgebraType, ST<:StateType, T<:Integer, C<:Number}
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
    _extract_compound_state_words(obj, cons) -> Vector{StateWord}

Extract compound StateWords (those with multiple expectations) from objective and constraints.

For state polynomial optimization, objectives may contain products of expectations like
`<A><B>`. To properly represent these in the moment matrix, the basis needs elements
of the form `<A>*I` and `<B>*I`. This function identifies all StateWords from the
objective and constraints that have multiple expectations (compound StateWords).

# Returns
A vector of StateWords that contain multiple expectations. These need to be passed
to `get_state_basis` to generate the extended basis elements.
"""
function _extract_compound_state_words(
    obj::NCStatePolynomial{C,ST,A,T},
    cons::Vector{NCStatePolynomial{C,ST,A,T}}
) where {C<:Number, ST<:StateType, A<:AlgebraType, T<:Integer}
    result = StateWord{ST,A,T}[]

    # Extract from objective
    for ncsw in monomials(obj)
        sw = ncsw.sw
        # A compound StateWord has more than one state symbol, OR has a non-identity
        # state symbol that differs from the identity (meaning it has actual expectations)
        if length(sw.state_syms) > 1 || (!isone(sw) && !isempty(sw.state_syms))
            push!(result, sw)
        end
    end

    # Extract from constraints
    for con in cons
        for ncsw in monomials(con)
            sw = ncsw.sw
            if length(sw.state_syms) > 1 || (!isone(sw) && !isempty(sw.state_syms))
                push!(result, sw)
            end
        end
    end

    return unique(result)
end

"""
    correlative_sparsity(pop::PolyOpt, order, elim_algo) -> CorrelativeSparsity

Decomposes a state polynomial optimization problem into a correlative sparsity pattern.

# Arguments
- `pop::PolyOpt{A,T,P}`: Polynomial optimization problem with NCStatePolynomial objective
- `order::Int`: Order of the moment relaxation
- `elim_algo::EliminationAlgorithm`: Algorithm for clique tree elimination

# Returns
- `CorrelativeSparsity{A,T,P,M,ST}`: Correlative sparsity structure for state polynomials
"""
function correlative_sparsity(
    pop::PolyOpt{A,T,P},
    order::Int,
    elim_algo::EliminationAlgorithm
) where {A<:AlgebraType,T<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,A,T}}
    registry = pop.registry
    all_cons = vcat(pop.eq_constraints, pop.ineq_constraints)

    # Build correlative graph and get index mappings
    G, sorted_indices, idx_to_node = get_state_correlative_graph(registry, pop.objective, all_cons)

    # Decompose graph into cliques
    clique_node_sets = clique_decomp(G, elim_algo)

    # Convert clique node positions back to variable indices
    node_to_idx = Dict{Int,T}(pos => idx for (idx, pos) in idx_to_node)
    cliques = map(clique_node_sets) do node_set
        sort([node_to_idx[node] for node in node_set])
    end

    # Assign constraints to cliques
    cliques_cons, global_cons = assign_state_constraint(cliques, all_cons, registry)

    # Generate moment matrix bases for each clique using subregistry + get_state_basis
    # Note: get_state_basis now automatically generates all (StateWord, Monomial) combinations
    # including compound forms like <M1><M2>*I needed for products of expectations
    M = NCStateWord{ST,A,T}
    cliques_moment_matrix_bases = Vector{Vector{M}}()

    for clique_indices in cliques
        # Create sub-registry for this clique
        sub_reg = subregistry(registry, clique_indices)
        # Generate NCStateWord basis - includes all compound forms automatically
        basis = get_state_basis(sub_reg, order; state_type=ST)
        # Sort and deduplicate
        push!(cliques_moment_matrix_bases, sorted_unique(basis))
    end

    # Compute degrees for localizing matrix basis truncation
    cliques_moment_matrix_bases_dg = [degree.(bs) for bs in cliques_moment_matrix_bases]

    # Generate localizing matrix bases (truncated based on constraint degree)
    cliques_localizing_bases = map(zip(eachindex(cliques), cliques_cons)) do (clique_idx, clique_cons_indices)
        if isempty(clique_cons_indices)
            return Vector{M}[]
        end

        # For each constraint, compute the reduced order basis
        cur_orders = order .- cld.(degree.(all_cons[clique_cons_indices]), 2)
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

    return CorrelativeSparsity{A,T,P,M,ST}(
        cliques,
        registry,
        all_cons,
        cliques_cons,
        global_cons,
        cliques_moment_matrix_bases,
        cliques_localizing_bases
    )
end

# =============================================================================
# State Term Sparsity
# =============================================================================

"""
    init_activated_supp(partial_obj::P, cons::Vector{P}, mom_mtx_bases::Vector{M})

Initialize the activated support for state polynomial term sparsity iteration.

This follows the NCTSSOS algorithm where the initial support includes:
1. Canonicalized objective monomials (symmetric_canon applied)
2. Constraint monomials

Note: Diagonal entries (monosquare terms) are NOT included by default, matching
NCTSSOS's `monosquare=false` default. This is crucial for proper term sparsity:
including all diagonal entries would create spurious edges in the term sparsity
graph because products of compound basis elements would match diagonals of other
basis elements, making the graph fully connected.

# Arguments
- `partial_obj::NCStatePolynomial`: Partial objective for this clique
- `cons::Vector{NCStatePolynomial}`: Constraint state polynomials
- `mom_mtx_bases::Vector{NCStateWord}`: Moment matrix basis NCStateWords (unused, kept for API)

# Returns
- `Vector{NCStateWord}`: Sorted union of canonicalized objective monomials and
  constraint monomials
"""
function init_activated_supp(
    partial_obj::P, cons::Vector{P}, _mom_mtx_bases::Vector{M}
) where {ST<:StateType, A<:AlgebraType, T<:Integer, C<:Number, P<:NCStatePolynomial{C,ST,A,T}, M<:NCStateWord{ST,A,T}}
    # Only include objective and constraint monomials (no monosquare/diagonal terms)
    # This matches NCTSSOS's default monosquare=false behavior
    return sorted_union(
        symmetric_canon.(monomials(partial_obj)),
        symmetric_canon.(mapreduce(monomials, vcat, cons; init=M[]))
    )
end

"""
    term_sparsities(initial_activated_supp, cons, mom_mtx_bases, localizing_mtx_bases, ts_algo)

Computes term sparsity structures for state polynomial optimization.

# Arguments
- `initial_activated_supp::Vector{NCStateWord}`: Initial activated support
- `cons::Vector{NCStatePolynomial}`: Constraint state polynomials
- `mom_mtx_bases::Vector{NCStateWord}`: Moment matrix basis NCStateWords
- `localizing_mtx_bases::Vector{Vector{NCStateWord}}`: Localizing matrix bases
- `ts_algo::EliminationAlgorithm`: Algorithm for term sparsity graphs

# Returns
- `Vector{TermSparsity{NCStateWord}}`: Term sparsity structures
"""
function term_sparsities(
    initial_activated_supp::Vector{M},
    cons::Vector{P},
    mom_mtx_bases::Vector{M},
    localizing_mtx_bases::Vector{Vector{M}},
    ts_algo::EliminationAlgorithm
) where {ST<:StateType, A<:AlgebraType, T<:Integer, C<:Number, P<:NCStatePolynomial{C,ST,A,T}, M<:NCStateWord{ST,A,T}}
    [
        [iterate_term_sparse_supp(initial_activated_supp, one(P), mom_mtx_bases, ts_algo)];
        map(zip(cons, localizing_mtx_bases)) do (poly, basis)
            iterate_term_sparse_supp(initial_activated_supp, poly, basis, ts_algo)
        end
    ]
end

"""
    get_term_sparsity_graph(cons_support, activated_supp, bases) for NCStateWord

Constructs a term sparsity graph for state polynomial constraints.

An edge is added between basis elements i and j if the canonicalized product
bi† * supp * bj (for any supp in cons_support) is present in the activated support.
The comparison uses `expval` to convert NCStateWords to StateWords, ensuring that
structurally different NCStateWords with the same expectation value are treated as equal.
"""
function get_term_sparsity_graph(
    cons_support::Vector{M}, activated_supp::Vector{M}, bases::Vector{M}
) where {ST<:StateType, A<:AlgebraType, T<:Integer, M<:NCStateWord{ST,A,T}}
    nterms = length(bases)
    G = SimpleGraph(nterms)
    # Convert activated support to StateWords via expval for proper comparison
    sorted_activated_supp = sort([symmetric_canon(expval(ncsw)) for ncsw in activated_supp])

    for i in 1:nterms, j in i+1:nterms
        for supp in cons_support
            # _neat_dot3 returns NCStateWord, simplify to get NCStatePolynomial
            connected_lr = simplify(_neat_dot3(bases[i], supp, bases[j]))
            connected_rl = simplify(_neat_dot3(bases[j], supp, bases[i]))
            # Check if any canonicalized StateWord is in activated support
            found = false
            for ncsw in monomials(connected_lr)
                canon_sw = symmetric_canon(expval(ncsw))
                if canon_sw in sorted_activated_supp
                    found = true
                    break
                end
            end
            if !found
                for ncsw in monomials(connected_rl)
                    canon_sw = symmetric_canon(expval(ncsw))
                    if canon_sw in sorted_activated_supp
                        found = true
                        break
                    end
                end
            end
            if found
                add_edge!(G, i, j)
                break
            end
        end
    end
    return G
end

"""
    iterate_term_sparse_supp(activated_supp, poly, basis, elim_algo) for NCStatePolynomial

Iteratively computes term sparsity support for a state polynomial.
"""
function iterate_term_sparse_supp(
    activated_supp::Vector{M}, poly::P, basis::Vector{M}, elim_algo::EliminationAlgorithm
) where {ST<:StateType, A<:AlgebraType, T<:Integer, C<:Number, P<:NCStatePolynomial{C,ST,A,T}, M<:NCStateWord{ST,A,T}}
    F = get_term_sparsity_graph(monomials(poly), activated_supp, basis)
    blocks = clique_decomp(F, elim_algo)
    map(block -> add_clique!(F, block), blocks)
    return TermSparsity(term_sparsity_graph_supp(F, basis, poly), map(x -> basis[x], blocks))
end

"""
    term_sparsity_graph_supp(G, basis, g) for NCStatePolynomial

Computes the support of a term sparsity graph for a state polynomial.
"""
function term_sparsity_graph_supp(
    G::SimpleGraph, basis::Vector{M}, g::P
) where {ST<:StateType, A<:AlgebraType, T<:Integer, C<:Number, P<:NCStatePolynomial{C,ST,A,T}, M<:NCStateWord{ST,A,T}}
    # _neat_dot3 returns NCStateWord, simplify to get NCStatePolynomial
    gsupp(a, b) = reduce(vcat, [monomials(simplify(_neat_dot3(a, g_supp, b))) for g_supp in monomials(g)])
    return union(
        [gsupp(basis[v], basis[v]) for v in vertices(G)]...,
        [gsupp(basis[e.src], basis[e.dst]) for e in edges(G)]...
    )
end

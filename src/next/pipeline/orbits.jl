# Signed actions through contract 1. No word representation is inspected here.
# In particular, the orbit table also accepts degree-changing
# actions; only the cyclic action used to split a block must preserve its rows.

struct CyclicRows{W,G}
    generator::G
    order::Int
    words::Vector{Vector{W}}
    phases::Vector{Vector{Phase}}
    stabilizers::Vector{Phase}
end

"Partition a closed row list into signed cyclic orbits. The order is explicit."
function CyclicRows(rows::AbstractVector{W}, g::AbstractAction, order::Integer) where {W<:AbstractWord}
    order > 0 || throw(ArgumentError("The cyclic order must be positive."))
    is_antiunitary(g) && throw(ArgumentError("The momentum generator must be unitary."))
    domain = Set(rows)
    length(domain) == length(rows) || throw(ArgumentError("The row list contains duplicates."))
    seen = Set{W}()
    orbits, phases, stabilizers = Vector{W}[], Vector{Phase}[], Phase[]
    for r in sort!(collect(domain))
        r in seen && continue
        ws, ps = W[], Phase[]
        w, p = r, one(Phase)
        while true
            w in domain || throw(ArgumentError("The cyclic generator does not preserve the row list; supply its orbit closure."))
            w in seen && throw(ArgumentError("The supplied action is not a permutation of the rows."))
            push!(seen, w); push!(ws, w); push!(ps, p)
            q, w = act(g, w); p = p * q
            w == r && break
            length(ws) < order || throw(ArgumentError("The declared cyclic order $order is too small."))
        end
        l = length(ws)
        order % l == 0 && mod(Int(p.k) * div(order, l), 4) == 0 ||
            throw(ArgumentError("The declared cyclic order $order does not close the signed row orbit."))
        push!(orbits, ws); push!(phases, ps); push!(stabilizers, p)
    end
    return CyclicRows(g, Int(order), orbits, phases, stabilizers)
end

struct OrbitFormSink{T}
    table::T
    terms::Dict{Int32,ComplexF64}
end
function push_term!(s::OrbitFormSink, c::Number, w::AbstractWord)
    m = orbit_coordinates!(s.table, w)
    m.re != 0 && (s.terms[m.re] = get(s.terms, m.re, 0.0im) + c*m.real_coefficient)
    m.im != 0 && (s.terms[m.im] = get(s.terms, m.im, 0.0im) + c*m.imag_coefficient)
    return s
end

function _close_orbit_form!(forms, terms; diagonal=false)
    for id in sort!(collect(keys(terms)))
        c = terms[id]
        if diagonal
            abs(imag(c)) <= 1e-10 * max(1, abs(c)) || error("Nonreal diagonal in orbit reduction.")
            c = complex(real(c))
        end
        abs(c) > 1e-12 && push_term!(forms, c, id)
    end
    close_form!(forms)
end

function _orbit_objective!(program, table, objective)
    sink = OrbitFormSink(table, Dict{Int32,ComplexF64}())
    for (c, w) in objective
        push_term!(sink, c, w)
    end
    for id in sort!(collect(keys(sink.terms)))
        c = sink.terms[id]
        abs(imag(c)) <= 1e-12 || throw(ArgumentError("The objective must be real."))
        iszero(real(c)) && continue
        push!(program.objective_ids, id); push!(program.objective_coefs, real(c))
    end
    program.nvars = length(table.labels)
    return validate(program), table
end

function _check_objective(objective, generators)
    isempty(objective) && return
    W = typeof(first(objective)[2])
    terms = Dict{W,ComplexF64}()
    for (c, w) in objective
        terms[w] = get(terms, w, 0im) + c
    end
    filter!(p -> !iszero(last(p)), terms)
    for g in generators
        image = Dict{W,ComplexF64}()
        for (w, c) in terms
            p, u = act(g, w)
            image[u] = get(image, u, 0im) + scale(p, is_antiunitary(g) ? conj(c) : c)
        end
        filter!(p -> !iszero(last(p)), image)
        image == terms || throw(ArgumentError("A supplied symmetry does not preserve the objective."))
    end
end

"""
    cyclic_program(rows, objective, generator, order; moment_generators=(), labels=())

Build momentum blocks directly from signed row orbits, including shorter
orbits. `labels` are commuting diagonal actions on these rows. Additional
`moment_generators` impose orbit relations, including degree-changing and
antiunitary maps; they do not claim a further block split. Returns program,
orbit table, and the (momentum, row-label) of each block.

This constructs no site-space moment matrix. A single sum over the second row orbit suffices:
translation invariance removes the sum over the first orbit.
"""
function cyclic_program(rows::AbstractVector{W}, objective, g::AbstractAction, order::Integer;
    moment_generators=(), labels=()) where {W<:AbstractWord}
    generators = (g, labels..., moment_generators...)
    _check_objective(objective, generators)
    orbits = CyclicRows(rows, g, order)
    table = OrbitMoments(W, generators)
    program = MomentProgram{ComplexF64}()
    sink = OrbitFormSink(table, Dict{Int32,ComplexF64}())
    rowlabels = map(orbits.words) do ws
        map(labels) do h
            is_antiunitary(h) && throw(ArgumentError("A row label must be unitary."))
            p, image = act(h, first(ws))
            image == first(ws) || throw(ArgumentError("A label action must fix each row up to phase."))
            for w in ws
                q, u = act(h, w)
                u == w && q == p || throw(ArgumentError("Row labels must be constant on each translation orbit."))
            end
            p.k
        end
    end
    sectors = Tuple{Int,Any}[]
    for k in 0:(order - 1), label in unique(rowlabels)
        active = [a for a in eachindex(orbits.words) if rowlabels[a] == label &&
            mod(4 * k * length(orbits.words[a]) - Int(orbits.stabilizers[a].k) * order, 4 * order) == 0]
        isempty(active) && continue
        for jj in eachindex(active), ii in jj:length(active)
            a, b = active[ii], active[jj]
            wa, wb = orbits.words[a], orbits.words[b]
            p, sa = star(first(wa))
            empty!(sink.terms)
            normalization = sqrt(length(wa) / length(wb))
            for t in eachindex(wb)
                c = normalization * cispi(-2 * k * (t - 1) / order)
                c = scale(p * orbits.phases[b][t], c)
                mul!(sink, sa, wb[t], c)
            end
            _close_orbit_form!(program.forms, sink.terms; diagonal=ii == jj)
        end
        push_block!(program, HPSD, length(active), length(sectors) + 1)
        push!(sectors, (k, label))
    end
    p, table = _orbit_objective!(program, table, objective)
    return _anti_real_program(p, generators), table, sectors
end

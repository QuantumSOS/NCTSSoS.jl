# Reference builders for the tests of the redesign. They use contract 1 only, so
# they work for each word type. They build the dense, unreduced relaxation: the
# reference that each reduced combination must match (acceptance levels A and B).
using NCTSSoS

if !@isdefined(Nx)
    const Nx = NCTSSoS.Next
end

"Moment ids for the reference builder: word key -> (id of real part, id of imaginary part); 0 means no variable."
struct ReferenceTable{W}
    ids::Dict{W,Tuple{Int32,Int32}}
    labels::Vector{Nx.MomentLabel{W}}
end
function ReferenceTable{W}() where {W}
    e = one(W)
    return ReferenceTable{W}(Dict(e => (Int32(1), Int32(0))), [Nx.MomentLabel{W}(e, 0x00)])
end
function moment_ids!(t::ReferenceTable{W}, key::W) where {W}
    return get!(t.ids, key) do
        kind = Nx.moment_reality(key)
        re = im = Int32(0)
        if kind != :imaginary
            push!(t.labels, Nx.MomentLabel{W}(key, 0x01)); re = Int32(length(t.labels))
        end
        if kind != :real
            push!(t.labels, Nx.MomentLabel{W}(key, 0x02)); im = Int32(length(t.labels))
        end
        (re, im)
    end
end

"A sink that writes the moment of each pushed term into the open form."
struct FormSink{W,C}
    table::ReferenceTable{W}
    forms::Nx.LinearForms{C}
end
function Nx.push_term!(s::FormSink{W,C}, c::Number, w::W) where {W,C}
    key, phase, conjugated = Nx.moment_key(w)
    re, im_id = moment_ids!(s.table, key)
    cc = Nx.scale(phase, convert(C, c))
    re != 0 && Nx.push_term!(s.forms, cc, re)
    im_id != 0 && Nx.push_term!(s.forms, cc * (conjugated ? -im : im), im_id)
    return s
end

"""
    dense_program(rows, objective; C = ComplexF64) -> (program, table)

The unreduced moment relaxation with one Hermitian block: entry `(i, j)` is the
moment of `star(rows[i]) * rows[j]`. `objective` is a vector of
`(coefficient, word)`; the program minimises the real part of its moment.
"""
function dense_program(rows::Vector{W}, objective; C::Type=ComplexF64) where {W<:Nx.AbstractWord}
    table = ReferenceTable{W}()
    program = Nx.MomentProgram{C}()
    sink = FormSink{W,C}(table, program.forms)
    n = length(rows)
    for j in 1:n, i in j:n
        p, s = Nx.star(rows[i])
        Nx.mul!(sink, s, rows[j], Nx.scale(p, one(C)))
        Nx.close_form!(program.forms)
    end
    Nx.push_block!(program, Nx.HPSD, n)
    weights = Dict{Int32,Float64}()
    for (c, w) in objective
        key, phase, conjugated = Nx.moment_key(w)
        re, im_id = moment_ids!(table, key)
        cc = Nx.scale(phase, complex(float(c)))
        re != 0 && (weights[re] = get(weights, re, 0.0) + real(cc))
        im_id != 0 && (weights[im_id] = get(weights, im_id, 0.0) + real(cc * (conjugated ? -im : im)))
    end
    for id in sort!(collect(keys(weights)))
        push!(program.objective_ids, id); push!(program.objective_coefs, weights[id])
    end
    program.nvars = length(table.labels)
    return Nx.validate(program), table
end

"""
    legacy_product_terms(A, T, a, b) -> Dict{Vector{Int},ComplexF64}

The product of two words in the current package, as legacy letters to
coefficient. `A` is the algebra type and `T` the index type of the registry.
"""
function legacy_product_terms(::Type{A}, ::Type{T}, a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer}) where {A,T}
    M = NormalMonomial{A,T}
    product = M(T.(a)) * M(T.(b))
    return Dict{Vector{Int},ComplexF64}(Int.(m.word) => ComplexF64(c) for (c, m) in terms(product))
end

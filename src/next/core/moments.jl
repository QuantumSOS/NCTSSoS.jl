"Open-addressed table of representative words. No aliases are inserted or canonicalized here."
mutable struct MomentTable{W<:AbstractWord}
    keys::Vector{W}
    ids::Vector{Tuple{Int32,Int32}}
    labels::Vector{MomentLabel{W}}
    count::Int
end
function MomentTable{W}(capacity::Integer=16) where {W<:AbstractWord}
    capacity >= 0 || throw(ArgumentError("Moment capacity must be nonnegative."))
    n = nextpow(2,max(8,2capacity))
    table = MomentTable{W}(Vector{W}(undef,n),fill((Int32(0),Int32(0)),n),MomentLabel{W}[],0)
    sizehint!(table.labels,max(1,2capacity))
    slot = _moment_slot(table,one(W))
    table.keys[slot] = one(W)
    table.ids[slot] = (Int32(1),Int32(0))
    push!(table.labels,MomentLabel(one(W),0x00))
    table.count = 1
    return table
end
Base.length(t::MomentTable) = t.count
nvars(t::MomentTable) = length(t.labels)
@inline function _moment_slot(t::MomentTable{W},key::W) where {W}
    mask = length(t.keys)-1
    slot = Int(hash(key,UInt(0)) & UInt(mask))+1
    @inbounds while t.ids[slot] != (Int32(0),Int32(0)) && t.keys[slot] != key
        slot = (slot & mask)+1
    end
    return slot
end
function reserve!(t::MomentTable{W},capacity::Integer) where {W}
    capacity >= 0 || throw(ArgumentError("Moment capacity must be nonnegative."))
    n = nextpow(2,max(8,2capacity))
    if n > length(t.keys)
        oldkeys,oldids = t.keys,t.ids
        t.keys = Vector{W}(undef,n)
        t.ids = fill((Int32(0),Int32(0)),n)
        for i in eachindex(oldids)
            if oldids[i] != (Int32(0),Int32(0))
                slot = _moment_slot(t,oldkeys[i])
                t.keys[slot],t.ids[slot] = oldkeys[i],oldids[i]
            end
        end
    end
    sizehint!(t.labels,max(length(t.labels),2capacity))
    return t
end
"Lookup only; (0,0) means absent. The caller supplies the orbit representative."
moment_ids(t::MomentTable{W},key::W) where {W} = t.ids[_moment_slot(t,key)]
"Insert a representative with its intrinsic or symmetry-imposed reality; repeated declarations must agree."
function moment_ids!(t::MomentTable{W},key::W,reality::Symbol=moment_reality(key)) where {W}
    reality in (:real,:imaginary,:complex) || throw(ArgumentError("Moment reality must be :real, :imaginary, or :complex."))
    slot = _moment_slot(t,key)
    existing = t.ids[slot]
    if existing != (Int32(0),Int32(0))
        stored = existing[1]==0 ? :imaginary : existing[2]==0 ? :real : :complex
        stored==reality || throw(ArgumentError("The representative already has reality $stored, not $reality."))
        return existing
    end
    if 2(t.count+1) > length(t.keys)
        reserve!(t,2(t.count+1))
        slot = _moment_slot(t,key)
    end
    length(t.labels) <= typemax(Int32)-2 || throw(ArgumentError("Moment ids exceed the Int32 limit."))
    re = im = Int32(0)
    if reality != :imaginary
        push!(t.labels,MomentLabel(key,0x01)); re = Int32(length(t.labels))
    end
    if reality != :real
        push!(t.labels,MomentLabel(key,0x02)); im = Int32(length(t.labels))
    end
    t.keys[slot],t.ids[slot] = key,(re,im)
    t.count += 1
    return re,im
end

struct MomentSink{W,C}
    table::MomentTable{W}
    forms::LinearForms{C}
end
function push_term!(sink::MomentSink{W,C},c::Number,w::W) where {W,C}
    iszero(c) && return sink
    key,p,conjugated = moment_key(w)
    re,im_id = moment_ids!(sink.table,key)
    C <: Real && im_id != 0 && throw(ArgumentError("These rows have complex moments; use C=ComplexF64 for the plain program."))
    cc = scale(p,convert(C,c))
    re != 0 && push_term!(sink.forms,cc,re)
    im_id != 0 && push_term!(sink.forms,cc*(conjugated ? -im : im),im_id)
    return sink
end

"All normal-form rows of degree at most order on one specified group of sites."
function basis_rows(::Type{W},algebra,group,order::Integer) where {W<:AbstractWord}
    order >= 0 || throw(ArgumentError("Row order must be nonnegative."))
    sites = sort!(collect(Int,group))
    allunique(sites) && all(>(0),sites) || throw(ArgumentError("Row sites must be positive and unique."))
    rows = W[one(W)]
    _extend_rows!(rows,algebra,sites,one(W),1,Int(order))
    sort!(rows)
    return rows
end
function _extend_rows!(rows::Vector{W},::Pauli,group,w::W,start,remaining) where {W<:PauliWord}
    remaining == 0 && return rows
    for j in start:length(group),axis in 0:2
        next = mul(w,pauli_letter(W,group[j],axis))[2]
        push!(rows,next)
        _extend_rows!(rows,Pauli(),group,next,j+1,remaining-1)
    end
    return rows
end
function _extend_rows!(rows::Vector{W},::Fermionic,group,w::W,start,remaining) where {W<:FermiWord}
    remaining == 0 && return rows
    n = length(group)
    for j in start:2n
        l = j<=n ? -group[n-j+1] : group[j-n]
        v = fermi_letter(W,l)
        next = FermiWord(w.c .| v.c,w.a .| v.a)
        push!(rows,next)
        _extend_rows!(rows,Fermionic(),group,next,j+1,remaining-1)
    end
    return rows
end
function _extend_rows!(rows::Vector{W},a::SiteAlgebra{A},group,w::W,start,remaining) where {A,W<:SiteWord{A}}
    remaining == 0 && return rows
    for j in start:length(group),label in 1:a.nletters
        v = site_letter(W,group[j],label)
        A !== Free && !isone(w) && _letter(w,degree(w))==v.bits && continue
        next = mul(w,v)[2]
        push!(rows,next)
        _extend_rows!(rows,a,group,next,j,remaining-1)
    end
    return rows
end

"Plain moment relaxation from rows and (coefficient,word) objective terms."
function plain_program(rows::Vector{W},objective; C::Type{<:Number}=ComplexF64,capacity::Integer=length(rows), word_string=string) where {W<:AbstractWord}
    !isempty(rows) && allunique(rows) || throw(ArgumentError("Rows must be nonempty and distinct."))
    one(W) in rows || throw(ArgumentError("Rows must contain the identity."))
    table = MomentTable{W}(capacity)
    program = MomentProgram{C}()
    sink = MomentSink(table,program.forms)
    n = length(rows)
    entries = n*(n+1)÷2
    sizehint!(program.forms.offsets,entries+1)
    sizehint!(program.forms.ids,entries)
    sizehint!(program.forms.coefs,entries)
    for j in 1:n,i in j:n
        p,s = star(rows[i])
        mul!(sink,s,rows[j],scale(p,one(C)))
        close_form!(program.forms)
    end
    push_block!(program,C <: Complex ? HPSD : PSD,n)
    weights = zeros(Float64,nvars(table))
    missing = nothing
    for (c,w) in objective
        iszero(c) && continue
        key,p,conjugated = moment_key(w)
        re,im_id = moment_ids(table,key)
        if (re,im_id)==(0,0)
            isnothing(missing) && (missing=W[])
            push!(missing,w)
            continue
        end
        cc = scale(p,complex(float(c)))
        re != 0 && (weights[re] += real(cc))
        im_id != 0 && (weights[im_id] += real(cc*(conjugated ? -im : im)))
    end
    if !isnothing(missing)
        required = maximum(w->cld(degree(w),2),missing)
        throw(ArgumentError("Objective words $(join(word_string.(unique(missing)), ", ")) are not covered by the rows. Use complete rows through degree $required on all objective sites, or add custom rows whose adjoint products cover these words."))
    end
    for id in eachindex(weights)
        if !iszero(weights[id])
            push!(program.objective_ids,Int32(id)); push!(program.objective_coefs,weights[id])
        end
    end
    program.nvars = nvars(table)
    return validate(program),table
end
function plain_program(r::Registry,order::Integer,objective; kwargs...)
    W = word_type(r,2order)
    return plain_program(basis_rows(W,r.algebra,1:r.nsites,order),objective;kwargs...)
end

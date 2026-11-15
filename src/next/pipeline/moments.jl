"A complex moment expressed on at most two real coordinates. Zero ids are omitted."
struct OrbitCoordinates
    re::Int32
    im::Int32
    # Phase relations give Gaussian integers with components in -1:1,
    # including ±(1±im) for a rotated real line. No floating basis is stored.
    real_coefficient::Complex{Int8}
    imag_coefficient::Complex{Int8}
end
const ZERO_ORBIT_COORDINATES = OrbitCoordinates(0,0,0,0)

mutable struct OrbitMoments{W,G,D,C}
    generators::G
    annihilators::D
    cyclic::C
    cache::Dict{W,OrbitCoordinates}
    representatives::MomentTable{W}
    labels::Vector{MomentLabel{W}}
    pending::Vector{Tuple{W,Phase,Bool,Bool}}
    relations::Dict{W,Tuple{Phase,Bool}}
    zero_orbit::Bool
    line::Int # -1: complex; otherwise z = im^line * conj(z)
    real_form::Bool
    covered::BitSet
end
function OrbitMoments(::Type{W}, generators; annihilators=(),real_form=false,cyclic=nothing) where {W<:AbstractWord}
    representatives = MomentTable{W}()
    t = OrbitMoments(generators,annihilators,cyclic,Dict{W,OrbitCoordinates}(),
        representatives,representatives.labels,Tuple{W,Phase,Bool,Bool}[],Dict{W,Tuple{Phase,Bool}}(),false,-1,real_form,BitSet([1]))
    m = orbit_coordinates!(t,one(W))
    m.re == 1 && m.im == 0 && m.real_coefficient == 1 ||
        throw(ArgumentError("A symmetry must fix the identity with phase +1."))
    t
end

# Relations have the convention L(v) = p * conj^b(L(seed)). A repeated
# vertex either forces zero or constrains the seed to a real line in C.
function _orbit_relation!(t, v, p, b, head)
    t.zero_orbit |= t.real_form && !isreal(p)
    relation=get(t.relations,v,nothing)
    if !isnothing(relation)
        q,c = relation
        if b == c
            t.zero_orbit |= p != q
        else
            line = Int((c ? q*conj(p) : conj(q)*p).k)
            t.zero_orbit |= t.line >= 0 && t.line != line
            t.line = line
        end
        return false
    else
        t.relations[v] = (p,b)
        push!(t.pending,(v,p,b,head))
        return true
    end
end

function _orbit_visit!(t,v,p,b)
    isnothing(t.cyclic) && return _orbit_relation!(t,v,p,b,false)
    # Close a translation orbit before queueing its normalizer actions.
    # A repeated vertex still checks its phase and conjugation relation.
    head=true
    while _orbit_relation!(t,v,p,b,head)
        q,v=act(t.cyclic.generator,v)
        p=conj(q)*p
        head=false
    end
    nothing
end

function orbit_coordinates!(t::OrbitMoments{W}, w::W) where {W}
    cached=get(t.cache,w,nothing)
    !isnothing(cached) && return cached
    for g in t.annihilators
        p,u = act(g,w)
        if u == w && p.k == 0x02 && (!is_antiunitary(g) || moment_reality(w) == :real)
            return ZERO_ORBIT_COORDINATES
        end
    end
    empty!(t.pending); empty!(t.relations)
    t.zero_orbit=false; t.line=t.real_form ? 0 : -1
    _orbit_visit!(t,w,one(Phase),false)
    rep=w; i=1
    while i <= length(t.pending)
        v,p,b,head=t.pending[i]
        isless(v,rep) && (rep=v)
        q,u=star(v)
        _orbit_visit!(t,u,conj(q)*conj(p),!b)
        u,q,c=moment_key(v)
        _orbit_visit!(t,u,(c ? q*conj(p) : conj(q)*p),xor(b,c))
        if !isnothing(t.cyclic) && head
            # g*T*g^-1 = T^±1: one image determines the entire cyclic orbit.
            for g in t.cyclic.normalizers
                q,u=act(g,v)
                _orbit_visit!(t,u,conj(q)*p,b)
            end
        end
        remaining=isnothing(t.cyclic) ? t.generators : t.cyclic.remaining
        for g in remaining
            q,u=act(g,v)
            anti=is_antiunitary(g)
            _orbit_visit!(t,u,conj(q)*(anti ? conj(p) : p),xor(b,anti))
        end
        i+=1
    end
    pr,br=t.relations[rep]
    line=t.line < 0 ? -1 : Int((pr*pr*(br ? conj(Phase(t.line)) : Phase(t.line))).k)
    re=im_id=Int32(0); alpha=beta=0.0im
    if !t.zero_orbit
        if isone(rep)
            line == 0 || throw(ArgumentError("Symmetry relations contradict the normalized identity."))
            re=Int32(1); alpha=1.0+0im
        else
            # On a rotated line the stored variable is still Re L(rep);
            # its imaginary part is supplied by alpha, not another table id.
            reality = line == -1 ? :complex : line == 2 ? :imaginary : :real
            re,im_id = moment_ids!(t.representatives,rep,reality)
            if line != 2
                alpha=line == 1 ? 1.0+1im : line == 3 ? 1.0-1im : 1.0+0im
            end
            if line in (-1,2)
                beta=1.0im
            end
        end
    end
    if isnothing(t.cyclic)
        for (v,p,b,_) in t.pending
            conjugated=xor(b,br)
            phase=p*(conjugated ? pr : conj(pr))
            t.cache[v]=OrbitCoordinates(re,im_id,scale(phase,conjugated ? conj(alpha) : alpha),
                scale(phase,conjugated ? conj(beta) : beta))
        end
    else
        # Reuse the eight coordinate maps across a lattice orbit's aliases.
        coordinates=let re=re,im_id=im_id,alpha=alpha,beta=beta,pr=pr
            ntuple(Val(8)) do i
                conjugated=i>4
                phase=Phase(i-1)*(conjugated ? pr : conj(pr))
                OrbitCoordinates(re,im_id,scale(phase,conjugated ? conj(alpha) : alpha),
                    scale(phase,conjugated ? conj(beta) : beta))
            end
        end
        for (v,p,b,_) in t.pending
            t.cache[v]=coordinates[Int(p.k)+1+4xor(b,br)]
        end
    end
    t.cache[w]
end

"Read orbit coordinates without registering words; an empty tuple is a symmetry-forced zero."
function moment_terms(t::OrbitMoments{W},w::AbstractWord) where {W}
    packed=try
        _pack_word(W,w)
    catch e
        e isa InexactError || rethrow()
        return nothing
    end
    m=get(t.cache,packed,nothing)
    if isnothing(m)
        # Every represented orbit is cached in full. An unseen orbit may
        # still be forced to zero, including by the annihilator shortcut.
        # Use a scratch table: querying must not grow the solved table.
        probe=OrbitMoments(W,t.generators;annihilators=t.annihilators,real_form=t.real_form,cyclic=t.cyclic)
        zero=orbit_coordinates!(probe,packed)
        return zero.re==0 && zero.im==0 ? () : nothing
    end
    m.re==0 && m.im==0 && return ()
    m.re==0 && return ((m.im,m.imag_coefficient),)
    m.im==0 && return ((m.re,m.real_coefficient),)
    return ((m.re,m.real_coefficient),(m.im,m.imag_coefficient))
end

moment_covered(t::OrbitMoments,id::Integer) = id in t.covered

public refresh_coverage!
"""
    refresh_coverage!(table::OrbitMoments, program::MomentProgram)

Replace the coverage snapshot with coordinates in the final forms or objective.
Builders call this after assembly and strengthening. After manually extending or
editing a program, call this again before solving and querying its moments.
Registration alone does not make a coordinate covered. ID 1 is always covered.
"""
function refresh_coverage!(t::OrbitMoments,program::MomentProgram)
    empty!(t.covered)
    union!(t.covered,(1,),program.forms.ids,program.objective_ids)
    return t
end

"Compatibility evaluation; polynomial coverage is checked by the model after cancellation."
function moment_value(t::OrbitMoments,w::AbstractWord,values)
    coordinates=moment_terms(t,w)
    isnothing(coordinates) && return nothing
    value=0.0+0.0im
    for (id,c) in coordinates
        value+=c*(id==1 ? 1.0 : values[id-1])
    end
    return value
end

"Compatibility accessor for real signed moments used by the M1 oracle tests."
function orbit_moment!(t::OrbitMoments,w)
    m=orbit_coordinates!(t,w)
    m.re == 0 && m.im == 0 && return (Int32(0),one(Phase))
    m.im == 0 || throw(ArgumentError("This moment needs complex coordinates; use orbit_coordinates!."))
    k=findfirst(==(m.real_coefficient),(1,im,-1,-im))
    isnothing(k) && throw(ArgumentError("This moment lies on a rotated real line; use orbit_coordinates!."))
    m.re,Phase(k-1)
end

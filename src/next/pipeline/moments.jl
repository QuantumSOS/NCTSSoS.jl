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

mutable struct OrbitMoments{W,G,D}
    generators::G
    annihilators::D
    cache::Dict{W,OrbitCoordinates}
    representatives::MomentTable{W}
    labels::Vector{MomentLabel{W}}
    pending::Vector{W}
    relations::Dict{W,Tuple{Phase,Bool}}
    zero_orbit::Bool
    line::Int # -1: complex; otherwise z = im^line * conj(z)
    real_form::Bool
end
function OrbitMoments(::Type{W}, generators; annihilators=(),real_form=false) where {W<:AbstractWord}
    representatives = MomentTable{W}()
    t = OrbitMoments(generators,annihilators,Dict{W,OrbitCoordinates}(),
        representatives,representatives.labels,W[],Dict{W,Tuple{Phase,Bool}}(),false,-1,real_form)
    m = orbit_coordinates!(t,one(W))
    m.re == 1 && m.im == 0 && m.real_coefficient == 1 ||
        throw(ArgumentError("A symmetry must fix the identity with phase +1."))
    t
end

# Relations have the convention L(v) = p * conj^b(L(seed)). A repeated
# vertex either forces zero or constrains the seed to a real line in C.
function _orbit_visit!(t, v, p, b)
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
    else
        t.relations[v] = (p,b)
        push!(t.pending,v)
    end
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
        v=t.pending[i]; p,b=t.relations[v]
        isless(v,rep) && (rep=v)
        q,u=star(v)
        _orbit_visit!(t,u,conj(q)*conj(p),!b)
        u,q,c=moment_key(v)
        _orbit_visit!(t,u,(c ? q*conj(p) : conj(q)*p),xor(b,c))
        for g in t.generators
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
    for v in t.pending
        p,b=t.relations[v]; conjugated=xor(b,br)
        phase=p*(conjugated ? pr : conj(pr))
        t.cache[v]=OrbitCoordinates(re,im_id,scale(phase,conjugated ? conj(alpha) : alpha),
            scale(phase,conjugated ? conj(beta) : beta))
    end
    t.cache[w]
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

# Only the square compiler's unsigned site maps and optional ordinary K
# use this table. Pauli words are Hermitian, and K kills precisely odd-Y
# orbits. There is no phase transport between distinct site aliases.
struct SquareOrbitMoments{W}
    L::Int
    points::Matrix{NTuple{2,Int}}
    conjugation::Bool
    representatives::NC.MomentTable{W}
    labels::Vector{NC.MomentLabel{W}}
    zero_words::Set{W}
    covered::BitSet
end

function SquareOrbitMoments(::Type{W}, L::Int, conjugation::Bool) where {W<:NC.PauliWord}
    L>=2 && W==NC.word_type(NC.Pauli(),L^2,0) || throw(ArgumentError("Square moment word width does not match the lattice."))
    table = NC.MomentTable{W}()
    maps = [(swap,sx,sy) for swap in (false,true) for sx in (1,-1) for sy in (1,-1)]
    points = [begin
        y,x = divrem(site-1,L)
        u,v = swap ? (y,x) : (x,y)
        (mod(sx*u,L),mod(sy*v,L))
    end for site in 1:L^2, (swap,sx,sy) in maps]
    SquareOrbitMoments(L,points,conjugation,table,table.labels,Set{W}(),BitSet([1]))
end

@inline _square_wrapped(i,L) = i >= L ? i-L : i

@inline function _square_remap(w::NC.PauliWord{K}, points, p, dx, dy, L) where {K}
    x = z = ntuple(_->UInt64(0),Val(K))
    for site in NC.sites(w)
        @inbounds u,v = points[site,p]
        dst = _square_wrapped(u+dx,L) + L*_square_wrapped(v+dy,L)
        a,b = ((site-1)>>>6)+1,(dst>>>6)+1
        srcbit,dstbit = UInt64(1)<<((site-1)&63),UInt64(1)<<(dst&63)
        @inbounds begin
            !iszero(w.x[a]&srcbit) && (x = Base.setindex(x,x[b]|dstbit,b))
            !iszero(w.z[a]&srcbit) && (z = Base.setindex(z,z[b]|dstbit,b))
        end
    end
    NC.PauliWord(x,z)
end

@inline _square_mask_less(a,b) = (a.x,a.z)<(b.x,b.z)

function _square_canonical(t::SquareOrbitMoments{NC.PauliWord{K}}, w::NC.PauliWord{K}) where {K}
    isone(w) && return w
    L = t.L
    best = w
    if L<=64 && (K==1 || 64%L==0)
        # Minimize the primary mask first (x unless the word is Z-only).
        # With whole rows in limbs, a minimizing y shift must put a primary
        # site in the first row of a limb. Otherwise moving up one row
        # reduces the first nonzero limb without a carry across limbs.
        primary = NC.MaskSites(all(iszero,w.x) ? w.z : w.x)
        rows_per_limb = div(64,L)
        for p in axes(t.points,2)
            shifts = UInt64(0)
            for site in primary
                @inbounds _,y = t.points[site,p]
                for boundary in 0:rows_per_limb:L-1
                    shifts |= UInt64(1)<<mod(boundary-y,L)
                end
            end
            while !iszero(shifts)
                dy = trailing_zeros(shifts)
                shifts &= shifts-UInt64(1)
                # Within the first nonzero limb, the highest occupied row
                # decides the first part of x-shift comparison. Its leftmost
                # bit must be at x=0 or another left shift improves the mask.
                limb,row = typemax(Int),-1
                for site in primary
                    @inbounds _,y = t.points[site,p]
                    y = _square_wrapped(y+dy,L)
                    a = (L*y)>>>6
                    if a<limb || (a==limb && y>row)
                        limb,row = a,y
                    end
                end
                for site in primary
                    @inbounds x,y = t.points[site,p]
                    _square_wrapped(y+dy,L)==row || continue
                    v = _square_remap(w,t.points,p,iszero(x) ? 0 : L-x,dy,L)
                    _square_mask_less(v,best) && (best=v)
                end
            end
        end
    else
        # Unaligned limbs can split a physical row. Keep the exact order
        # without asserting the whole-row anchor rule for those widths.
        for p in axes(t.points,2), dy in 0:L-1, dx in 0:L-1
            v = _square_remap(w,t.points,p,dx,dy,L)
            _square_mask_less(v,best) && (best=v)
        end
    end
    best
end

_square_k_zero(t,w) = t.conjugation && isodd(sum(count_ones(x&z) for (x,z) in zip(w.x,w.z)))
_square_in_bounds(t,w) = iszero(t.L^2&63) || iszero((w.x[end]|w.z[end])>>(t.L^2&63))

function NC.orbit_coordinates!(t::SquareOrbitMoments{W}, w::W) where {W}
    _square_in_bounds(t,w) || throw(ArgumentError("Moment word is outside the square registry."))
    _square_register!(t,_square_canonical(t,w))
end

function _square_register!(t::SquareOrbitMoments, rep)
    if _square_k_zero(t,rep)
        # Ward rows on a zero word can contain nonzero terms, so keep the
        # canonical word even though it owns no moment variable.
        push!(t.zero_words,rep)
        return NC.ZERO_ORBIT_COORDINATES
    end
    id,_ = NC.moment_ids!(t.representatives,rep,:real)
    NC.OrbitCoordinates(id,0,1,0)
end

function _square_moment_terms(t::SquareOrbitMoments{W},w,allowed) where {W}
    w isa W && _square_in_bounds(t,w) || return nothing
    _square_k_zero(t,w) && return ()
    allowed || return nothing
    id,_ = NC.moment_ids(t.representatives,_square_canonical(t,w))
    iszero(id) ? nothing : ((id,1),)
end

NC.moment_terms(t::SquareOrbitMoments,w::NC.AbstractWord) = _square_moment_terms(t,w,true)
NC.moment_covered(t::SquareOrbitMoments,id::Integer) = id in t.covered
function NC.refresh_coverage!(t::SquareOrbitMoments,program::NC.MomentProgram)
    empty!(t.covered)
    union!(t.covered,(1,),program.forms.ids,program.objective_ids)
    t
end

_pricing_existing_terms(t::SquareOrbitMoments,w,allowed::Bool) = _square_moment_terms(t,w,allowed)
_pricing_native_table(::SquareOrbitMoments) = true

# Unsigned site maps commute with full collective spin. Ordinary K
# conjugates D_x and D_z to their negatives and commutes with D_y.
function NC._ward_seed_words(t::SquareOrbitMoments{W},derivations) where {W}
    all(derivations) do d
        d isa NC.SpinRotation && d.support == ntuple(length(d.support)) do i
            n = clamp(t.L^2-64(i-1),0,64)
            typemax(UInt64) >> (64-n)
        end
    end || throw(ArgumentError("Compact square Ward seeds require full collective spin rotations."))
    W[[l.word for l in t.labels];collect(t.zero_words)]
end
NC._ward_key(t::SquareOrbitMoments{W},w::W) where {W} = _square_canonical(t,w)
# Unsigned site maps preserve odd-Y parity and carry no moment phase.
NC._ward_coordinates!(t::SquareOrbitMoments{W},w::W,key::W) where {W} = _square_register!(t,key)

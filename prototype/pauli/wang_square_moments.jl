# Product-moment lookup for the checked square B4 group, not a row or cone
# reduction. Multiplication phases stay in the caller's OrbitFormSink.
struct WangSquareMoments{W}
    L::Int
    points::Matrix{NTuple{2,Int}}
    representatives::NC.MomentTable{W}
    labels::Vector{NC.MomentLabel{W}}
    cache::Dict{W,Int32}
    covered::BitSet
end

# Parity already proves each zero from the packed word, without an orbit walk.
NC._with_pso_zero_memo(f, ::WangSquareMoments, order) = f()

function _wang_square_table(H,sym,layout)
    n = H.registry.nsites
    L = isqrt(n)
    W = eltype(H.data.words)
    !isnothing(layout) && layout.shape === :square && layout.extent == L ||
        throw(ArgumentError("Compact B4 moments require a discovered square chart."))
    # Not every graph automorphism of C4 × C4 belongs to spatial D4.
    points = _square_points(L)
    any(axes(points,2)) do p
        all(eachindex(layout.grid)) do i
            x,y = points[i,p]
            layout.grid[i] == 1+x+L*y
        end
    end || throw(ArgumentError("The discovered chart is not a physical square chart."))
    spatial = [_symmetry_action(W,n,p) for p in _wang_spatial_maps(layout)]
    signs = _wang_uniform_parities(H,[NC._sector_terms(H,W).terms])
    K = NC.CliffordMap(NC.AntiunitaryPauli(one(W)))
    Set(sym.generators) == Set([spatial;signs;_workflow_axes(W,n)]) && sym.moment_only == (K,) ||
        throw(ArgumentError("Compact B4 requires precisely translations, D4, signed axes and conjugation."))
    reps = NC.MomentTable{W}()
    WangSquareMoments(L,points,reps,reps.labels,Dict(one(W)=>Int32(1)),BitSet([1]))
end

# Uniform X and Z conjugations require equal axis-count parities. Ordinary
# conjugation kills odd Y on a Hermitian word. Thus only all-even survives,
# and every uniform signed axis permutation then has phase +1.
function _wang_even_axes(w)
    x,z,y = UInt64(0),UInt64(0),UInt64(0)
    for i in eachindex(w.x)
        x ⊻= w.x[i]; z ⊻= w.z[i]; y ⊻= w.x[i] & w.z[i]
    end
    iseven(count_ones(x) | count_ones(z) | count_ones(y))
end

function _wang_square_key(t,w)
    # All nonempty orbits have an image containing the origin. Enumerating
    # precisely those images gives the same set from every orbit member;
    # sorting disjoint masks also removes global axis permutations.
    best = nothing
    for p in axes(t.points,2), site in NC.sites(w)
        x,y = t.points[site,p]
        v = _square_remap(w,t.points,p,t.L-x,t.L-y,t.L)
        a,b,c = v.x .& .~v.z, v.x .& v.z, v.z .& .~v.x
        b < a && ((a,b) = (b,a))
        c < b && ((b,c) = (c,b))
        b < a && ((a,b) = (b,a))
        image = NC.PauliWord(a .| b,b .| c)
        (isnothing(best) || _square_mask_less(image,best)) && (best=image)
    end
    best
end

# -1 is absent; 0 is symmetry-forced zero. In particular an LSO probe must
# not register a moment or mark it covered. Cache queried aliases only.
function _wang_square_id!(t,w,register)
    _square_in_bounds(t,w) || throw(ArgumentError("Moment outside the square registry."))
    _wang_even_axes(w) || return Int32(0)
    id = get(t.cache,w,Int32(-1))
    id > 0 && return id
    key = _wang_square_key(t,w)
    id = first(register ? NC.moment_ids!(t.representatives,key,:real) : NC.moment_ids(t.representatives,key))
    iszero(id) && return Int32(-1)
    return t.cache[w] = id
end

function NC.orbit_coordinates!(t::WangSquareMoments{W},w::W) where {W}
    id = _wang_square_id!(t,w,true)
    iszero(id) ? NC.ZERO_ORBIT_COORDINATES : NC.OrbitCoordinates(id,0,1,0)
end

function NC.moment_terms(t::WangSquareMoments{W},w::NC.AbstractWord) where {W}
    w isa W && _square_in_bounds(t,w) || return nothing
    id = _wang_square_id!(t,w,false)
    id < 0 ? nothing : iszero(id) ? () : ((id,Complex{Int8}(1)),)
end

NC.moment_covered(t::WangSquareMoments,id::Integer) = id in t.covered
function NC.refresh_coverage!(t::WangSquareMoments,program::NC.MomentProgram)
    empty!(t.covered)
    union!(t.covered,(1,),program.forms.ids,program.objective_ids)
    t
end

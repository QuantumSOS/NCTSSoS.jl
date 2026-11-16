"""
    NCTSSoS.operator_name(registry, letter) -> Symbol

Look up a degree-one operator word in registry. Detached words are interpreted
in the supplied registry; use NCTSSoS.variables(registry)[name] for named lookup.
For a scalar symbol use NCTSSoS.support_letters to inspect its operator dependencies.
"""
function operator_name(r::Registry,letter::AbstractWord)
    degree(letter)==1 && (!(letter isa FunctionalWord) || _ftail(letter)==1) ||
        throw(ArgumentError("operator_name needs one operator letter; use word_letters or support_letters first."))
    packed=try
        _pack_word(_model_word_type(r),letter)
    catch err
        err isa Union{MethodError,InexactError} || rethrow()
        throw(ArgumentError("operator_name needs a letter of the registry's algebra and site range."))
    end
    location=packed isa FunctionalWord ? _letter_location(packed) : only(_locations(packed))
    name=findfirst(==(location),r.names)
    isnothing(name) && throw(ArgumentError("No variable is named at site $(location[1]), label $(location[2]) in this registry."))
    return name
end

function _polynomial(r::Registry,w::W) where {W<:AbstractWord}
    data=PolynomialBuffer{W,ComplexF64}(1)
    push_term!(data,1,w)
    return _polynomial(r,data)
end

"""
    NCTSSoS.basis(registry, degree)

Return a sized iterator of coefficient-one polynomials through degree, including
the identity. Rows belong to registry and can be used in Relaxation(rows=...).
Normal words are generated directly in the registry's modeling word capacity;
polynomial wrappers are created on iteration. Negative or excessive degrees error.
"""
function basis(r::Registry,d::Integer)
    word_type(r,d) # Check the algebra's supported degree before allocating rows.
    W=_model_word_type(r)
    W<:FunctionalWord && d>_fcap(W) && throw(ArgumentError(
        "This registry's modeling words support degree at most $(_fcap(W)); reduce the requested basis degree."))
    return (_polynomial(r,w) for w in basis_rows(W,r.algebra,1:r.nsites,d))
end

"""
    NCTSSoS.block_rows(built::NCTSSoS.BuiltRelaxation)

Iterate the closed moment/localizing family as vectors of words, before symmetry
reduction. These are not the final cones: strengthening cones are excluded and
symmetry-reduced cones may use linear combinations of words. Use built.sizes.blocks
for the final cone count. Vectors on the block route are borrowed read-only data.
The plain route regenerates its rows on demand; this inspection may allocate,
but does not change the program, moment table or prepared solver cache.
"""
function block_rows(built::BuiltRelaxation)
    closure=built.sizes.closure
    isnothing(closure) || return (block.rows for block in closure.family)
    p,r=built.problem,built.relaxation
    W=eltype(built.table.keys)
    result=[_rows(p,r,W)]
    order=isnothing(r.order) ? maximum(degree,r.rows) : r.order
    registry=p.objective.registry
    for c in p.constraints
        c isa OperatorConstraint || continue
        push!(result,basis_rows(W,registry.algebra,1:registry.nsites,_localizing_order(order,c.polynomial)))
    end
    return result
end

"""
    NCTSSoS.block_sides(built::NCTSSoS.BuiltRelaxation) -> Vector{Int}

Return the side of each PSD or HPSD cone in the final built program, in program
order, after symmetry reduction and including PSD/HPSD strengthening cones.
Exclude ZERO, NONNEGATIVE and SOC blocks. An HPSD side is its complex matrix
dimension, before any backend realification; these are not native task sizes.
Unlike block_rows, this describes final cones, not their unreduced word families.
"""
block_sides(built::BuiltRelaxation) = [b.side for b in built.program.blocks if b.kind in (PSD,HPSD)]

_show_algebra(io,a) = print(io,nameof(typeof(a)))
_show_algebra(io,::SiteAlgebra{A}) where {A} = print(io,nameof(A))
function _show_algebra(io,a::FunctionalAlgebra{F}) where {F}
    print(io,nameof(F),"("); _show_algebra(io,a.base); print(io,")")
end
function Base.show(io::IO,r::Registry)
    print(io,"Registry("); _show_algebra(io,r.algebra)
    print(io,", ",r.nsites," sites, ",length(r)," variables)")
end
function Base.show(io::IO,::MIME"text/plain",r::Registry)
    show(io,r)
    isempty(r.names) && return
    print(io,"\nNames: ")
    names=sort!(collect(keys(r.names));by=name->r.names[name])
    for (i,name) in enumerate(names)
        i>12 && get(io,:limit,true) && return print(io,", …")
        i>1 && print(io,", ")
        print(io,name)
    end
end

function Base.show(io::IO,values::Union{AbstractVector{<:Polynomial},AbstractVector{<:AbstractVector{<:Polynomial}}})
    print(io,"[")
    for (i,value) in enumerate(values)
        if i>8 && get(io,:limit,true)
            print(io,", … (",length(values)," entries)")
            break
        end
        i>1 && print(io,", ")
        show(io,value)
    end
    print(io,"]")
end
Base.show(io::IO,::MIME"text/plain",values::Union{AbstractVector{<:Polynomial},AbstractVector{<:AbstractVector{<:Polynomial}}}) = show(io,values)

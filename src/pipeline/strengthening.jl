abstract type StrengtheningKind end

"Positivity of the normalized density matrix on these Pauli sites; valid for each state."
struct DensityMatrix <: StrengtheningKind
    sites::Vector{Int}
    DensityMatrix(sites) = new(_strengthening_sites(sites))
end

"""
Named Kramers-invariant cone; r is an odd-Y Pauli label on 3 or 4 sites;
valid for each state. For k sites, the sum of squared expectations of the
kept words is at most 2^(k−1)−1.
"""
struct InvariantCone <: StrengtheningKind
    sites::Vector{Int}
    r::String
    function InvariantCone(sites,r::AbstractString)
        ss=_strengthening_sites(sites)
        length(ss) in (3,4) && length(r)==length(ss) && all(c->c in "IXYZ",r) && isodd(count(==('Y'),r)) ||
            throw(ArgumentError("InvariantCone needs 3 or 4 sites and a matching IXYZ label with odd Y parity."))
        new(ss,String(r))
    end
end

"Commutator equalities for the moment words of the selected block family; valid for eigenstates of the objective."
struct Stationarity <: StrengtheningKind end

"Positivity of the objective's second variation on each selected row block; valid for the ground space of the objective."
struct StateOptimality <: StrengtheningKind end

_validity(::Union{DensityMatrix,InvariantCone}) = "valid for each state"
_validity(::Stationarity) = "valid for eigenstates of the objective"
_validity(::StateOptimality) = "valid for the ground space of the objective"
_strengthening_width(k::Union{DensityMatrix,InvariantCone},problem,order) = length(k.sites)
_strengthening_width(::Union{Stationarity,StateOptimality},problem,order) = 2order+degree(problem.objective)

function _strengthening_sites(sites)
    ss=Int[sites...]
    !isempty(ss) && issorted(ss) && allunique(ss) && all(>(0),ss) ||
        throw(ArgumentError("Strengthening sites must be nonempty, positive, sorted and distinct."))
    ss
end

struct StrengthenedStates{V,K}
    base::V
    kinds::K
end
function Base.show(io::IO,v::StrengthenedStates)
    show(io,v.base)
    for kind in v.kinds
        print(io,"; ",nameof(typeof(kind)),": ",_validity(kind))
    end
end
Base.show(io::IO,::MIME"text/plain",v::StrengthenedStates)=show(io,v)

# An operator-algebra cone, not a site-space moment matrix. Matrices live
# only on the requested density subsystem and are shared across its orbit.
struct PauliExtra{W}
    words::Vector{W}
    phases::Vector{Phase}
    matrices::Vector{Matrix{ComplexF64}}
    radius::Float64
    kind::ConeKind
    source::Int
end
_extra_key(e::PauliExtra)=(e.kind,e.radius,sort(e.words))
function _extra_image(e::PauliExtra,g)
    images=[act(g,w) for w in e.words]
    PauliExtra(last.(images),e.phases.*first.(images),e.matrices,e.radius,e.kind,e.source)
end

function _pauli_extra(kind,W,source)
    k=length(kind.sites)
    words=W[]; matrices=Matrix{ComplexF64}[]
    paulis=([1 0;0 1],[0 1;1 0],[0 -im;im 0],[1 0;0 -1])
    for code in 0:4^k-1
        digits=digits!(zeros(Int,k),code;base=4)
        word=foldl((w,sa)->last(mul(w,pauli_letter(W,sa[1],sa[2]-1))),
            ((site,axis) for (site,axis) in zip(kind.sites,digits) if axis!=0);init=one(W))
        if kind isa InvariantCone
            code==0 && continue
            # K contributes one minus sign per Y; conjugation by P_r
            # contributes one per unequal nonidentity pair.
            parity=count(==(2),digits)+count(zip(digits,kind.r)) do (a,b)
                a!=0 && b!='I' && "IXYZ"[a+1]!=b
            end
            isodd(parity) && continue
        else
            push!(matrices,foldl(kron,(ComplexF64.(paulis[a+1]) for a in digits))/2^k)
        end
        push!(words,word)
    end
    PauliExtra(words,fill(one(Phase),length(words)),matrices,
        kind isa InvariantCone ? sqrt(2.0^(k-1)-1) : 1.0,kind isa InvariantCone ? SOC : HPSD,source)
end

# Close the operator supports, not the group. Density positivity depends on
# the complete matrix algebra, so signed changes of its Pauli basis describe
# the same cone. SOC signs also leave its Euclidean norm unchanged.
function _extra_orbits(initial,generators)
    extras=empty(initial); index=Dict{Any,Int}(); representatives=Int[]
    for root in initial
        haskey(index,_extra_key(root)) && continue
        push!(extras,root); index[_extra_key(root)]=length(extras)
        push!(representatives,length(extras)); cursor=length(extras)
        while cursor<=length(extras)
            for g in generators
                image=_extra_image(extras[cursor],g)
                key=_extra_key(image)
                haskey(index,key) && continue
                push!(extras,image); index[key]=length(extras)
            end
            cursor+=1
        end
    end
    extras,representatives
end

function _append_extra!(program,table,extra;face=nothing)
    if extra.kind==SOC
        push_term!(program.forms,extra.radius,1); close_form!(program.forms)
        sink=OrbitFormSink(table,Dict{Int32,ComplexF64}())
        for (w,p) in zip(extra.words,extra.phases)
            empty!(sink.terms); push_term!(sink,scale(p,1.0+0im),w)
            _close_orbit_form!(program.forms,sink.terms;diagonal=true)
        end
        push_block!(program,SOC,length(extra.words)+1,extra.source)
        return 0.0
    end
    side=size(first(extra.matrices),1)
    coefficients=Dict{Int32,Matrix{ComplexF64}}()
    for (w,p,A) in zip(extra.words,extra.phases,extra.matrices)
        m=orbit_coordinates!(table,w)
        for (id,c) in ((m.re,m.real_coefficient),(m.im,m.imag_coefficient))
            id==0 && continue
            B=get!(() -> zeros(ComplexF64,side,side),coefficients,id)
            B .+= scale(p,c).*A
        end
    end
    ids=sort!(collect(keys(coefficients)))
    matrices=[coefficients[id] for id in ids]
    if !isnothing(face) && size(face.basis,2)<side
        _density_face_equations!(program,ids,matrices,face,extra.source)
        B=ComplexF64.(face.basis)
        matrices=[B'*A*B for A in matrices]
    end
    rows,error=_coefficient_sectors(matrices)
    for U in rows
        parts=[U'*A*U for A in matrices]
        real_block=all(A->maximum(abs,imag.(A))<=1e-11,parts)
        n=size(U,2)
        for j in 1:n,i in j:n
            for (id,A) in zip(ids,parts)
                c=real_block || i==j ? complex(real(A[i,j])) : A[i,j]
                abs(c)>1e-12 && push_term!(program.forms,c,id)
            end
            close_form!(program.forms)
        end
        push_block!(program,real_block ? PSD : HPSD,n,extra.source)
    end
    error
end

function _strengthen!(program,table,problem,kinds,closure;sectors=())
    isempty(kinds) && return nothing,nothing,1:0
    site_kinds=filter(k->k isa Union{DensityMatrix,InvariantCone},kinds)
    isempty(site_kinds) || problem.objective.registry.algebra isa Pauli ||
        throw(ArgumentError("DensityMatrix and InvariantCone require Pauli variables; for bosonic variables use Sector, and for fermionic variables use FermionicConditions or Sector."))
    all(k->maximum(k.sites)<=problem.objective.registry.nsites,site_kinds) ||
        throw(ArgumentError("A strengthening site is outside the registry; select existing sites."))
    W=eltype(table.representatives.keys)
    initial=PauliExtra{W}[_pauli_extra(kind,W,length(problem.constraints)+i)
        for (i,kind) in enumerate(kinds) if kind isa Union{DensityMatrix,InvariantCone}]
    extras,representatives=_extra_orbits(initial,closure.block_generators)
    start=length(program.blocks)
    faces=[_density_face(extras[i],sectors) for i in representatives]
    errors=[_append_extra!(program,table,extras[i];face) for (i,face) in zip(representatives,faces)]
    density_faces=[(;source=extras[i].source,before=size(face.basis,1),after=size(face.basis,2),sectors=face.sources)
        for (i,face) in zip(representatives,faces) if !isnothing(face)]
    hamiltonian=_hamiltonian_strengthen!(program,table,problem,kinds,closure)
    fermionic=_fermionic_strengthen!(program,table,problem,kinds,closure)
    sector_start=length(program.blocks)+1
    equations=_sector_equalities!(program,table,closure,sectors)
    sector_blocks=sector_start:length(program.blocks)
    program.nvars=length(table.labels)
    validate(program)
    noptimal=count(k->k isa StateOptimality,kinds)
    cones=isnothing(hamiltonian) ? extras[representatives] : [extras[representatives];hamiltonian.representatives]
    isnothing(fermionic) || (cones=[cones;fermionic.representatives])
    report=(initial_cones=length(initial)+noptimal*closure.initial_blocks+(isnothing(fermionic) ? 0 : fermionic.initial_cones),
        closed_cones=length(extras)+noptimal*closure.closed_blocks+(isnothing(fermionic) ? 0 : fermionic.closed_cones),representatives=cones,
        added_blocks=count(b->b.kind!=ZERO,program.blocks[start+1:end]),
        added_equalities=sum((b.side for b in program.blocks[start+1:end] if b.kind==ZERO);init=0),first_block=start+1,
        reconstruction_error=maximum(errors;init=0.0),kinds=kinds,hamiltonian=hamiltonian,fermionic=fermionic,density_faces=density_faces)
    return report,equations,sector_blocks
end

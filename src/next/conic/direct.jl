abstract type DirectBackend <: AbstractSolver end

"SCS native moment-form writer. Load SCS before solving; keywords are SCS settings."
struct SCSBackend{O<:NamedTuple} <: DirectBackend
    options::O
end
SCSBackend(;kwargs...) = SCSBackend((;kwargs...))
"Mosek affine-conic moment-form writer. Keywords use MSK_IPAR_/DPAR_/SPAR_ parameter names."
struct MosekBackend{O<:NamedTuple} <: DirectBackend
    options::O
end
MosekBackend(;kwargs...) = MosekBackend((;kwargs...))
"BPSDP native SOS writer. Load BPSDP before solving; keywords are BPSDP.Settings fields."
struct BPSDPBackend{O<:NamedTuple} <: DirectBackend
    options::O
end
BPSDPBackend(;kwargs...) = BPSDPBackend((;kwargs...))

"Construct a backend problem without solving it (and without needing a Mosek license)."
function prepare_standard(form,backend::DirectBackend)
    throw(ArgumentError("The extension for $(typeof(backend)) is not loaded; load its solver package with using SCS, Mosek, or BPSDP before constructing the problem."))
end
solve_standard(form,backend::DirectBackend) = solve_standard(prepare_standard(form,backend))

function _check_append_only(old::MomentStandardForm,form::MomentStandardForm)
    m,k = length(old.b),length(old.cones)
    unchanged = form.c==old.c && form.offset==old.offset && size(form.A,2)==size(old.A,2) &&
        length(form.b)>=m && length(form.cones)>=k && form.b[1:m]==old.b &&
        form.A[1:m,:]==old.A && form.cones[1:k]==old.cones
    unchanged || throw(ArgumentError("An open solver handle only accepts appended constraints with the same objective and variables; start a fresh solve for other changes."))
    return nothing
end

"Realify HPSD blocks as [Re(M) -Im(M); Im(M) Re(M)], preserving all other cones."
function real_psd_form(form::MomentStandardForm)
    any(c->c.kind==:hpsd,form.cones) || return form
    I,J,V = Int[],Int[],Float64[]
    cones = DataCone[]
    row = 0
    for cone in form.cones
        first = row+1
        if cone.kind!=:hpsd
            for k in cone.first:cone.first+cone.length-1
                row+=1; push!(I,row); push!(J,k); push!(V,1.0)
            end
            push!(cones,DataCone(cone.kind,cone.side,first,cone.length,cone.label))
            continue
        end
        n = cone.side
        # Position of a real/imaginary lower-triangle component in the original svec.
        coordinate(i,j) = cone.first+(j-1)*(2n-j+1)+2(i-j)
        for j in 1:2n,i in j:2n
            row+=1
            a,b = mod1(i,n),mod1(j,n)
            if (i>n)==(j>n)
                k = coordinate(a,b)-(a==b ? 0 : 1)
                push!(I,row); push!(J,k); push!(V,1.0)
            elseif a!=b
                k = coordinate(max(a,b),min(a,b))
                push!(I,row); push!(J,k); push!(V,a>b ? 1.0 : -1.0)
            end
        end
        push!(cones,DataCone(:psd,2n,first,n*(2n+1),cone.label))
    end
    transform = sparse(I,J,V,row,length(form.b))
    return MomentStandardForm(transform*form.A,transform*form.b,form.c,form.offset,cones)
end

"""
    real_program(program)

Convert the reduced Hermitian blocks to real PSD cones. Already real blocks
keep their side. A genuinely complex block uses [Re(A) -Im(A); Im(A) Re(A)].
This fallback doubles that reduced block; it does not claim an antiunitary
fixed-basis reduction. No unreduced site-space matrix is constructed.
"""
function real_program(p::MomentProgram)
    out = MomentProgram{Float64}(p.nvars)
    append!(out.objective_ids, p.objective_ids)
    append!(out.objective_coefs, p.objective_coefs)
    for b in p.blocks
        b.kind in (PSD, HPSD) || throw(ArgumentError("The M1 real conversion expects only PSD blocks."))
        is_real = all(isreal(p.forms.coefs[t]) for f in form_range(b) for t in term_range(p.forms,f))
        n, side = b.side, is_real ? b.side : 2b.side
        for j in 1:side, i in j:side
            row, col = mod1(i,n), mod1(j,n)
            imaginary = i > n && j <= n
            conjugated = row < col
            f = b.first_form + tri_offset(n, max(row,col), min(row,col))
            for t in term_range(p.forms,f)
                c = p.forms.coefs[t]
                coefficient = imaginary ? (conjugated ? -imag(c) : imag(c)) : real(c)
                iszero(coefficient) || push_term!(out.forms, coefficient, p.forms.ids[t])
            end
            close_form!(out.forms)
        end
        push_block!(out, PSD, side, b.label)
    end
    validate(out)
end

_anti_real_program(p, generators) = any(is_antiunitary, generators) ? real_program(p) : p

_conjugate_action(g,h) = compose(compose(g,h),inv(g))
_action_image(g,w,::Nothing) = act(g,w)
function _action_image(g,w,before)
    p,u = act(before,w)
    q,v = act(g,u)
    q*(is_antiunitary(g) ? conj(p) : p),v
end

# Q has orthonormal columns; image is their anti-linear image. A preserving
# involution admits a fixed real basis of the same complex dimension.
function _fixed_real_columns(Q,image)
    B=Q'*image
    LinearAlgebra.norm(image-Q*B)<=1e-9 || return nothing
    LinearAlgebra.norm(B*conj(B)-LinearAlgebra.I)<=1e-9 || return nothing
    fixed=Vector{ComplexF64}[]
    for j in axes(Q,2),odd in (false,true)
        v=odd ? (Q[:,j]-image[:,j])/2im : (Q[:,j]+image[:,j])/2
        for u in fixed
            c=LinearAlgebra.dot(u,v)
            abs(c)>1e-12 && (v .-= real(c).*u)
        end
        norm=LinearAlgebra.norm(v)
        norm>1e-9 || continue
        push!(fixed,v/norm)
    end
    length(fixed)==size(Q,2) || error("Antiunitary fixed-space rank mismatch.")
    fixed
end

function _real_finite_sectors(rows::AbstractVector{W},sectors,group,anti) where {W}
    result=[Vector{Tuple{ComplexF64,W}}[] for _ in sectors]
    index=Dict(w=>i for (i,w) in enumerate(rows))
    for (s,sector) in enumerate(sectors)
        isempty(sector) && continue
        Q=zeros(ComplexF64,length(rows),length(sector))
        for (j,row) in enumerate(sector),(c,w) in row
            Q[index[w],j]=c
        end
        found=false
        for g in group.elements
            image=zeros(ComplexF64,size(Q)); closed=true
            for (i,w) in enumerate(rows)
                p,u=_action_image(g,w,anti)
                if !haskey(index,u)
                    closed=false; break
                end
                for j in axes(Q,2)
                    image[index[u],j]+=scale(p,conj(Q[i,j]))
                end
            end
            closed || continue
            fixed=_fixed_real_columns(Q,image)
            isnothing(fixed) && continue
            result[s]=[[(c,rows[i]) for (i,c) in enumerate(v) if abs(c)>1e-12] for v in fixed]
            found=true; break
        end
        found || throw(ArgumentError("No same-side real basis was certified: no supplied antiunitary coset element preserves this projected sector as an involution. A degree-changing antiunitary needs a closed row family; exchanged or quaternionic sectors need a different decomposition. No doubled block was substituted."))
    end
    result
end

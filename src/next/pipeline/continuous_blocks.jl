# Numerical discovery never licenses a split: the coefficientwise certificate
# below checks the complete emitted decomposition, including repeated copies.
function _spectral_copies(matrices)
    seed=sum((sin(k*sqrt(2.0))/max(1.0,norm(A)))*A for (k,A) in enumerate(matrices))
    spectrum=eigen(Hermitian(seed))
    Q=spectrum.vectors
    groups=_eigenvalue_groups(spectrum.values,1e-8*max(1.0,norm(seed)))
    transformed=[Q'*A*Q for A in matrices]
    scales=[max(1.0,norm(A)) for A in transformed]
    transport=Dict{Int,Matrix{ComplexF64}}()
    components=Vector{Int}[]
    for root in eachindex(groups)
        haskey(transport,root) && continue
        d=length(groups[root])
        transport[root]=Matrix{ComplexF64}(I,d,d)
        nodes=[root]; cursor=1
        while cursor<=length(nodes)
            i=nodes[cursor]; cursor+=1
            for j in eachindex(groups)
                haskey(transport,j) && continue
                length(groups[j])==d || continue
                for (A,coefficient_scale) in zip(transformed,scales)
                    edge=A[groups[i],groups[j]]
                    scale=norm(edge)/sqrt(d)
                    scale>1e-8*coefficient_scale || continue
                    norm(edge*edge'/scale^2-I,Inf)<=1e-7 || continue
                    # The polar factor avoids accumulating a transport's
                    # singular-value error along a spanning tree.
                    F=svd(edge)
                    transport[j]=F.V*F.U'*transport[i]
                    push!(nodes,j)
                    break
                end
            end
        end
        push!(components,nodes)
    end
    return [[hcat([(Q[:,groups[i]]*transport[i])[:,a] for i in nodes]...)
        for a in 1:length(groups[first(nodes)])] for nodes in components]
end

"Check U's unitarity and A_k = Σ_i Σ_a U_ia B_ki U_ia' for every coefficient."
function _certify_coefficient_copies(matrices,copies;atol=1e-9)
    n=size(first(matrices),1)
    U=hcat((V for component in copies for V in component)...)
    size(U)==(n,n) || return (accepted=false,reconstruction_error=Inf,orthogonality_error=Inf,reduced=nothing)
    orthogonality=norm(U'*U-I,Inf)
    # B is retained as complex even when all entries happen to be real.
    reduced=[[first(component)'*A*first(component) for A in matrices] for component in copies]
    # Certify the Hermitian matrices actually emitted, including real diagonals.
    for component in reduced,B in component
        B .= (B+B')/2
    end
    error=0.0
    for (k,A) in enumerate(matrices)
        reconstructed=zeros(ComplexF64,n,n)
        for (component,coefficients) in zip(copies,reduced),V in component
            reconstructed .+= V*coefficients[k]*V'
        end
        error=max(error,norm(A-reconstructed,Inf)/max(1.0,norm(A,Inf)))
    end
    return (accepted=max(orthogonality,error)<=atol,reconstruction_error=error,
        orthogonality_error=orthogonality,reduced=reduced)
end

function _block_coefficients(program,block,ids)
    matrices=[zeros(ComplexF64,block.side,block.side) for _ in ids]
    index=Dict(id=>i for (i,id) in enumerate(ids))
    for j in 1:block.side,i in j:block.side
        f=block.first_form+tri_offset(block.side,i,j)
        for t in term_range(program.forms,f)
            A=matrices[index[program.forms.ids[t]]]
            c=program.forms.coefs[t]
            A[i,j]+=c
            i!=j && (A[j,i]+=conj(c))
        end
    end
    return matrices
end

function _copy_cone!(out,program,block)
    for f in form_range(block)
        for t in term_range(program.forms,f)
            push_term!(out.forms,program.forms.coefs[t],program.forms.ids[t])
        end
        close_form!(out.forms)
    end
    push_block!(out,block.kind,block.side,block.label)
end

"""
Split only family PSD/HPSD blocks, after finite reduction and exact continuous
moment elimination. No strengthening cone is changed. Dense storage is bounded
by both side and 16*K*n² bytes for the complex coefficient matrices; discovery
and certification also need working matrices. Every unsuccessful proposal keeps
the original sparse block. Tolerances use maximum-entry norms: absolute for
U'U-I and relative to max(1, maxabs(A_k)) for coefficient reconstruction.
"""
function _continuous_blocks(program;family_blocks=length(program.blocks),max_side=256,
    max_bytes=128*1024^2,atol=1e-9)
    out=MomentProgram{ComplexF64}(program.nvars)
    append!(out.objective_ids,program.objective_ids)
    append!(out.objective_coefs,program.objective_coefs)
    records=NamedTuple[]
    for (input,block) in enumerate(program.blocks)
        if input>family_blocks
            _copy_cone!(out,program,block)
            continue
        end
        ids=sort!(unique([program.forms.ids[t] for f in form_range(block) for t in term_range(program.forms,f)]))
        bytes=16*length(ids)*Int128(block.side)^2
        reason=block.kind ∉ (PSD,HPSD) ? :not_psd : block.side>max_side ? :side_cap :
            bytes>max_bytes ? :memory_cap : isempty(ids) ? :zero_block : :candidate
        first_output=length(out.blocks)+1
        reconstruction=orthogonality=0.0
        multiplicities=[1]
        if reason==:candidate
            matrices=_block_coefficients(program,block,ids)
            copies=_spectral_copies(matrices)
            if length(copies)==1 && length(only(copies))==1
                reason=:no_reduction
            else
                certificate=_certify_coefficient_copies(matrices,copies;atol)
                reconstruction=certificate.reconstruction_error
                orthogonality=certificate.orthogonality_error
                reason=certificate.accepted ? :split : :certificate_failed
                if certificate.accepted
                    multiplicities=length.(copies)
                    for (component,matrices) in zip(copies,certificate.reduced)
                        side=size(first(component),2)
                        # A real cone needs an already certified real route and
                        # a real change of basis; real-looking coefficients alone
                        # do not license a real representation of a complex group.
                        real_block=block.kind==PSD && all(V->all(isreal,V),component) && all(A->all(isreal,A),matrices)
                        for j in 1:side,i in j:side
                            for (id,A) in zip(ids,matrices)
                                c=i==j ? complex(real(A[i,j])) : A[i,j]
                                iszero(c) || push_term!(out.forms,c,id)
                            end
                            close_form!(out.forms)
                        end
                        push_block!(out,real_block ? PSD : HPSD,side,block.label)
                    end
                end
            end
        end
        reason==:split || _copy_cone!(out,program,block)
        outputs=collect(first_output:length(out.blocks))
        push!(records,(input=input,outputs=outputs,reason=reason,before=block.side,
            sides=[out.blocks[i].side for i in outputs],multiplicities=multiplicities,
            coefficient_bytes=bytes,reconstruction_error=reconstruction,orthogonality_error=orthogonality))
    end
    psd_records=filter(r->program.blocks[r.input].kind in (PSD,HPSD),records)
    kept=Dict(reason=>count(r->r.reason==reason,records) for reason in unique(r.reason for r in records) if reason!=:split)
    report=(method=:coefficient_spectral_transport,blocks=records,split=count(r->r.reason==:split,records),kept=kept,
        largest_before=maximum((r.before for r in psd_records);init=0),
        largest_after=maximum((maximum(r.sides) for r in psd_records);init=0),
        reconstruction_error=maximum((r.reconstruction_error for r in records);init=0.0),
        orthogonality_error=maximum((r.orthogonality_error for r in records);init=0.0),
        max_coefficient_bytes=maximum((r.coefficient_bytes for r in records);init=Int128(0)),
        side_cap=max_side,memory_cap=max_bytes,tolerance=atol)
    return validate(out),report
end

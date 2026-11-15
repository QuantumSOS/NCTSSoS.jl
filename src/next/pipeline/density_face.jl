# A density algebra can contain a sector even after a Clifford map changes
# its site support. Test membership in its complete operator basis, not sites.
function _density_face(extra,sectors)
    extra.kind==HPSD && !isempty(sectors) || return nothing
    index=Dict(w=>i for (i,w) in enumerate(extra.words))
    contained=filter(s->all(w->haskey(index,w),keys(last(s).terms)),sectors)
    isempty(contained) && return nothing
    C=Complex{Rational{BigInt}}
    n=size(first(extra.matrices),1)
    operators=map(contained) do (source,s)
        A=zeros(C,n,n)
        for (w,c) in s.terms
            k=index[w]
            # extra.matrices represent rho, so remove their normalization.
            factor=scale(conj(extra.phases[k]),c)*n
            P=extra.matrices[k]
            for j in 1:n,i in 1:n
                iszero(P[i,j]) || (A[i,j]+=factor*C(P[i,j]))
            end
        end
        source=>A
    end
    _density_kernel(operators)
end

function _density_kernel(operators)
    C=eltype(last(first(operators)))
    n=size(last(first(operators)),1)
    pivots=Dict{Int,Dict{Int,C}}()
    for (source,A) in operators
        for i in 1:n
            row=Dict(j=>A[i,j] for j in 1:n if !iszero(A[i,j]))
            _insert_generator_equation!(pivots,row;constant=false)
        end
        length(pivots)==n && throw(ArgumentError(
            "Sector $source has an empty common kernel for DensityMatrix: the selected value is not an eigenvalue of the operator on these sites within the preceding sectors; choose a compatible sector value."))
    end
    free=setdiff(1:n,keys(pivots))
    B=zeros(C,n,length(free))
    for (k,j) in enumerate(free); B[j,k]=1; end
    for i in sort!(collect(keys(pivots));rev=true),k in eachindex(free)
        for (j,c) in pivots[i]
            j==i || iszero(B[j,k]) || (B[i,k]-=c*B[j,k])
        end
    end
    # These independent rows span all (Q-q) rows. R*rho=0 includes
    # off-diagonal entries; keeping only B'*rho*B would be a relaxation.
    R=zeros(C,length(pivots),n)
    for (k,i) in enumerate(sort!(collect(keys(pivots)))),(j,c) in pivots[i]
        R[k,j]=c
    end
    (;basis=B,annihilator=R,sources=first.(operators))
end

function _density_face_equations!(program,ids,matrices,face,source)
    C=eltype(face.basis)
    R=sparse(face.annihilator)
    products=[R*C.(sparse(A)) for A in matrices]
    seen=Set{Vector{Pair{Int32,Rational{BigInt}}}}()
    count=0
    for j in axes(face.annihilator,2),i in axes(face.annihilator,1),part in (real,imag)
        row=[id=>part(A[i,j]) for (id,A) in zip(ids,products) if !iszero(part(A[i,j]))]
        isempty(row) && continue
        pivot=last(first(row))
        key=[id=>c/pivot for (id,c) in row]
        key in seen && continue
        push!(seen,key)
        for (id,c) in key; push_term!(program.forms,c,id); end
        close_form!(program.forms)
        count+=1
    end
    count==0 || push_block!(program,ZERO,count,source)
    nothing
end

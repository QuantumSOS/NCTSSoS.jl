# Common Hermitian commutant, restricted first to the eigenspaces of one
# generic coefficient combination. This is small-subsystem linear algebra,
# not continuous-group reduction. Check every coefficient after splitting.
using LinearAlgebra: I, norm, eigen, Hermitian, svd

function _coefficient_sectors(matrices;atol=1e-9)
    n=size(first(matrices),1)
    eye=Matrix{ComplexF64}(I,n,n)
    all(A->norm(A-A',Inf)<=atol,matrices) || error("Non-Hermitian strengthening coefficient")
    seed=sum((sin(j*sqrt(2.0))/max(1.0,norm(A)))*A for (j,A) in enumerate(matrices))
    spectrum=eigen(Hermitian(seed))
    groups=_eigenvalue_groups(spectrum.values,1e-8*max(1.0,norm(seed)))
    candidates=Matrix{ComplexF64}[]
    real_family=all(A->all(isreal,A),matrices)
    for group in groups,(j,q) in enumerate(group),p in group[1:j]
        u,v=spectrum.vectors[:,p],spectrum.vectors[:,q]
        push!(candidates,p==q ? u*u' : (u*v'+v*u')/sqrt(2.0))
        (p==q || real_family) || push!(candidates,im*(u*v'-v*u')/sqrt(2.0))
    end
    for A in matrices
        length(candidates)==1 && break
        A=A/max(1.0,norm(A))
        commutators=[A*X-X*A for X in candidates]
        all(C->norm(C,Inf)<=atol,commutators) && continue
        C=reduce(hcat,vec.(commutators))
        decomposition=svd([real.(C);imag.(C)])
        null=findall(s->s<=1e-10*max(1.0,maximum(decomposition.S)),decomposition.S)
        isempty(null) && error("Strengthening commutant lost the identity")
        basis=reduce(hcat,vec.(candidates))*decomposition.V[:,null]
        candidates=[reshape(basis[:,j],n,n) for j in axes(basis,2)]
    end
    splitter=sum(sin(j*sqrt(2.0))*X for (j,X) in enumerate(candidates))
    spectrum=eigen(Hermitian(splitter))
    groups=_eigenvalue_groups(spectrum.values,1e-7*max(1.0,norm(splitter)))
    length(groups)==1 && return [eye],0.0
    rows=[spectrum.vectors[:,group] for group in groups]
    residual=maximum(norm(A-sum(U*(U'*A*U)*U' for U in rows),Inf)/max(1.0,norm(A,Inf)) for A in matrices)
    residual<=atol || return [eye],0.0 # Keep the original cone if numerical splitting is uncertified.
    rows,residual
end

function _eigenvalue_groups(values,tolerance)
    groups=Vector{Int}[]
    for i in eachindex(values)
        if isempty(groups) || abs(values[i]-values[first(last(groups))])>tolerance
            push!(groups,[i])
        else
            push!(last(groups),i)
        end
    end
    groups
end

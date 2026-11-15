using LinearAlgebra
BLAS.set_num_threads(4)

function ring_ground(n)
    H=zeros(1<<n,1<<n)
    for state in 0:(1<<n)-1, i in 0:n-1
        j=mod(i+1,n)
        different=((state>>i)&1)!=((state>>j)&1)
        H[state+1,state+1]+=different ? -0.25 : 0.25
        if different
            target=state ⊻ (1<<i) ⊻ (1<<j)
            H[target+1,state+1]+=0.5
        end
    end
    e=eigen(Symmetric(H))
    energy=first(e.values); psi=e.vectors[:,1]
    println("RING n=",n," energy=",energy," scaled_to_20=",20energy/n,
        " eigen_residual=",norm(H*psi-energy*psi))
    energy,psi
end

for n in (4,6,8)
    ring_ground(n)
end

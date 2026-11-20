using LinearAlgebra

function bonds(lattice, L)
    iseven(L) && L >= 4 || error("expected even extent ≥ 4 (no doubled periodic bonds)")
    lattice == "chain" && return [(i, mod1(i + 1, L)) for i in 1:L]
    lattice == "square" || error("unknown lattice $lattice")
    return vcat([(x + L*(y-1), mod1(x+1,L) + L*(y-1)) for x in 1:L for y in 1:L],
                [(x + L*(y-1), x + L*(mod1(y+1,L)-1)) for x in 1:L for y in 1:L])
end

# Each spin pair contributes ±1/4 on the diagonal and 1/2 on a spin exchange.
# With an even number of spins and SU(2) invariance every multiplet has Sz=0.
function exact_energy(lattice, L)
    edges = bonds(lattice, L)
    n = lattice == "chain" ? L : L^2
    n <= 16 || error("ED is only supported for at most 16 sites")
    states = [s for s in 0:(1<<n)-1 if count_ones(s) == n÷2]
    indices = zeros(Int, 1<<n)
    indices[states .+ 1] = eachindex(states)
    h = zeros(length(states), length(states))
    for (col, s) in enumerate(states), (i, j) in edges
        if ((s>>(i-1)) ⊻ (s>>(j-1))) & 1 == 0
            h[col, col] += 0.25
        else
            h[col, col] -= 0.25
            row = indices[(s ⊻ (1<<(i-1)) ⊻ (1<<(j-1))) + 1]
            h[row, col] += 0.5
        end
    end
    return only(eigvals!(Symmetric(h), 1:1)) / n
end

# Periodic spin-1/2 XXX ring only. Ground-state quantum numbers and logarithmic
# Bethe equations: https://pure.uva.nl/ws/files/29914613/PhysRevB..93.pdf, (4–5).
# The +1/4 per site restores the constant subtracted in that paper's H.
function bethe_ring(n)
    iseven(n) && n >= 4 || error("expected even periodic ring")
    m = n ÷ 2
    q = (1:m) .- (m+1)/2
    u = asinh.(tanpi.(2q/n)) ./ pi
    for _ in 1:30
        d = u .- transpose(u)
        f = n .* atan.(2u) .- pi .* q .- vec(sum(atan.(d); dims=2))
        if norm(f,Inf) < 1e-12 && all(>(0),diff(u)) && norm(u+reverse(u),Inf) < 1e-12
            return (; energy=0.25-sum(inv.(u.^2 .+ 0.25))/(2n), residual=norm(f,Inf))
        end
        k = inv.(1 .+ d.^2)
        u -= (Diagonal(2n ./ (1 .+ 4u.^2) .- vec(sum(k; dims=2))) + k) \ f
    end
    error("ground-state Bethe equations did not converge")
end

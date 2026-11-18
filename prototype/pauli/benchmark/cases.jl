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

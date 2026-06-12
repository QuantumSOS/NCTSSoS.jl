#!/usr/bin/env julia
# Micro-profile of NormalMonomial hot paths + end-to-end attribution.
# Run with: julia --project benchmark/monomial_profile.jl

using NCTSSoS
using NCTSSoS: simplify!, simplify, neat_dot, _neat_dot3, _mul_term,
    _unchecked_monomial, _process_terms!, NormalMonomial, PauliAlgebra,
    NonCommutativeAlgebra, FermionicAlgebra
using Printf
using Profile
using Random

# ─── tiny harness: ns/op and bytes/op ────────────────────────────────────────
function measure(f::F; N::Int=100_000) where {F}
    f()  # warmup / compile
    GC.gc()
    t0 = time_ns()
    allocs = @allocated begin
        for _ in 1:N
            f()
        end
    end
    t1 = time_ns()
    return (t1 - t0) / N, allocs / N
end

results = Tuple{String,Float64,Float64}[]
function report!(name, f; N=100_000)
    ns, bytes = measure(f; N)
    push!(results, (name, ns, bytes))
    @printf("%-52s %10.1f ns/op %10.1f B/op\n", name, ns, bytes)
end

println("Julia ", VERSION, " | NCTSSoS monomial hot-path profile")
println("="^88)

function run_micros()
# ─── Pauli algebra ───────────────────────────────────────────────────────────
reg, (sx, sy, sz) = create_pauli_variables(1:10)
Tp = eltype(sx[1].word)
println("\n[Pauli]  word eltype = $Tp")

m_a = monomials(sx[1] * sy[3] * sz[5])[1]   # degree-3 canonical monomial
m_b = monomials(sy[2] * sz[3] * sx[7])[1]
wa, wm, wb = m_a.word, m_b.word, monomials(sz[4] * sx[6])[1].word

report!("pauli: _mul_term (deg3 x deg3)", () -> _mul_term(m_a, m_b))
report!("pauli: full * -> Polynomial", () -> m_a * m_b)
report!("pauli: simplify!(neat_dot3) moment kernel",
    () -> simplify!(PauliAlgebra, _neat_dot3(wa, wm, wb)))

# pure simplify! on a pre-allocated scratch word (steady-state cost, no alloc expected)
scratch = zeros(Tp, 9)
master = vcat(reverse(wa), wm, wb)
report!("pauli: simplify! only (reused buffer, deg9 raw)",
    () -> begin
        resize!(scratch, length(master))
        copyto!(scratch, master)
        simplify!(PauliAlgebra, scratch)
    end)

report!("pauli: hash(m)", () -> hash(m_a))
report!("pauli: m == m", () -> m_a == m_b)
report!("pauli: isless", () -> isless(m_a, m_b))
report!("pauli: adjoint(m)", () -> adjoint(m_a))

# Dict workload: the moment-matrix dedup pattern
basis_monos = unique([monomials(sx[i] * sy[j])[1] for i in 1:10 for j in 1:10 if i != j])
d = Dict(m => i for (i, m) in enumerate(basis_monos))
report!("pauli: Dict lookup (basis dedup pattern)",
    () -> d[basis_monos[37]])

# sort! of term tuples — the _process_terms! pattern
C = ComplexF64
master_terms = [(rand(C), m) for m in basis_monos]
buf = copy(master_terms)
report!("pauli: _process_terms! (90 terms)",
    () -> begin
        resize!(buf, length(master_terms))
        copyto!(buf, master_terms)
        _process_terms!(buf, C)
    end; N=10_000)

# polynomial product: Heisenberg-style H * H
H = sum(0.25 * (sx[i] * sx[i+1] + sy[i] * sy[i+1] + sz[i] * sz[i+1]) for i in 1:9)
report!("pauli: H * H (27-term x 27-term poly)", () -> H * H; N=1_000)

# ─── NonCommutative algebra ──────────────────────────────────────────────────
regn, (x,) = create_noncommutative_variables([("x", 1:10)])
Tn = eltype(x[1].word)
println("\n[NonCommutative]  word eltype = $Tn")
nm_a = monomials(x[1] * x[2] * x[3])[1]
nm_b = monomials(x[2] * x[3] * x[4])[1]
report!("nc: _mul_term (deg3 x deg3)", () -> _mul_term(nm_a, nm_b))
report!("nc: full * -> Polynomial", () -> nm_a * nm_b)
report!("nc: hash(m)", () -> hash(nm_a))

# ─── Fermionic algebra (PBW path) ────────────────────────────────────────────
regf, (a, ad) = create_fermionic_variables(1:6)
println("\n[Fermionic]  word eltype = $(eltype(a[1].word))")
fm_a = monomials(ad[1] * a[2])[1]
fm_b = monomials(ad[2] * a[3])[1]
report!("fermi: full * -> Polynomial (reordering path)", () -> fm_a * fm_b)
report!("fermi: adjoint via poly", () -> adjoint(Polynomial(fm_a)); N=10_000)
return reg, sx, sy, sz
end

reg, sx, sy, sz = run_micros()

# ─── end-to-end attribution ──────────────────────────────────────────────────
println("\n[End-to-end] compute_sparsity + moment_relax, Heisenberg n=10, order=3")
H8 = sum(0.25 * (sx[i] * sx[i+1] + sy[i] * sy[i+1] + sz[i] * sz[i+1]) for i in 1:9)
pop = polyopt(H8, reg)
cfg = SolverConfig(optimizer=nothing, order=3)

build() = begin
    sp = compute_sparsity(pop, cfg)
    NCTSSoS.moment_relax(pop, sp.corr_sparsity, sp.cliques_term_sparsities)
end
build()  # compile
GC.gc()
stats = @timed build()
@printf("wall: %.2f s | alloc: %.2f GiB | gc: %.1f%%\n",
    stats.time, stats.bytes / 2^30, 100 * stats.gctime / stats.time)

Profile.clear()
@profile for _ in 1:3; build(); end
io = IOBuffer()
Profile.print(io; format=:flat, sortedby=:count, mincount=10)
lines = split(String(take!(io)), '\n')
srclines = filter(l -> occursin("src/", l) && occursin("NCTSSoS", l), lines)
println("\nTop NCTSSoS src frames (flat, sorted by count):")
foreach(println, first(srclines, 30))

# allocation attribution
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 0.01 build()
allocres = Profile.Allocs.fetch()
sitecounts = Dict{String,Int}()
for a in allocres.allocs
    for fr in a.stacktrace
        file = String(fr.file)
        if occursin("NCTSSoS/src", file) || occursin("NCTSSoS.jl/src", file)
            key = string(basename(file), ":", fr.line, " ", fr.func)
            sitecounts[key] = get(sitecounts, key, 0) + 1
            break
        end
    end
end
println("\nTop allocation sites (sampled, NCTSSoS src only):")
for (k, v) in sort(collect(sitecounts); by=last, rev=true)[1:min(end, 20)]
    @printf("%6d  %s\n", v, k)
end

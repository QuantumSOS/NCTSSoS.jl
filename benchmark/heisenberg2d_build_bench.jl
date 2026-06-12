# 5x5 2D nearest-neighbor Heisenberg: build-path benchmark
# (compute_sparsity + moment_relax), no solver.
using NCTSSoS, Printf

const L = 5
site(i, j) = (i - 1) * L + j           # 1..25
reg, (sx, sy, sz) = create_pauli_variables(1:L*L)

edges = Tuple{Int,Int}[]
for i in 1:L, j in 1:L
    i < L && push!(edges, (site(i, j), site(i + 1, j)))
    j < L && push!(edges, (site(i, j), site(i, j + 1)))
end
H = sum(0.25 * (sx[a] * sx[b] + sy[a] * sy[b] + sz[a] * sz[b]) for (a, b) in edges)
println("5x5 Heisenberg: $(L*L) sites, $(length(edges)) bonds, $(length(NCTSSoS.terms(H))) terms")

pop = polyopt(H, reg)
order = parse(Int, get(ENV, "BENCH_ORDER", "2"))
cfg = SolverConfig(optimizer=nothing, order=order)

build() = begin
    sp = compute_sparsity(pop, cfg)
    mp = NCTSSoS.moment_relax(pop, sp.corr_sparsity, sp.cliques_term_sparsities)
    (sp, mp)
end

println("order = $order | warmup/compile run...")
t_first = @elapsed build()
@printf("first call (incl. compile): %.1f s\n", t_first)

runs = parse(Int, get(ENV, "BENCH_RUNS", "3"))
for r in 1:runs
    GC.gc()
    stats = @timed build()
    @printf("run %d: wall %.2f s | alloc %.2f GiB | gc %.1f%%\n",
        r, stats.time, stats.bytes / 2^30, 100 * stats.gctime / stats.time)
end

# ─── structural fingerprint (A/B equivalence proof) ──────────────────────────
# SHA-256 over cliques, term-sparsity supports, block bases, total basis,
# n_unique, and every constraint-matrix polynomial (coeffs + words).
using SHA
let
    sp = compute_sparsity(pop, cfg)
    mp = NCTSSoS.moment_relax(pop, sp.corr_sparsity, sp.cliques_term_sparsities)

    ctx = SHA.SHA256_CTX()
    feed(x::AbstractString) = update!(ctx, codeunits(x))
    feedmono(m) = feed(string("w[", join(Int.(m.word), ","), "]"))
    feedpoly(p) = begin
        feed("poly:")
        for (c, m) in NCTSSoS.terms(p)
            feed(string(c)); feedmono(m)
        end
    end

    feed("cliques:")
    foreach(c -> feed(string(Int.(c))), sp.corr_sparsity.cliques)
    feed("term_sparsities:")
    for tss in sp.cliques_term_sparsities, ts in tss
        feed("supp:"); foreach(feedmono, ts.term_sparse_graph_supp)
        feed("blocks:")
        for bb in ts.block_bases
            feed("b:"); foreach(feedmono, bb)
        end
    end
    feed("total_basis:")
    foreach(feedmono, mp.total_basis)
    feed(string("n_unique:", mp.n_unique_moment_matrix_elements))
    feed("constraints:")
    for (cone, mat) in mp.constraints
        feed(string(cone, size(mat)))
        foreach(feedpoly, mat)
    end
    feedpoly(mp.objective)

    @printf("cliques: %d (sizes %s) | total_basis: %d | n_unique: %d\n",
        length(sp.corr_sparsity.cliques), string(sort(length.(sp.corr_sparsity.cliques))),
        length(mp.total_basis), mp.n_unique_moment_matrix_elements)
    println("SHA-256 fingerprint: ", bytes2hex(digest!(ctx)))
end

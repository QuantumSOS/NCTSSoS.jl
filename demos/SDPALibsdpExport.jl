module SDPALibsdpExport

# Export a NCTSSoS `MomentProblem` to a libsdp-friendly complex SDPA-sparse
# format that the `complex_sdp_solver` (BPSDP backend) can ingest.
#
# Math, in one paragraph:
#
# A `MomentProblem` is "minimize ⟨H⟩ over moments y_α s.t. M(y) ⪰ 0 and
# Σ ℓ y = 0".  The moment matrix M is block-diagonal in our PQG demo, with each
# block a complex Hermitian PSD constraint produced by NCTSSoS as a
# `Matrix{Polynomial}` whose (r,c) entry is a polynomial in the moments.
# libsdp's primal variable X is exactly that block-diagonal Hermitian PSD
# matrix.  Each canonical moment α appears in some HPSD block at a "pivot"
# position pivot[α]; we *define* y_α := X[pivot[α]].  Every other HPSD entry
# (b,r,c) gets a linear equality that says "X[(b,r,c)] equals the polynomial
# expansion in the y_α's".  Every entry of every `:Zero` matrix in the
# `MomentProblem` becomes a linear equality among the y_α's.  The identity
# moment is pinned to 1.
#
# libsdp computes constraints as `Re ⟨F_k, X⟩ = b_k` (where `⟨A,B⟩ = tr(A^† B)`
# in the Hermitian sense, and the `Re` is taken because b is real-valued).
# A complex linear functional on X therefore splits into TWO libsdp rows:
# one for the real part of the equation, one for the imaginary part.  We
# encode `Σ v_α X[pivot[α]] = b` as
#
#   real-part row:  F[pivot[α]] = conj(v_α),       b_real = Re(b)
#   imag-part row:  F[pivot[α]] = i·conj(v_α),     b_imag = Im(b)
#
# This is sound because Re ⟨F, X⟩ = Σ Re(conj(F[α]) · X[α]) — see
# `ComplexSDPHelper::evaluate_Au` in libsdp/src/sdp_helper.cc.
#
# Output: a single text file in the libsdp-complex SDPA-sparse format
# defined here:
#
#   <m>                                   # number of equality constraints
#   <nblocks>                             # number of PSD blocks
#   <n_1> <n_2> ... <n_nblocks>           # block dimensions (positive ints)
#   <b_1> <b_2> ... <b_m>                 # right-hand sides (real doubles)
#   <k> <blk> <r> <c> <re> <im>           # one entry per line; k=0 is F_0 (objective),
#   ...                                   # k=1..m is F_1..F_m
#
# Indices are 1-based.  Entries are written verbatim — there is NO automatic
# Hermitian mirroring (libsdp's evaluate_Au sums user-supplied entries
# directly).  A companion text metadata sidecar records the pivot map and
# constraint provenance so a downstream user can map the libsdp solution back
# to moments.

using NCTSSoS
using NCTSSoS: NormalMonomial, Polynomial, MomentProblem, AlgebraType,
    FermionicAlgebra, symmetric_canon, expval, monomials, coefficients
using Printf

export export_libsdp, LibsdpExportSummary

# A 1-indexed (block, row, col) location inside the block-diagonal X.
const Pivot = Tuple{Int,Int,Int}

# A pivot together with the multiplicative phase φ such that
#     X[pivot] = φ · y_α     ⇔     y_α = X[pivot] / φ.
# For fermionic problems |φ| = 1 (typically ±1, occasionally ±i for some algebras).
struct PhasedPivot
    pos::Pivot
    phase::ComplexF64
end

# (k, blk, r, c, re, im) — one row in the libsdp-complex sparse file.
struct SparseEntry
    k::Int
    blk::Int
    r::Int
    c::Int
    re::Float64
    im::Float64
end

struct LibsdpExportSummary
    n_blocks::Int
    n_hpsd_blocks::Int
    n_aux_blocks::Int
    block_dims::Vector{Int}
    n_constraints::Int
    n_objective_entries::Int
    n_constraint_entries::Int
    n_unique_moments::Int
    n_hpsd_pivots::Int
    n_aux_pivots::Int
    n_zero_constraint_components::Int
    n_dropped_orphan_terms::Int
    sense::Symbol
end

# -----------------------------------------------------------------------------
# Pivot discovery
#
# Walk every :HPSD block.  For each (b, r, c) entry M[r,c] that is a single-term
# polynomial `coef * mono` (with |coef| close to 1 — the typical bare-moment-
# matrix case for fermionic problems), record (b, r, c) as the pivot for
# canon(mono) if no pivot is registered yet.  We prefer entries with coef == 1
# (purely "X[pivot] = y_α") so the moment-binding equalities remain trivial
# at pivots.

const PIVOT_COEF_TOL = 1e-10

# Accept any unit-magnitude complex coefficient as a pivot phase.
# Fermionic anticommutation routinely produces coef = -1; complex algebras may
# also produce ±i.  We do *not* accept coefficients with magnitude far from 1
# because then the entry is a sum or scaled term, not a clean pivot candidate.
function _is_unit_coef(coef)
    return abs(abs(coef) - 1.0) < PIVOT_COEF_TOL
end

function _prefer_real_unit(existing_phase::ComplexF64, candidate_phase::ComplexF64)
    # Prefer +1 > -1 > +i > -i > anything else, so the resulting equalities
    # have the simplest possible coefficients.
    score(c) = (abs(imag(c)) < PIVOT_COEF_TOL ? 0 : 1, -real(c), -imag(c))
    return score(candidate_phase) < score(existing_phase)
end

function _canonical_of(mono)
    return symmetric_canon(expval(mono))
end

function _terms(poly::Polynomial)
    # NCTSSoS doesn't expose `terms` on Polynomial uniformly; rebuild from
    # `coefficients` + `monomials`.
    return zip(coefficients(poly), monomials(poly))
end

function _polynomial_terms_with_canon(poly::Polynomial)
    out = Tuple{ComplexF64, Any}[]
    for (coef, mono) in _terms(poly)
        iszero(coef) && continue
        canon = _canonical_of(mono)
        push!(out, (ComplexF64(coef), canon))
    end
    return out
end

"""
    discover_pivots(mp; orphans_per_block=32) ->
        (pivots, hpsd_block_sizes, hpsd_block_indices, n_hpsd_blocks, n_aux_blocks)

Build the pivot map.

First pass: scan `:HPSD` block entries, register single-term unit-coef
entries as pivots for the corresponding canonical moment.  Second pass:
for every canonical moment that appears anywhere in the problem (objective,
any `:HPSD` block, any `:Zero` matrix) but lacks a pivot, allocate an
auxiliary 2x2-style Hermitian PSD block so the orphan moment lives at a
clean off-diagonal entry with two free positive slacks anchoring it.

Multiple orphans share an auxiliary block of dimension `orphans_per_block + 1`:
row/column 1 holds a shared positive anchor `s` on the diagonal, the off-
diagonal entries `(1, k+1)` host the orphan moments `y_k`, and the diagonal
positions `(k+1, k+1)` hold per-orphan slacks `t_k`.  Schur complement says
the block is PSD iff every `t_k >= 0` and `s >= sum_k |y_k|^2 / t_k`, so
the orphan moments are effectively free.

Returns:

* `pivots :: Dict{Canon, PhasedPivot}` — includes both HPSD-derived pivots
  and aux-block pivots.
* `block_dims :: Vector{Int}` — dimensions of all blocks (HPSD then aux),
  in libsdp ordering.
* `hpsd_block_indices :: Vector{Int}` — index in `mp.constraints` of each
  HPSD entry (used to map the 2nd pass back to constraint provenance).
* `n_hpsd_blocks :: Int`, `n_aux_blocks :: Int` — split counts.
"""
function discover_pivots(mp; orphans_per_block::Int = 32)
    orphans_per_block >= 1 ||
        throw(ArgumentError("orphans_per_block must be >= 1"))

    pivots = Dict{Any, PhasedPivot}()
    block_dims = Int[]
    hpsd_block_indices = Int[]

    # ---- Pass 1: HPSD single-term unit-coef pivots ----------------------
    for (k, (cone, mat)) in enumerate(mp.constraints)
        cone === :HPSD || continue
        n = size(mat, 1)
        push!(block_dims, n)
        push!(hpsd_block_indices, k)
        b = length(block_dims)

        for j in 1:n, i in 1:n
            poly = mat[i, j]
            iszero(poly) && continue

            n_terms = 0
            cand_canon = nothing
            cand_coef = ComplexF64(0)
            for (coef, mono) in _terms(poly)
                iszero(coef) && continue
                n_terms += 1
                if n_terms == 1
                    cand_canon = _canonical_of(mono)
                    cand_coef = ComplexF64(coef)
                else
                    break
                end
            end
            n_terms == 1 || continue
            _is_unit_coef(cand_coef) || continue

            existing = get(pivots, cand_canon, nothing)
            if existing === nothing ||
                _prefer_real_unit(existing.phase, cand_coef)
                pivots[cand_canon] = PhasedPivot((b, i, j), cand_coef)
            end
        end
    end
    n_hpsd_blocks = length(block_dims)

    # ---- Pass 2: collect orphans in deterministic order -----------------
    seen_orphans = Set{Any}()
    orphan_order = Any[]

    function maybe_record!(canon)
        if !haskey(pivots, canon) && !(canon in seen_orphans)
            push!(seen_orphans, canon)
            push!(orphan_order, canon)
        end
    end

    for (coef, mono) in _terms(mp.objective)
        iszero(coef) && continue
        maybe_record!(_canonical_of(mono))
    end
    for (cone, mat) in mp.constraints
        for entry in mat
            iszero(entry) && continue
            for (coef, mono) in _terms(entry)
                iszero(coef) && continue
                maybe_record!(_canonical_of(mono))
            end
        end
    end

    # ---- Pass 3: allocate auxiliary blocks for orphans ------------------
    n_aux_blocks = 0
    if !isempty(orphan_order)
        idx_in_pack = 0
        b_aux = 0
        for canon in orphan_order
            if idx_in_pack == 0
                # Provisionally allocate a full-size aux block; we trim the
                # last one below.
                push!(block_dims, orphans_per_block + 1)
                n_aux_blocks += 1
                b_aux = length(block_dims)
            end
            idx_in_pack += 1
            pivots[canon] = PhasedPivot((b_aux, 1, 1 + idx_in_pack), ComplexF64(1))
            idx_in_pack >= orphans_per_block && (idx_in_pack = 0)
        end
        # Trim the trailing aux block to its actual content.
        # If last_idx_in_pack is k > 0, the block holds k orphans + 1 anchor.
        last_used = mod(length(orphan_order) - 1, orphans_per_block) + 1
        block_dims[end] = last_used + 1
    end

    return pivots, block_dims, hpsd_block_indices, n_hpsd_blocks, n_aux_blocks
end

# -----------------------------------------------------------------------------
# Polynomial → linear functional in pivot positions
#
# For a polynomial p(y) = Σ c_α y_α, expand to f(X) = Σ c_α X[pivot[α]].
# Each canonical α whose pivot is unknown becomes an "orphan" — collected
# and reported.  In a well-formed PQG MomentProblem there should be no
# orphans because `_validate_polynomial_relaxation_support` guarantees it,
# but we surface any leak to keep the export deterministic.

function _accumulate_lin!(
    lin::Dict{Pivot, ComplexF64},
    pivots::Dict{<:Any, PhasedPivot},
    poly::Polynomial,
    multiplier::ComplexF64,
    orphans::Dict{Any, Int},
)
    for (coef, mono) in _terms(poly)
        iszero(coef) && continue
        canon = _canonical_of(mono)
        ph = get(pivots, canon, nothing)
        if ph === nothing
            orphans[canon] = get(orphans, canon, 0) + 1
            continue
        end
        # X[pivot] = phase · y_α  ⇒  coef · y_α = (coef / phase) · X[pivot]
        contribution = multiplier * ComplexF64(coef) / ph.phase
        lin[ph.pos] = get(lin, ph.pos, ComplexF64(0)) + contribution
    end
    return lin
end

# -----------------------------------------------------------------------------
# Constraint emission
#
# For a complex linear equation Σ_p v_p · X[p] = rhs (rhs ∈ ℂ) split into the
# real-part row and the imag-part row.  Drop a row if all its coefficients are
# numerically zero AND the rhs is zero.

const COEF_TOL = 1e-13

function _emit_complex_equation!(
    entries::Vector{SparseEntry},
    rhs_vec::Vector{Float64},
    constraint_log::Vector{NamedTuple},
    next_k::Ref{Int},
    lin::Dict{Pivot, ComplexF64},
    rhs::ComplexF64,
    provenance::AbstractString,
)
    # Real-part row:  F[p] = conj(v_p),  rhs = Re(rhs).
    # Sort sparse positions so regenerated files do not churn gratuitously.
    lin_pairs = sort!(collect(lin); by = first)
    re_entries = SparseEntry[]
    re_nonzero = false
    for (pos, v) in lin_pairs
        c = conj(v)
        if abs(c) > COEF_TOL
            re_nonzero = true
            push!(re_entries, SparseEntry(0, pos[1], pos[2], pos[3],
                                          real(c), imag(c)))
        end
    end
    if re_nonzero || abs(real(rhs)) > COEF_TOL
        k = next_k[]
        next_k[] += 1
        for e in re_entries
            push!(entries, SparseEntry(k, e.blk, e.r, e.c, e.re, e.im))
        end
        push!(rhs_vec, real(rhs))
        push!(constraint_log, (; k = k, kind = :real, provenance = String(provenance)))
    end

    # Imag-part row:  F[p] = i·conj(v_p),  rhs = Im(rhs).
    im_entries = SparseEntry[]
    im_nonzero = false
    for (pos, v) in lin_pairs
        c = im * conj(v)
        if abs(c) > COEF_TOL
            im_nonzero = true
            push!(im_entries, SparseEntry(0, pos[1], pos[2], pos[3],
                                          real(c), imag(c)))
        end
    end
    if im_nonzero || abs(imag(rhs)) > COEF_TOL
        k = next_k[]
        next_k[] += 1
        for e in im_entries
            push!(entries, SparseEntry(k, e.blk, e.r, e.c, e.re, e.im))
        end
        push!(rhs_vec, imag(rhs))
        push!(constraint_log, (; k = k, kind = :imag, provenance = String(provenance)))
    end
end

# -----------------------------------------------------------------------------
# Top-level export

"""
    export_libsdp(mp; outdir, basename="problem", sense=:min) -> LibsdpExportSummary

Walk a `MomentProblem` and write a libsdp-complex SDPA-sparse file plus a
metadata sidecar.

# Files written

* `<outdir>/<basename>.dat-c`     — the SDPA-sparse problem (see header for format).
* `<outdir>/<basename>_meta.txt`  — pivot map, block layout, sense, constraint
                                    provenance counts.

# Arguments

* `mp::MomentProblem` — built by `moment_relax(...)`.
* `outdir::AbstractString` — directory; created if missing.
* `basename::AbstractString` — file stem; default "problem".
* `sense::Symbol` — `:min` (default) or `:max`.

The objective is `Re ⟨F_0, X⟩`.  The libsdp Python driver should call
`sdp_solver.solve(...)` with `procedure = "minimize"` (for `:min`) or
`procedure = "maximize"` (for `:max`).

# Returns

A `LibsdpExportSummary` describing the exported problem (block count, total
entries, etc.).
"""
function export_libsdp(mp;
                      outdir::AbstractString,
                      basename::AbstractString = "problem",
                      sense::Symbol = :min,
                      orphans_per_block::Int = 32)
    sense in (:min, :max) || throw(ArgumentError("sense must be :min or :max"))
    mkpath(outdir)

    pivots, block_dims, hpsd_block_indices, n_hpsd_blocks, n_aux_blocks =
        discover_pivots(mp; orphans_per_block = orphans_per_block)
    n_blocks = length(block_dims)
    n_aux_pivots = length(pivots) - count(pp -> pp.pos[1] <= n_hpsd_blocks, values(pivots))
    n_hpsd_pivots = length(pivots) - n_aux_pivots

    @info("discovered blocks",
          n_hpsd_blocks,
          n_aux_blocks,
          hpsd_block_dims = block_dims[1:n_hpsd_blocks],
          aux_block_dims_min_max = isempty(block_dims[(n_hpsd_blocks+1):end]) ?
              "—" :
              "$(minimum(block_dims[(n_hpsd_blocks+1):end])) / $(maximum(block_dims[(n_hpsd_blocks+1):end]))",
          n_hpsd_pivots,
          n_aux_pivots)

    # Map :HPSD entry index in mp.constraints back to its 1-indexed block id.
    hpsd_block_id = Dict{Int, Int}()
    for (b, k) in enumerate(hpsd_block_indices)
        hpsd_block_id[k] = b
    end

    # Buffers we will build up.
    entries = SparseEntry[]
    rhs_vec = Float64[]
    constraint_log = NamedTuple[]
    orphans = Dict{Any, Int}()
    next_k = Ref{Int}(1)  # 0 is reserved for the objective

    # ---- Objective F_0 ---------------------------------------------------
    #
    # We want Re ⟨F_0, X⟩ = Re Σ_α h_α y_α = ⟨H⟩ (real, since H is Hermitian).
    # Encode F_0 with entry F_0[pivot[α]] = conj(h_α).
    obj_lin = Dict{Pivot, ComplexF64}()
    obj_orphans = Dict{Any, Int}()
    _accumulate_lin!(obj_lin, pivots, mp.objective, ComplexF64(1), obj_orphans)
    n_objective_entries = 0
    for (pos, v) in sort!(collect(obj_lin); by = first)
        c = conj(v)  # libsdp computes Re(conj(F[α]) · X[α]) = Re(v · X[α]) = Re(v · y_α)
        abs(c) > COEF_TOL || continue
        push!(entries, SparseEntry(0, pos[1], pos[2], pos[3], real(c), imag(c)))
        n_objective_entries += 1
    end
    for (k, v) in obj_orphans
        orphans[k] = get(orphans, k, 0) + v
    end

    # ---- Normalization: y_one = 1 ----------------------------------------
    one_canon = _canonical_of(one(eltype(mp.total_basis)))
    one_ph = get(pivots, one_canon, nothing)
    one_ph === nothing &&
        error("Identity moment has no pivot.  Is the basis missing the identity entry?")
    # Want y_one = 1, i.e. X[pivot] = phase · 1 = phase.
    norm_lin = Dict{Pivot, ComplexF64}(one_ph.pos => ComplexF64(1))
    _emit_complex_equation!(entries, rhs_vec, constraint_log, next_k,
                            norm_lin, one_ph.phase,
                            "normalization (identity moment = 1)")

    # ---- HPSD entry-binding equalities -----------------------------------
    n_zero_components = 0
    for (kc, (cone, mat)) in enumerate(mp.constraints)
        if cone === :HPSD
            b = hpsd_block_id[kc]
            n = size(mat, 1)
            for j in 1:n, i in 1:n
                pos = (b, i, j)
                poly = mat[i, j]
                # Skip if this entry IS the pivot of the (single) canonical
                # mono with coef = 1 — the equation would be trivial 0 = 0.
                if !iszero(poly)
                    canon_terms = _polynomial_terms_with_canon(poly)
                    if length(canon_terms) == 1 && _is_unit_coef(canon_terms[1][1])
                        ph = get(pivots, canon_terms[1][2], nothing)
                        if ph !== nothing && ph.pos === pos &&
                            abs(ph.phase - canon_terms[1][1]) < PIVOT_COEF_TOL
                            continue  # X[pos] := phase · y_α — definition row is implicit.
                        end
                    end
                end

                # Build LHS = X[pos] - poly_in_pivots.
                lin = Dict{Pivot, ComplexF64}(pos => ComplexF64(1))
                _accumulate_lin!(lin, pivots, poly, ComplexF64(-1), orphans)

                _emit_complex_equation!(entries, rhs_vec, constraint_log, next_k,
                                        lin, ComplexF64(0.0, 0.0),
                                        "HPSD-bind block=$(b), entry=($i,$j)")
            end
        elseif cone === :Zero
            n_zero_components += 1
            n = size(mat, 1)
            # Hermitian Zero matrices: only the upper triangle (i ≤ j) is
            # informative.  The (j, i) entry is the conjugate of the (i, j)
            # entry, so encoding both would be redundant.
            for j in 1:n, i in 1:j
                poly = mat[i, j]
                iszero(poly) && continue
                lin = Dict{Pivot, ComplexF64}()
                _accumulate_lin!(lin, pivots, poly, ComplexF64(1), orphans)
                _emit_complex_equation!(entries, rhs_vec, constraint_log, next_k,
                                        lin, ComplexF64(0.0, 0.0),
                                        "Zero-row constraint=$(kc), entry=($i,$j)")
            end
        else
            error("Unsupported cone $(cone) at constraint index $kc")
        end
    end

    n_dropped_orphan_terms = isempty(orphans) ? 0 : sum(values(orphans))
    if !isempty(orphans)
        # We allocated aux blocks for every canonical mono observed during
        # pivot discovery, so any remaining orphans imply the relaxation
        # basis is internally inconsistent (e.g., a polynomial gained terms
        # after pivot discovery).  Surface loudly.
        @warn("$(length(orphans)) canonical moments with no pivot — these terms were dropped (implicit ⟨·⟩ = 0).  Total occurrences: $n_dropped_orphan_terms.  This usually means the polynomial set changed between pivot discovery and constraint emission.")
    end

    n_constraints = length(rhs_vec)
    n_constraint_entries = length(entries) - n_objective_entries

    # ---- Write the .dat-c file -------------------------------------------
    dat_path = joinpath(outdir, basename * ".dat-c")
    open(dat_path, "w") do io
        @printf(io, "* libsdp-complex SDPA-sparse export of NCTSSoS MomentProblem\n")
        @printf(io, "* sense=%s\n", string(sense))
        @printf(io, "* pivots=%d (hpsd=%d, aux=%d), blocks=%d (hpsd=%d, aux=%d)\n",
                length(pivots), n_hpsd_pivots, n_aux_pivots,
                n_blocks, n_hpsd_blocks, n_aux_blocks)
        @printf(io, "* constraints=%d, F_entries=%d (obj=%d)\n",
                n_constraints, length(entries), n_objective_entries)
        @printf(io, "%d\n", n_constraints)
        @printf(io, "%d\n", n_blocks)
        for n in block_dims
            @printf(io, "%d ", n)
        end
        @printf(io, "\n")
        for b in rhs_vec
            @printf(io, "%.17g ", b)
        end
        @printf(io, "\n")
        for e in entries
            @printf(io, "%d %d %d %d %.17g %.17g\n",
                    e.k, e.blk, e.r, e.c, e.re, e.im)
        end
    end

    # ---- Write the metadata sidecar --------------------------------------
    meta_path = joinpath(outdir, basename * "_meta.txt")
    open(meta_path, "w") do io
        @printf(io, "# libsdp-complex SDPA-sparse export metadata\n")
        @printf(io, "sense                       = %s\n", sense)
        @printf(io, "n_blocks                    = %d\n", n_blocks)
        @printf(io, "n_hpsd_blocks               = %d\n", n_hpsd_blocks)
        @printf(io, "n_aux_blocks                = %d\n", n_aux_blocks)
        @printf(io, "hpsd_block_dims             = %s\n",
                join(block_dims[1:n_hpsd_blocks], " "))
        if n_aux_blocks > 0
            aux_dims = block_dims[(n_hpsd_blocks+1):end]
            @printf(io, "aux_block_dim_first_last    = %d / %d\n",
                    first(aux_dims), last(aux_dims))
        end
        @printf(io, "n_constraints               = %d\n", n_constraints)
        @printf(io, "n_objective_entries         = %d\n", n_objective_entries)
        @printf(io, "n_constraint_entries        = %d\n", n_constraint_entries)
        @printf(io, "n_unique_canonical_moments  = %d\n", length(pivots))
        @printf(io, "n_hpsd_pivots               = %d\n", n_hpsd_pivots)
        @printf(io, "n_aux_pivots                = %d\n", n_aux_pivots)
        @printf(io, "n_zero_constraint_matrices  = %d\n", n_zero_components)
        @printf(io, "n_dropped_orphan_terms      = %d\n", n_dropped_orphan_terms)

        kind_counts = Dict{Symbol, Int}()
        for log in constraint_log
            kind_counts[log.kind] = get(kind_counts, log.kind, 0) + 1
        end
        @printf(io, "real_part_rows              = %d\n", get(kind_counts, :real, 0))
        @printf(io, "imag_part_rows              = %d\n", get(kind_counts, :imag, 0))

        provenance_counts = Dict{String, Int}()
        for log in constraint_log
            family = first(split(log.provenance, " "))
            provenance_counts[family] = get(provenance_counts, family, 0) + 1
        end
        @printf(io, "\n# constraint provenance (family => count)\n")
        for (family, cnt) in sort!(collect(provenance_counts); by = first)
            @printf(io, "%-20s %d\n", family, cnt)
        end

        if !isempty(orphans)
            @printf(io, "\n# dropped orphan canonical moments (canon => occurrences)\n")
            shown = 0
            for (canon, cnt) in sort!(collect(orphans); by = x -> string(x[1]))
                @printf(io, "%s  %d\n", string(canon), cnt)
                shown += 1
                shown >= 32 && break
            end
            shown < length(orphans) &&
                @printf(io, "... (%d more elided)\n", length(orphans) - shown)
        end
    end

    return LibsdpExportSummary(n_blocks,
                               n_hpsd_blocks,
                               n_aux_blocks,
                               block_dims,
                               n_constraints,
                               n_objective_entries,
                               n_constraint_entries,
                               length(pivots),
                               n_hpsd_pivots,
                               n_aux_pivots,
                               n_zero_components,
                               n_dropped_orphan_terms,
                               sense)
end

end # module

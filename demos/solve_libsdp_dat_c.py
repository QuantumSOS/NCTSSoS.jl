#!/usr/bin/env python3
"""Read a libsdp-complex SDPA-sparse file (produced by SDPALibsdpExport.jl)
and solve it with libsdp's complex Hermitian BPSDP backend.

The complex BPSDP solver computes the primal "max/min Re(<F0, X>) s.t.
Re(<Fk, X>) = bk, X >= 0 (block-diagonal Hermitian)" exactly as our exporter
emits.  Per libsdp/src/sdp_helper.cc, the only valid `sdp_algorithm` for the
complex solver is "bpsdp" — that is what we want here.

Usage:
    python demos/solve_libsdp_dat_c.py --problem path/to/h4_periodic_moment_sos.dat-c
    python demos/solve_libsdp_dat_c.py --problem path/to/problem.dat-c --maxiter 200000 --tol 1e-7 --verbose

If your environment lacks libsdp, the script will still parse the file and
print summary statistics.  Pass --no-solve to force structure-only mode.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Iterable, Optional


def _parse_dat_c(path: Path):
    """Parse the libsdp-complex SDPA-sparse format.

    Returns:
        b           — list[float] (length m)
        Fi          — list of complex_sdp_matrix-like dicts of length m + 1.
                      Fi[0] is the objective; Fi[1:] are F_1..F_m.
        block_dims  — list[int] (length nblocks)
    """
    with path.open("r") as f:
        # Skip header lines beginning with "*" or quote characters.
        m = nblocks = None
        block_dims = None
        b = None
        for line in f:
            stripped = line.strip()
            if not stripped or stripped[0] in ('*', '"'):
                continue
            m = int(stripped)
            break

        for line in f:
            stripped = line.strip()
            if not stripped or stripped[0] in ('*', '"'):
                continue
            nblocks = int(stripped)
            break

        for line in f:
            stripped = line.strip()
            if not stripped or stripped[0] in ('*', '"'):
                continue
            block_dims = [abs(int(x)) for x in stripped.split()]
            break

        for line in f:
            stripped = line.strip()
            if not stripped or stripped[0] in ('*', '"'):
                continue
            b = [float(x) for x in stripped.split()]
            break

        if m is None or nblocks is None or block_dims is None or b is None:
            raise ValueError(f"Malformed header in {path}")

        if len(block_dims) != nblocks:
            raise ValueError(f"block dims length {len(block_dims)} != nblocks {nblocks}")
        if len(b) != m:
            raise ValueError(f"b length {len(b)} != m {m}")

        # Collect entries grouped by k.
        entries_by_k = [[] for _ in range(m + 1)]
        for line in f:
            stripped = line.strip()
            if not stripped or stripped[0] in ('*', '"'):
                continue
            tok = stripped.split()
            if len(tok) < 6:
                raise ValueError(f"Bad entry line: {stripped!r}")
            k, blk, r, c = (int(x) for x in tok[:4])
            re = float(tok[4])
            im = float(tok[5])
            if not (0 <= k <= m):
                raise ValueError(f"Constraint id {k} out of range [0, {m}]")
            entries_by_k[k].append((blk, r, c, complex(re, im)))

    return b, entries_by_k, block_dims


def _summarize(b, entries_by_k, block_dims, sense=None):
    print("== Parsed problem ==")
    print(f"  blocks                 : {len(block_dims)}")
    print(f"    sizes (min/max/sum)  : {min(block_dims)} / {max(block_dims)} "
          f"/ {sum(block_dims)}")
    print(f"    sum of n_b^2 (primal): {sum(d * d for d in block_dims)}")
    print(f"  m (constraints)        : {len(b)}")
    print(f"  objective entries      : {len(entries_by_k[0])}")
    constraint_entries = sum(len(e) for e in entries_by_k[1:])
    print(f"  constraint entries     : {constraint_entries}")
    if sense is not None:
        print(f"  sense                  : {sense}")

    # Constraint-RHS distribution (most are zero; spot-check non-zero ones).
    nonzero_b = [(k, bk) for k, bk in enumerate(b, start=1) if abs(bk) > 1e-15]
    print(f"  constraints with nonzero RHS: {len(nonzero_b)}")
    if nonzero_b[:8]:
        for k, bk in nonzero_b[:8]:
            print(f"    k={k:6d}  b={bk:+.6g}  entries={len(entries_by_k[k])}")
        if len(nonzero_b) > 8:
            print(f"    ... ({len(nonzero_b) - 8} more)")

    # Empty constraints (0 = 0) are libsdp-fatal; flag them up front.
    empties = [k for k, e in enumerate(entries_by_k[1:], start=1)
               if not e and abs(b[k - 1]) < 1e-15]
    print(f"  trivially empty constraints (0 = 0): {len(empties)}")


def _build_libsdp_inputs(b, entries_by_k, block_dims):
    """Build the (b, Fi, dimensions) triple for libsdp's complex_sdp_solver."""
    from libsdp.sdp_helper import complex_sdp_matrix

    Fi = []
    for k, entries in enumerate(entries_by_k):
        Fk = complex_sdp_matrix()
        for (blk, r, c, val) in entries:
            Fk.block_number.append(blk)
            Fk.row.append(r)
            Fk.column.append(c)
            Fk.value.append(val)
        Fi.append(Fk)

    return list(b), Fi, list(block_dims)


def _solve(b, entries_by_k, block_dims, *, maxiter, tol, verbose, procedure):
    import numpy as np
    from libsdp import sdp_options
    from libsdp.sdp_helper import complex_sdp_solver
    from libsdp.sdpa_file_io import clean_complex_sdpa_problem

    b_clean, Fi_clean = clean_complex_sdpa_problem(b, _build_libsdp_inputs(
        b, entries_by_k, block_dims)[1])

    options = sdp_options()
    options.sdp_algorithm = "bpsdp"
    options.procedure = procedure  # "minimize" or "maximize"
    options.guess_type = "zero"
    options.sdp_error_convergence = tol
    options.sdp_objective_convergence = tol
    options.cg_convergence = max(1e-12, tol * 1e-3)
    options.print_level = 1 if verbose else 0

    print("\n== Calling libsdp complex BPSDP ==")
    print(f"  procedure              : {procedure}")
    print(f"  sdp_error_convergence  : {tol:.1e}")
    print(f"  cg_convergence         : {options.cg_convergence:.1e}")
    print(f"  maxiter                : {maxiter}")

    sdp = complex_sdp_solver(options, Fi_clean, list(block_dims))
    t0 = time.time()
    x = sdp.solve(b_clean, maxiter)
    t1 = time.time()
    x = np.array(x)

    c = np.array(sdp.get_c())
    Ax = np.array(sdp.get_Au(x))
    primal_obj = float(np.real(np.dot(np.conj(c), x)))
    primal_err = float(np.linalg.norm(np.asarray(Ax) - np.asarray(b_clean)))

    print("\n== libsdp BPSDP results ==")
    print(f"  walltime               : {t1 - t0:.2f} s")
    print(f"  primal objective       : {primal_obj:+.12f}")
    print(f"  ||Ax - b||             : {primal_err:.3e}")

    try:
        y = np.array(sdp.get_y())
        dual_obj = float(np.dot(b_clean, y))
        print(f"  dual objective         : {dual_obj:+.12f}")
        print(f"  duality gap            : {abs(primal_obj - dual_obj):.3e}")
    except Exception:
        pass

    return x


def _detect_sense(path: Path) -> str:
    """Best-effort: read the header comments to detect sense=min|max."""
    try:
        with path.open("r") as f:
            for _ in range(8):
                line = f.readline()
                if not line:
                    break
                if "sense=" in line:
                    return line.split("sense=", 1)[1].strip().split()[0]
    except OSError:
        pass
    return "min"


def main(argv: Optional[Iterable[str]] = None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--problem", required=True, type=Path,
                        help="Path to the .dat-c file produced by SDPALibsdpExport.jl.")
    parser.add_argument("--maxiter", type=int, default=200_000,
                        help="Maximum BPSDP iterations (default: 200000).")
    parser.add_argument("--tol", type=float, default=1e-7,
                        help="sdp_error_convergence and sdp_objective_convergence (default: 1e-7).")
    parser.add_argument("--procedure", choices=("minimize", "maximize"), default=None,
                        help="Override the sense detected from the file header.")
    parser.add_argument("--verbose", action="store_true",
                        help="Pass print_level=1 to libsdp.")
    parser.add_argument("--no-solve", action="store_true",
                        help="Skip the libsdp call; only parse and print stats.")
    args = parser.parse_args(argv)

    if not args.problem.exists():
        print(f"error: {args.problem} not found", file=sys.stderr)
        return 2

    sense = args.procedure or ("minimize" if _detect_sense(args.problem) == "min" else "maximize")

    print(f"Reading {args.problem} ...")
    t0 = time.time()
    b, entries_by_k, block_dims = _parse_dat_c(args.problem)
    print(f"  parse walltime         : {time.time() - t0:.2f} s")

    _summarize(b, entries_by_k, block_dims, sense=sense)

    if args.no_solve:
        return 0

    try:
        import libsdp  # noqa: F401
    except Exception as exc:
        print(f"\nlibsdp not importable ({exc}); rerun with --no-solve to skip the solve call.",
              file=sys.stderr)
        return 3

    _solve(b, entries_by_k, block_dims,
           maxiter=args.maxiter, tol=args.tol,
           verbose=args.verbose, procedure=sense)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

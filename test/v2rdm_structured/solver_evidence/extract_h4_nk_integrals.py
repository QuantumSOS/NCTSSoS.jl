#!/usr/bin/env python3
"""Extract physical periodic H4 active-space integrals for arbitrary Nk.

This is the PySCF side of the H4/Nk sweep.  It uses the Figure-1-matching
1D-chain setup used by the earlier Nk=2 run:

  * H4 unit cell at x = 0,1,2,3 Å, lattice a = 4 Å
  * gth-tzvp basis, no explicit pseudopotential
  * cell.dimension = 1, low_dim_ft_type = inf_vacuum
  * KRHF with GDF density fitting (PySCF RSDF rejects 1D low-dimensional cells)
  * first 8 spatial MOs per k-point as the active space

Outputs are intentionally boring:

  <outdir>/h4_chain_nk<N>_integrals_drop<TOL>.txt
      Unnormalized h1e/eri text dump consumed by the Julia NCTSSoS exporter.
  <outdir>/h4_chain_nk<N>.npz
      Full NumPy archive for audit/reuse.
  <outdir>/h4_chain_nk<N>_meta.json
      HF, active-HF, and energy shift metadata.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


def build_cell(args: argparse.Namespace):
    from pyscf.pbc import gto

    cell = gto.Cell()
    cell.atom = """
H  0.0  0.0  0.0
H  1.0  0.0  0.0
H  2.0  0.0  0.0
H  3.0  0.0  0.0
"""
    cell.a = [[4.0, 0.0, 0.0], [0.0, args.vacuum, 0.0], [0.0, 0.0, args.vacuum]]
    cell.basis = args.basis
    if args.pseudo.lower() != "none":
        cell.pseudo = args.pseudo
    cell.dimension = args.dimension
    if args.low_dim_ft_type.lower() != "none":
        cell.low_dim_ft_type = args.low_dim_ft_type
    cell.unit = "angstrom"
    cell.verbose = args.verbose
    cell.build()
    return cell


def find_conserving_k4(kpts: np.ndarray, lattice_vectors: np.ndarray, ik1: int, ik2: int, ik3: int) -> int | None:
    target = kpts[ik1] + kpts[ik2] - kpts[ik3]
    for ik4, k4 in enumerate(kpts):
        diff = target - k4
        frac = np.dot(diff, lattice_vectors.T) / (2.0 * np.pi)
        if np.allclose(frac, np.round(frac), atol=1e-6):
            return ik4
    return None


def active_hf_energy(h1e: dict[int, np.ndarray], eri: dict[tuple[int, int, int, int], np.ndarray], *, nk: int, n_active_elec: int) -> float:
    nocc = n_active_elec // 2

    e1 = 0.0 + 0.0j
    for k in range(nk):
        for i in range(nocc):
            e1 += 2.0 * h1e[k][i, i] / nk

    e2 = 0.0 + 0.0j
    for (k1, k2, k3, k4), raw in eri.items():
        v = raw / (nk * nk)
        for i in range(nocc):
            for j in range(nocc):
                if k1 == k3 and k2 == k4:
                    e2 += 2.0 * v[i, i, j, j]
                if k1 == k4 and k2 == k3:
                    e2 -= v[i, j, j, i]
    return float(np.real(e1 + e2))


def drop_label(drop_tol: float) -> str:
    if drop_tol == 0:
        return "0"
    exponent = int(round(math.log10(drop_tol)))
    if math.isclose(drop_tol, 10.0**exponent):
        return f"1e{exponent}"
    return f"{drop_tol:g}".replace("+", "")


def write_text_integrals(path: Path, h1e: dict[int, np.ndarray], eri: dict[tuple[int, int, int, int], np.ndarray], meta: dict[str, Any], drop_tol: float) -> tuple[int, int]:
    n_h1e = 0
    n_eri = 0
    with path.open("w", encoding="utf-8") as f:
        f.write("# Physical periodic H4 active-space integrals from PySCF\n")
        f.write("# format: h1e k p q Re Im; eri k1 k2 k3 k4 p r q s Re Im\n")
        f.write("# entries are unnormalized; Julia exporter divides by Nk and Nk^2\n")
        f.write(f"# metadata: {json.dumps(meta, sort_keys=True)}\n\n")

        nk = int(meta["nk"])
        norb = int(meta["n_active_orb"])
        for k in range(nk):
            block = h1e[k]
            for p in range(norb):
                for q in range(norb):
                    value = complex(block[p, q])
                    if abs(value) <= drop_tol:
                        continue
                    f.write(f"h1e {k:2d} {p:4d} {q:4d} {value.real: .12e} {value.imag: .12e}\n")
                    n_h1e += 1

        for key in sorted(eri):
            block = eri[key]
            k1, k2, k3, k4 = key
            for p in range(norb):
                for r in range(norb):
                    for q in range(norb):
                        for s in range(norb):
                            value = complex(block[p, r, q, s])
                            if abs(value) <= drop_tol:
                                continue
                            f.write(
                                f"eri {k1:2d} {k2:2d} {k3:2d} {k4:2d} "
                                f"{p:4d} {r:4d} {q:4d} {s:4d} "
                                f"{value.real: .12e} {value.imag: .12e}\n"
                            )
                            n_eri += 1
    return n_h1e, n_eri


def jsonable(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def extract(args: argparse.Namespace) -> None:
    from pyscf.pbc import df, scf

    t_all = time.time()
    cell = build_cell(args)
    kpts = cell.make_kpts([args.nk, 1, 1])

    print("== physical H4 periodic integral extraction ==")
    print(
        f"nk={args.nk} active_orb={args.n_active_orb} active_elec={args.n_active_elec} "
        f"basis={args.basis} pseudo={args.pseudo} vacuum={args.vacuum:g}A "
        f"dimension={cell.dimension} low_dim_ft_type={getattr(cell, 'low_dim_ft_type', None)} "
        f"density_fit={args.density_fit}"
    )
    print(f"AOs/cell={cell.nao_nr()} electrons/cell={cell.nelectron}")
    for ik, kpt in enumerate(kpts):
        print(f"k[{ik}] = {kpt}")

    mf = scf.KRHF(cell, kpts)
    if args.density_fit == "gdf":
        mf.with_df = df.GDF(cell, kpts)
    elif args.density_fit == "rsdf":
        mf.with_df = df.RSDF(cell, kpts)
    else:  # argparse should catch this
        raise ValueError(args.density_fit)
    mf.conv_tol = args.scf_tol

    t0 = time.time()
    ehf = float(mf.kernel())
    print(f"HF/cell = {ehf: .12f} Ha ({time.time() - t0:.1f}s)")

    norb = args.n_active_orb
    for ik, coeff in enumerate(mf.mo_coeff):
        if coeff.shape[1] < norb:
            raise ValueError(f"k[{ik}] has only {coeff.shape[1]} MOs; need {norb}")
        print(f"MO energies k[{ik}] first {norb + 2}: {mf.mo_energy[ik][:norb + 2]}")

    hcore_ao = mf.get_hcore()
    h1e: dict[int, np.ndarray] = {}
    for ik in range(args.nk):
        c = mf.mo_coeff[ik][:, :norb]
        h1e[ik] = c.conj().T @ hcore_ao[ik] @ c

    lat_vecs = cell.lattice_vectors()
    eri: dict[tuple[int, int, int, int], np.ndarray] = {}
    t0 = time.time()
    for ik1 in range(args.nk):
        for ik2 in range(args.nk):
            for ik3 in range(args.nk):
                ik4 = find_conserving_k4(kpts, lat_vecs, ik1, ik2, ik3)
                if ik4 is None:
                    continue
                c1 = mf.mo_coeff[ik1][:, :norb]
                c2 = mf.mo_coeff[ik2][:, :norb]
                c3 = mf.mo_coeff[ik3][:, :norb]
                c4 = mf.mo_coeff[ik4][:, :norb]
                # PySCF returns chemist order for this argument order:
                # (p k1, r k3 | q k2, s k4), stored [p,r,q,s].
                block = mf.with_df.ao2mo(
                    [c1, c3, c2, c4],
                    [kpts[ik1], kpts[ik3], kpts[ik2], kpts[ik4]],
                    compact=False,
                ).reshape(norb, norb, norb, norb)
                eri[(ik1, ik2, ik3, ik4)] = block
                print(f"eri {(ik1, ik2, ik3, ik4)} max|V|={np.max(np.abs(block)):.6e}")
    print(f"ERI extraction walltime = {time.time() - t0:.1f}s")

    active_hf = active_hf_energy(h1e, eri, nk=args.nk, n_active_elec=args.n_active_elec)
    energy_shift = ehf - active_hf
    print(f"HF active/cell      = {active_hf: .12f} Ha")
    print(f"constant shift/cell = {energy_shift: .12f} Ha")

    args.outdir.mkdir(parents=True, exist_ok=True)
    stem = f"h4_chain_nk{args.nk}"
    npz_path = args.outdir / f"{stem}.npz"
    txt_path = args.outdir / f"{stem}_integrals_drop{drop_label(args.drop_tol)}.txt"
    meta_path = args.outdir / f"{stem}_meta.json"

    save: dict[str, Any] = {
        "nk": args.nk,
        "n_active_orb": norb,
        "n_active_elec": args.n_active_elec,
        "kpoints": kpts,
        "ehf": ehf,
        "active_hf": active_hf,
        "energy_shift": energy_shift,
        "basis": args.basis,
        "pseudo": args.pseudo,
        "vacuum_angstrom": args.vacuum,
        "dimension": cell.dimension,
        "low_dim_ft_type": getattr(cell, "low_dim_ft_type", None) or "none",
        "density_fit": args.density_fit,
        "scf_tol": args.scf_tol,
        "drop_tol": args.drop_tol,
    }
    for ik in range(args.nk):
        save[f"h1e_k{ik}"] = h1e[ik]
        save[f"h1e_norm_k{ik}"] = h1e[ik] / args.nk
        save[f"mo_energy_k{ik}"] = mf.mo_energy[ik]
        save[f"mo_coeff_k{ik}"] = mf.mo_coeff[ik]
    for key, block in eri.items():
        k1, k2, k3, k4 = key
        save[f"eri_{k1}_{k2}_{k3}_{k4}"] = block
        save[f"eri_norm_{k1}_{k2}_{k3}_{k4}"] = block / (args.nk * args.nk)
    np.savez(npz_path, **save)

    meta = {key: jsonable(value) for key, value in save.items() if not key.startswith(("h1e", "eri", "mo_"))}
    meta["script"] = Path(__file__).name
    meta["wall_seconds"] = time.time() - t_all
    n_h1e, n_eri = write_text_integrals(txt_path, h1e, eri, meta, args.drop_tol)
    meta["text_h1e_entries"] = n_h1e
    meta["text_eri_entries"] = n_eri
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"wrote npz  {npz_path}")
    print(f"wrote text {txt_path} ({n_h1e} h1e, {n_eri} eri entries)")
    print(f"wrote meta {meta_path}")
    print(f"total walltime = {time.time() - t_all:.1f}s")


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nk", type=int, required=True)
    p.add_argument("--n-active-orb", type=int, default=8)
    p.add_argument("--n-active-elec", type=int, default=4)
    p.add_argument("--basis", default="gth-tzvp")
    p.add_argument("--pseudo", default="none")
    p.add_argument("--vacuum", type=float, default=10.0)
    p.add_argument("--dimension", type=int, default=1)
    p.add_argument("--low-dim-ft-type", default="inf_vacuum")
    p.add_argument("--density-fit", choices=("gdf", "rsdf"), default="gdf")
    p.add_argument("--scf-tol", type=float, default=1e-10)
    p.add_argument("--drop-tol", type=float, default=1e-12)
    p.add_argument("--verbose", type=int, default=3)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args(argv)
    if args.nk < 1:
        p.error("--nk must be positive")
    if args.n_active_orb < 1:
        p.error("--n-active-orb must be positive")
    if args.n_active_elec < 1:
        p.error("--n-active-elec must be positive")
    return args


def main(argv: list[str] | None = None) -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    extract(parse_args(sys.argv[1:] if argv is None else argv))


if __name__ == "__main__":
    main()

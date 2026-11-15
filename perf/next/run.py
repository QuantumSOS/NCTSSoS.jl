#!/usr/bin/env python3
"""One bounded Julia process per row; stdlib only. Input checkout is read-only."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time
import tomllib

ROOT = Path(__file__).resolve().parents[2]


def integral_input(checkout, n, output):
    source = checkout / f"problems/h{n}_sto3g/integrals/h{n}_R1.0_sto3g.json"
    data = json.loads(source.read_text())
    path = output / f"h{n}.toml"
    # JSON arrays/numbers are valid TOML values; retain only the numeric inputs.
    keys = ("h_ij", "eri", "nuclear_repulsion", "n_spatial_orbitals", "n_electrons", "n_alpha", "n_beta")
    path.write_text("\n".join(f"{key} = {json.dumps(data[key])}" for key in keys))
    return path


def run(args, row):
    target, route, n, order, variant = row
    if args.save_program and route != "next":
        raise ValueError("--save-program requires route next")
    name = "-".join(map(str, row))
    if args.mode != "build":
        name += "-" + args.mode
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    path = output / f"{name}.toml"
    if path.exists() and not args.replace:
        raise FileExistsError(f"{path}: use --replace to rerun explicitly")
    path.unlink(missing_ok=True)
    data = integral_input(args.checkout, n, output) if target == "hydrogen" else "-"
    cmd = ["julia", f"--project={args.project}", "--startup-file=no", str(ROOT / "perf/next/probe.jl"),
           target, route, str(n), str(order), variant, str(data), str(path), args.mode]
    if args.save_program:
        cmd.append("save_program")
    env = dict(os.environ, JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT / "perf/next").glob("*.jl")}
    begin = time.monotonic()
    peak = 0
    reason = "exit"
    print(f"START {name}: estimate 1–10 min; hard limits {args.seconds}s / {args.gb} GB", flush=True)
    with (output / f"{name}.log").open("w") as log:
        log.write(f"COMMAND {' '.join(cmd)}\n"); log.flush()
        proc = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        while proc.poll() is None:
            try:
                status = Path(f"/proc/{proc.pid}/status").read_text()
                rss = int(next(line.split()[1] for line in status.splitlines() if line.startswith("VmRSS:"))) * 1024
                peak = max(peak, rss)
            except (FileNotFoundError, StopIteration):
                pass
            if peak > args.gb * 10**9:
                reason = "memory_limit"
            elif time.monotonic() - begin > args.seconds:
                reason = "time_limit"
            if reason != "exit":
                proc.kill(); proc.wait(); break
            time.sleep(0.1)
    elapsed = time.monotonic() - begin
    if reason == "exit" and proc.returncode:
        reason = "error"
    result = tomllib.loads(path.read_text()) if path.exists() else {}
    result.update(exit_code=proc.returncode, stop_reason=reason, wall_seconds=elapsed,
                  sampled_peak_rss_bytes=peak, command=cmd, git_revision=revision, worker_sha256=hashes)
    peak = max(peak, result.get("peak_rss_bytes", 0))
    result["peak_rss_bytes"] = peak
    (output / f"{name}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"END {name}: {elapsed:.2f}s; {reason}; stage={result.get('stage')}; peak={peak / 1e9:.3f} GB", flush=True)
    print("| case | warm build s | bytes | peak GB | rows | moments | vars | blocks | side | lower s | lower rows/cols/nnz |")
    median = lambda key: f"{statistics.median(result[key]):.6g}" if result.get(key) else "—"
    print(f"| {name} | {median('build_seconds')} | {median('build_bytes')} | {peak/1e9:.3f} | "
          + " | ".join(str(result.get(k, "—")) for k in ("rows", "moments", "variables", "blocks", "side"))
          + f" | {median('lower_seconds')} | " + "/".join(str(result.get(k, "—")) for k in ("lower_rows", "lower_cols", "nonzeros")) + " |", flush=True)
    return 1 if reason == "error" else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("chain", "torus", "hydrogen", "bell"))
    parser.add_argument("route", choices=("next", "legacy"))
    parser.add_argument("n", type=int)
    parser.add_argument("order", type=int)
    parser.add_argument("variant")
    parser.add_argument("--project", default=".")
    parser.add_argument("--checkout", type=Path, default=ROOT.parent / "repos/NCTSSOSBenchmark")
    parser.add_argument("--output", type=Path, default=ROOT / ".amp/in/orb-backup/scale")
    parser.add_argument("--seconds", type=float, default=600)
    parser.add_argument("--gb", type=float, default=20)
    parser.add_argument("--mode", choices=("build", "profile", "solve"), default="build")
    parser.add_argument("--save-program", action="store_true",
                        help="Serialize the first Next program outside the build timer for equivalence checks")
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()
    raise SystemExit(run(args, (args.target, args.route, args.n, args.order, args.variant)))

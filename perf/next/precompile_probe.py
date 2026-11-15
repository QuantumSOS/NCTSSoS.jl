#!/usr/bin/env python3
"""Measure isolated workload blocks in a disposable source copy, never this checkout.

The environment must develop --source and contain SCS and Mosek. A private depot retains
dependency caches; only this package's and its solver extensions' images are removed
before each repetition. Initial dependency compilation is setup, not workload cost.
P1/P2/... select one package block, E1/E2/... one SCS block, M1 a Mosek block; none and all
provide the controls. Isolated savings overlap and must not be added.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
BLOCK = re.compile(r"^@compile_workload let\n.*?^end\n", re.M | re.S)
CASES = ("xxz6", "xxz12front", "hubbard", "ball", "CHSH")
EXTRA_CASES = ("heisenberg8continuous", "fermionic4number", "heisenberg8detect", "xxz6density")


def execute(command, log, env):
    start = time.monotonic()
    with log.open("w") as stream:
        stream.write("COMMAND " + " ".join(map(str, command)) + "\n")
        stream.flush()
        subprocess.run(command, env=env, cwd=ROOT, stdout=stream,
                       stderr=subprocess.STDOUT, check=True, timeout=600)
    return time.monotonic() - start


def main(args):
    source = args.source.resolve()
    if source == ROOT or not source.is_relative_to(ROOT / ".amp/in/orb-backup"):
        raise ValueError("Use a disposable source copy under .amp/in/orb-backup, not the checkout")
    paths = {"P": source / "src/next/precompile.jl", "E": source / "ext/NCTSSoSSCSExt.jl",
             "M": source / "ext/NCTSSoSMosekExt.jl"}
    if args.backend == "mosek" and any(case not in EXTRA_CASES for case in args.cases):
        raise ValueError("Mosek probes require the extra cases")
    probes = [(case, "scs" if args.backend == "both" else args.backend, case) for case in args.cases]
    if args.backend == "both":
        probes += [(case, "mosek", "mosek:" + case) for case in args.cases if case in EXTRA_CASES]
    originals = {key: path.read_text() for key, path in paths.items()}
    blocks = {key: list(BLOCK.finditer(text)) for key, text in originals.items()}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    depot = (args.depot or output / "depot").resolve()
    if not depot.is_relative_to(ROOT / ".amp/in/orb-backup"):
        raise ValueError("Use a private benchmark depot under .amp/in/orb-backup")
    depot.mkdir(parents=True, exist_ok=True)
    probe = ROOT / "perf/next/first_call.jl"
    extra_probe = ROOT / "perf/next/first_call_extra.jl"
    fingerprints = {key: hashlib.sha256(text.encode()).hexdigest() for key, text in originals.items()}
    fingerprints["probe"] = hashlib.sha256(probe.read_bytes()).hexdigest()
    if any(case in EXTRA_CASES for case in args.cases):
        fingerprints["extra_probe"] = hashlib.sha256(extra_probe.read_bytes()).hexdigest()
    records = []
    try:
        for variant in args.variants:
            if variant not in ("none", "all"):
                if variant[0] not in paths or not 1 <= int(variant[1:]) <= len(blocks[variant[0]]):
                    raise ValueError(f"No workload {variant}")
            for key, path in paths.items():
                counter = iter(range(1, len(blocks[key]) + 1))
                path.write_text(BLOCK.sub(lambda m: m[0] if variant in ("all", key + str(next(counter))) else "", originals[key]))
            for rep in range(1, args.replicates + 1):
                name = f"{variant}-{rep}"
                result_path = output / (name + ".json")
                if result_path.exists():
                    result = json.loads(result_path.read_text())
                    assert result["source_hashes"] == fingerprints
                    assert result["backend"] == args.backend
                    assert set(result["cases"]) == {key for _, _, key in probes}
                    records.append(result)
                    continue
                for module in ("NCTSSoS", "NCTSSoSSCSExt", "NCTSSoSMosekExt"):
                    for cache in depot.glob("compiled/v*/" + module):
                        shutil.rmtree(cache)
                env = dict(os.environ, JULIA_DEPOT_PATH=f"{depot}:{Path.home() / '.julia'}",
                           JULIA_NUM_THREADS="1", JULIA_NUM_PRECOMPILE_TASKS="1",
                           JULIA_PKG_PRECOMPILE_AUTO="0", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
                julia = ["julia", f"--project={args.environment.resolve()}", "--startup-file=no"]
                print(f"START {name}: estimate 2–5 min including {len(probes)} fresh-process probes", flush=True)
                begin = time.monotonic()
                log = output / (name + "-precompile.log")
                compile_code = ('using Pkg; '
                                '@assert realpath(Pkg.dependencies()[Pkg.project().dependencies["NCTSSoS"]].source) == realpath(ARGS[1]); '
                                'Pkg.precompile(["NCTSSoS", "SCS", "Mosek"])')
                elapsed = execute(julia + ["-e", compile_code, str(source)], log, env)
                times = {}
                for line in log.read_text().splitlines():
                    match = re.search(r"([\d.]+) ms\s+✓ NCTSSoS(.*)$", line)
                    if match:
                        key = ("extension" if "NCTSSoSSCSExt" in match[2] else
                               "mosek_extension" if "NCTSSoSMosekExt" in match[2] else "package")
                        times[key] = float(match[1]) / 1000
                if not {"package", "extension"}.issubset(times):
                    raise RuntimeError(f"Expected both compile timings in {log}; found {times}")
                result = dict(variant=variant, repetition=rep, backend=args.backend, source_hashes=fingerprints,
                              precompile_seconds=times, precompile_process_seconds=elapsed, cases={})
                for case, backend, key in probes:
                    log = output / (name + "-" + key.replace(":", "-") + ".log")
                    script = extra_probe if case in EXTRA_CASES else probe
                    backend_args = [backend] if case in EXTRA_CASES else []
                    wall = execute(julia + [str(script), case] + backend_args, log, env)
                    if "Precompiling" in log.read_text():
                        raise RuntimeError(f"Unexpected precompilation in first-call probe {log}")
                    line = next(line for line in log.read_text().splitlines() if line.startswith("PROBE "))
                    values = dict(field.split("=", 1) for field in line.split()[1:])
                    result["cases"][key] = {field: float(value) if field not in ("case", "status") else value
                                             for field, value in values.items()}
                    result["cases"][key]["process_seconds"] = wall
                    if values["status"] != "OPTIMAL":
                        raise RuntimeError(f"Probe did not reach OPTIMAL: {log}")
                result["wall_seconds"] = time.monotonic() - begin
                result_path.write_text(json.dumps(result, indent=2) + "\n")
                records.append(result)
                print(f"END {name}: {result['wall_seconds']:.2f}s; package {times['package']:.3f}s; extension {times['extension']:.3f}s", flush=True)
    finally:
        for key, path in paths.items():
            path.write_text(originals[key])
    print("| variant | repetition | package s | extension s | " + " | ".join(args.cases) + " |")
    print("|---|---:|---:|---:|" + "---:|" * len(args.cases))
    for result in records:
        times = result["precompile_seconds"]
        print(f"| {result['variant']} | {result['repetition']} | {times['package']:.3f} | {times['extension']:.3f} | "
              + " | ".join(str(result['cases'][case]['first_call']) for case in args.cases) + " |")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--environment", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--depot", type=Path)
    parser.add_argument("--variants", nargs="+", default=["none", "P1", "P2", "P3", "E1", "E2", "E3", "all"])
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--cases", nargs="+", default=list(CASES))
    parser.add_argument("--backend", choices=("scs", "mosek", "both"), default="scs")
    main(parser.parse_args())

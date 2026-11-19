#!/usr/bin/env python3
"""Serial, fresh-process Heisenberg benchmark runner (Linux; Python 3.11+)."""

import argparse
import dataclasses
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import time
import tomllib

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BASE = "e59863669f0fb30497c84f219240c1cf1bf38973"
REF = {"chain": "0b152dfcd9bcde3403aa4d967642ce5261c80712",
       "square": "9a3ae67ff2bb6da1c04b72318ecbfcec71db5495"}
GIB = 1 << 30
ENERGY_ATOL = 1e-6


@dataclasses.dataclass(frozen=True)
class Case:
    lattice: str
    extent: int
    reference_extra: int
    exact_control: bool = False
    strengthened: bool = False

    @property
    def id(self):
        return f"{self.lattice}-{self.extent}" + ("-rdm" if self.strengthened else "")


CASES = [Case("chain", 4, 1, True), Case("chain", 6, 2, True),
         Case("chain", 8, 3, True), Case("chain", 8, 3, True, True),
         Case("square", 4, 0, True), Case("square", 4, 0, True, True),
         *(Case("chain", L, 3, False, strengthened) for L in (60, 80) for strengthened in (False, True)),
         Case("square", 6, 0), Case("square", 6, 0, False, True)]


def run_process(cmd, log, timeout, env):
    # Wait on Julia itself, not a shell or a pipeline; Linux wait4 ru_maxrss is KiB.
    # Use a separate group so a timeout cannot leave solver descendants behind.
    start = time.monotonic()
    expired = False
    with log.open("wb") as output:
        proc = subprocess.Popen(cmd, env=env, stdout=output, stderr=subprocess.STDOUT,
                                start_new_session=True)
        try:
            while True:
                pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
                if pid:
                    break
                if time.monotonic() - start >= timeout:
                    expired = True
                    os.killpg(proc.pid, signal.SIGKILL)
                    _, status, usage = os.wait4(proc.pid, 0)
                    break
                time.sleep(0.02)
        except BaseException:
            os.killpg(proc.pid, signal.SIGKILL)
            _, status, _ = os.wait4(proc.pid, 0)
            raise
        finally:
            proc.returncode = os.waitstatus_to_exitcode(status)
    return {"wall_seconds": time.monotonic() - start, "peak_rss_kib": usage.ru_maxrss,
            "exit_code": proc.returncode, "timed_out": expired}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def cgroup_guard():
    row = next(line for line in Path("/proc/self/cgroup").read_text().splitlines()
               if line.startswith("0::"))
    group = Path("/sys/fs/cgroup") / row[3:].lstrip("/")
    cap = (group / "memory.max").read_text().strip()
    swap = (group / "memory.swap.max").read_text().strip()
    if cap == "max" or int(cap) > 100 * GIB or swap != "0":
        raise SystemExit("Use a dedicated cgroup v2 with MemoryMax ≤ 100 GiB and MemorySwapMax=0")
    return group, int(cap)


def oom_kills(group):
    return int(dict(line.split() for line in (group / "memory.events").read_text().splitlines())["oom_kill"])


def stock_diagnostics(row, log):
    # GSB does not return its model. On this pin it prints a termination
    # line for EVERY nonoptimal result, then prints the returned objective.
    # Read only the native final iteration, never an earlier good iterate.
    row["status_source"] = "stock GSB output and Mosek final iteration"
    for key, pattern in {
        "basis_seconds": r"Obtained the monomial basis in (\S+) seconds",
        "block_seconds": r"Finished block-diagonalization in (\S+) seconds",
        "lso_seconds": r"Added linear state optimality constraints in (\S+) seconds",
        "pso_seconds": r"Added PSD state optimality constraints in (\S+) seconds",
        "rdm_seconds": r"Added the rdm constraint in (\S+) seconds",
        "solve_seconds": r"SDP solving time: (\S+) seconds",
    }.items():
        match = re.search(pattern, log)
        if match:
            value = float(match[1])
            if math.isfinite(value):
                row[f"stock_{key}"] = value
    status = re.search(r"^termination status: (\S+)$", log, re.M)
    objective = re.findall(r"^optimum = (\S+)$", log, re.M)
    if status or len(objective) != 1 or float(objective[0]) != row.get("energy_per_site"):
        if status:
            row["termination_status"] = status[1]
        return
    row["termination_status"] = "OPTIMAL"
    table = re.search(r"^ITE\s+PFEAS\s+DFEAS\s+GFEAS\s+PRSTATUS\s+POBJ\s+DOBJ\s+MU\s+TIME[^\n]*\n"
                      r"((?:^\d+[^\n]*\n)+)^Optimizer terminated\.", log, re.M)
    if not table:
        return
    values = table[1].splitlines()[-1].split()
    if len(values) != 9:
        return
    pf, df, _, _, p, d, _, _ = map(float, values[1:])
    if not all(map(math.isfinite, (pf, df, p, d))):
        return
    row.update(stock_pfeas=pf, stock_dfeas=df, primal_objective_per_site=p,
               dual_objective_per_site=d, solver_objective_relative_gap=abs(p-d)/max(1, abs(p), abs(d)))
    row["stock_mosek_final_iteration_seconds"] = float(values[-1])


def validate_bound(row, case, exact):
    e = row.get("energy_per_site")
    p = row.get("primal_objective_per_site")
    feasible = (max(row.get("stock_pfeas", math.inf), row.get("stock_dfeas", math.inf)) <= 1e-8
                if row.get("method") == "reference" else
                row.get("primal_status") == row.get("dual_status") == "FEASIBLE_POINT")
    if (row["outcome"] != "ok" or row.get("termination_status") != "OPTIMAL" or not feasible or
            row.get("solver_objective_relative_gap", math.inf) > 1e-7 or
            e is None or not math.isfinite(e) or p is None or abs(e - p) > ENERGY_ATOL):
        if row["outcome"] == "ok":
            row["outcome"] = "solver_rejected"
        return
    if case.strengthened and row.get("method") != "reference" and (
            row.get("rdm_psd_blocks", 0) <= 0 or (case.lattice == "chain" and row.get("pso_input_rows", 0) <= 0)):
        row["outcome"] = "missing_strengthening"
        return
    bonds_per_site = 1 if case.lattice == "chain" else 2
    floor = -0.75 * bonds_per_site   # Each bond has minimum eigenvalue -3/4.
    neel = -0.25 * bonds_per_site    # A staggered Z product state on either even lattice.
    if exact is not None:
        row["signed_error_per_site"] = e - exact
        row["absolute_error_per_site"] = abs(e - exact)
        row["exact_gap_per_site"] = exact - e
        row["relative_exact_gap"] = (exact - e) / abs(exact)
    if not floor - ENERGY_ATOL <= e <= (exact if exact is not None else neel) + ENERGY_ATOL:
        row["outcome"] = "energy_check_failed"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true", help="ED/topology only; no solver or package setup")
    parser.add_argument("--env", type=Path, help="shared Julia environment from setup.jl")
    parser.add_argument("--out", required=True, type=Path, help="new result directory outside the checkout")
    parser.add_argument("--reviewed-cvm-commit", help="full reviewed integrated HEAD; attest this is the Tencent CVM")
    parser.add_argument("--cases", nargs="+", choices=[c.id for c in CASES])
    parser.add_argument("--supports", nargs="+", choices=["graph_distance", "local_products", "full"],
                        default=["graph_distance", "local_products"])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--timeout", type=int, default=1800, help="seconds per child; same for both methods")
    args = parser.parse_args()
    if args.repeats < 3 or args.timeout <= 0:
        parser.error("use at least 3 repeats and a positive timeout")
    if platform.system() != "Linux":
        parser.error("Linux required for wait4 peak RSS in KiB")
    out = args.out.resolve()
    env = dict(os.environ, JULIA_NUM_THREADS="1,0", JULIA_NUM_GC_THREADS="1,0",
               JULIA_LOAD_PATH="@:@stdlib", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", BLIS_NUM_THREADS="1")
    julia = ["julia", "--startup-file=no", "--threads=1,0", "--gcthreads=1,0"]
    # Measure the actual VM executable, not a juliaup or other launcher process.
    julia[0] = subprocess.check_output(julia + ["-e", "print(joinpath(Sys.BINDIR, Base.julia_exename()))"],
                                       env=env, text=True, timeout=args.timeout)
    worker = str(HERE / "worker.jl")
    if args.smoke:
        out.mkdir(parents=True)
        row = run_process(julia + [worker, "smoke"], out / "smoke.log", args.timeout, env)
        print(json.dumps(row))
        return int(row["exit_code"] != 0)

    if not args.env or out.is_relative_to(REPO):
        parser.error("measurement needs --env and a new --out outside the checkout")
    head = git("rev-parse", "HEAD")
    if args.reviewed_cvm_commit != head or git("status", "--porcelain"):
        parser.error("measurement needs a clean, reviewed integrated HEAD (--reviewed-cvm-commit)")
    subprocess.run(["git", "merge-base", "--is-ancestor", BASE, head], cwd=REPO, check=True)
    license_path = Path(env.get("MOSEKLM_LICENSE_FILE", "~/mosek/mosek.lic")).expanduser().resolve()
    if not license_path.is_file():
        parser.error("one readable local Mosek license file is required; no solver fallback")
    env["MOSEKLM_LICENSE_FILE"] = str(license_path)
    mem_kib = int(Path("/proc/meminfo").read_text().splitlines()[0].split()[1])
    if not 110 * GIB <= mem_kib * 1024 <= 130 * GIB:
        parser.error("measurement requires the agreed 128 GB Tencent CVM")
    group, memory_cap = cgroup_guard()
    root = args.env.resolve()
    projects = {lattice: root / lattice for lattice in REF}
    locks = {lattice: {name: (project / name).read_bytes() for name in ("Project.toml", "Manifest.toml")}
             for lattice, project in projects.items()}
    out.mkdir(parents=True, exist_ok=False)
    for lattice, files in locks.items():
        (out / lattice).mkdir()
        for name, data in files.items():
            (out / lattice / name).write_bytes(data)
    selected = [c for c in CASES if c.id in (args.cases or [c.id for c in CASES[:7]])]
    metadata = {"prototype_revision": head, "reference_revisions": REF,
                "lock_sha256": {lat: {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
                                for lat, files in locks.items()}, "supports": args.supports,
                "julia_executable": julia[0],
                "uname": list(platform.uname()), "mem_total_kib": mem_kib,
                "cpu_affinity": sorted(os.sched_getaffinity(0)), "memory_max_bytes": memory_cap,
                "cpu": subprocess.check_output(["lscpu"], text=True), "started_utc": time.strftime("%FT%TZ", time.gmtime()),
                "threads": 1, "repeats": args.repeats, "child_timeout_seconds": args.timeout,
                "energy_atol_per_site": ENERGY_ATOL, "cases": [dataclasses.asdict(c) for c in selected]}
    (out / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def invoke(method, label, case):
        data = out / f"{label}.toml"
        cmd = julia + [f"--project={projects[case.lattice]}", worker, method, str(data), case.lattice, str(case.extent), str(case.reference_extra),
                       str(case.strengthened).lower()]
        before = oom_kills(group)
        log = out / f"{label}.log"
        row = run_process(cmd, log, args.timeout, env)
        row.update(method=method, case=case.id, label=label)
        if data.exists():
            row.update(tomllib.loads(data.read_text()))
        row["outcome"] = ("timeout" if row["timed_out"] else "oom" if oom_kills(group) > before
                          else "worker_error" if row["exit_code"] != 0 or not data.exists() else "ok")
        if method == "reference" and row["outcome"] == "ok":
            stock_diagnostics(row, log.read_text())
        return row

    failed = False
    with (out / "results.jsonl").open("w", buffering=1) as results:
        def record(row):
            results.write(json.dumps(row, allow_nan=False) + "\n")

        for lattice in sorted({case.lattice for case in selected}):
            label = f"preflight-{lattice}"
            preflight = invoke("preflight", label, next(c for c in selected if c.lattice == lattice))
            record(preflight)
            if (preflight["outcome"] != "ok" or preflight.get("mosek_version") != "11.2.8" or
                    preflight.get("julia_version") != "1.12.7"):
                raise SystemExit(f"Mosek/Julia preflight failed; see {out / (label + '.log')}")
        oracles, audits = {}, {}
        for case in selected:
            exact = None
            if case.exact_control or case.lattice == "chain":
                method = "ed" if case.exact_control else "bethe"
                key = (case.lattice, case.extent, method)
                if key not in oracles:
                    oracles[key] = invoke(method, f"{case.id}-{method}", case)
                    record(oracles[key])
                oracle = oracles[key]
                if oracle["outcome"] != "ok":
                    failed = True
                    continue
                exact = oracle["exact_energy_per_site"]
            methods = ["reference"]
            for support in dict.fromkeys(args.supports):
                key = (case.lattice, case.extent, support)
                if key not in audits:
                    audit = invoke(f"audit_{support}", f"{case.id}-audit-{support}", case)
                    audits[key] = audit
                    record(audit)
                if audits[key]["outcome"] == "ok":
                    methods.append(support)
                else:
                    failed = True
            for i in range(args.repeats):
                rows = {}
                for method in (methods if i % 2 == 0 else methods[::-1]):
                    row = invoke(method, f"{case.id}-{i+1}-{method}", case)
                    if method != "reference":
                        audit = audits[case.lattice, case.extent, method]
                        row["comparison_lane"] = audit["comparison_lane"]
                        if row["outcome"] == "ok" and (row.get("prototype_row_counts") != [audit["selected_rows"]] or
                            row.get("translation_orders") != audit["retained_translation_orders"] or
                            row.get("point_order") != audit["retained_point_order"]):
                            row["outcome"] = "audit_mismatch"
                    validate_bound(row, case, exact)
                    rows[method] = row
                    record(row)
                reference = rows.get("reference")
                for method, row in rows.items():
                    if method == "reference" or not reference or reference["outcome"] != "ok" or row["outcome"] != "ok":
                        continue
                    audit = audits[case.lattice, case.extent, method]
                    record({"case": case.id, "replicate": i + 1, "candidate": method,
                            "comparison_lane": audit["comparison_lane"], "row_span_match": audit["row_span_match"],
                            "verified_symmetry_match": audit["verified_symmetry_match"], "strength_equivalence": "not_established",
                            "prototype_minus_reference_per_site": row["energy_per_site"] - reference["energy_per_site"]})
                rejected = {m for m, row in rows.items() if row["outcome"] != "ok"}
                failed |= bool(rejected)
                methods = [m for m in methods if m not in rejected]
                if not methods:
                    break
    if (git("rev-parse", "HEAD") != head or git("status", "--porcelain") or
            any((projects[lat] / name).read_bytes() != data for lat, files in locks.items() for name, data in files.items())):
        raise SystemExit("source or manifest changed during the run; discard these results")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())

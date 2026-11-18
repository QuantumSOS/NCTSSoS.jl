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
import signal
import subprocess
import time
import tomllib

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BASE = "e59863669f0fb30497c84f219240c1cf1bf38973"
REF = "0b152dfcd9bcde3403aa4d967642ce5261c80712"
GIB = 1 << 30
ENERGY_ATOL = 1e-6


@dataclasses.dataclass(frozen=True)
class Case:
    lattice: str
    extent: int
    reference_extra: int
    prototype_pair_range: int
    exact_control: bool = False
    strengthened: bool = False

    @property
    def id(self):
        return f"{self.lattice}-{self.extent}" + ("-rdm" if self.strengthened else "")


# Frozen order-2 moment bases; non-equivalent optional strengthenings.
# 'extra' and 'pair_range' are deliberately distinct, not equivalent knobs.
CASES = [Case("chain", 4, 1, 2, True), Case("chain", 6, 2, 3, True),
         Case("square", 4, 0, 4, True),
         *(Case("chain", L, 3, 4) for L in (24, 40, 64)),
         Case("chain", 8, 3, 4, True, True), Case("chain", 40, 3, 4, False, True)]


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


def validate_bound(row, case, exact):
    e = row.get("energy_per_site")
    p = row.get("primal_objective_per_site")
    if (row["outcome"] != "ok" or row.get("termination_status") != "OPTIMAL" or
            row.get("primal_status") != "FEASIBLE_POINT" or row.get("dual_status") != "FEASIBLE_POINT" or
            row.get("solver_objective_relative_gap", math.inf) > 1e-7 or
            e is None or not math.isfinite(e) or p is None or abs(e - p) > ENERGY_ATOL):
        if row["outcome"] == "ok":
            row["outcome"] = "solver_rejected"
        return
    if case.strengthened and (row.get("rdm_psd_blocks", 0) <= 0 or row.get("stationarity_equalities", 0) <= 0):
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
    project = args.env.resolve()
    manifest = project / "Manifest.toml"
    fingerprint = hashlib.sha256(manifest.read_bytes()).hexdigest()
    out.mkdir(parents=True)
    (out / "Manifest.toml").write_bytes(manifest.read_bytes())
    (out / "Project.toml").write_bytes((project / "Project.toml").read_bytes())
    julia += [f"--project={project}", worker]
    selected = [c for c in CASES if c.id in (args.cases or ["chain-4", "chain-6", "square-4", "chain-40"])]
    metadata = {"prototype_revision": head, "reference_revision": REF, "manifest_sha256": fingerprint,
                "julia_executable": julia[0],
                "uname": list(platform.uname()), "mem_total_kib": mem_kib,
                "cpu_affinity": sorted(os.sched_getaffinity(0)), "memory_max_bytes": memory_cap,
                "cpu": subprocess.check_output(["lscpu"], text=True), "started_utc": time.strftime("%FT%TZ", time.gmtime()),
                "threads": 1, "repeats": args.repeats, "child_timeout_seconds": args.timeout,
                "energy_atol_per_site": ENERGY_ATOL, "cases": [dataclasses.asdict(c) for c in selected]}
    (out / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def invoke(method, label, case):
        data = out / f"{label}.toml"
        cmd = julia + [method, str(data), case.lattice, str(case.extent), str(case.reference_extra),
                       str(case.prototype_pair_range), str(case.exact_control).lower(), str(case.strengthened).lower()]
        before = oom_kills(group)
        log = out / f"{label}.log"
        row = run_process(cmd, log, args.timeout, env)
        row.update(method=method, case=case.id, label=label)
        if data.exists():
            row.update(tomllib.loads(data.read_text()))
        row["outcome"] = ("timeout" if row["timed_out"] else "oom" if oom_kills(group) > before
                          else "worker_error" if row["exit_code"] != 0 or not data.exists() else "ok")
        if (row["outcome"] == "worker_error" and method == "reference" and case.lattice == "square" and
                "UndefVarError: `resort` not defined in `QMBCertify`" in log.read_text(errors="replace")):
            row["outcome"] = "upstream_failure"
        return row

    failed = False
    with (out / "results.jsonl").open("w", buffering=1) as results:
        def record(row):
            results.write(json.dumps(row, allow_nan=False) + "\n")

        preflight = invoke("preflight", "preflight", CASES[0])
        record(preflight)
        if preflight["outcome"] != "ok":
            raise SystemExit(f"Mosek preflight failed; see {out / 'preflight.log'}")
        for case in selected:
            exact = None
            if case.exact_control:
                oracle = invoke("ed", f"{case.id}-ed", case)
                record(oracle)
                if oracle["outcome"] != "ok":
                    failed = True
                    continue
                exact = oracle["exact_energy_per_site"]
            for i in range(args.repeats):
                rows = {}
                methods = ("reference", "prototype") if i % 2 == 0 else ("prototype", "reference")
                for method in methods:
                    row = invoke(method, f"{case.id}-{i+1}-{method}", case)
                    validate_bound(row, case, exact)
                    rows[method] = row
                    record(row)
                pair = {"case": case.id, "replicate": i + 1, "strength_equivalence": "not_established"}
                if all(r["outcome"] == "ok" for r in rows.values()):
                    pair["prototype_minus_reference_per_site"] = (
                        rows["prototype"]["energy_per_site"] - rows["reference"]["energy_per_site"])
                else:
                    failed = True
                record(pair)
                # Do not pay for the same failure again; preserve any successful mate.
                if "prototype_minus_reference_per_site" not in pair:
                    break
    if git("rev-parse", "HEAD") != head or git("status", "--porcelain") or hashlib.sha256(manifest.read_bytes()).hexdigest() != fingerprint:
        raise SystemExit("source or manifest changed during the run; discard these results")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())

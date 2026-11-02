"""Print the complete scale table from retained records, including failures."""
import json
import statistics
import tomllib
from pathlib import Path

root = Path(__file__).parent / "results" / "dev382_scale"
print("| Record | Warm samples | Build s | Solver s | Process s | Peak GiB | Bound | Residual l1 | Gap | Outcome |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
for path in sorted(root.glob("*.json")):
    process = json.loads(path.read_text())
    data_path = path.with_suffix(".toml")
    data = tomllib.loads(data_path.read_text()) if data_path.exists() else {}
    rows = data.get("samples", [])
    warm = [s for s in rows if not s["cold"] and s.get("stage") == "complete"]
    measured = warm or rows

    def metric(key):
        values = [s[key] for s in measured if key in s]
        if not values:
            return "—"
        if len(values) == 1:
            return f"{values[0]:.3f}†"
        return f"{statistics.median(values):.3f} [{min(values):.3f}, {max(values):.3f}]"

    complete = [s for s in rows if s.get("stage") == "complete"]
    bound = f'{complete[-1]["bound"]:.10f}' if complete else "—"
    residual = f'{max(s["residual_l1"] for s in complete):.2e}' if complete else "—"
    gap = f'{max(s["gap"] for s in complete):.2e}' if complete else "—"
    if process["reason"] != "exited":
        outcome = process["reason"]
    elif process["returncode"] != 0:
        outcome = f'exit {process["returncode"]}'
    elif rows:
        outcome = rows[-1].get("status", rows[-1]["stage"])
        if rows[-1].get("valid"):
            outcome += "; validated"
    else:
        outcome = "no sample"
    peak = max(process["peak_rss_bytes"], data.get("peak_rss_bytes", 0))/1024**3
    print(f"| [{path.stem}](results/dev382_scale/{path.stem}.toml) "
          f"([command/limits](results/dev382_scale/{path.stem}.json)) | {len(warm)} | "
          f"{metric('construction_s')} | {metric('solver_s')} | {process['wall_s']:.2f} | "
          f"{peak:.3f} | {bound} | {residual} | {gap} | {outcome} |")
print("\nTimes are median [minimum, maximum] for warmed samples; † means one cold sample.")
print("Process time and peak resident memory include startup, compilation, all samples, and validation.")
print("Failures have no validated bound; an absent bound is not replaced by an infeasible solver iterate.")

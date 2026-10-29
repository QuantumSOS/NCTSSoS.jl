"""Run one benchmark with the predeclared identical wall/RSS budgets."""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

prefix = Path(sys.argv[1])
command = sys.argv[3:]
assert sys.argv[2] == "--" and command
started = time.monotonic()
peak = 0
reason = "exited"
with prefix.with_suffix(".log").open("w") as log:
    process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    while process.poll() is None:
        try:
            fields = Path(f"/proc/{process.pid}/status").read_text().splitlines()
            rss = max(int(line.split()[1]) * 1024 for line in fields
                      if line.startswith(("VmRSS:", "VmHWM:")))
            peak = max(peak, rss)
        except (FileNotFoundError, ValueError, ProcessLookupError):
            pass
        if peak > 24 * 1024**3 or time.monotonic() - started > 900:
            reason = "rss_limit" if peak > 24 * 1024**3 else "wall_timeout"
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
            break
        time.sleep(0.1)
    code = process.wait()
record = dict(command=command, wall_s=time.monotonic()-started,
              peak_rss_bytes=peak, reason=reason, returncode=code,
              wall_limit_s=900, rss_limit_bytes=24*1024**3)
prefix.with_suffix(".json").write_text(json.dumps(record, indent=2)+"\n")
print(json.dumps(record), flush=True)

"""Run a validation command with a wall-time and process-RSS safety limit."""
import os
import subprocess
import sys
import time

seconds = float(os.environ.get("RUN_SECONDS", "600"))
p = subprocess.Popen(sys.argv[1:])
started = time.monotonic()
peak = 0
while p.poll() is None:
    try:
        with open(f"/proc/{p.pid}/statm") as f:
            rss = int(f.read().split()[1]) * os.sysconf("SC_PAGE_SIZE")
        peak = max(peak, rss)
    except FileNotFoundError:
        continue
    if rss > 16 * 2**30 or time.monotonic() - started > seconds:
        p.terminate()
        try:
            p.wait(timeout=5)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
        print(f"LIMIT_EXCEEDED: 16 GiB RSS or {seconds:g} s; "
              f"elapsed_seconds={time.monotonic()-started:.3f} "
              f"peak_RSS_GiB={peak/2**30:.3f}", flush=True)
        sys.exit(124)
    time.sleep(0.5)
print(f"elapsed_seconds={time.monotonic()-started:.3f} peak_RSS_GiB={peak/2**30:.3f}", flush=True)
sys.exit(p.returncode)

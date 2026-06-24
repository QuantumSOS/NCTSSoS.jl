#!/usr/bin/env python3
"""Fetch a live AutoDL Pro SSH command, then submit the N=100 run.

Requires an AutoDL developer token in AUTODL_TOKEN unless --token-env is used.
The script deliberately does not print root_password or the token.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


API_HOST = "https://api.autodl.com"


def request(token: str, method: str, path: str, body: dict[str, Any] | None = None) -> dict[str, Any]:
    data = None if body is None else json.dumps(body).encode("utf-8")
    req = urllib.request.Request(
        API_HOST + path,
        data=data,
        method=method,
        headers={
            "Authorization": token,
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            payload = resp.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"AutoDL API HTTP {exc.code}: {detail}") from exc
    except urllib.error.URLError as exc:
        raise SystemExit(f"AutoDL API request failed: {exc}") from exc

    try:
        parsed = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"AutoDL API returned non-JSON payload: {payload[:500]}") from exc
    if parsed.get("code") != "Success":
        msg = parsed.get("msg") or parsed
        raise SystemExit(f"AutoDL API error: {msg}")
    return parsed


def list_instances(token: str, page_size: int = 50) -> list[dict[str, Any]]:
    page = 1
    out: list[dict[str, Any]] = []
    while True:
        parsed = request(
            token,
            "POST",
            "/api/v1/dev/instance/pro/list",
            {"page_index": page, "page_size": page_size},
        )
        data = parsed.get("data") or {}
        out.extend(data.get("list") or [])
        max_page = int(data.get("max_page") or page)
        if page >= max_page:
            return out
        page += 1


def snapshot(token: str, instance_uuid: str) -> dict[str, Any]:
    parsed = request(
        token,
        "GET",
        "/api/v1/dev/instance/pro/snapshot",
        {"instance_uuid": instance_uuid},
    )
    data = parsed.get("data")
    if not isinstance(data, dict):
        raise SystemExit(f"snapshot for {instance_uuid} had no data")
    return data


def power_on(token: str, instance_uuid: str) -> None:
    request(
        token,
        "POST",
        "/api/v1/dev/instance/pro/power_on",
        {"instance_uuid": instance_uuid, "payload": "gpu"},
    )


def choose_instance(instances: list[dict[str, Any]], instance_uuid: str | None) -> dict[str, Any]:
    if instance_uuid:
        for item in instances:
            if item.get("uuid") == instance_uuid:
                return item
        raise SystemExit(f"instance {instance_uuid} not found in AutoDL Pro instance list")

    running = [item for item in instances if item.get("status") == "running"]
    if len(running) == 1:
        return running[0]
    if len(running) > 1:
        print("Multiple running AutoDL Pro instances found; pass --instance:", file=sys.stderr)
        for item in running:
            print_instance(item, file=sys.stderr)
        raise SystemExit(2)

    print("No running AutoDL Pro instance found.", file=sys.stderr)
    if instances:
        print("Known instances:", file=sys.stderr)
        for item in instances[:20]:
            print_instance(item, file=sys.stderr)
    raise SystemExit(2)


def print_instance(item: dict[str, Any], *, file: Any = sys.stdout) -> None:
    print(
        f"- {item.get('uuid')} status={item.get('status')} "
        f"name={item.get('name')!r} gpu={item.get('gpu_spec_uuid')} region={item.get('region_sign')}",
        file=file,
    )


def run_submit(repo: Path, ssh_command: str, submit: bool, check_only: bool) -> int:
    script = repo / "perf" / "submit_heisenberg_n100_autodl_from_ssh_command.sh"
    cmd = [str(script)]
    if check_only:
        cmd.append("--check-only")
    cmd.append(ssh_command)
    if not submit and not check_only:
        print(ssh_command)
        return 0
    return subprocess.call(cmd, cwd=repo)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--instance", help="AutoDL Pro instance UUID. Required if more than one instance is running.")
    parser.add_argument("--token-env", default="AUTODL_TOKEN", help="Environment variable containing the developer token.")
    parser.add_argument("--list", action="store_true", help="List Pro instances and exit.")
    parser.add_argument("--power-on", action="store_true", help="Power on --instance before fetching the SSH command.")
    parser.add_argument("--wait", type=int, default=0, help="Seconds to wait after --power-on before snapshot.")
    parser.add_argument("--check-only", action="store_true", help="Fetch SSH command and run readiness checks only.")
    parser.add_argument("--print-ssh", action="store_true", help="Print the fetched SSH command instead of submitting.")
    args = parser.parse_args()

    token = os.environ.get(args.token_env, "").strip()
    if not token:
        raise SystemExit(f"missing AutoDL developer token: set {args.token_env}")

    repo = Path(__file__).resolve().parents[1]
    instances = list_instances(token)
    if args.list:
        for item in instances:
            print_instance(item)
        return 0

    item = choose_instance(instances, args.instance)
    instance_uuid = str(item["uuid"])
    if args.power_on:
        if not args.instance:
            raise SystemExit("--power-on requires --instance to avoid starting the wrong paid instance")
        print(f"Powering on AutoDL instance {instance_uuid}...")
        power_on(token, instance_uuid)
        if args.wait:
            time.sleep(args.wait)

    snap = snapshot(token, instance_uuid)
    ssh_command = str(snap.get("ssh_command") or "").strip()
    if not ssh_command:
        proxy_host = snap.get("proxy_host")
        ssh_port = snap.get("ssh_port")
        if proxy_host and ssh_port:
            ssh_command = f"ssh -p {ssh_port} root@{proxy_host}"
    if not ssh_command:
        raise SystemExit(f"AutoDL snapshot for {instance_uuid} did not include an SSH command")

    print(f"Using AutoDL instance {instance_uuid}; status={item.get('status')}; ssh_command={ssh_command}")
    return run_submit(repo, ssh_command, submit=not args.print_ssh, check_only=args.check_only)


if __name__ == "__main__":
    raise SystemExit(main())

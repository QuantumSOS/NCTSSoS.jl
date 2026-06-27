#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  perf/submit_heisenberg_n100_autodl_from_ssh_command.sh [--check-only] "ssh -p PORT root@HOST"

This creates a temporary SSH config with Host autodl, then uses the existing
.easy-ssh.conf remote_dir. It does not edit ~/.ssh/config.

Environment:
  NCTS_AUTODL_IDENTITY   SSH key, default: ~/.ssh/id_ed25519
  NCTS_REMOTE_CMD        Remote command, default: bash perf/run_heisenberg_n100_mosek_remote.sh
  NCTS_MONITOR           1 to monitor after submit, 0 to return after submit. Default: 1
EOF
}

check_only=0
if [[ "${1:-}" == "--check-only" ]]; then
    check_only=1
    shift
fi

if [[ $# -ne 1 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 2
fi

repo_root="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$repo_root"

if [[ ! -f .easy-ssh.conf ]]; then
    echo "error: .easy-ssh.conf is missing" >&2
    exit 1
fi

real_home="${HOME:?}"
identity="${NCTS_AUTODL_IDENTITY:-$real_home/.ssh/id_ed25519}"
if [[ ! -f "$identity" ]]; then
    echo "error: SSH identity not found: $identity" >&2
    exit 1
fi
real_ssh="$(command -v ssh)"

parsed="$(
    python3 - "$1" <<'PY'
import shlex
import sys

tokens = shlex.split(sys.argv[1])
if not tokens or tokens[0] != "ssh":
    raise SystemExit("expected an ssh command")

port = None
user = None
host = None
i = 1
while i < len(tokens):
    tok = tokens[i]
    if tok == "-p" and i + 1 < len(tokens):
        port = tokens[i + 1]
        i += 2
        continue
    if tok.startswith("-p") and len(tok) > 2:
        port = tok[2:]
        i += 1
        continue
    if tok == "-l" and i + 1 < len(tokens):
        user = tokens[i + 1]
        i += 2
        continue
    if tok.startswith("-"):
        # Skip a conservative set of options with values. AutoDL's documented
        # command is simple, but this keeps harmless OpenSSH flags from being
        # misread as the host.
        if tok in {"-i", "-F", "-J", "-o", "-b", "-c", "-m"} and i + 1 < len(tokens):
            i += 2
        else:
            i += 1
        continue
    host = tok
    i += 1

if host is None:
    raise SystemExit("could not parse SSH host")
if "@" in host:
    parsed_user, parsed_host = host.rsplit("@", 1)
    user = user or parsed_user
    host = parsed_host
if not port:
    raise SystemExit("could not parse SSH port; expected -p PORT")
if not user:
    user = "root"
if not port.isdigit():
    raise SystemExit(f"invalid SSH port: {port}")

print(user)
print(host)
print(port)
PY
)"

user="$(printf '%s\n' "$parsed" | sed -n '1p')"
host="$(printf '%s\n' "$parsed" | sed -n '2p')"
port="$(printf '%s\n' "$parsed" | sed -n '3p')"

tmp_home="$(mktemp -d "${TMPDIR:-/tmp}/nctssos-autodl-home.XXXXXX")"
trap 'rm -rf "$tmp_home"' EXIT
mkdir -m 700 "$tmp_home/.ssh"
mkdir -m 700 "$tmp_home/bin"
touch "$tmp_home/.ssh/known_hosts"
chmod 600 "$tmp_home/.ssh/known_hosts"

cat >"$tmp_home/.ssh/config" <<EOF
Host autodl
  HostName $host
  Port $port
  User $user
  IdentityFile $identity
  IdentitiesOnly yes
  StrictHostKeyChecking accept-new
  UserKnownHostsFile $tmp_home/.ssh/known_hosts
  ServerAliveInterval 30
  ServerAliveCountMax 3
EOF
chmod 600 "$tmp_home/.ssh/config"

cat >"$tmp_home/bin/ssh" <<EOF
#!/usr/bin/env bash
exec "$real_ssh" -F "$tmp_home/.ssh/config" "\$@"
EOF
chmod 700 "$tmp_home/bin/ssh"

run_with_autodl_ssh() {
    PATH="$tmp_home/bin:$PATH" RSYNC_RSH="$tmp_home/bin/ssh" "$@"
}

echo "AutoDL endpoint: $user@$host:$port"
echo "Remote directory from .easy-ssh.conf:"
sed -n "s/^remote_dir=//p" .easy-ssh.conf
echo

echo "Checking SSH/easy-ssh readiness..."
run_with_autodl_ssh env NCTS_AUTODL_CHECK_TIMEOUT="${NCTS_AUTODL_CHECK_TIMEOUT:-8}" EASY_SSH_CONNECT_TIMEOUT="${EASY_SSH_CONNECT_TIMEOUT:-30}" bash perf/check_autodl_ssh.sh
run_with_autodl_ssh env EASY_SSH_CONNECT_TIMEOUT="${EASY_SSH_CONNECT_TIMEOUT:-30}" easy-ssh status

if [[ "$check_only" == "1" ]]; then
    echo "check-only: endpoint is ready"
    exit 0
fi

remote_cmd="${NCTS_REMOTE_CMD:-bash perf/run_heisenberg_n100_mosek_remote.sh}"
echo
echo "Submitting N=100 sparse-degree-4 Heisenberg MOSEK run."
echo "Estimated elapsed after sync: 5-15 minutes if Julia/MOSEK packages are already warm; first setup can take longer."
echo "Remote command: $remote_cmd"
run_with_autodl_ssh env EASY_SSH_CONNECT_TIMEOUT="${EASY_SSH_CONNECT_TIMEOUT:-30}" easy-ssh submit "$remote_cmd"

if [[ "${NCTS_MONITOR:-1}" == "1" ]]; then
    run_with_autodl_ssh easy-ssh monitor
fi

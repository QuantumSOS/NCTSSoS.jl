#!/usr/bin/env bash
set -euo pipefail

host="${1:-autodl}"
timeout_s="${NCTS_AUTODL_CHECK_TIMEOUT:-20}"
rc=0
banner_ok=0

section() {
    printf '\n## %s\n' "$1"
}

section "ssh -G $host"
ssh -G "$host" 2>/dev/null | awk '
    /^(user|hostname|port|proxycommand|identityfile|connecttimeout) / { print }
'

hostname="$(ssh -G "$host" 2>/dev/null | awk '$1 == "hostname" { print $2; exit }')"
port="$(ssh -G "$host" 2>/dev/null | awk '$1 == "port" { print $2; exit }')"

section "DNS"
doh_json=""
doh_ip=""
if command -v dig >/dev/null 2>&1; then
    printf 'local dig: '
    dig +short "$hostname" || true
fi
if command -v curl >/dev/null 2>&1; then
    printf 'DoH A: '
    doh_json="$(curl -fsSL --max-time "$timeout_s" \
        -H 'accept: application/dns-json' \
        "https://cloudflare-dns.com/dns-query?name=${hostname}&type=A" || true)"
    printf '%s' "$doh_json" | tr -d '\n'
    printf '\n'
    doh_ip="$(printf '%s' "$doh_json" |
        perl -ne 'while (/"data":"([0-9]+(?:\.[0-9]+){3})"/g) { print $1; exit }')"
fi

read_banner() {
    local target="$1"
perl -MIO::Socket::INET -e '
    my ($h, $p, $t) = @ARGV;
    my $s = IO::Socket::INET->new(PeerHost => $h, PeerPort => $p, Proto => "tcp", Timeout => $t)
        or die "connect failed\n";
    eval {
        local $SIG{ALRM} = sub { die "timeout\n" };
        alarm $t;
        my $n = read($s, my $buf, 80);
        alarm 0;
        if (defined($n) && $n > 0) {
            $buf =~ s/\r?\n\z//;
            print "$buf\n";
            exit($buf =~ /^SSH-/ ? 0 : 2);
        } else {
            print "no bytes before close\n";
            exit 2;
        }
    };
    if ($@) {
        chomp $@;
        print "$@\n";
        exit 2;
    }
' "$target" "$port" "$timeout_s"
}

section "TCP"
if command -v nc >/dev/null 2>&1; then
    nc -vz -w "$timeout_s" "$hostname" "$port" || true
fi

section "SSH banner"
if read_banner "$hostname"; then
    banner_ok=1
else
    rc=1
fi

if [[ -n "$doh_ip" && "$doh_ip" != "$hostname" ]]; then
    section "DoH IP TCP"
    if command -v nc >/dev/null 2>&1; then
        nc -vz -w "$timeout_s" "$doh_ip" "$port" || true
    fi

    section "DoH IP SSH banner"
    if read_banner "$doh_ip"; then
        banner_ok=1
    fi
fi

if (( ! banner_ok )); then
    rc=1
fi

section "raw ssh"
ssh -o BatchMode=yes -o ConnectTimeout="$timeout_s" "$host" echo ok || rc=1

section "easy-ssh status"
if command -v easy-ssh >/dev/null 2>&1; then
    EASY_SSH_CONNECT_TIMEOUT="$timeout_s" easy-ssh status || rc=1
else
    echo "easy-ssh not found"
    rc=1
fi

exit "$rc"

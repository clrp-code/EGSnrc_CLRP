#!/usr/bin/env bash
# Quick regression checks for eb-setup.sh — run before pushing feature/eb-setup.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

pass=0
fail=0

ok()   { printf 'ok: %s\n' "$*"; pass=$((pass + 1)); }
bad()  { printf 'FAIL: %s\n' "$*" >&2; fail=$((fail + 1)); }

# 1) configure launcher dry-run (no doubled my_dir path)
out=$(unset EGS_CONFIG HEN_HOUSE EGS_HOME; ./eb-setup.sh install --dry-run 2>&1 || true)
if echo "$out" | grep -q 'HEN_HOUSE/scripts && unset HEN_HOUSE EGS_HOME EGS_CONFIG && ./configure'; then
    ok 'configure dry-run uses cd + ./configure'
else
    bad 'configure dry-run missing expected cd ./configure pattern'
fi

# 2) subshell ./configure must not double scripts path in my_dir
if [[ -f HEN_HOUSE/scripts/configure ]]; then
    my_dir=$(
        unset HEN_HOUSE EGS_HOME EGS_CONFIG
        cd HEN_HOUSE/scripts
        sh -c 'my_dir=$(pwd)/$(dirname $0); echo "$my_dir"' ./configure
    )
    if [[ "$my_dir" == *"/scripts//"* ]] || [[ "$my_dir" == *"//Users/"* ]]; then
        bad "configure my_dir doubled: $my_dir"
    else
        ok "configure my_dir sane: $my_dir"
    fi
else
    bad 'HEN_HOUSE/scripts/configure not found'
fi

# 3) path-length guard rejects long REPO_ROOT on install
guard_out=$(bash -c '
    source "'"$ROOT"'/scripts/lib/common.sh"
    REPO_ROOT="/Users/marc/Developer/scratch/EGSnrc_CLRP-test"
    EB_CMD=install
    check_mortran_path_lengths
' 2>&1 || true)
if echo "$guard_out" | grep -q 'path too long'; then
    ok 'Mortran path guard rejects long clone path on install'
else
    bad 'Mortran path guard did not reject long path'
fi

# 4) expand_user_path
if [[ "$(bash -c 'source "'"$ROOT"'/scripts/lib/common.sh"; expand_user_path "~/foo/bar"')" == "${HOME}/foo/bar" ]]; then
    ok 'expand_user_path handles ~/'
else
    bad 'expand_user_path'
fi

echo "---"
printf '%d passed, %d failed\n' "$pass" "$fail"
(( fail == 0 ))

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

# 5) tarball install-dir arg (must capture $2; ~ expanded)
tar_out=$(bash -c '
    source "'"$ROOT"'/scripts/lib/common.sh"
    source "'"$ROOT"'/scripts/lib/tarball.sh"
    DRY_RUN=1
    tarball_extract_install_dir "/fake/archive.tar.gz" "~/tarball-test/eb-release"
' 2>&1)
if echo "$tar_out" | grep -q "${HOME}/tarball-test/eb-release"; then
    ok 'tarball_extract_install_dir expands install-dir'
else
    bad "tarball_extract_install_dir: $tar_out"
fi

# 6) configure log EGS_HOME wins over a conflicting --egs-home
resolve_out=$(bash -c '
    source "'"$ROOT"'/scripts/lib/common.sh"
    export EGS_CONFIG="'"$HOME"'/scratch/tarball-test/eb-release/HEN_HOUSE/specs/tarball.conf"
    EGS_HOME_OVERRIDE="/tmp/eb-setup-wrong-egs_home"
    resolve_paths
    echo "$EGS_HOME_RESOLVED"
' 2>&1)
if [[ -f "$HOME/scratch/tarball-test/eb-release/HEN_HOUSE/log/configure-tarball-marc.log" ]]; then
    log_home=$(grep -E '^EGS_HOME:[[:space:]]+' \
        "$HOME/scratch/tarball-test/eb-release/HEN_HOUSE/log/configure-tarball-marc.log" \
        | head -1 | sed -E 's/^EGS_HOME:[[:space:]]+//;s/[[:space:]]+$//')
    if [[ -n "$log_home" && "$resolve_out" == "${log_home%/}/" ]]; then
        ok 'resolve_paths prefers finalize log over --egs-home'
    else
        bad "resolve_paths should use configure log ($log_home), got: $resolve_out"
    fi
else
    ok 'resolve_paths configure-log test skipped (no tarball install)'
fi

# 7) pickup after configure (uses scratch install if present)
if [[ -f "$HOME/scratch/eb/HEN_HOUSE/specs/test.conf" ]]; then
    pickup_out=$(bash -c '
        source "'"$ROOT"'/scripts/lib/common.sh"
        HEN_HOUSE="'"$HOME"'/scratch/eb/HEN_HOUSE"
        REPO_ROOT="'"$HOME"'/scratch/eb"
        pickup_egsnrc_env_after_configure
        echo "$EGS_CONFIG_RESOLVED|$EGS_HOME_RESOLVED"
    ' 2>&1)
    if echo "$pickup_out" | grep -q 'specs/test.conf' && echo "$pickup_out" | grep -q 'egs_home'; then
        ok 'pickup_egsnrc_env_after_configure finds test.conf + egs_home'
    else
        bad "pickup failed: $pickup_out"
    fi
fi

# 8) eb-env.sh is EGSnrc-only (no CLRP bashrc additions)
env_out=$(bash -c '
    tmp=$(mktemp -d)
    mkdir -p "$tmp/HEN_HOUSE/scripts" "$tmp/HEN_HOUSE/specs"
    printf "my_machine = linux\nHEN_HOUSE = %s/HEN_HOUSE\n" "$tmp" > "$tmp/HEN_HOUSE/specs/test.conf"
    source "'"$ROOT"'/scripts/lib/common.sh"
    REPO_ROOT="$tmp"
    HEN_HOUSE="$tmp/HEN_HOUSE"
    EGS_CONFIG_RESOLVED="$tmp/HEN_HOUSE/specs/test.conf"
    EGS_HOME_RESOLVED="'"$HOME"'/scratch/egs_home/"
    emit_shell_setup >/dev/null
    cat "$tmp/eb-env.sh"
' 2>&1)
if echo "$env_out" | grep -q 'egsnrc_bashrc_additions' \
    && ! echo "$env_out" | grep -q 'clrp_bashrc_additions'; then
    ok 'eb-env.sh sources egsnrc_bashrc_additions only'
else
    bad "eb-env.sh should not reference clrp_bashrc_additions: $env_out"
fi

# 9) release tarball injects release.mk with EGS_RELEASE + both SHAs
if [[ -f "$ROOT/HEN_HOUSE/user_codes/egs_brachy/egs_brachy/Makefile" ]]; then
    rm_out="$("$ROOT/scripts/build-release-tarball.sh" 9.9.9-test 2>&1)" || true
    rel_mk="$ROOT/dist/EGSnrc_CLRP-egs_brachy-9.9.9-test/HEN_HOUSE/specs/release.mk"
    if [[ -f "$rel_mk" ]] \
        && grep -q 'EGS_RELEASE.*9.9.9-test' "$rel_mk" \
        && grep -q 'EGS_CLRP_HASH' "$rel_mk" \
        && grep -q 'EGS_BRACHY_HASH' "$rel_mk"; then
        ok 'build-release-tarball.sh writes release.mk with dual SHAs'
        rm -rf "$ROOT/dist/EGSnrc_CLRP-egs_brachy-9.9.9-test" \
               "$ROOT/dist/EGSnrc_CLRP-egs_brachy-9.9.9-test.tar.gz"
    else
        bad "build-release-tarball.sh release.mk: $rm_out"
    fi
else
    ok 'release.mk test skipped (egs_brachy submodule not initialized)'
fi

# 10) write_release_mk on git install tree
if [[ -d "$HOME/scratch/eb/.git" && -f "$HOME/scratch/eb/HEN_HOUSE/specs/test.conf" ]]; then
    wmk_out=$(bash -c '
        source "'"$ROOT"'/scripts/lib/common.sh"
        export EGS_CONFIG="'"$HOME"'/scratch/eb/HEN_HOUSE/specs/test.conf"
        resolve_paths
        write_release_mk
        cat "'"$HOME"'/scratch/eb/HEN_HOUSE/specs/release.mk"
    ' 2>&1)
    if echo "$wmk_out" | grep -q 'EGS_CLRP_HASH' && echo "$wmk_out" | grep -q 'EGS_BRACHY_HASH'; then
        ok 'write_release_mk writes dual SHAs on git install'
    else
        bad "write_release_mk: $wmk_out"
    fi
fi

echo "---"
printf '%d passed, %d failed\n' "$pass" "$fail"
(( fail == 0 ))

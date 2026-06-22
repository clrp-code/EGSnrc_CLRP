#!/usr/bin/env bash
# Scratch cwd smoke test (NRC #1399 I/O + #934 include paths)
set -euo pipefail

SCRATCH="${SCRATCH:-$HOME/Developer/scratch}"
EGS_HOME="${EGS_HOME:-$SCRATCH/egs_home/}"
HEN_HOUSE="${HEN_HOUSE:-/Users/marc/Developer/clrp/EGSnrc-eb/HEN_HOUSE/}"
MY_MACHINE="${MY_MACHINE:-eb-dev}"

export EGS_HOME HEN_HOUSE EGS_CONFIG="${EGS_CONFIG:-$HEN_HOUSE/specs/eb-dev.conf}"
export PATH="$EGS_HOME/bin/$MY_MACHINE:$PATH"

EB="$EGS_HOME/bin/$MY_MACHINE/egs_brachy"
TEMPLATE="${TEMPLATE:-$(dirname "$0")/smoke_ra.egsinp}"
PASS=0
FAIL=0

pass() { printf 'PASS: %s\n' "$*"; PASS=$((PASS + 1)); }
fail() { printf 'FAIL: %s\n' "$*" >&2; FAIL=$((FAIL + 1)); }

cleanup() {
    rm -rf "$SCRATCH"/egsrun_*smoke_ra* "$SCRATCH"/smoke_ra.egslst "$SCRATCH"/smoke_ra.mederr \
           "$SCRATCH"/smoke_ra.phantom.3ddose "$SCRATCH"/smoke_ra.phantom.edep.3ddose \
           "$SCRATCH"/smoke_ra.egsdat 2>/dev/null || true
    rm -rf "$EGS_HOME/egs_brachy"/egsrun_*smoke_ra* "$EGS_HOME/egs_brachy"/smoke_ra.egslst \
           "$EGS_HOME/egs_brachy"/smoke_ra.mederr "$EGS_HOME/egs_brachy"/smoke_ra.phantom.3ddose \
           "$EGS_HOME/egs_brachy"/smoke_ra.phantom.edep.3ddose 2>/dev/null || true
    rm -f "$EGS_HOME/egs_brachy/smoke_ra.egsinp" 2>/dev/null || true
}

run_smoke() {
    "$EB" -i smoke_ra -s >/dev/null
}

[[ -x "$EB" ]] || { echo "egs_brachy not found: $EB"; exit 1; }
[[ -f "$TEMPLATE" ]] || { echo "missing $TEMPLATE"; exit 1; }

cleanup
[[ -f "$SCRATCH/smoke_ra.egsinp" ]] || cp "$TEMPLATE" "$SCRATCH/smoke_ra.egsinp"

echo "=== Test 1: input in cwd → outputs in cwd ==="
cd "$SCRATCH"
run_smoke
if [[ -f "$SCRATCH/smoke_ra.egslst" || -f "$SCRATCH/smoke_ra.mederr" || -f "$SCRATCH/smoke_ra.phantom.3ddose" ]]; then
    pass "outputs in $SCRATCH"
else
    fail "no smoke_ra outputs in $SCRATCH"
fi
if [[ -f "$EGS_HOME/egs_brachy/smoke_ra.egslst" ]]; then
    fail "outputs leaked to EGS_HOME/egs_brachy"
else
    pass "EGS_HOME/egs_brachy clean"
fi

echo "=== Test 2: bare name → resolve from EGS_HOME/egs_brachy ==="
cleanup
cp "$TEMPLATE" "$EGS_HOME/egs_brachy/smoke_ra.egsinp"
rm -f "$SCRATCH/smoke_ra.egsinp" "$SCRATCH"/smoke_ra.egslst "$SCRATCH"/smoke_ra.mederr 2>/dev/null || true
cd "$SCRATCH"
[[ -f smoke_ra.egsinp ]] && fail "scratch input should be absent" || pass "no input in cwd"
run_smoke
if [[ -f "$EGS_HOME/egs_brachy/smoke_ra.egslst" || -f "$EGS_HOME/egs_brachy/smoke_ra.mederr" ]]; then
    pass "EGS_HOME fallback works; outputs next to resolved input"
else
    fail "bare name smoke_ra did not produce outputs in EGS_HOME/egs_brachy"
fi
if [[ -f "$SCRATCH/smoke_ra.egslst" ]]; then
    fail "outputs incorrectly written to cwd when input was only in EGS_HOME"
else
    pass "cwd clean when input resolved from EGS_HOME"
fi

echo "=== Test 3: \$EGS_HOME in include file and media/muen paths (#934) ==="
grep -q 'include file = \$EGS_HOME' "$TEMPLATE" && pass 'include file uses $EGS_HOME' || fail 'include file missing $EGS_HOME'
grep -q 'material data file = \$EGS_HOME' "$TEMPLATE" && pass 'material data file uses $EGS_HOME' || fail 'material data file missing $EGS_HOME'
grep -q 'muen file = \$EGS_HOME' "$TEMPLATE" && pass 'muen file uses $EGS_HOME' || fail 'muen file missing $EGS_HOME'

cleanup
[[ -f "$SCRATCH/smoke_ra.egsinp" ]] || cp "$TEMPLATE" "$SCRATCH/smoke_ra.egsinp"
rm -f "$EGS_HOME/egs_brachy/smoke_ra.egsinp"

echo "=== Summary: $PASS passed, $FAIL failed ==="
[[ "$FAIL" -eq 0 ]]

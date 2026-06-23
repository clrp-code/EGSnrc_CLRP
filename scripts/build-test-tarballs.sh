#!/usr/bin/env bash
# Build two sample release tarballs for local install + update testing.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${ROOT}/scripts/build-release-tarball.sh"
STAMP_FILE="${ROOT}/HEN_HOUSE/user_codes/egs_brachy/egs_brachy/.eb-setup-tarball-test-stamp"

[[ -x "$BUILD" ]] || { echo "error: missing $BUILD" >&2; exit 1; }

echo "=== Building tarball pair for eb-setup testing ==="
"$BUILD" "1.0.0-alpha.1"
"${ROOT}/scripts/build-eb-setup-tarball.sh" "1.0.0-alpha.1"

# Second tarball: same tree + a tiny marker so update rsync has something to merge.
echo "test-alpha.2" > "$STAMP_FILE"
"$BUILD" "1.0.0-alpha.2"
"${ROOT}/scripts/build-eb-setup-tarball.sh" "1.0.0-alpha.2"
rm -f "$STAMP_FILE"

echo
echo "Test tarballs:"
ls -lh "${ROOT}/dist/"EGSnrc_CLRP-egs_brachy-1.0.0-alpha.*.tar.gz \
       "${ROOT}/dist/"EGSnrc_CLRP-eb-setup-1.0.0-alpha.*.tar.gz 2>/dev/null || true
echo
echo "Next: see docs/eb-setup-testing.md § Tarball testing walkthrough"

#!/usr/bin/env bash
# Build the slim eb-setup bootstrap tarball (installer only — no HEN_HOUSE).
# Same layout as CI .github/workflows/release-tarball.yml.
#
# Usage:
#   ./scripts/build-eb-setup-tarball.sh 1.0.0-alpha.1
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

VERSION="${1:-}"
OUT_DIR="${ROOT}/dist"

usage() {
    cat <<EOF
Usage: $(basename "$0") VERSION

Build dist/EGSnrc_CLRP-eb-setup-VERSION.tar.gz (bootstrap installer only).

Contains eb-setup.sh and scripts/ — user runs install to download the full release.

Example:
  $(basename "$0") 1.0.0-alpha.1
EOF
}

[[ -n "$VERSION" ]] || { usage; exit 1; }
[[ -f "$ROOT/eb-setup.sh" ]] || { echo "error: run from EGSnrc_CLRP repo root" >&2; exit 1; }
[[ -d "$ROOT/scripts/lib" ]] || { echo "error: scripts/lib missing" >&2; exit 1; }

DIR="EGSnrc_CLRP-eb-setup-${VERSION}"
ARCHIVE="${OUT_DIR}/${DIR}.tar.gz"
STAGING="${OUT_DIR}/${DIR}"

rm -rf "$STAGING"
mkdir -p "$STAGING"
cp "$ROOT/eb-setup.sh" "$STAGING/"
cp -R "$ROOT/scripts" "$STAGING/"

tar -czf "$ARCHIVE" -C "$OUT_DIR" "$DIR"
echo "Wrote $ARCHIVE"
echo "  $(du -h "$ARCHIVE" | cut -f1)  $(tar -tzf "$ARCHIVE" | wc -l | tr -d ' ') files"
echo "  extract anywhere, then:  ./eb-setup.sh install"

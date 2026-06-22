#!/usr/bin/env bash
# Build a release tarball locally (same layout as .github/workflows/release-tarball.yml).
# Use for eb-setup tarball testing without pushing tags to GitHub.
#
# Usage:
#   ./scripts/build-test-tarballs.sh              # builds alpha.1 and alpha.2 pair
#   ./scripts/build-release-tarball.sh 1.0.0-alpha.1
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

VERSION="${1:-}"
OUT_DIR="${ROOT}/dist"

usage() {
    cat <<EOF
Usage: $(basename "$0") VERSION

Build dist/EGSnrc_CLRP-egs_brachy-VERSION.tar.gz (end-user, no .git).

Examples:
  $(basename "$0") 1.0.0-alpha.1
  ./scripts/build-test-tarballs.sh   # alpha.1 + alpha.2 for install/update tests

Run from feature/eb-setup with submodules initialized:
  git submodule update --init --recursive
EOF
}

[[ -n "$VERSION" ]] || { usage; exit 1; }
[[ -f "$ROOT/eb-setup.sh" ]] || { echo "error: run from EGSnrc_CLRP repo root" >&2; exit 1; }
[[ -f "$ROOT/HEN_HOUSE/specs/unix.spec" ]] || { echo "error: HEN_HOUSE missing" >&2; exit 1; }
[[ -f "$ROOT/HEN_HOUSE/user_codes/egs_brachy/egs_brachy/Makefile" ]] || {
    echo "error: egs_brachy submodule not initialized (git submodule update --init --recursive)" >&2
    exit 1
}

DIR="EGSnrc_CLRP-egs_brachy-${VERSION}"
ARCHIVE="${OUT_DIR}/EGSnrc_CLRP-egs_brachy-${VERSION}.tar.gz"

mkdir -p "${OUT_DIR}/${DIR}"
rsync -a \
    --exclude='.git/' \
    --exclude='.github/' \
    --exclude='.cursor/' \
    --exclude='dist/' \
    --exclude='egs_home/' \
    --exclude='eb-env.sh' \
    --exclude='HEN_HOUSE/lib/' \
    --exclude='HEN_HOUSE/bin/' \
    --exclude='HEN_HOUSE/log/' \
    --exclude='HEN_HOUSE/egs++/dso/' \
    ./ "${OUT_DIR}/${DIR}/"

tar -czf "$ARCHIVE" -C "$OUT_DIR" "$DIR"
echo "Wrote $ARCHIVE"
echo "  $(du -h "$ARCHIVE" | cut -f1)  $(tar -tzf "$ARCHIVE" | wc -l | tr -d ' ') files"

#!/usr/bin/env bash
# eb-setup.sh — install, update, sync, and diagnose egs_brachy
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/lib/common.sh
source "$ROOT/scripts/lib/common.sh"
source "$ROOT/scripts/lib/detect.sh"
source "$ROOT/scripts/lib/dirty.sh"
source "$ROOT/scripts/lib/sync.sh"
source "$ROOT/scripts/lib/tarball.sh"
source "$ROOT/scripts/lib/git_ops.sh"

main() {
    parse_args "$@"
    resolve_paths
    case "$EB_CMD" in
        check)   preflight_tools; print_check_report ;;
        sync)    cmd_sync ;;
        update)  cmd_update ;;
        install) cmd_install ;;
        env)     cmd_env ;;
        help)    usage ;;
        *)       die "unknown command: $EB_CMD (try: help)" ;;
    esac
}

main "$@"

# shellcheck shell=bash
# Shared helpers for eb-setup.sh

EB_SETUP_VERSION="0.1.0-dev"

EB_CMD=""
FROM_TARBALL=0
INSTALL_DIR=""
TARBALL_PATH=""
TARBALL_TMPDIR=""
TARBALL_PAYLOAD_ROOT=""
EGS_HOME_OVERRIDE=""
NON_INTERACTIVE=0
EB_YES=0
DRY_RUN=0
STRICT=0
STASH=0

EB_ROOT=""
REPO_ROOT=""
HEN_HOUSE=""
EGS_CONFIG_RESOLVED=""
EGS_HOME_RESOLVED=""
MY_MACHINE=""
EB_SUBMODULE=""

# Mortran embeds HEN_HOUSE/EGS_CONFIG paths in machine.macros (%C80). Long clone
# paths make pegs4 fail with "FATAL STRING OR STATEMENT TOO LONG" (configure.log).
EB_MAX_REPO_ROOT_LEN=42

log()  { printf 'eb-setup: %s\n' "$*"; }
warn() { printf 'eb-setup: warning: %s\n' "$*" >&2; }
die()  { printf 'eb-setup: error: %s\n' "$*" >&2; exit "${2:-1}"; }

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

run() {
    if (( DRY_RUN )); then
        printf 'eb-setup: [dry-run] '
        printf '%q ' "$@"
        printf '\n'
        return 0
    fi
    "$@"
}

eb_setup_root() {
    local src="${BASH_SOURCE[1]:-${BASH_SOURCE[0]}}"
    while [[ -L "$src" ]]; do src="$(readlink "$src")"; done
    cd "$(dirname "$src")/../.." && pwd
}

parse_args() {
    EB_CMD="${1:-check}"
    shift || true
    if [[ "$EB_CMD" == "-h" || "$EB_CMD" == "--help" || "$EB_CMD" == "help" ]]; then
        EB_CMD="help"
        return 0
    fi
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --from-tarball)
                [[ -n "${2:-}" && "$2" != -* ]] || die "--from-tarball requires PATH to release .tar.gz"
                FROM_TARBALL=1
                TARBALL_PATH="$2"
                shift 2
                ;;
            --tarball)
                FROM_TARBALL=1
                TARBALL_PATH="${2:?}"
                shift 2
                ;;
            --install-dir)    INSTALL_DIR="${2:?}"; shift 2 ;;
            --egs-home)       EGS_HOME_OVERRIDE="${2:?}"; shift 2 ;;
            --non-interactive) NON_INTERACTIVE=1; shift ;;
            --yes)            EB_YES=1; shift ;;
            --dry-run)        DRY_RUN=1; shift ;;
            --strict)         STRICT=1; shift ;;
            --stash)          STASH=1; shift ;;
            -h|--help)        EB_CMD="help"; shift ;;
            *) die "unknown option: $1" ;;
        esac
    done
    if (( STRICT && STASH )); then die "--strict and --stash are mutually exclusive"; fi
}

usage() {
    cat <<EOF
eb-setup ${EB_SETUP_VERSION} — install, update, sync, and diagnose egs_brachy

Usage: eb-setup.sh <command> [options]

Commands: check | install | update | sync | env | help

Options: --from-tarball PATH --install-dir PATH --egs-home PATH
         --non-interactive --yes --dry-run --strict --stash

Tarball:  update --from-tarball release.tar.gz
          install --from-tarball release.tar.gz --install-dir PATH
          (--tarball PATH is an alias for --from-tarball PATH)

Testing: branch feature/eb-setup; scratch ~/Developer/scratch
Docs:    docs/eb-setup-testing.md  docs/branch-layout.md
EOF
}

_config_value() {
    local key="$1" file="$2"
    grep -E "^${key}[[:space:]]*=" "$file" 2>/dev/null | head -1 | sed -E "s/^${key}[[:space:]]*=[[:space:]]*//" | tr -d '\r'
}

# Expand leading ~ (bash does not expand ~ inside quoted --egs-home values).
expand_user_path() {
    local p="$1"
    if [[ "$p" == "~" ]]; then
        echo "$HOME"
    elif [[ "$p" == "~/"* ]]; then
        echo "${HOME}/${p:2}"
    else
        echo "$p"
    fi
}

# Read EGS_HOME recorded by configure/finalize (authoritative after install configure).
_egs_home_from_configure_logs() {
    local eh="" log log_dir="${HEN_HOUSE%/}/log"
    [[ -d "$log_dir" ]] || return 0
    for log in "$log_dir"/configure-*-"${USER}.log" "$log_dir"/configure-*.log; do
        [[ -f "$log" ]] || continue
        eh=$(grep -E '^EGS_HOME:[[:space:]]+' "$log" | head -1 \
            | sed -E 's/^EGS_HOME:[[:space:]]+//;s/[[:space:]]+$//')
        [[ -n "$eh" ]] && break
        eh=$(grep -E 'EGS_HOME[[:space:]]*=' "$log" | tail -1 \
            | sed -E 's/.*EGS_HOME[[:space:]]*=[[:space:]]*//;s/[[:space:]]+$//')
        [[ -n "$eh" ]] && break
    done
    [[ -n "$eh" ]] && expand_user_path "$eh"
}

resolve_paths() {
    EB_ROOT="$(eb_setup_root)"
    REPO_ROOT="$EB_ROOT"
    if [[ -n "${EGS_CONFIG:-}" && -f "${EGS_CONFIG}" ]]; then
        EGS_CONFIG_RESOLVED="$EGS_CONFIG"
        local hh="$(_config_value HEN_HOUSE "$EGS_CONFIG")"
        hh="${hh%/}"
        if [[ -n "$hh" && -d "$hh" ]]; then
            HEN_HOUSE="$hh"
            REPO_ROOT="$(cd "$HEN_HOUSE/.." && pwd)"
        fi
        MY_MACHINE="$(_config_value my_machine "$EGS_CONFIG")"
    elif [[ -d "$EB_ROOT/HEN_HOUSE" ]]; then
        HEN_HOUSE="$EB_ROOT/HEN_HOUSE"
    fi
    if [[ -n "$HEN_HOUSE" ]]; then EB_SUBMODULE="$HEN_HOUSE/user_codes/egs_brachy"; fi
    # After configure, finalize's log wins over --egs-home (user may have accepted a different path).
    local eh=""
    if [[ -n "$HEN_HOUSE" ]]; then
        eh="$(_egs_home_from_configure_logs)"
    fi
    if [[ -n "$eh" ]]; then
        EGS_HOME_RESOLVED="$eh"
    elif [[ -n "$EGS_HOME_OVERRIDE" ]]; then
        EGS_HOME_RESOLVED="$(expand_user_path "$EGS_HOME_OVERRIDE")"
    elif [[ -n "${EGS_HOME:-}" ]]; then
        EGS_HOME_RESOLVED="$(expand_user_path "$EGS_HOME")"
    elif [[ -f "$HOME/.egsnrcrc" ]]; then
        EGS_HOME_RESOLVED="$(expand_user_path "$(grep 'EGS_HOME' "$HOME/.egsnrcrc" | head -1 | sed -E 's/.*EGS_HOME[^=]*=[[:space:]]*//; s/[[:space:]]*$//')")"
    elif [[ -d "${REPO_ROOT}/egs_home" ]]; then
        EGS_HOME_RESOLVED="${REPO_ROOT}/egs_home/"
    fi
    normalize_egs_home
}

# EGSnrc makefiles use $(EGS_HOME)$(BIN_SUBDIR) — EGS_HOME must end with /
normalize_egs_home() {
    [[ -n "$EGS_HOME_RESOLVED" ]] || return 0
    EGS_HOME_RESOLVED="${EGS_HOME_RESOLVED%/}/"
}

eb_source_path() { echo "${HEN_HOUSE}/user_codes/egs_brachy/egs_brachy"; }
eb_dest_path()    { echo "${EGS_HOME_RESOLVED}egs_brachy"; }
eb_executable()   { echo "${EGS_HOME_RESOLVED}bin/${MY_MACHINE}/egs_brachy"; }

# Suggest shell profile for env snippet (matches finalize_egs_foruser logic).
shell_profile_hint() {
    local sh
    sh=$(basename "${SHELL:-bash}")
    case "$sh" in
        zsh)  echo "~/.zshrc" ;;
        bash) echo "~/.bashrc" ;;
        tcsh|csh) echo "~/.cshrc" ;;
        *)    echo "~/.profile or your shell startup file" ;;
    esac
}

# Print and save EGSnrc env lines (configure only prints these — does not edit the profile).
emit_shell_setup() {
    local profile snippet
    profile=$(shell_profile_hint)
    snippet="${REPO_ROOT:-.}/eb-env.sh"
    resolve_paths
    [[ -n "$EGS_CONFIG_RESOLVED" && -n "$EGS_HOME_RESOLVED" ]] || return 0

    cat >"$snippet" <<EOF
# egs_brachy environment — generated by eb-setup.sh
# Add to $(shell_profile_hint), or run:  source "$snippet"
export EGS_CONFIG="$EGS_CONFIG_RESOLVED"
export EGS_HOME="$EGS_HOME_RESOLVED"
source "${HEN_HOUSE%/}/scripts/egsnrc_bashrc_additions"
EOF

    echo
    echo "================================================================================"
    echo "  Shell setup (required for new terminals)"
    echo "================================================================================"
    echo
    echo "  configure/finalize do not edit your profile. Add these lines to $profile"
    echo "  (or source the saved file below):"
    echo
    sed 's/^/    /' "$snippet"
    echo
    echo "  Saved: $snippet"
    echo "  Re-print anytime:  ./eb-setup.sh env"
    echo
    echo "================================================================================"
    echo
}

export_egs_env() {
    normalize_egs_home
    [[ -n "$EGS_CONFIG_RESOLVED" ]] && export EGS_CONFIG="$EGS_CONFIG_RESOLVED"
    [[ -n "$EGS_HOME_RESOLVED" ]]    && export EGS_HOME="$EGS_HOME_RESOLVED"
    [[ -n "$HEN_HOUSE" ]]            && export HEN_HOUSE="${HEN_HOUSE%/}/"
}

preflight_tools() {
    need_cmd bash; need_cmd make; need_cmd rsync
    (( FROM_TARBALL )) || need_cmd git
}

# Write HEN_HOUSE/specs/release.mk from git SHAs (git installs; skipped on tarball trees without .git).
write_release_mk() {
    [[ -n "$HEN_HOUSE" ]] || return 0
    (( DRY_RUN )) && return 0
    local spec_file="${HEN_HOUSE%/}/specs/release.mk"
    local clrp="" eb="" rel="" have_git=0

    if [[ -n "$REPO_ROOT" ]] && git -C "$REPO_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        have_git=1
        clrp="$(git -C "$REPO_ROOT" rev-parse --short=7 HEAD)"
        rel="$(git -C "$REPO_ROOT" describe --tags --abbrev=0 2>/dev/null \
            | sed -e 's/^v//' -e 's/^egs_brachy-//' || true)"
    fi
    if [[ -n "${EB_SUBMODULE:-}" ]] && git -C "$EB_SUBMODULE" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        have_git=1
        eb="$(git -C "$EB_SUBMODULE" rev-parse --short=7 HEAD)"
    fi
    if (( ! have_git )); then
        [[ -f "$spec_file" ]] && return 0
        return 0
    fi
    mkdir -p "${HEN_HOUSE%/}/specs"
    {
        echo "# Generated by eb-setup.sh — do not edit."
        [[ -n "$rel" ]]  && printf 'EGS_RELEASE = -DEGS_RELEASE="\\"%s\\""\n' "$rel"
        [[ -n "$clrp" ]] && printf 'EGS_CLRP_HASH = -DEGS_CLRP_HASH="\\"%s\\""\n' "$clrp"
        [[ -n "$eb" ]]   && printf 'EGS_BRACHY_HASH = -DEGS_BRACHY_HASH="\\"%s\\""\n' "$eb"
    } >"$spec_file"
    log "wrote $spec_file (CLRP=${clrp:-?} egs_brachy=${eb:-?})"
}

# configure/finalize run in a subshell — pick up paths they created for this session.
pickup_egsnrc_env_after_configure() {
    local spec_dir="${HEN_HOUSE%/}/specs" f m newest="" newest_m=0
    for f in "$spec_dir"/*.conf; do
        [[ -f "$f" ]] || continue
        grep -q '^my_machine[[:space:]]*=' "$f" || continue
        m=$(stat -f %m "$f" 2>/dev/null || stat -c %Y "$f")
        if (( m > newest_m )); then newest_m=$m; newest=$f; fi
    done
    [[ -n "$newest" ]] || die "no EGSnrc config found in $spec_dir after configure"
    EGS_CONFIG_RESOLVED="$newest"
    MY_MACHINE="$(_config_value my_machine "$newest")"

    local eh=""
    eh="$(_egs_home_from_configure_logs)"
    if [[ -n "$eh" ]]; then
        EGS_HOME_RESOLVED="$eh"
    elif [[ -n "$EGS_HOME_OVERRIDE" ]]; then
        EGS_HOME_RESOLVED="$(expand_user_path "$EGS_HOME_OVERRIDE")"
    elif [[ -d "${REPO_ROOT}/egs_home" ]]; then
        EGS_HOME_RESOLVED="${REPO_ROOT}/egs_home/"
    fi
    normalize_egs_home
    [[ -n "$EGS_HOME_RESOLVED" ]] || die "could not determine EGS_HOME - use --egs-home PATH"
    log "detected EGS_CONFIG=$EGS_CONFIG_RESOLVED"
    log "detected EGS_HOME=$EGS_HOME_RESOLVED"
}

# Pegs4/mortran fails when absolute paths baked into machine.macros are too long.
check_mortran_path_lengths() {
    [[ -n "$REPO_ROOT" ]] || return 0
    local n=${#REPO_ROOT}
    if (( n > EB_MAX_REPO_ROOT_LEN )); then
        warn "repo path is ${n} chars (limit ~${EB_MAX_REPO_ROOT_LEN} for Mortran):"
        warn "  $REPO_ROOT"
        warn "Use a shorter clone dir, e.g.:"
        warn "  --install-dir \"\$HOME/scratch/eb\""
        warn "  ln -s \"\$PWD\" \"\$HOME/scratch/eb\" && cd \"\$HOME/scratch/eb\""
        if [[ "$EB_CMD" == "install" ]]; then
            die "path too long for EGSnrc configure (see configure.log pegs4 / Mortran stop 12)"
        fi
    fi
}

# shellcheck shell=bash
# Release tarball install/update (no git). See .github/workflows/release-tarball.yml.

# Slim bootstrap tree: eb-setup.sh + scripts only (no HEN_HOUSE).
is_eb_setup_bootstrap_tree() {
    [[ -n "${REPO_ROOT:-}" ]] || return 1
    [[ -f "$REPO_ROOT/eb-setup.sh" ]] || return 1
    [[ -d "$REPO_ROOT/scripts/lib" ]] || return 1
    [[ ! -f "$REPO_ROOT/HEN_HOUSE/specs/unix.spec" ]] || return 1
}

# End-user release tree: extracted tarball with eb-setup.sh, no .git at repo root.
is_release_install_tree() {
    [[ -n "${REPO_ROOT:-}" ]] || return 1
    [[ -f "$REPO_ROOT/eb-setup.sh" ]] || return 1
    [[ -f "$REPO_ROOT/HEN_HOUSE/specs/unix.spec" ]] || return 1
    [[ ! -d "$REPO_ROOT/.git" ]] || return 1
}

# Where bootstrap install placed (or will place) the full release tree.
resolve_install_root() {
    local f ir default
    default="$(eb_default_install_dir)"
    f="${EB_ROOT}/.eb-setup/install-root"
    if [[ -f "$f" ]]; then
        ir="$(expand_user_path "$(head -1 "$f")")"
        [[ -d "$ir/HEN_HOUSE" ]] && { echo "$ir"; return 0; }
    fi
    if [[ -d "$default/HEN_HOUSE" ]]; then
        echo "$default"
        return 0
    fi
    return 1
}

write_install_root() {
    local root="$1"
    mkdir -p "${EB_ROOT}/.eb-setup"
    printf '%s\n' "$(expand_user_path "$root")" > "${EB_ROOT}/.eb-setup/install-root"
}

install_dest_taken() {
    local dest="$1"
    [[ -f "$dest/eb-setup.sh" || -d "$dest/HEN_HOUSE" ]]
}

_tarball_cache_root() {
    local root=""
    root="$(resolve_install_root 2>/dev/null || true)"
    [[ -n "$root" ]] && { echo "$root"; return 0; }
    if is_eb_setup_bootstrap_tree; then
        echo "${EB_ROOT}"
        return 0
    fi
    echo "${REPO_ROOT:-${HOME}/.cache/eb-setup}"
}

_read_installed_release_version() {
    local spec="${HEN_HOUSE%/}/specs/release.mk"
    [[ -f "$spec" ]] || return 0
    grep -E '^EGS_RELEASE[[:space:]]*=' "$spec" 2>/dev/null | head -1 \
        | sed -E 's/^EGS_RELEASE[[:space:]]*=[[:space:]]*-DEGS_RELEASE="\\"([^"\\]+)\\""/\1/'
}

# Pick end-user tarball asset URL from GitHub releases JSON (read stdin).
_tarball_pick_asset_url() {
    local repo="$1" tag="${2:-}" json
    json="$(cat)"
    need_cmd python3
    JSON_PAYLOAD="$json" python3 - "$repo" "$tag" <<'PY'
import json, os, re, sys

repo = sys.argv[1]
tag = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] else ""
data = json.loads(os.environ["JSON_PAYLOAD"])
pattern = re.compile(r"^EGSnrc_CLRP-egs_brachy-(.+)\.tar\.gz$")

def asset_url(release):
    for asset in release.get("assets") or []:
        name = asset.get("name") or ""
        if "-src" in name:
            continue
        m = pattern.match(name)
        if m:
            return asset.get("browser_download_url") or "", m.group(1)
    return "", ""

if tag:
    want = tag if tag.startswith("egs_brachy-") else f"egs_brachy-{tag}"
    for rel in data if isinstance(data, list) else [data]:
        if rel.get("draft"):
            continue
        if rel.get("tag_name") == want:
            url, ver = asset_url(rel)
            if url:
                print(url)
                print(ver)
                sys.exit(0)
    sys.exit(1)

for rel in data if isinstance(data, list) else []:
    if rel.get("draft"):
        continue
    url, ver = asset_url(rel)
    if url:
        print(url)
        print(ver)
        sys.exit(0)
sys.exit(1)
PY
}

tarball_fetch_release_info() {
    local repo="${EB_RELEASE_REPO:-clrp-code/EGSnrc_CLRP}"
    local tag="${EB_RELEASE_TAG:-}"
    need_cmd curl
    local api_out pick_out url version
    if [[ -n "$tag" ]]; then
        local want="$tag"
        [[ "$want" == egs_brachy-* ]] || want="egs_brachy-$want"
        api_out="$(curl -fsSL "https://api.github.com/repos/${repo}/releases/tags/${want}")" \
            || die "GitHub release not found: ${want} (set EB_RELEASE_TAG or use --from-tarball)"
    else
        api_out="$(curl -fsSL "https://api.github.com/repos/${repo}/releases?per_page=30")" \
            || die "could not fetch GitHub releases for ${repo} (network? use --from-tarball PATH)"
    fi
    pick_out="$(printf '%s' "$api_out" | _tarball_pick_asset_url "$repo" "$tag")" \
        || die "no end-user tarball asset on GitHub release (expected EGSnrc_CLRP-egs_brachy-*.tar.gz)"
    url="$(sed -n '1p' <<<"$pick_out")"
    version="$(sed -n '2p' <<<"$pick_out")"
    printf '%s\n%s\n' "$url" "$version"
}

# Download latest (or EB_RELEASE_TAG / EB_RELEASE_URL) end-user tarball. Returns 1 if already current.
tarball_download_release() {
    local url version dest cache_dir fname installed info
    if [[ -n "${EB_RELEASE_URL:-}" ]]; then
        url="$EB_RELEASE_URL"
        version="$(basename "$url" .tar.gz | sed 's/^EGSnrc_CLRP-egs_brachy-//')"
    else
        info="$(tarball_fetch_release_info)"
        url="$(sed -n '1p' <<<"$info")"
        version="$(sed -n '2p' <<<"$info")"
    fi
    installed="$(_read_installed_release_version)"
    if [[ -n "$installed" && "$installed" == "$version" ]]; then
        log "already at release $version (HEN_HOUSE/specs/release.mk)"
        return 1
    fi
    cache_dir="$(_tarball_cache_root)/.eb-setup/cache"
    mkdir -p "$cache_dir"
    fname="EGSnrc_CLRP-egs_brachy-${version}.tar.gz"
    dest="$cache_dir/$fname"
    log "downloading release ${version}..."
    if (( DRY_RUN )); then
        printf 'eb-setup: [dry-run] curl -fsSL -o %q %q\n' "$dest" "$url"
        TARBALL_PATH="$dest"
        return 0
    fi
    run curl -fsSL -o "$dest" "$url"
    TARBALL_PATH="$dest"
    log "saved $dest"
}

tarball_require_path() {
    [[ -n "$TARBALL_PATH" ]] || die "--from-tarball PATH is required (path to release .tar.gz)"
    TARBALL_PATH="$(expand_user_path "$TARBALL_PATH")"
    [[ -f "$TARBALL_PATH" ]] || (( DRY_RUN )) || die "tarball not found: $TARBALL_PATH"
    case "$TARBALL_PATH" in
        *.tar.gz|*.tgz) ;;
        *) die "tarball must be .tar.gz or .tgz: $TARBALL_PATH" ;;
    esac
}

# Find top-level extract dir containing HEN_HOUSE/eb-setup.sh.
tarball_find_payload_root() {
    local extract_dir="$1" d
    if [[ -f "$extract_dir/HEN_HOUSE/specs/unix.spec" && -f "$extract_dir/eb-setup.sh" ]]; then
        echo "$extract_dir"
        return 0
    fi
    for d in "$extract_dir"/*; do
        [[ -d "$d" ]] || continue
        if [[ -f "$d/HEN_HOUSE/specs/unix.spec" && -f "$d/eb-setup.sh" ]]; then
            echo "$d"
            return 0
        fi
    done
    return 1
}

tarball_extract_to_temp() {
    local archive="$1" tmp
    need_cmd tar
    tmp="$(mktemp -d "${TMPDIR:-/tmp}/eb-setup-tarball.XXXXXX")"
    TARBALL_TMPDIR="$tmp"
    TARBALL_PAYLOAD_ROOT=""
    log "extracting $(basename "$archive")..."
    if (( DRY_RUN )); then
        printf 'eb-setup: [dry-run] tar -xzf %q -C %q\n' "$archive" "$tmp"
        TARBALL_PAYLOAD_ROOT="$tmp/${EB_INSTALL_TREE_NAME:-EGSnrc_CLRP}"
        return 0
    fi
    run tar -xzf "$archive" -C "$tmp"
    TARBALL_PAYLOAD_ROOT="$(tarball_find_payload_root "$tmp")" \
        || die "tarball does not look like an EGSnrc_CLRP release (missing HEN_HOUSE/specs/unix.spec)"
}

tarball_cleanup_temp() {
    [[ -n "${TARBALL_TMPDIR:-}" && -d "$TARBALL_TMPDIR" ]] || return 0
    if (( DRY_RUN )); then return 0; fi
    rm -rf "$TARBALL_TMPDIR"
    TARBALL_TMPDIR=""
}

# Rsync release tree into an existing install, keeping configure/build artifacts.
tarball_apply_to_repo() {
    local src="${TARBALL_PAYLOAD_ROOT:?}/" dst="${REPO_ROOT:?}/"
    local -a opts=(-a)
    (( DRY_RUN )) && opts+=(-n -i)
    log "merging release into $REPO_ROOT..."
    run rsync "${opts[@]}" \
        --exclude='egs_home/' \
        --exclude='eb-env.sh' \
        --exclude='HEN_HOUSE/lib/' \
        --exclude='HEN_HOUSE/bin/' \
        --exclude='HEN_HOUSE/log/' \
        --exclude='HEN_HOUSE/egs++/dso/' \
        --exclude='HEN_HOUSE/specs/' \
        "$src" "$dst"
    # Add new spec templates only; never overwrite user *.conf from configure.
    if [[ -d "$src/HEN_HOUSE/specs" ]]; then
        run rsync "${opts[@]}" --ignore-existing \
            "$src/HEN_HOUSE/specs/" "${dst}HEN_HOUSE/specs/"
    fi
}

tarball_extract_install_dir() {
    local archive="$1" dest="$2"
    dest="$(expand_user_path "$dest")"
    tarball_extract_to_temp "$archive"
    if (( DRY_RUN )); then
        printf 'eb-setup: [dry-run] rsync release to %q\n' "$dest"
        return 0
    fi
    mkdir -p "$dest"
    run rsync -a "${TARBALL_PAYLOAD_ROOT}/" "$dest/"
    REPO_ROOT="$dest"
    HEN_HOUSE="${dest%/}/HEN_HOUSE"
    EB_SUBMODULE="${HEN_HOUSE}/user_codes/egs_brachy"
}

cmd_update_from_tarball() {
    tarball_require_path
    [[ -n "$EGS_CONFIG_RESOLVED" && -f "$EGS_CONFIG_RESOLVED" ]] \
        || die "EGS_CONFIG must be set for tarball update"
    [[ -n "$HEN_HOUSE" && -d "$HEN_HOUSE" ]] || die "HEN_HOUSE not found"
    if [[ -d "$REPO_ROOT/.git" ]]; then
        die "git checkout detected — use: eb-setup.sh update (without --from-tarball)"
    fi
    if (( STASH || STRICT )); then
        warn "--stash/--strict apply to git updates only; ignored for --from-tarball"
    fi
    tarball_extract_to_temp "$TARBALL_PATH"
    trap tarball_cleanup_temp EXIT
    tarball_apply_to_repo
    if (( DRY_RUN )); then
        log "dry-run complete (skipped egs++ rebuild and sync)"
        return 0
    fi
    if [[ -f "${TARBALL_PAYLOAD_ROOT}/HEN_HOUSE/specs/release.mk" ]]; then
        mkdir -p "${HEN_HOUSE%/}/specs"
        run cp "${TARBALL_PAYLOAD_ROOT}/HEN_HOUSE/specs/release.mk" \
            "${HEN_HOUSE%/}/specs/release.mk"
        log "updated HEN_HOUSE/specs/release.mk from release tarball"
    fi
    log "rebuilding egs++..."
    run make -C "$HEN_HOUSE/egs++"
    sync_egs_brachy_from_henhouse
    ensure_egs_home_bin
    export_egs_env
    run make -C "$(eb_dest_path)"
    log "tarball update complete"
    emit_shell_setup
}

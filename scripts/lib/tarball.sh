# shellcheck shell=bash
# Release tarball install/update (no git). See .github/workflows/release-tarball.yml.

tarball_require_path() {
    [[ -n "$TARBALL_PATH" ]] || die "--from-tarball PATH is required (path to release .tar.gz)"
    TARBALL_PATH="$(expand_user_path "$TARBALL_PATH")"
    [[ -f "$TARBALL_PATH" ]] || die "tarball not found: $TARBALL_PATH"
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
        TARBALL_PAYLOAD_ROOT="$tmp/EGSnrc_CLRP-egs_brachy"
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
    local archive="$1" dest
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
    log "rebuilding egs++..."
    run make -C "$HEN_HOUSE/egs++"
    sync_egs_brachy_from_henhouse
    ensure_egs_home_bin
    export_egs_env
    run make -C "$(eb_dest_path)"
    log "tarball update complete"
    emit_shell_setup
}

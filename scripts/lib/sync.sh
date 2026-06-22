# shellcheck shell=bash
rsync_preview_legend() {
    echo "(rsync preview — >f..t.... = file will be copied/updated from submodule → EGS_HOME)"
}

confirm_overwrites() {
    local src dst diff n
    src="$(eb_source_path)/"; dst="$(eb_dest_path)/"
    [[ -d "$dst" ]] || return 0
    diff="$(rsync -ani --exclude='eb_tests/**/results/' "$src" "$dst" 2>/dev/null | grep -E '^[<>*]' || true)"
    [[ -z "$diff" ]] && return 0
    n=$(echo "$diff" | wc -l | tr -d ' ')
    rsync_preview_legend
    echo "$n file(s) would change (showing up to 30):"
    echo "$diff" | head -30 | sed -E 's/^>[fh][.cpstoguax]*[[:space:]]+/  → /'
    if (( EB_YES )); then return 0; fi
    read -r -p "Proceed with sync? [y/N] " ans
    [[ "$ans" == [yY] || "$ans" == [yY][eE][sS] ]] || die "sync cancelled"
}

sync_egs_brachy_from_henhouse() {
    local src dst rsync_opts=(-a --backup --suffix=".pre-sync-$(date +%Y%m%d)")
    [[ -n "$HEN_HOUSE" && -n "$EGS_HOME_RESOLVED" ]] || die "HEN_HOUSE and EGS_HOME must be set"
    src="$(eb_source_path)/"; dst="$(eb_dest_path)/"
    [[ -d "$src" ]] || die "source not found: $src"
    mkdir -p "$dst"
    if (( DRY_RUN )); then
        run rsync -ani --exclude='eb_tests/**/results/' "$src" "$dst"
        return 0
    fi
    confirm_overwrites
    rsync_opts+=(--exclude='eb_tests/**/results/')
    if (( EB_YES )); then rsync_opts+=(--delete-after); fi
    log "syncing egs_brachy to EGS_HOME..."
    run rsync "${rsync_opts[@]}" "$src" "$dst"
}

cmd_sync() {
    resolve_paths; preflight_tools
    sync_egs_brachy_from_henhouse
    if (( DRY_RUN )); then
        log "dry-run complete (skipped make)"
        return 0
    fi
    ensure_egs_home_bin
    export_egs_env
    log "compiling egs_brachy..."
    run make -C "$(eb_dest_path)"
    log "sync complete"
}

ensure_egs_home_bin() {
    [[ -n "$EGS_HOME_RESOLVED" && -n "$MY_MACHINE" ]] || return 0
    mkdir -p "${EGS_HOME_RESOLVED}bin/${MY_MACHINE}"
}

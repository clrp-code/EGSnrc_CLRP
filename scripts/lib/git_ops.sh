# shellcheck shell=bash
run_egsnrc_configure() {
    local cfg_dir="${HEN_HOUSE%/}/scripts"
    local script=configure
    (( NON_INTERACTIVE )) && script=configure.expect
    [[ -f "$cfg_dir/$script" ]] || die "not found: $cfg_dir/$script"
    # configure sets my_dir=$(pwd)/$(dirname $0). Must invoke as ./configure
    # from scripts/ in this shell — bash -c 'exec ./configure' leaves $0 absolute on macOS.
    if (( DRY_RUN )); then
        printf 'eb-setup: [dry-run] (cd %q && unset HEN_HOUSE EGS_HOME EGS_CONFIG && ./%s)\n' \
            "$cfg_dir" "$script"
        return 0
    fi
    local cfg_status=0
    (
        unset HEN_HOUSE EGS_HOME EGS_CONFIG
        cd "$cfg_dir" || exit 1
        ./"$script"
    ) || cfg_status=$?
    (( cfg_status == 0 )) || die "EGSnrc configure failed (see HEN_HOUSE/log/configure*.log)"
}

cmd_update() {
    resolve_paths; preflight_tools
    if (( FROM_TARBALL )); then
        cmd_update_from_tarball
        return
    fi
    [[ -n "$REPO_ROOT" ]] || die "cannot find EGSnrc repo root"
    if (( STASH )); then stash_git_repos; fi
    enforce_tier1_git_policy
    log "pulling EGSnrc (egs_brachy branch)..."
    run git -C "$REPO_ROOT" pull origin egs_brachy || \
        run git -C "$REPO_ROOT" pull clrp egs_brachy || \
        warn "could not pull egs_brachy — on branch $(_git_branch "$REPO_ROOT")?"
    if [[ -d "$EB_SUBMODULE/.git" || -f "$EB_SUBMODULE/.git" ]]; then
        log "pulling egs_brachy submodule..."
        run git -C "$EB_SUBMODULE" pull origin main
    fi
    log "rebuilding egs++..."
    run make -C "$HEN_HOUSE/egs++"
    sync_egs_brachy_from_henhouse
    ensure_egs_home_bin
    export_egs_env
    run make -C "$(eb_dest_path)"
    if (( STASH )); then pop_git_stashes; fi
    log "update complete"
    emit_shell_setup
}

cmd_install() {
    resolve_paths; preflight_tools; check_mortran_path_lengths
    if (( FROM_TARBALL )); then
        tarball_require_path
        local dest
        dest="$(expand_user_path "${INSTALL_DIR:-$HOME/scratch/eb}")"
        if [[ -f "$dest/eb-setup.sh" || -d "$dest/HEN_HOUSE" ]]; then
            die "install dir already exists: $dest (use: update --from-tarball PATH)"
        fi
        tarball_extract_install_dir "$TARBALL_PATH" "$dest"
        trap tarball_cleanup_temp EXIT
    else
        if [[ ! -d "$REPO_ROOT/.git" ]]; then
            local dest
            dest="$(expand_user_path "${INSTALL_DIR:-$HOME/scratch/eb}")"
            log "cloning CLRP fork to $dest..."
            run git clone https://github.com/clrp-code/EGSnrc_CLRP.git "$dest"
            REPO_ROOT="$dest"; HEN_HOUSE="$dest/HEN_HOUSE"; EB_SUBMODULE="$HEN_HOUSE/user_codes/egs_brachy"
        fi
        local cur_branch="$(_git_branch "$REPO_ROOT")"
        if [[ "$cur_branch" == "egs_brachy" ]]; then
            log "on egs_brachy branch"
        elif [[ -f "$REPO_ROOT/eb-setup.sh" ]]; then
            warn "staying on branch $cur_branch (eb-setup testing branch)"
        else
            log "checking out egs_brachy branch..."
            run git -C "$REPO_ROOT" checkout egs_brachy
        fi
        log "initializing egs_brachy submodule..."
        run git -C "$REPO_ROOT" submodule update --init --recursive
    fi
    if [[ -z "$EGS_CONFIG_RESOLVED" ]]; then
        log "configure EGSnrc (interactive)..."
        run_egsnrc_configure
        pickup_egsnrc_env_after_configure
        export_egs_env
        log "configure complete — continuing with sync"
    fi
    cmd_sync
    emit_shell_setup
    log "install complete"
}

cmd_env() {
    resolve_paths
    emit_shell_setup
}

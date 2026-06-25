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
        # --egs-home becomes finalize's default suggestion (still overridable interactively).
        if [[ -n "$EGS_HOME_OVERRIDE" ]]; then
            export EGS_HOME="$(expand_user_path "$EGS_HOME_OVERRIDE")"
        fi
        cd "$cfg_dir" || exit 1
        ./"$script"
    ) || cfg_status=$?
    (( cfg_status == 0 )) || die "EGSnrc configure failed (see HEN_HOUSE/log/configure*.log)"
}

_git_checkout_and_submodule() {
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
}

_cmd_install_from_tarball_file() {
    local dest
    dest="$(eb_install_dest)"
    if install_dest_taken "$dest"; then
        die "install dir already exists: $dest (use: eb-setup.sh update)"
    fi
    check_mortran_path_lengths_for "$dest"
    tarball_require_path
    tarball_extract_install_dir "$TARBALL_PATH" "$dest"
    trap tarball_cleanup_temp EXIT
    if is_eb_setup_bootstrap_tree; then
        write_install_root "$dest"
    fi
}

_cmd_install_bootstrap() {
    local dest bootstrap_root="$EB_ROOT"
    dest="$(eb_install_dest)"
    if install_dest_taken "$dest"; then
        die "install dir already exists: $dest (use: eb-setup.sh update)"
    fi
    check_mortran_path_lengths_for "$dest"
    if (( INSTALL_GIT )); then
        log "cloning CLRP fork to $dest..."
        run git clone https://github.com/clrp-code/EGSnrc_CLRP.git "$dest"
        REPO_ROOT="$dest"
        HEN_HOUSE="${dest%/}/HEN_HOUSE"
        EB_SUBMODULE="${HEN_HOUSE}/user_codes/egs_brachy"
        _git_checkout_and_submodule
    else
        log "bootstrap installer — fetching release to $dest"
        tarball_download_release \
            || die "could not download release (network? use: install --from-tarball PATH)"
        tarball_extract_install_dir "$TARBALL_PATH" "$dest"
        trap tarball_cleanup_temp EXIT
    fi
    write_install_root "$dest"
    REPO_ROOT="$dest"
    HEN_HOUSE="${dest%/}/HEN_HOUSE"
    EB_SUBMODULE="${HEN_HOUSE}/user_codes/egs_brachy"
    if [[ "$bootstrap_root" != "$REPO_ROOT" ]]; then
        log "for day-to-day use:  cd $REPO_ROOT && source ./eb-env.sh"
    fi
}

_cmd_install_release_inplace() {
    local dest repo
    dest="$(eb_install_dest)"
    repo="$(cd "$REPO_ROOT" && pwd -P)"
    dest="$(expand_user_path "$dest")"
    mkdir -p "$dest" 2>/dev/null || true
    dest="$(cd "$dest" && pwd -P)"
    if [[ "$repo" == "$dest" ]]; then
        log "release tree detected — in-place install"
    else
        if install_dest_taken "$dest"; then
            die "install dir already exists: $dest (use: eb-setup.sh update)"
        fi
        log "installing release tree to $dest"
        check_mortran_path_lengths_for "$dest"
        run rsync -a "${REPO_ROOT%/}/" "${dest%/}/"
        REPO_ROOT="$dest"
        HEN_HOUSE="${dest%/}/HEN_HOUSE"
        EB_SUBMODULE="${HEN_HOUSE}/user_codes/egs_brachy"
        log "for day-to-day use:  cd $REPO_ROOT && source ./eb-env.sh"
    fi
    check_mortran_path_lengths
}

_cmd_install_git_clone() {
    local dest
    dest="$(eb_install_dest)"
    if install_dest_taken "$dest"; then
        die "install dir already exists: $dest"
    fi
    check_mortran_path_lengths_for "$dest"
    log "cloning CLRP fork to $dest..."
    run git clone https://github.com/clrp-code/EGSnrc_CLRP.git "$dest"
    REPO_ROOT="$dest"
    HEN_HOUSE="${dest%/}/HEN_HOUSE"
    EB_SUBMODULE="${HEN_HOUSE}/user_codes/egs_brachy"
    _git_checkout_and_submodule
}

_cmd_install_git_tree() {
    check_mortran_path_lengths
    _git_checkout_and_submodule
}

_cmd_install_finish() {
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

cmd_update() {
    resolve_paths; preflight_tools
    if is_eb_setup_bootstrap_tree && [[ -z "${HEN_HOUSE:-}" || ! -d "$HEN_HOUSE" ]]; then
        die "no installation found — run: eb-setup.sh install (default: $(eb_default_install_dir))"
    fi
    if (( FROM_TARBALL )); then
        cmd_update_from_tarball
        return
    fi
    if is_release_install_tree; then
        if tarball_download_release; then
            FROM_TARBALL=1
            cmd_update_from_tarball
        else
            log "update complete (no newer release)"
        fi
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
    resolve_paths; preflight_tools

    if (( FROM_TARBALL )); then
        _cmd_install_from_tarball_file
    elif is_eb_setup_bootstrap_tree; then
        _cmd_install_bootstrap
    elif is_release_install_tree; then
        _cmd_install_release_inplace
    elif [[ -d "$REPO_ROOT/.git" ]]; then
        _cmd_install_git_tree
    elif (( INSTALL_GIT )); then
        _cmd_install_git_clone
    else
        die "cannot install from this location (extract eb-setup tarball, or use install --git)"
    fi

    _cmd_install_finish
}

cmd_env() {
    resolve_paths
    emit_shell_setup
}

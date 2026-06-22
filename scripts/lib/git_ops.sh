# shellcheck shell=bash
cmd_update() {
    resolve_paths; preflight_tools
    (( FROM_TARBALL )) && die "update requires git checkout (not --from-tarball)"
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
}

cmd_install() {
    resolve_paths; preflight_tools
    if (( ! FROM_TARBALL )); then
        if [[ ! -d "$REPO_ROOT/.git" ]]; then
            local dest="${INSTALL_DIR:-$HOME/Developer/scratch/EGSnrc_CLRP}"
            log "cloning CLRP fork to $dest..."
            run git clone https://github.com/clrp-code/EGSnrc_CLRP.git "$dest"
            REPO_ROOT="$dest"; HEN_HOUSE="$dest/HEN_HOUSE"; EB_SUBMODULE="$HEN_HOUSE/user_codes/egs_brachy"
        fi
        log "checking out egs_brachy branch and submodule..."
        run git -C "$REPO_ROOT" checkout egs_brachy
        run git -C "$REPO_ROOT" submodule update --init --recursive
    fi
    if [[ -z "$EGS_CONFIG_RESOLVED" ]]; then
        log "configure EGSnrc (interactive)..."
        warn "unset HEN_HOUSE EGS_HOME EGS_CONFIG if switching installs"
        if (( NON_INTERACTIVE )); then
            run "$HEN_HOUSE/scripts/configure.expect"
        else
            run "$HEN_HOUSE/scripts/configure"
        fi
        die "after configure, set EGS_CONFIG/EGS_HOME and re-run: eb-setup.sh sync"
    fi
    cmd_sync
    cmd_env
}

cmd_env() {
    resolve_paths
    echo "# Add to ~/.bashrc or ~/.zshrc:"
    echo "# EGS_HOME must end with / (EGSnrc uses \$(EGS_HOME)bin/...)"
    [[ -n "$EGS_CONFIG_RESOLVED" ]] && echo "export EGS_CONFIG=\"$EGS_CONFIG_RESOLVED\""
    [[ -n "$EGS_HOME_RESOLVED" ]]    && echo "export EGS_HOME=\"$EGS_HOME_RESOLVED\""
    echo "# optional CLRP aliases (exeb, cdeb):"
    [[ -n "$HEN_HOUSE" ]] && echo "source \"${HEN_HOUSE%/}/scripts/clrp_bashrc_additions\""
}

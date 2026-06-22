# shellcheck shell=bash
run_egsnrc_configure() {
    local cfg_dir="${HEN_HOUSE%/}/scripts"
    local script=configure
    (( NON_INTERACTIVE )) && script=configure.expect
    [[ -f "$cfg_dir/$script" ]] || die "not found: $cfg_dir/$script"
    # EGSnrc configure sets my_dir=$(pwd)/$(dirname $0); run ./configure from scripts/.
    run env -u HEN_HOUSE -u EGS_HOME -u EGS_CONFIG \
        bash -c 'cd "$1" && exec ./"$2"' _ "$cfg_dir" "$script"
}

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
            local dest
            dest="$(expand_user_path "${INSTALL_DIR:-$HOME/Developer/scratch/EGSnrc_CLRP}")"
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

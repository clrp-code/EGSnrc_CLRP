# shellcheck shell=bash
SCENARIO=""; SCENARIO_DESC=""

_git_branch() { git -C "$1" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown; }
_git_remote_url() {
    git -C "$1" remote get-url origin 2>/dev/null || git -C "$1" remote get-url clrp 2>/dev/null || echo unknown
}
_submodule_ok() { [[ -d "$EB_SUBMODULE/egs_brachy" && -f "$EB_SUBMODULE/egs_brachy/Makefile" ]]; }
_is_clrp_remote() { [[ "$1" == *EGSnrc_CLRP* || "$1" == *clrp-code* ]]; }
_is_nrc_remote() { [[ "$1" == *nrc-cnrc/EGSnrc* ]]; }

detect_scenario() {
    if [[ -z "$HEN_HOUSE" || ! -d "$HEN_HOUSE" ]]; then
        if is_eb_setup_bootstrap_tree 2>/dev/null; then
            SCENARIO="A0"
            SCENARIO_DESC="Bootstrap installer — run install to fetch release"
            return
        fi
        SCENARIO="A"; SCENARIO_DESC="Greenfield — no EGSnrc detected"; return; fi
    if [[ -n "$REPO_ROOT" ]] && git -C "$REPO_ROOT" rev-parse --git-dir >/dev/null 2>&1; then
        local remote="$(_git_remote_url "$REPO_ROOT")"
        if _is_nrc_remote "$remote" && ! _is_clrp_remote "$remote"; then
            SCENARIO="B"; SCENARIO_DESC="NRC EGSnrc detected (not CLRP fork)"; return; fi
        if ! _submodule_ok; then
            SCENARIO="C"; SCENARIO_DESC="CLRP fork present but egs_brachy submodule incomplete"; return; fi
    elif (( FROM_TARBALL )) && ! _submodule_ok; then
        SCENARIO="C"; SCENARIO_DESC="Tarball missing egs_brachy sources"; return; fi
    if [[ -z "$EGS_CONFIG_RESOLVED" || ! -f "$EGS_CONFIG_RESOLVED" ]]; then
        SCENARIO="D"; SCENARIO_DESC="Sources ready but EGSnrc not configured"; return; fi
    if [[ ! -f "$(eb_dest_path)/Makefile" ]]; then
        SCENARIO="E"; SCENARIO_DESC="Configured but egs_brachy not in EGS_HOME"; return; fi
    if [[ -x "$(eb_executable)" ]]; then
        SCENARIO="F"; SCENARIO_DESC="Fully installed"; return; fi
    SCENARIO="E"; SCENARIO_DESC="egs_brachy in EGS_HOME but not built"
}

print_check_report() {
    log "checking egs_brachy installation..."; echo
    echo "Environment:"
    echo "  EGS_CONFIG   ${EGS_CONFIG_RESOLVED:-not set}"
    echo "  HEN_HOUSE    ${HEN_HOUSE:-not found}"
    echo "  EGS_HOME     ${EGS_HOME_RESOLVED:-not set}"
    echo "  my_machine   ${MY_MACHINE:-unknown}"; echo
    if is_eb_setup_bootstrap_tree 2>/dev/null && [[ -z "${HEN_HOUSE:-}" || ! -d "$HEN_HOUSE" ]]; then
        echo "Install type:   eb-setup bootstrap (installer only)"
        echo "  bootstrap      ${EB_ROOT:-$REPO_ROOT}"
        echo "  install target $(eb_default_install_dir) (or --install-dir)"
        echo
    elif [[ -n "$REPO_ROOT" ]] && git -C "$REPO_ROOT" rev-parse --git-dir >/dev/null 2>&1; then
        echo "Git (Tier 1):"
        echo "  EGSnrc repo    $(_git_remote_url "$REPO_ROOT") @ $(_git_branch "$REPO_ROOT")"
        git -C "$REPO_ROOT" status --short | sed 's/^/    /' || echo "    (clean)"
        if [[ -d "$EB_SUBMODULE/.git" || -f "$EB_SUBMODULE/.git" ]]; then
            echo "  egs_brachy sub @ $(_git_branch "$EB_SUBMODULE")"
            git -C "$EB_SUBMODULE" status --short | sed 's/^/    /' || echo "    (clean)"
        fi; echo
    elif [[ -n "$REPO_ROOT" && -d "$HEN_HOUSE" ]]; then
        echo "Install type:   release tarball (no git)"
        if [[ -n "${EB_ROOT:-}" && "$EB_ROOT" != "$REPO_ROOT" && -f "$EB_ROOT/eb-setup.sh" \
            && ! -f "$EB_ROOT/HEN_HOUSE/specs/unix.spec" ]]; then
            echo "  bootstrap      $EB_ROOT"
            echo "  install tree   $REPO_ROOT"
        fi
        echo;
    fi
    echo "EGS_HOME (Tier 2):"
    if [[ -n "$EGS_HOME_RESOLVED" && -d "$(eb_dest_path)" ]]; then
        echo "  egs_brachy/    present"
        if [[ -x "$(eb_executable)" ]]; then echo "  executable     $(eb_executable)"
        else echo "  executable     not found"; fi
    else echo "  egs_brachy/    not present"; fi; echo
    detect_scenario
    echo "Scenario: ${SCENARIO} — ${SCENARIO_DESC}"; echo
    case "$SCENARIO" in
        A0) echo "Next: eb-setup.sh install  (default install: ~/EGSnrc_CLRP)" ;;
        A|C|D) echo "Next: eb-setup.sh install" ;;
        B) echo "Next: fresh install to ~/Developer/scratch (see docs/eb-setup-testing.md)" ;;
        E) echo "Next: eb-setup.sh sync" ;;
        F) echo "Next: eb-setup.sh update"; echo "       (offline: update --from-tarball PATH)"; echo "       cd \"$(eb_dest_path)\" && make test" ;;
    esac
}

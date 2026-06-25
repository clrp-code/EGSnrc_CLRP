# shellcheck shell=bash
declare -a EB_STASH_REFS=()

_git_has_tracked_changes() { git -C "$1" diff --quiet 2>/dev/null && git -C "$1" diff --cached --quiet 2>/dev/null; }
_git_has_untracked() { [[ -n "$(git -C "$1" ls-files --others --exclude-standard 2>/dev/null)" ]]; }

_git_dirty_kind() {
    local repo="$1" tracked=0
    _git_has_tracked_changes "$repo" || tracked=1
    if (( tracked )); then echo tracked
    elif _git_has_untracked "$repo"; then echo untracked
    else echo clean; fi
}

enforce_tier1_git_policy() {
    local repo label kind
    for entry in "EGSnrc:$REPO_ROOT" "egs_brachy:$EB_SUBMODULE"; do
        label="${entry%%:*}"; repo="${entry#*:}"
        [[ -n "$repo" ]] || continue
        git -C "$repo" rev-parse --git-dir >/dev/null 2>&1 || continue
        kind="$(_git_dirty_kind "$repo")"
        [[ "$kind" == clean ]] && continue
        if [[ "$kind" == tracked ]]; then
            if (( STASH )); then continue; fi
            git -C "$repo" status --short | sed 's/^/  /'
            die "aborting: ${label} has tracked changes (try --stash)" 2
        fi
        if [[ "$kind" == untracked ]]; then
            if (( STRICT )); then
                git -C "$repo" ls-files --others --exclude-standard | sed 's/^/  /'
                die "aborting: ${label} has untracked files (--strict)" 2
            fi
            warn "${label} has untracked files only (continuing)"
        fi
    done
}

stash_git_repos() {
    local repo msg="eb-setup-$(date +%Y%m%d-%H%M%S)"
    EB_STASH_REFS=()
    for repo in "$REPO_ROOT" "$EB_SUBMODULE"; do
        [[ -n "$repo" ]] || continue
        git -C "$repo" rev-parse --git-dir >/dev/null 2>&1 || continue
        [[ "$(_git_dirty_kind "$repo")" == clean ]] && continue
        log "stashing $(basename "$repo")..."
        if (( DRY_RUN )); then run git -C "$repo" stash push -u -m "$msg"; continue; fi
        if git -C "$repo" stash push -u -m "$msg" >/dev/null; then
            EB_STASH_REFS+=("$repo")
        fi
    done
}

pop_git_stashes() {
    local repo
    for (( i=${#EB_STASH_REFS[@]}-1; i>=0; i-- )); do
        repo="${EB_STASH_REFS[$i]}"
        log "restoring stash in $(basename "$repo")..."
        run git -C "$repo" stash pop || die "stash pop failed in $repo — see: git -C \"$repo\" stash list"
    done
}

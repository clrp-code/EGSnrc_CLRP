# nrc-cnrc upstream integration (EGSnrc-eb)

This fork tracks [clrp-code/EGSnrc_CLRP](https://github.com/clrp-code/EGSnrc_CLRP) on branch `egs_brachy`, with `nrc` pointing at [nrc-cnrc/EGSnrc](https://github.com/nrc-cnrc/EGSnrc).

Some NRC changes land on `develop` before we merge them into CLRP. For features still in open PRs, we carry them on dedicated **`nrc/*` branches** so day-to-day CLRP work stays clean.

## Remotes

| Remote | Repository | Use |
|--------|------------|-----|
| `clrp` | clrp-code/EGSnrc_CLRP | CLRP fork; push target |
| `nrc` | nrc-cnrc/EGSnrc | Official upstream |
| `mchamberland-pr` | mchamberland/EGSnrc | Open NRC PR branches (optional) |

Add the PR author remote once:

```bash
git remote add mchamberland-pr https://github.com/mchamberland/EGSnrc.git
git fetch mchamberland-pr
```

## Branch layout

| Branch | Purpose |
|--------|---------|
| `egs_brachy` | Main CLRP line; tracks `clrp/egs_brachy` |
| `nrc/run-anywhere-pr1399` | Pre-merge [nrc-cnrc/EGSnrc#1399](https://github.com/nrc-cnrc/EGSnrc/pull/1399) (run-anywhere I/O) |
| `feature/*` | CLRP feature work (may update the egs_brachy submodule) |

**egs_brachy** itself is a submodule at `HEN_HOUSE/user_codes/egs_brachy` with its own branches.

## nrc/run-anywhere-pr1399

[PR #1399](https://github.com/nrc-cnrc/EGSnrc/pull/1399) adds run-anywhere behaviour:

- Resolve bare input names from the **current directory first**, then `$EGS_HOME/<user_code>/`
- Write run artifacts next to the input (or under `-d` / `--output-dir`)
- Wrapper and Qt GUI updates; macOS `@loader_path` for egs++ dylibs

This branch contains **only the five `#1212` commits** cherry-picked onto `egs_brachy` (not a full merge of `mchamberland-pr/run-anywhere`, which would pull unrelated `develop` history and conflict heavily).

Cherry-picked commits (oldest first):

1. `5c7b7db3` — Mortran core
2. `84a672ee` — egs++
3. `c8159a9d` — wrapper scripts
4. `613eaa0f` — egs_gui
5. `a1d5a4b8` — `@loader_path` dylibs

### Local development

Use this branch as the EGSnrc-eb base when you want run-anywhere for testing (e.g. eb_tests from any directory):

```bash
cd /path/to/EGSnrc-eb
git checkout nrc/run-anywhere-pr1399
# configure / build when ready — see main README
```

Combine with egs_brachy feature work by checking out the submodule branch you need; the submodule pointer is independent of the EGSnrc core branch.

Example: source-coordinate-transform work

```bash
# EGSnrc core with run-anywhere
git checkout nrc/run-anywhere-pr1399

# egs_brachy feature (submodule)
cd HEN_HOUSE/user_codes/egs_brachy
git checkout feature/source-coordinate-transform
```

### When upstream merges PR #1399

1. Fetch official upstream:
   ```bash
   git fetch nrc develop
   ```
2. Merge into CLRP main line:
   ```bash
   git checkout egs_brachy
   git merge nrc/develop
   ```
3. Reconcile `nrc/run-anywhere-pr1399`:
   - If the merged patches match, delete the branch — it is redundant.
   - If NRC squash-rebased, you may get conflicts in the same files; resolve by **preferring upstream** (`nrc/develop`).
   - Optionally revert the cherry-pick range on `egs_brachy` before merging `nrc/develop` if git reports duplicate changes.

This does **not** corrupt history permanently. Worst case: one-time merge conflicts in overlapping files (`egs_application.cpp`, scripts, etc.).

### Do not

- Merge `mchamberland-pr/run-anywhere` wholesale into `egs_brachy` (large conflict surface).
- Mix nrc-cnrc cherry-picks and unrelated CLRP features in a single commit.
- Assume PR branches stay rebased — re-fetch before adding more commits.

## Refreshing an nrc/* branch

```bash
git fetch mchamberland-pr
git checkout nrc/run-anywhere-pr1399
git rebase egs_brachy   # after CLRP moves forward
# If the PR was force-pushed, re-cherry-pick the five commits onto a fresh branch from egs_brachy
```

## Adding another pre-merge NRC PR

1. Branch from `egs_brachy`: `git checkout -b nrc/<short-name>-pr<number> egs_brachy`
2. Cherry-pick **only** the PR’s functional commits (inspect with `git log nrc/develop..<pr-branch>`).
3. Document the branch in this file.
4. When merged upstream, merge `nrc/develop` into `egs_brachy` and drop the `nrc/*` branch.

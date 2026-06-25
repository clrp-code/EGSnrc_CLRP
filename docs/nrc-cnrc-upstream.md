# nrc-cnrc upstream integration (EGSnrc-eb)

> **Run-anywhere adopted early.** [nrc-cnrc/EGSnrc#1399](https://github.com/nrc-cnrc/EGSnrc/pull/1399)
> (run-anywhere I/O) is cherry-picked onto `egs_brachy` before upstream merge.
> When NRC merges #1399 into `develop`, reconcile via `git merge nrc/develop` and
> **prefer upstream** on conflicts. See [When upstream merges PR #1399](#when-upstream-merges-pr-1399) below.

This fork tracks [clrp-code/EGSnrc_CLRP](https://github.com/clrp-code/EGSnrc_CLRP) on branch `egs_brachy`, with `nrc` pointing at [nrc-cnrc/EGSnrc](https://github.com/nrc-cnrc/EGSnrc).

Some NRC changes land on `develop` before we merge them into CLRP. For features still in open PRs, we may carry them on dedicated **`nrc/*` branches** for tracking; do not merge those branches wholesale into `egs_brachy`.

## Remotes

| Remote | Repository | Use |
|--------|------------|-----|
| `clrp` | clrp-code/EGSnrc_CLRP | CLRP fork; push target |
| `nrc` | nrc-cnrc/EGSnrc | Official upstream |
| `fork` | mchamberland/EGSnrc | Open NRC PR branches (optional) |

Add the PR author remote once:

```bash
git remote add fork https://github.com/mchamberland/EGSnrc.git
git fetch fork
```

## Run-anywhere on egs_brachy

[PR #1399](https://github.com/nrc-cnrc/EGSnrc/pull/1399) adds run-anywhere behaviour:

- Resolve bare input names from the **current directory first**, then `$EGS_HOME/<user_code>/`
- Write run artifacts next to the input (or under `-d` / `--output-dir`)
- Wrapper and Qt GUI updates; macOS `@rpath` for egs++ dylibs

Six functional `#1212` commits cherry-picked onto `egs_brachy` (oldest first; `@loader_path` commit omitted — superseded by `@rpath`):

1. Mortran core
2. egs++
3. Wrapper scripts
4. egs_gui
5. `@rpath` install names + embedded DSO rpath (macOS)
6. JCF `.lock` / uniform RCO paths via `getOutputDir()` (matches upstream `3ce3e860`)

**egs_brachy** user code is a submodule at `HEN_HOUSE/user_codes/egs_brachy` with its own branches.

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
3. On conflicts in run-anywhere files: resolve by **preferring upstream** (`nrc/develop`).
4. If git reports duplicate changes from the early cherry-picks, optionally revert the cherry-pick range before merging:
   ```bash
   git revert <oldest>^..<newest>
   git merge nrc/develop
   ```
5. Delete obsolete `nrc/run-anywhere-pr1399` branch if it still exists.

This does **not** corrupt history permanently. Worst case: one-time merge conflicts in overlapping files (`egs_application.cpp`, scripts, etc.).

### Do not

- Merge `fork/run-anywhere` wholesale into `egs_brachy` (large conflict surface).
- Mix nrc-cnrc cherry-picks and unrelated CLRP features in a single commit.
- Assume PR branches stay rebased — re-fetch before adding more commits.

## Adding another pre-merge NRC PR

1. Branch from `egs_brachy`: `git checkout -b nrc/<short-name>-pr<number> egs_brachy`
2. Cherry-pick **only** the PR's functional commits (inspect with `git log nrc/develop..<pr-branch>`).
3. Document the branch in this file.
4. When merged upstream, merge `nrc/develop` into `egs_brachy` and drop the `nrc/*` branch.

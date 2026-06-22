# eb-setup testing guide

## Branch layout

See [branch-layout.md](branch-layout.md) for the full policy. Short version:

| Branch | Purpose |
|--------|---------|
| `egs_brachy` | User/release line |
| `dev` | General CLRP integration |
| `feature/eb-setup` | **Current** — `eb-setup.sh` and install tooling |
| `feature/*` | Other focused tasks |

While working on `eb-setup.sh`, stay on **`feature/eb-setup`**. That branch includes NRC [#1399](https://github.com/nrc-cnrc/EGSnrc/pull/1399) (run from any directory) and [#934](https://github.com/nrc-cnrc/EGSnrc/issues/934) (`$EGS_HOME` in `include file`).

| Repo | Branch | Purpose |
|------|--------|---------|
| EGSnrc-eb | `feature/eb-setup` | eb-setup.sh development |
| egs_brachy submodule | `feature/eb-setup` (when mirrored) | standalone copy of tool |

## Scratch directory

Use `~/Developer/scratch` for isolated install tests — never the main dev tree.

```bash
# Example fresh install test
cd ~/Developer/scratch
rm -rf EGSnrc_CLRP-test
git clone https://github.com/clrp-code/EGSnrc_CLRP.git EGSnrc_CLRP-test
cd EGSnrc_CLRP-test
git checkout egs_brachy
git submodule update --init --recursive
# copy or use eb-setup.sh from feature branch
./eb-setup.sh install
```

## Quick checks (main dev tree)

```bash
cd /Users/marc/Developer/clrp/EGSnrc-eb
export EGS_CONFIG="$PWD/HEN_HOUSE/specs/eb-dev.conf"
export EGS_HOME=~/Developer/scratch/egs_home   # or your test EGS_HOME
./eb-setup.sh check
./eb-setup.sh sync --dry-run
```

## Scratch smoke test

Automated script (run after `eb-setup.sh sync` and rebuild):

```bash
cd /Users/marc/Developer/clrp/EGSnrc-eb   # on feature/eb-setup
export EGS_CONFIG="$PWD/HEN_HOUSE/specs/eb-dev.conf"
export EGS_HOME=~/Developer/scratch/egs_home/
./docs/scratch-smoke-test.sh
```

**Verified 2026-06-21** — all 6 checks pass:

| Test | Behaviour |
|------|-----------|
| Input in cwd | `egs_brachy -i smoke_ra -s` from `~/Developer/scratch` writes outputs in scratch |
| Bare name fallback | `smoke_ra.egsinp` only in `$EGS_HOME/egs_brachy/` runs from scratch; outputs next to resolved input |
| `$EGS_HOME` include | `include file = $EGS_HOME/egs_brachy/lib/transport/...` works (#934) |

Notes:

- Use **`-s`** (simple run control) when cwd ≠ `$EGS_HOME/egs_brachy` — default JCF still opens lock files under `$EGS_HOME/egs_brachy/` (known gap).
- **`material data file`** / **`muen file`** — expand `$EGS_HOME` after this patch (`replace_env` in pegsless Fortran; `egsExpandPath` for egs_brachy `muen file`). Use `$VAR/` at the start of the path (Fortran does not support `%VAR%`).
- Outputs are written **next to the resolved input file**, not always cwd (PR #1399).

Manual run:

```bash
cd ~/Developer/scratch
export EGS_CONFIG=/Users/marc/Developer/clrp/EGSnrc-eb/HEN_HOUSE/specs/eb-dev.conf
export EGS_HOME=~/Developer/scratch/egs_home/
export PATH="$EGS_HOME/bin/eb-dev:$PATH"
source "$EGS_CONFIG" 2>/dev/null || true
# optional: source $HEN_HOUSE/scripts/clrp_bashrc_additions  # defines exeb
egs_brachy -i smoke_ra -s
# or: exeb smoke_ra -s
```

Input template: [`docs/smoke_ra.egsinp`](smoke_ra.egsinp)


## Dirty tree tests

```bash
# Tier 1 abort (tracked change in HEN_HOUSE)
touch HEN_HOUSE/egs++/dummy.txt
./eb-setup.sh update          # should abort
./eb-setup.sh update --stash  # should stash, update, pop

# Tier 2 (EGS_HOME) — should not block update
echo test > $EGS_HOME/egs_brachy/my_test.egsinp
./eb-setup.sh update          # should proceed
```

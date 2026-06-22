# eb-setup testing guide

> **Branch policy:** [branch-layout.md](branch-layout.md) · **Status checklist:** [branch-tracker.md](branch-tracker.md)

This branch (`feature/eb-setup`) is **eb-setup only** — no run-anywhere (#1399). For run-anywhere scratch tests, checkout **`dev`**.

## Quick checks (main dev tree)

```bash
cd /Users/marc/Developer/clrp/EGSnrc-eb
git checkout feature/eb-setup
export EGS_CONFIG="$PWD/HEN_HOUSE/specs/eb-dev.conf"   # or your config
export EGS_HOME=~/Developer/scratch/egs_home/
./eb-setup.sh check
./eb-setup.sh sync --dry-run
```

## Scratch directory (fresh install test)

Use `~/Developer/scratch` — never the main dev tree.

Use a **clean shell** so your main install does not leak in (`unset` is enough — do not edit `~/.zshrc`):

```bash
cd ~/Developer/scratch
rm -rf EGSnrc_CLRP-test
git clone https://github.com/clrp-code/EGSnrc_CLRP.git EGSnrc_CLRP-test
cd EGSnrc_CLRP-test
git checkout feature/eb-setup
git submodule update --init --recursive

unset EGS_CONFIG HEN_HOUSE EGS_HOME
# use $HOME or unquoted ~ — quoted "~/..." is not expanded by the shell
./eb-setup.sh install --egs-home "$HOME/Developer/scratch/egs_home_test/"
```

`install` runs `HEN_HOUSE/scripts/configure` as `./configure` from that directory (required by EGSnrc).

After `configure`, set `EGS_CONFIG` / `EGS_HOME` and run:

```bash
./eb-setup.sh sync
./eb-setup.sh check
```

## Run-anywhere smoke test (on `dev`, not this branch)

```bash
git checkout dev
./docs/scratch-smoke-test.sh
```

See [eb-setup-testing.md on `dev`](https://github.com/clrp-code/EGSnrc_CLRP/blob/dev/docs/eb-setup-testing.md) for full scratch I/O tests (#1399 + #934).

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

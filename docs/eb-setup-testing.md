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

Use a **short clone path** — Mortran fails pegs4 if `HEN_HOUSE` paths are too long (`FATAL STRING OR STATEMENT TOO LONG` in `configure.log`). Prefer `$HOME/scratch/eb` (not `EGSnrc_CLRP-test`).

Use a **clean shell** so your main install does not leak in (`unset` is enough — do not edit `~/.zshrc`):

```bash
mkdir -p ~/scratch
cd ~/scratch
rm -rf eb
git clone https://github.com/clrp-code/EGSnrc_CLRP.git eb
cd eb
git checkout feature/eb-setup
git pull
git submodule update --init --recursive

unset EGS_CONFIG HEN_HOUSE EGS_HOME
./eb-setup.sh install --egs-home "$HOME/scratch/egs_home/"
```

Before pushing `feature/eb-setup`, run `./scripts/test-eb-setup.sh` from the repo root.

`install` runs `HEN_HOUSE/scripts/configure` as `./configure` from that directory (required by EGSnrc).

After `configure`, `install` auto-detects paths, runs `sync`, then prints a **shell setup** block and writes `eb-env.sh` in the repo root. Add those lines to your profile (or `source ./eb-env.sh` in the current shell only).

```bash
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

## Tarball testing

See **[release-tarballs.md](release-tarballs.md)** for building local sample tarballs and the full install/update walkthrough.

Quick start:

```bash
./scripts/build-test-tarballs.sh    # dist/EGSnrc_CLRP-egs_brachy-1.0.0-alpha.{1,2}.tar.gz
```

```bash
# Fresh install from release (no git):
./eb-setup.sh install --from-tarball ~/path/to/dist/EGSnrc_CLRP-egs_brachy-1.0.0-alpha.1.tar.gz \
  --install-dir ~/scratch/tarball-test/eb-release

# Update existing tarball install (EGS_CONFIG/EGS_HOME must be set):
cd ~/scratch/tarball-test/eb-release && source ./eb-env.sh
./eb-setup.sh update --from-tarball ~/path/to/dist/EGSnrc_CLRP-egs_brachy-1.0.0-alpha.2.tar.gz
```

Preserves on update: `HEN_HOUSE/specs/*.conf`, `lib/`, `bin/`, `log/`, `egs++/dso/`, and all of `$EGS_HOME`.

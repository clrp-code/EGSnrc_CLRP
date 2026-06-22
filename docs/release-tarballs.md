# Release tarballs

## Local test builds (start here)

Do **not** need GitHub Actions or tags for initial eb-setup tarball testing.

```bash
cd /path/to/EGSnrc-eb   # feature/eb-setup
git submodule update --init --recursive
./scripts/build-test-tarballs.sh
```

Produces in `dist/`:

| File | Purpose |
|------|---------|
| `EGSnrc_CLRP-egs_brachy-1.0.0-alpha.1.tar.gz` | fresh `install --from-tarball` |
| `EGSnrc_CLRP-egs_brachy-1.0.0-alpha.2.tar.gz` | `update --from-tarball` (includes a small marker file) |

Single version: `./scripts/build-release-tarball.sh 1.0.0-alpha.1`

Layout matches CI (`.github/workflows/release-tarball.yml`): no `.git`, vendored `egs_brachy` submodule contents, no local `lib/`/`bin/`/`dso/` build artifacts.

## Tarball testing walkthrough

Use a **short install path** (Mortran path limits), separate from your git scratch tree:

```bash
mkdir -p ~/scratch/tarball-test
cd ~/scratch/tarball-test

# 1) Install from alpha.1
rm -rf eb-release
/path/to/EGSnrc-eb/eb-setup.sh install \
  --from-tarball /path/to/EGSnrc-eb/dist/EGSnrc_CLRP-egs_brachy-1.0.0-alpha.1.tar.gz \
  --install-dir ~/scratch/tarball-test/eb-release \
  --egs-home "$HOME/scratch/tarball-egs_home/"

# configure interactively when prompted; install continues to sync + eb-env.sh

# 2) New shell / source env
cd ~/scratch/tarball-test/eb-release
source ./eb-env.sh
./eb-setup.sh check    # expect scenario F after configure+sync

# 3) Dry-run update
./eb-setup.sh update --from-tarball \
  /path/to/EGSnrc-eb/dist/EGSnrc_CLRP-egs_brachy-1.0.0-alpha.2.tar.gz \
  --dry-run

# 4) Real update + test
./eb-setup.sh update --from-tarball \
  /path/to/EGSnrc-eb/dist/EGSnrc_CLRP-egs_brachy-1.0.0-alpha.2.tar.gz
cd "$EGS_HOME/egs_brachy" && make test

# 5) Marker file from alpha.2 should appear after update:
test -f "$EGS_HOME/egs_brachy/.eb-setup-tarball-test-stamp" && cat "$EGS_HOME/egs_brachy/.eb-setup-tarball-test-stamp"
```

Checklist:

- [ ] `install --from-tarball` → configure → sync → scenario F
- [ ] `update --from-tarball --dry-run` shows rsync preview
- [ ] `update --from-tarball` preserves `HEN_HOUSE/specs/*.conf`
- [ ] `make test` passes after update
- [ ] `check` shows **release tarball (no git)**
- [ ] egs log banner shows `EGSnrc <version> for …` (from `release.mk` in tarball)

## Version banner (`EGS_RELEASE` / `GIT_HASH`)

Git checkouts pick up release label and commit hash at **compile time** via make (`HEN_HOUSE/specs/all_common.spec`):

- **`EGS_RELEASE`** — nearest annotated git tag (`git describe --tags --abbrev=0`), with leading `v` or `egs_brachy-` stripped
- **`GIT_HASH`** — short commit hash (empty if not in a git repo)

End-user tarballs have **no `.git`**, so `scripts/build-release-tarball.sh` writes `HEN_HOUSE/specs/release.mk` into the tarball:

```makefile
EGS_RELEASE = -DEGS_RELEASE="\"1.0.0-alpha.1\""
GIT_HASH = -DGIT_HASH="\"cab17a98\""
```

After install/configure/sync, recompiling shows in the egs log, e.g.:

```
EGSnrc 1.0.0-alpha.1 for arm-apple-darwin…
```

If `EGS_RELEASE` is unset, the banner is `EGSnrc for …` (no fake “version 4”).

## GitHub Actions

Workflow: `.github/workflows/release-tarball.yml`

| Trigger | What happens |
|---------|----------------|
| **`workflow_dispatch`** (manual, version input) | Builds end-user tarball only; uploads to **Actions artifacts** — no GitHub Release |
| **Tag push** `egs_brachy-*` | Builds end-user + source tarballs; publishes a **GitHub Release** |

### Draft and prerelease (alpha/beta tags)

When you push a tag whose name contains `alpha` or `beta` (e.g. `egs_brachy-1.0.0-alpha.1`), the workflow sets:

```yaml
draft: true       # if tag contains alpha or beta
prerelease: true
```

| Setting | What users see |
|---------|----------------|
| **Draft** | Release hidden from the public Releases list; only people with the direct URL can download assets. Good for smoke-testing CI output before announcing. |
| **Prerelease** | Release is visible but marked “Pre-release”; GitHub does **not** show it as **Latest**. Appropriate for alpha/beta builds users may try voluntarily. |
| **Neither** (e.g. tag `egs_brachy-1.0.0`) | Normal public release; GitHub may mark it **Latest**. |

To publish a tested alpha: open the draft release on GitHub → Edit → uncheck **Draft** (keep **Pre-release** checked until GA).

Suggested rollout:

1. **Local** — `./scripts/build-test-tarballs.sh`
2. **CI dry run** — Actions → Release tarball → Run workflow → download artifact
3. **Internal GH test** — `git tag egs_brachy-1.0.0-alpha.1 && git push clrp egs_brachy-1.0.0-alpha.1` → draft + prerelease release
4. **Ship** — tag `egs_brachy-1.0.0` (no alpha/beta) → full public release

## End-user vs source tarball

| Asset | CI name | eb-setup uses |
|-------|---------|---------------|
| End-user (no `.git`) | `EGSnrc_CLRP-egs_brachy-VERSION.tar.gz` | **yes** — `install` / `update --from-tarball` |
| Source (with `.git`) | `EGSnrc_CLRP-egs_brachy-VERSION-src.tar.gz` | no — for developers who want git history |

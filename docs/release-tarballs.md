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

## GitHub Actions (when ready)

Workflow: `.github/workflows/release-tarball.yml` — runs on tag push `egs_brachy-*`.

### Keeping builds off the public “latest” download

There is no hidden “do not download” flag. Use release **visibility** instead:

| Mechanism | Effect |
|-----------|--------|
| **Draft release** | Not listed publicly; only people with the URL see assets. Best for CI smoke tests. |
| **Prerelease** | Visible but marked pre-release; GitHub won't treat it as Latest. Use for alpha/beta. |
| **Tag naming** | e.g. `egs_brachy-1.0.0-alpha.1` — signals intent; pair with prerelease. |
| **Manual workflow only** | Change trigger from `push: tags` to `workflow_dispatch` until you're ready for tag-driven releases. |

Suggested rollout:

1. **Now** — local `scripts/build-test-tarballs.sh` only.
2. **CI dry run** — add `workflow_dispatch` job that uploads artifacts to the **Actions run** (not a Release); download from the workflow page.
3. **Internal GH test** — push tag, workflow creates a **draft** + **prerelease** release.
4. **Ship** — remove draft/prerelease; tag `egs_brachy-1.0.0` as full release.

Example workflow addition (later):

```yaml
on:
  workflow_dispatch:
  push:
    tags:
      - 'egs_brachy-*'

# in softprops/action-gh-release@v2:
with:
  draft: ${{ contains(github.ref_name, 'alpha') || contains(github.ref_name, 'beta') }}
  prerelease: ${{ contains(github.ref_name, 'alpha') || contains(github.ref_name, 'beta') }}
```

## End-user vs source tarball

| Asset | CI name | eb-setup uses |
|-------|---------|---------------|
| End-user (no `.git`) | `EGSnrc_CLRP-egs_brachy-VERSION.tar.gz` | **yes** — `install` / `update --from-tarball` |
| Source (with `.git`) | `EGSnrc_CLRP-egs_brachy-VERSION-src.tar.gz` | no — for developers who want git history |

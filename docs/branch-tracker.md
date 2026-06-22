# Branch tracker (CLRP EGSnrc-eb)

Last updated: 2026-06-22

## At a glance

| Branch | Run-anywhere? | eb-setup? | Ships to users? |
|--------|---------------|-----------|-----------------|
| **`egs_brachy`** | No | After merge | **Yes** |
| **`dev`** | **Yes** | No | No (testing only) |
| **`feature/eb-setup`** | **No** | **In progress** | → merges to `egs_brachy` |

## `feature/eb-setup` (you are here)

- Base: `egs_brachy` (#934 only)
- Adding: `eb-setup.sh`, `scripts/lib/`, release tarball CI, install docs
- **Not included:** run-anywhere, scratch smoke test (those are on `dev`)

## Next steps (eb-setup)

- [x] Split branch from `egs_brachy` (no run-anywhere)
- [ ] Commit `eb-setup.sh` bundle
- [x] Scratch install + `make test` end-to-end
- [x] **Tarball** install/update walk-through — see [release-tarballs.md](release-tarballs.md)
- [ ] Merge `feature/eb-setup` → `egs_brachy`

## Run-anywhere (use `dev`)

```bash
git checkout dev
./docs/scratch-smoke-test.sh
```

Full tracker on `dev` branch: see `docs/branch-tracker.md` there.

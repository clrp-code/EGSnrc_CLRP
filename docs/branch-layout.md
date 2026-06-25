# EGSnrc_CLRP branch layout

How we use branches in [clrp-code/EGSnrc_CLRP](https://github.com/clrp-code/EGSnrc_CLRP).

**Living checklist:** see [branch-tracker.md](branch-tracker.md).

## Everyday branches

| Branch | When to use | Merge target |
|--------|-------------|--------------|
| **`egs_brachy`** | User/release line | (users) |
| **`dev`** | Run-anywhere (#1399) testing | `egs_brachy` after NRC merges |
| **`feature/eb-setup`** | **`eb-setup.sh`** install tooling | **`egs_brachy`** when ready |

## Policy

- **`feature/eb-setup`** is based on **`egs_brachy`** and contains **no run-anywhere**.
- **`dev`** holds run-anywhere until [PR #1399](https://github.com/nrc-cnrc/EGSnrc/pull/1399) lands upstream.
- Do **not** merge `dev` into `egs_brachy` wholesale before NRC merges run-anywhere.

See [nrc-cnrc-upstream.md](nrc-cnrc-upstream.md) on the `dev` branch for upstream PR tracking.

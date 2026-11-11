# Lean setup

The root Lake project uses the `physlib` git submodule as a local dependency.
Its `lean-toolchain` is kept in sync with Physlib (`leanprover/lean4:v4.35.0`),
and Physlib pins Mathlib to the matching release.

Initialize submodules doing anything Lean-related:

```sh
git submodule update --init --recursive
```

`.devcontainer/setup-lean.sh` download's Physlib's Lake cache, builds Physlib, and then builds the
root Lean project.

Add Lean modules under `lean/`.

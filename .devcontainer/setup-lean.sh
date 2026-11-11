#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
physlib_dir="${repo_root}/physlib"

if [[ ! -f "${physlib_dir}/lakefile.toml" || ! -f "${physlib_dir}/lean-toolchain" ]]; then
  echo "physlib submodule is missing or not initialized: ${physlib_dir}" >&2
  echo "Run 'git submodule update --init --recursive' and rebuild the dev container." >&2
  exit 1
fi

cd "${physlib_dir}"
if ! lake exe get_cache; then
  echo "Warning: Physlib cache download failed; continuing with a regular Lake build." >&2
fi
lake build

cd "${repo_root}"
lake build

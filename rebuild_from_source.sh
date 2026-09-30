#!/usr/bin/env bash
# Rebuild the installed-package files in this folder from source/RiboseQC.
#
# R does not run a package from its .R source: it loads the built files (R/RiboseQC.rdb and .rdx).
# So after a change in source/RiboseQC, run this script with the R that will use the package.
# It builds into a temporary folder, and only if that succeeds replaces the built files here.
# Everything else in this folder (.git, source/, README.md, this script) is left alone.
#
# Do not run it while jobs that use RiboseQC are running.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build="$(mktemp -d)"
trap 'rm -rf "$build"' EXIT

echo "Building with $(R --version | head -1)"
R CMD INSTALL --library="$build" "$here/source/RiboseQC"

for item in "$build"/RiboseQC/*; do
    name="$(basename "$item")"
    rm -rf "${here:?}/$name"
    cp -R "$item" "$here/$name"
done

echo "Done: $(grep '^Version' "$here/DESCRIPTION"), $(grep '^Built' "$here/DESCRIPTION")"

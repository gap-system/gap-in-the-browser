#!/usr/bin/env bash
#
# Build a GitHub-Pages-sized "GAP in the browser" site from a clean GAP
# clone. The result lands in build/gap/web-example/; deploy.sh publishes
# it to the gh-pages branch.
#
# Requirements: git, curl, docker or podman, and a working native GAP
# build toolchain (autotools, C compiler) -- the native build only exists
# to run "make html", since a git checkout contains no built manuals and
# the in-browser help would otherwise come up empty.
#
# Until the emscripten polish work is merged upstream, we build from the
# ChrisJefferson fork. Repoint at gap-system/gap (branch master) once it
# is up to date:
GAP_REPO="${GAP_REPO:-https://github.com/ChrisJefferson/gap.git}"
GAP_REF="${GAP_REF:-new-polish}"
PKG_URL="${PKG_URL:-https://github.com/gap-system/PackageDistro/releases/download/latest/packages.tar.gz}"

set -euo pipefail
cd "$(dirname "$0")"

# Packages removed from the shipped site. GitHub Pages caps a published
# site at 1 GB; the full package distribution plus GAP itself is ~1.1 GB,
# and much of that can never work in the browser anyway. A package is on
# this list because either
#   (a) its AvailabilityTest refuses to load without a compiled kernel
#       module or external binary (neither exists under wasm), or
#   (b) it needs (transitively) a package from (a), so it can never load
#       either, or
#   (c) it loads but is useless in a browser (X11 GUIs, interfaces to
#       external systems, the package manager).
# Notable casualties of (b): hap (via nq), simpcomp (via io), rcwa (via
# fr -> io), unitlib (via io), agt (via digraphs). If GAP's wasm build
# ever gains kernel modules for io/orb/digraphs, revisit this list --
# statically linking just io and orb would revive ~150 MB of pure-GAP
# packages.
PRUNE_PACKAGES=(
    4ti2interface ace agt anupq browse caratinterface cddinterface
    classicalmaximals cohomolo crypting curlinterface cvec
    datastructures deepthought digraphs ferret fining float fplsa fr
    francy fwtree genss hap hapcryst help images io io_forhomalg itc
    json jupyterkernel jupyterviz kan kbmag majoranaalgebras matgrp
    nconvex normalizinterface nq openmath orb origami packagemaker
    packagemanager polymaking profiling rcwa recog scscp semigroups
    sgpviz simpcomp singular unitlib walrus xgap xmod xmodalg
    zeromqinterface
)

SRC=build/gap
mkdir -p build

# The GAP package distribution (~540 MB) is cached across runs; delete
# build/packages.tar.gz to pick up a newer release.
if [[ ! -f build/packages.tar.gz ]]; then
    echo ">> Downloading GAP package distribution"
    curl -fL -o build/packages.tar.gz.part "$PKG_URL"
    mv build/packages.tar.gz.part build/packages.tar.gz
fi

echo ">> Fresh clone of $GAP_REPO ($GAP_REF)"
rm -rf "$SRC"
git clone --depth 1 --branch "$GAP_REF" "$GAP_REPO" "$SRC"

echo ">> Extracting packages"
mkdir "$SRC/pkg"
tar -xzf build/packages.tar.gz -C "$SRC/pkg"

echo ">> Native GAP build (only needed to build the manuals)"
# The manuals must be built BEFORE pruning: doc/make_doc fails the build
# on any unresolved reference, and the GAP manuals cross-reference the
# manuals of packages we are about to delete (resolved via each package's
# shipped manual.six). "make html" skips the PDF versions, which we would
# delete anyway, so no TeX is needed.
(
    cd "$SRC"
    ./autogen.sh
    ./configure
    make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
    make html
)

echo ">> Pruning packages that cannot work in the browser"
for p in "${PRUNE_PACKAGES[@]}"; do
    # A missing directory means the distribution renamed or dropped the
    # package: the list above is stale and needs a human look.
    if [[ ! -d "$SRC/pkg/$p" ]]; then
        echo "Error: pkg/$p not found in the package distribution." >&2
        echo "Update PRUNE_PACKAGES in $0." >&2
        exit 1
    fi
    rm -rf "${SRC:?}/pkg/$p"
done
# grape works without nauty (only the isomorphism functions need it);
# its bundled nauty sources are 8 MB we can never compile or run.
# Likewise guava's Leon binaries.
rm -rf "$SRC"/pkg/grape/nauty* "$SRC"/pkg/guava/src "$SRC"/pkg/guava/bin

# The in-browser help reads the txt/HTML manuals; the package PDFs are
# ~60 MB of dead weight. Delete before the wasm build, so gap-fs.json
# never lists them.
echo ">> Deleting PDFs"
find "$SRC/pkg" "$SRC/doc" -name '*.pdf' -delete

echo ">> Building wasm GAP and assembling the site"
(cd "$SRC" && etc/emscripten/build-in-docker.sh)

# Strip every .gitignore from the assembled site. Many GAP packages ship
# a developer .gitignore that excludes their *built* doc artifacts
# (doc/manual.six, doc/chap*.html, ...) -- the very files we just built
# and want to serve. Left in place, deploy.sh's "git add" honours those
# rules and silently drops the docs from the published site. They are
# dead files in a static site regardless.
find "$SRC/web-example" -name .gitignore -delete

# Without this GitHub Pages runs the site through Jekyll, which drops
# files starting with "_" and chokes on a 45000-file tree.
touch "$SRC/web-example/.nojekyll"

SIZE=$(du -sh "$SRC/web-example" | cut -f1)
cat <<EOF

Site built: $SRC/web-example ($SIZE; the GitHub Pages limit is 1 GB)
Test it:    (cd $SRC/web-example && ../etc/emscripten/serve.py)
            then open http://localhost:8080/
Publish:    ./deploy.sh
EOF

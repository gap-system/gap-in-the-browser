# GAP in the browser

Build and publish a self-contained [GAP](https://www.gap-system.org)
website: the GAP kernel compiled to WebAssembly, plus the library and as
much of the package distribution as fits, served as a static GitHub
Pages site. The build machinery itself lives in the GAP repository under
`etc/emscripten/`; this repository just drives it from a clean clone and
holds the published site.

## Usage

```sh
./build.sh     # ~30-60 minutes: clean clone, native build (for the
               # manuals), wasm build in a pinned emsdk container
./deploy.sh    # force-push the site to the gh-pages branch
```

Requirements: git, curl, docker or podman, and a native GAP build
toolchain (autotools, C/C++ compiler). The native build exists only to
run `make doc` -- a git checkout contains no built manuals, and without
them the in-browser help is empty.

Knobs (environment variables): `GAP_REPO` and `GAP_REF` select the GAP
clone (currently defaulting to the ChrisJefferson fork until the
emscripten polish is merged upstream), `PKG_URL` the package
distribution tarball. The tarball is cached as `build/packages.tar.gz`;
delete it to pick up a newer distribution release.

## What gets left out

GitHub Pages caps a published site at 1 GB and GAP plus its full package
distribution is ~1.1 GB, so `build.sh` prunes:

- packages requiring a compiled kernel module or external binary --
  neither can exist in the wasm build, so these can never load
  (semigroups, digraphs, io, orb, browse, anupq, ...);
- packages needing one of those, transitively (hap via nq, simpcomp and
  unitlib via io, rcwa via fr, agt via digraphs, ...);
- packages that load but cannot do anything useful in a browser (X11
  GUIs, interfaces to external systems, the package manager);
- all PDF manuals (the help system reads the text/HTML versions) and
  grape's bundled nauty sources.

The full list, with reasoning, is `PRUNE_PACKAGES` in `build.sh`. The
result is ~110 packages and a site of roughly 650 MB.

## Hosting notes

The terminal needs `SharedArrayBuffer`, which browsers only enable on
cross-origin-isolated pages. GitHub Pages cannot send the required
COOP/COEP headers, so the site includes `coi-serviceworker.js` as a
workaround. GitHub Pages must be configured (repository settings ->
Pages) to serve the `gh-pages` branch, root directory.

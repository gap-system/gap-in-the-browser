#!/usr/bin/env bash
#
# Publish build/gap/web-example/ to the gh-pages branch as a single
# parentless commit and force-push it. Rebuilding and redeploying never
# grows the repository: each deployment replaces the whole branch, so
# dead wasm blobs from old deployments become unreachable and GitHub
# garbage-collects them.

set -euo pipefail
cd "$(dirname "$0")"

SITE=build/gap/web-example
if [[ ! -f "$SITE/gap.wasm" || ! -f "$SITE/.nojekyll" ]]; then
    echo "Error: no built site at $SITE -- run ./build.sh first." >&2
    exit 1
fi

# Stage the site with a throwaway index so the real working tree and
# index are never touched.
GIT_INDEX_FILE="$(mktemp -u)"
export GIT_INDEX_FILE
trap 'rm -f "$GIT_INDEX_FILE"' EXIT

git --work-tree="$SITE" add -A
TREE=$(git write-tree)
COMMIT=$(git commit-tree "$TREE" -m "Deploy $(cat "$SITE/build-id")")
git update-ref refs/heads/gh-pages "$COMMIT"
git push -f origin gh-pages

echo "Deployed $(cat "$SITE/build-id") to gh-pages."

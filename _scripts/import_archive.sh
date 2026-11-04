#!/bin/sh
# One-shot copy of the MkDocs-era versions (2.0 … 3.8.2) from the gh-pages branch into R2.
#
#   _scripts/import_archive.sh [destination [version ...]]
#
# destination defaults to r2:nextflow-training-docs/training. Versions default to every entry of
# gh-pages/versions.json except 0.dev. Each goes to <destination>/<version>/.
# latest/ and versions.json (without 0.dev) are written only if absent, so a run after
# the first Astro release leaves them alone. Re-runs skip files already uploaded.
#
# Needs git, jq and an rclone remote named r2, for example:
#   export RCLONE_CONFIG_R2_TYPE=s3 RCLONE_CONFIG_R2_PROVIDER=Cloudflare RCLONE_CONFIG_R2_REGION=auto \
#     RCLONE_CONFIG_R2_NO_CHECK_BUCKET=true RCLONE_CONFIG_R2_ACCESS_KEY_ID=… \
#     RCLONE_CONFIG_R2_SECRET_ACCESS_KEY=… RCLONE_CONFIG_R2_ENDPOINT=https://<account>.r2.cloudflarestorage.com
set -eu

dest=${1:-r2:nextflow-training-docs/training}
[ $# -gt 0 ] && shift
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# gh-pages is about 9 GB, so check out one version at a time.
git clone --quiet --depth 1 --branch gh-pages --filter=blob:none --sparse \
    https://github.com/nextflow-io/training.git "$work/gh-pages"
cd "$work/gh-pages"
git sparse-checkout set --no-cone /versions.json

versions=${*:-$(jq -r '.[].version | select(. != "0.dev")' versions.json)}
# Not `rclone cat`: on R2 it exits 0 for a missing object.
exists() { [ -n "$(rclone lsf "$dest/$1" 2>/dev/null)" ]; }
if ! exists latest/index.html; then versions="$versions latest"; fi

for v in $versions; do
    echo "== $v"
    git sparse-checkout set --no-cone /versions.json "/$v/"
    [ -d "$v" ] || { echo "$v is not on gh-pages" >&2; exit 1; }
    rclone copy "$v" "$dest/$v" --checksum --fast-list --transfers 32 --checkers 64 --stats-one-line --stats 30s
done

if exists versions.json; then
    echo "$dest/versions.json exists; left as is."
else
    jq 'map(select(.version != "0.dev"))' versions.json | rclone rcat "$dest/versions.json"
fi

#!/usr/bin/env bash
# Adapted from https://github.com/shz9/1KG_HGDP_ld/blob/main/scripts/prepare_for_zenodo.sh
# Same package layout: ANCESTRY.tar.gz containing ANCESTRY/chr_N/<Zarr keys>.
set -euo pipefail
cd "$(dirname "$0")/.."
INPUT_DIR=${INPUT_DIR:-${1:-results}}
OUTPUT_DIR=${OUTPUT_DIR:-${2:-results_zenodo}}
[[ -d "$INPUT_DIR" ]] || { echo "Missing LD directory: $INPUT_DIR" >&2; exit 1; }
for cmd in unzip tar sha256sum; do command -v "$cmd" >/dev/null; done
mkdir -p "$OUTPUT_DIR"
shopt -s nullglob
archives=("$INPUT_DIR"/*/chr_*.zip)
[[ ${#archives[@]} -gt 0 ]] || { echo "No ancestry/chr_*.zip files in $INPUT_DIR" >&2; exit 1; }
tmp_dir=
partial=
cleanup() {
    [[ -z "$tmp_dir" ]] || rm -rf -- "$tmp_dir"
    [[ -z "$partial" ]] || rm -f -- "$partial"
    return 0
}
trap cleanup EXIT
for ancestry_dir in "$INPUT_DIR"/*/; do
    zip_files=("$ancestry_dir"/chr_*.zip)
    [[ ${#zip_files[@]} -gt 0 ]] || continue
    ancestry=$(basename "$ancestry_dir")
    tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/cag-zenodo.XXXXXX")
    for zip_file in "${zip_files[@]}"; do
        chr=$(basename "$zip_file" .zip)
        unzip -tq "$zip_file" >/dev/null
        mkdir -p "$tmp_dir/$ancestry/$chr"
        unzip -q "$zip_file" -d "$tmp_dir/$ancestry/$chr"
    done
    partial="$OUTPUT_DIR/.${ancestry}.tar.gz.part"
    tar -czf "$partial" -C "$tmp_dir" "$ancestry"
    mv -- "$partial" "$OUTPUT_DIR/${ancestry}.tar.gz"
    partial=
    (cd "$OUTPUT_DIR" && sha256sum "${ancestry}.tar.gz" > "${ancestry}.tar.gz.sha256")
    rm -rf -- "$tmp_dir"
    tmp_dir=
    echo "Created $OUTPUT_DIR/${ancestry}.tar.gz"
done

#!/usr/bin/env bash
set -euo pipefail

src="$1"
dst="$2"

mkdir -p "$(dirname "$dst")"

case "$src" in
    *.gz)
        # Stream-copy gzip FASTQ so the workflow owns a normalized path.
        cp -f "$src" "$dst"
        ;;
    *)
        gzip -c "$src" > "$dst"
        ;;
esac

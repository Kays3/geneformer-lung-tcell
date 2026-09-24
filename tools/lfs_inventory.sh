#!/usr/bin/env bash
# Report which git-LFS objects in a clone are present and which are still
# un-pulled pointers.
#
#   bash tools/lfs_inventory.sh <clone> [<clone> ...]
#
# An un-pulled LFS file is a ~130-byte text stub, not the content. Nothing
# fails at checkout; it fails later, in whatever tries to read it, with an
# error about the format rather than about the file being absent --
# `SafetensorError: header too large` for a model, an unpickling error for a
# dictionary. The cost is paid by whoever is furthest from the cause.
#
# Exit status: 0 = no pointers found, 1 = pointers found, 2 = usage error.
# Read-only: this never writes to the clones it inspects and never fetches.
#
# WHY EVERY FILE IS SCANNED, AND NOT JUST THE SMALL ONES
# An LFS pointer is small, so filtering to small files first looks like an
# obvious speed-up. Do not. A prefilter that is even slightly wrong removes
# the cases before the check runs and the result is a confident "none found".
# Reading 45 bytes of every file is cheap enough that the prefilter buys
# nothing worth that risk.
#
# Portability: bash 3.2 (macOS) and 5.x; no `declare -A`, no `find -printf`,
# no `stat` flags, since those differ between GNU and BSD.

set -u

if [ "$#" -eq 0 ]; then
    echo "usage: $(basename "$0") <clone-dir> [<clone-dir> ...]" >&2
    echo "  reports un-pulled git-LFS pointers; read-only, never fetches" >&2
    exit 2
fi

found_any=0

for repo in "$@"; do
    if [ ! -d "$repo" ]; then
        printf '== %s\n   ABSENT\n' "$repo"
        continue
    fi

    total=0
    pointers=0
    listing=""

    while IFS= read -r f; do
        total=$((total + 1))
        # An LFS pointer's first line is exactly this. Match the prefix only;
        # the oid and size follow on later lines.
        if head -c 45 "$f" 2>/dev/null | grep -q '^version https://git-lfs'; then
            pointers=$((pointers + 1))
            want=$(grep -a '^size ' "$f" 2>/dev/null | awk '{print $2}')
            listing="$listing   ${f#"$repo"/}  (expects ${want:-unknown} bytes)
"
        fi
    done <<EOF
$(find "$repo" -path "$repo/.git" -prune -o -type f -print 2>/dev/null)
EOF

    # How many objects this clone actually holds, for context: a clone with
    # pointers and a populated store has fetched some objects and not others,
    # which is a different situation from one that fetched nothing.
    objects=0
    if [ -d "$repo/.git/lfs/objects" ]; then
        objects=$(find "$repo/.git/lfs/objects" -type f 2>/dev/null | wc -l | tr -d ' ')
    fi

    printf '== %s\n' "$repo"
    printf '   files scanned: %s   un-pulled pointers: %s   lfs objects held: %s\n' \
           "$total" "$pointers" "$objects"
    if [ "$pointers" -gt 0 ]; then
        found_any=1
        printf '%s' "$listing"
    fi
done

if [ "$found_any" -ne 0 ]; then
    printf '\nUn-pulled pointers found. To fix, fetch the named paths rather than\n'
    printf 'everything -- `git lfs pull` will also bring down anything else the\n'
    printf 'repository tracks -- then verify each result against the oid recorded\n'
    printf 'in its own pointer file:\n\n'
    printf '  git lfs pull --include "<path>"\n'
    printf '  sha256sum <path>        # must equal the oid the pointer stated\n'
    exit 1
fi

exit 0

#!/usr/bin/env bash
# Sync this repository between the laptop, thinkstation2, and GitHub.
#
# The failure this prevents: committing on the laptop and on thinkstation2
# without pulling in between, which forks history and needs a manual rebase to
# untangle. Run this before you start work on a machine and after you finish.
#
#   ./tools/sync.sh            # rebase onto origin, push, then fast-forward ts2
#   ./tools/sync.sh --check    # report divergence and exit, change nothing
#   ./tools/sync.sh --no-ts2   # GitHub only (use when ts2 is unreachable)
#
# Direct laptop -> ts2 pushes work because ts2 sets
# `receive.denyCurrentBranch=updateInstead`, which updates its working tree in
# place. That requires ts2's tree to be clean; the script checks first and tells
# you rather than failing halfway.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

TS2_HOST="${TS2_HOST:-ts2}"
TS2_PATH="${TS2_PATH:-workspace/geneformer-lung-tcell}"
BRANCH="${BRANCH:-main}"
CHECK_ONLY=0
SKIP_TS2=0

for arg in "$@"; do
    case "$arg" in
        --check) CHECK_ONLY=1 ;;
        --no-ts2) SKIP_TS2=1 ;;
        -h|--help) sed -n '2,18p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) printf 'unknown option: %s\n' "$arg" >&2; exit 2 ;;
    esac
done

say() { printf '\n\033[1m%s\033[0m\n' "$*"; }
warn() { printf '  ! %s\n' "$*" >&2; }

ts2_online() {
    [[ "$SKIP_TS2" -eq 1 ]] && return 1
    ssh -o ConnectTimeout=8 -o BatchMode=yes "$TS2_HOST" true 2>/dev/null
}

say "Local state"
if [[ -n "$(git status --porcelain)" ]]; then
    warn "working tree is dirty; commit or stash before syncing"
    git status --short | sed 's/^/    /'
    exit 1
fi
printf '  %s at %s\n' "$BRANCH" "$(git rev-parse --short HEAD)"

say "Fetching origin"
git fetch --quiet origin
behind="$(git rev-list --count "HEAD..origin/$BRANCH")"
ahead="$(git rev-list --count "origin/$BRANCH..HEAD")"
printf '  behind origin: %s   ahead of origin: %s\n' "$behind" "$ahead"

TS2_UP=0
if ts2_online; then
    TS2_UP=1
    ts2_head="$(ssh -o BatchMode=yes "$TS2_HOST" "cd $TS2_PATH && git rev-parse HEAD" 2>/dev/null || echo unknown)"
    ts2_dirty="$(ssh -o BatchMode=yes "$TS2_HOST" "cd $TS2_PATH && git status --porcelain | head -1" 2>/dev/null || true)"
    printf '  ts2 at %s%s\n' "${ts2_head:0:7}" "$([[ -n "$ts2_dirty" ]] && echo ' (DIRTY)')"
else
    warn "ts2 unreachable or skipped; GitHub-only sync"
fi

if [[ "$CHECK_ONLY" -eq 1 ]]; then
    say "Check only; nothing changed"
    exit 0
fi

if [[ "$behind" -gt 0 ]]; then
    say "Rebasing onto origin/$BRANCH"
    git rebase "origin/$BRANCH"
fi

if [[ "$(git rev-list --count "origin/$BRANCH..HEAD")" -gt 0 ]]; then
    say "Pushing to origin"
    git push origin "$BRANCH"
else
    say "origin already up to date"
fi

if [[ "$TS2_UP" -eq 1 ]]; then
    if [[ -n "${ts2_dirty:-}" ]]; then
        warn "ts2 working tree is dirty; skipping direct push so nothing there is clobbered"
        warn "commit or stash on ts2, then re-run"
    else
        say "Pushing to ts2 (updates its working tree in place)"
        git push "$TS2_HOST" "$BRANCH"
    fi
fi

say "Final state"
printf '  local  %s\n' "$(git rev-parse HEAD)"
printf '  origin %s\n' "$(git rev-parse "origin/$BRANCH")"
if [[ "$TS2_UP" -eq 1 ]]; then
    printf '  ts2    %s\n' "$(ssh -o BatchMode=yes "$TS2_HOST" "cd $TS2_PATH && git rev-parse HEAD")"
fi

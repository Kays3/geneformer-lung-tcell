#!/usr/bin/env bash
# Move the lab's Geneformer data out of the retired per-machine home directory
# into /srv/lab, group-owned by labusers.
#
# Why: /home/thinkstation{1,2} is mode 750 owned by an account that is no longer
# used, so no current user can read it. Copying it into a personal home would
# work today and break again the next time accounts change. A group-owned path
# with setgid survives account churn and lets every labusers member read it.
#
# /srv/lab is LOCAL to each node. Run this on each node separately; the two
# machines hold different arms of the perturbation run, so neither copy is
# redundant and there is no merge between them.
#
#   sudo bash tools/migrate_to_srv_lab.sh --inventory   # read-only, shows what would move
#   sudo bash tools/migrate_to_srv_lab.sh --migrate     # copy, then set group + setgid
#   sudo bash tools/migrate_to_srv_lab.sh --verify      # re-check counts and permissions
#
# The source is never deleted. Remove it yourself once you have verified the
# copy and no longer need the old account.
set -euo pipefail

LAB_ROOT="${LAB_ROOT:-/srv/lab}"
LAB_GROUP="${LAB_GROUP:-labusers}"
# Defaults to the retired account matching this machine's short hostname.
SRC_USER="${SRC_USER:-}"
MODE=""

for arg in "$@"; do
    case "$arg" in
        --inventory) MODE=inventory ;;
        --migrate)   MODE=migrate ;;
        --verify)    MODE=verify ;;
        --src-user=*) SRC_USER="${arg#*=}" ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) printf 'unknown option: %s\n' "$arg" >&2; exit 2 ;;
    esac
done
[[ -z "$MODE" ]] && { printf 'pick one of --inventory | --migrate | --verify\n' >&2; exit 2; }

if [[ -z "$SRC_USER" ]]; then
    for candidate in thinkstation1 thinkstation2; do
        [[ -d "/home/$candidate/workspace" ]] && SRC_USER="$candidate"
    done
fi
[[ -z "$SRC_USER" ]] && { printf 'no retired workspace found; pass --src-user=NAME\n' >&2; exit 1; }

SRC="/home/$SRC_USER/workspace"
say()  { printf '\n\033[1m%s\033[0m\n' "$*"; }
ok()   { printf '  %s\n' "$*"; }
warn() { printf '  ! %s\n' "$*" >&2; }

[[ $EUID -eq 0 ]] || { warn "must run as root (sudo)"; exit 1; }
[[ -d "$SRC" ]] || { warn "source not found: $SRC"; exit 1; }

# Everything the analysis actually needs: the experiment workspace, the model
# weights and venv-adjacent assets, and the spatial raw data that lives inside
# the old repo checkout rather than in KD.
SUBDIRS=(
    "KD"
    "geneformer-uv-starter"
)
# A virtualenv's compiled packages and shebangs embed absolute paths from the
# machine that built it, so copying one to a new prefix produces an interpreter
# that half-works and fails obscurely. Model weights and the Geneformer checkout
# inside the same tree ARE data and must come across, so exclude only the venv
# and rebuild it with geneformer_uv_setup/scripts/bootstrap_workspace.sh.
RSYNC_EXCLUDES=(--exclude '.venv/' --exclude '__pycache__/' --exclude '*.pyc')
EXTRA_PATHS=(
    "/home/$SRC_USER/workspace/geneformer-lung-tcell/sclc_validation/audit/source_metadata/GSE263196_RAW"
)

say "Source: $SRC   (retired account: $SRC_USER)"
say "Target: $LAB_ROOT   group: $LAB_GROUP"

case "$MODE" in
inventory)
    say "Sizes (this is what will be copied)"
    total=0
    for d in "${SUBDIRS[@]}"; do
        if [[ -d "$SRC/$d" ]]; then
            sz="$(du -sb "$SRC/$d" 2>/dev/null | cut -f1)"
            total=$((total + sz))
            printf '  %-28s %8s  %s files\n' "$d" "$(numfmt --to=iec "$sz")" \
                   "$(find "$SRC/$d" -type f 2>/dev/null | wc -l)"
        else
            warn "$d: absent"
        fi
    done
    for p in "${EXTRA_PATHS[@]}"; do
        if [[ -d "$p" ]]; then
            sz="$(du -sb "$p" | cut -f1)"; total=$((total + sz))
            printf '  %-28s %8s\n' "$(basename "$p")" "$(numfmt --to=iec "$sz")"
        fi
    done
    printf '\n  TOTAL %s\n' "$(numfmt --to=iec "$total")"
    say "Free space on $(df -h --output=target "$(dirname "$LAB_ROOT")" | tail -1)"
    df -h "$(dirname "$LAB_ROOT")" | tail -1
    avail="$(df -B1 --output=avail "$(dirname "$LAB_ROOT")" | tail -1)"
    if (( avail < total * 11 / 10 )); then
        warn "less than 110% of the payload free -- do not migrate yet"
    else
        ok "sufficient free space"
    fi
    say "Nothing was changed. Re-run with --migrate to copy."
    ;;

migrate)
    say "Creating $LAB_ROOT"
    install -d -m 2775 -g "$LAB_GROUP" "$LAB_ROOT"
    ok "$(stat -c '%n  %U:%G %a' "$LAB_ROOT")"

    for d in "${SUBDIRS[@]}"; do
        [[ -d "$SRC/$d" ]] || { warn "$d: absent, skipped"; continue; }
        say "Copying $d"
        # -H preserves hardlinks, -A/-X keep ACLs and xattrs, no --delete so the
        # run is idempotent and can never remove anything already in the target.
        rsync -aHAX "${RSYNC_EXCLUDES[@]}" --info=progress2 "$SRC/$d/" "$LAB_ROOT/$d/"
        ok "done: $LAB_ROOT/$d"
    done

    for p in "${EXTRA_PATHS[@]}"; do
        [[ -d "$p" ]] || { warn "$(basename "$p"): absent, skipped"; continue; }
        say "Copying $(basename "$p")"
        install -d -m 2775 -g "$LAB_GROUP" "$LAB_ROOT/spatial_raw"
        rsync -aHAX --info=progress2 "$p/" "$LAB_ROOT/spatial_raw/$(basename "$p")/"
        ok "done: $LAB_ROOT/spatial_raw/$(basename "$p")"
    done

    say "Applying group ownership and setgid"
    chgrp -R "$LAB_GROUP" "$LAB_ROOT"
    chmod -R g+rX "$LAB_ROOT"
    # setgid on directories so anything created later inherits labusers; this is
    # what stops the permission problem recurring.
    find "$LAB_ROOT" -type d -exec chmod g+s {} +
    ok "group=$LAB_GROUP, dirs setgid, group-readable"
    say "Next: rebuild the Python environment"
    ok "the .venv was deliberately NOT copied; rebuild it with"
    ok "  bash geneformer_uv_setup/scripts/bootstrap_workspace.sh"
    ok "then point PYTHON_BIN at the new interpreter in ~/.config/geneformer-lung-tcell/paths.env"
    say "Source left untouched at $SRC -- delete it yourself after verifying"
    ;;

verify)
    say "Verifying $LAB_ROOT"
    fail=0
    for d in "${SUBDIRS[@]}"; do
        [[ -d "$SRC/$d" ]] || continue
        a="$(find "$SRC/$d" -type f -not -path '*/.venv/*' -not -name '*.pyc' 2>/dev/null | wc -l)"
        b="$(find "$LAB_ROOT/$d" -type f 2>/dev/null | wc -l)"
        if [[ "$a" == "$b" ]]; then
            ok "$d: $b files match"
        else
            warn "$d: source $a vs target $b files"; fail=1
        fi
    done
    bad_group="$(find "$LAB_ROOT" ! -group "$LAB_GROUP" -print -quit 2>/dev/null || true)"
    [[ -n "$bad_group" ]] && { warn "not group $LAB_GROUP: $bad_group"; fail=1; } || ok "group ownership uniform"
    nosetgid="$(find "$LAB_ROOT" -type d ! -perm -g+s -print -quit 2>/dev/null || true)"
    [[ -n "$nosetgid" ]] && { warn "directory without setgid: $nosetgid"; fail=1; } || ok "setgid set on all directories"
    say "$([[ $fail -eq 0 ]] && echo 'VERIFIED' || echo 'PROBLEMS FOUND - do not delete the source')"
    exit $fail
    ;;
esac

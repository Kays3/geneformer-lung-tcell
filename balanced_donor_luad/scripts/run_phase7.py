"""Phase 7 launcher: refuses to read any Phase 6 output until every entry in
registration/required_gates.json has its result on disk with the required value
(Amendment 3g.5). A skipped gate is a failed launch, not a later discovery."""
import argparse
import json
import os
import sys


def _get(obj, path):
    for k in path:
        if not isinstance(obj, dict) or k not in obj:
            raise KeyError("/".join(map(str, path)))
        obj = obj[k]
    return obj


def check_gates(base, manifest_path):
    """Return a list of failures (empty = all satisfied)."""
    manifest = json.load(open(manifest_path))
    failures = []
    for g in manifest["gates"]:
        f = os.path.join(base, g["file"])
        if not os.path.exists(f):
            failures.append(f"{g['id']}: missing {g['file']}"); continue
        try:
            if g["type"] == "text_contains":
                if g["text"] not in open(f).read():
                    failures.append(f"{g['id']}: text not found in {g['file']}")
                continue
            data = json.load(open(f))
            if g["type"] == "json_field":
                v = _get(data, g["path"])
                if v != g["equals"]:
                    failures.append(f"{g['id']}: {'/'.join(g['path'])} = {v!r}, required {g['equals']!r}")
            elif g["type"] == "json_field_exists":
                _get(data, g["path"])
            elif g["type"] == "json_equals":
                if data != g["value"]:
                    failures.append(f"{g['id']}: {g['file']} differs from the registered value")
            else:
                failures.append(f"{g['id']}: unknown gate type {g['type']!r}")
        except KeyError as e:
            failures.append(f"{g['id']}: field {e} absent in {g['file']}")
        except json.JSONDecodeError as e:
            failures.append(f"{g['id']}: unparseable {g['file']}: {e}")
    return failures


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", required=True, help="balanced_donor_luad/ directory of the repo")
    p.add_argument("--phase6-root", required=True)
    p.add_argument("--design", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--ovx-index", default=None,
                   help="overexpress position index (Amendment 3h); default phase7_prep/ovx_index.json")
    a = p.parse_args()
    failures = check_gates(a.base, os.path.join(a.base, "registration", "required_gates.json"))
    if failures:
        print("PHASE 7 LAUNCH REFUSED:\n  " + "\n  ".join(failures), file=sys.stderr)
        sys.exit(2)
    end = json.load(open(os.path.join(a.base, "phase7_prep", "probe_end_result.json")))
    host_drift = not (end["delete"]["registered_PASS"] and end["overexpress"]["registered_PASS"])
    sys.argv = [sys.argv[0], "--phase6-root", a.phase6_root, "--design", a.design,
                "--rules", os.path.join(a.base, "registration", "phase7_rules.json"), "--out", a.out,
                "--host-drift", "true" if host_drift else "false",
                "--ovx-index", a.ovx_index or os.path.join(a.base, "phase7_prep", "ovx_index.json")]
    import analyse
    analyse.main()


if __name__ == "__main__":
    main()

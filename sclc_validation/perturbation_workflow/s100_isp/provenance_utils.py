"""Shared SHA-256 helper for the S100 ISP analysis layer.

Same principle as run_targeted_panel.py's runner_provenance()/
input_file_sha256() (added 2026-09-23 after three stale-pinned-hash
incidents): every file this analysis layer reads gets hashed from its own
bytes on disk at the moment it is read, not trusted from a value recorded
earlier. Kept as a tiny standalone helper here rather than imported from
run_targeted_panel.py, because importing that module executes argparse
against sys.argv as a side effect of import (ARGS = parse_args() at module
scope) -- not a safe dependency for a pure analysis library to carry.
"""
from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(path) -> str:
    """SHA-256 of a file's bytes on disk right now. Never raises: a hash
    failure (missing file, permission error) is recorded as an unresolved
    marker so a provenance field can never crash a run that only wanted a
    diagnostic value -- same failure mode as runner_provenance()."""
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except Exception as exc:  # pragma: no cover - diagnostic only
        return f"<unresolved: {exc}>"

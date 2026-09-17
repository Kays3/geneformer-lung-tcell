"""Patch-application idempotence test for apply_bf16_patch.sh.

Runs anywhere (no GPU, no torch, no geneformer needed): builds a throwaway
git repo containing verbatim copies of the two real upstream files this
patch touches (fetched from the pinned f45a6c7 checkout; see
fixtures/README.md), then exercises the apply script exactly the way it
will be used on ts1 -- apply once (should patch), apply again (should
no-op), and confirm a wrong pin is rejected.

pytest sclc_validation/bf16_bench/tests/test_patch_idempotence.py -v
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

BENCH_DIR = Path(__file__).resolve().parent.parent
APPLY_SCRIPT = BENCH_DIR / "apply_bf16_patch.sh"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _make_fake_checkout(tmp_path: Path) -> tuple[Path, str]:
    """A git repo containing verbatim (unpatched) copies of the two real
    files, at HEAD, so apply_bf16_patch.sh's patch/pin logic runs against
    exactly the same content it will see on ts1."""
    repo = tmp_path / "fake_geneformer"
    (repo / "geneformer").mkdir(parents=True)
    shutil.copy(FIXTURES_DIR / "emb_extractor.orig.py", repo / "geneformer" / "emb_extractor.py")
    shutil.copy(FIXTURES_DIR / "perturber_utils.orig.py", repo / "geneformer" / "perturber_utils.py")
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "add", "-A"], cwd=repo, check=True)
    subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-q", "-m", "sim pin"],
        cwd=repo, check=True,
    )
    sim_commit = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], cwd=repo, check=True, capture_output=True, text=True,
    ).stdout.strip()
    return repo, sim_commit


def _run_apply(repo: Path, sim_commit: str) -> subprocess.CompletedProcess:
    # `patch` needs a writable TMPDIR for its scratch file; inherit ours
    # rather than stripping the environment (an `env={...}`-only override
    # left BSD patch falling back to bare /tmp, which this sandbox blocks --
    # not a bug in apply_bf16_patch.sh, just an artifact of over-stripping
    # the test's own subprocess environment).
    env = {**os.environ, "EXPECTED_GENEFORMER_COMMIT": sim_commit}
    return subprocess.run(
        ["bash", str(APPLY_SCRIPT), str(repo)],
        env=env,
        capture_output=True, text=True,
    )


@pytest.fixture
def fake_checkout(tmp_path):
    return _make_fake_checkout(tmp_path)


def test_first_apply_succeeds(fake_checkout):
    repo, sim_commit = fake_checkout
    result = _run_apply(repo, sim_commit)
    assert result.returncode == 0, result.stderr
    assert "Applied" in result.stdout
    patched = (repo / "geneformer" / "perturber_utils.py").read_text()
    assert "torch.cat(tensor_list).float().cpu().numpy()" in patched


def test_second_apply_is_noop(fake_checkout):
    repo, sim_commit = fake_checkout
    first = _run_apply(repo, sim_commit)
    assert first.returncode == 0, first.stderr
    before = (repo / "geneformer" / "perturber_utils.py").read_text()

    second = _run_apply(repo, sim_commit)
    assert second.returncode == 0, second.stderr
    assert "already applied" in second.stdout

    after = (repo / "geneformer" / "perturber_utils.py").read_text()
    assert before == after  # not double-patched


def test_wrong_pin_is_rejected(fake_checkout):
    repo, _sim_commit = fake_checkout
    result = _run_apply(repo, "0000000")
    assert result.returncode != 0
    assert "expected 0000000" in result.stderr

#!/usr/bin/env python3
"""CPU-only planning and validation tests for the T4 GPU runner."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import run_t4_overexpression as runner


class T4RunnerTests(unittest.TestCase):
    def test_cuda_safe_default_nproc_is_one(self):
        self.assertEqual(runner.NPROC, 1)

    def test_nonserial_nproc_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "require --nproc 1"):
            runner.validate_nproc(2)

    def test_stale_raw_outputs_are_cleared(self):
        with tempfile.TemporaryDirectory() as directory:
            raw_dir = Path(directory)
            prefix = "t4_program_exhaustion_normal"
            stale = raw_dir / f"in_silico_overexpress_{prefix}_old_raw.pickle"
            stale.write_bytes(b"partial")
            runner.clear_raw_outputs(raw_dir, prefix)
            self.assertFalse(stale.exists())

    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = runner.load_manifest(runner.MANIFEST)

    def test_full_plan_has_273_units(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            work = runner.select_work(
                self.manifest, Path(directory), "all", list(runner.STATE_NAMES)
            )
        self.assertEqual(len(work), 273)

    def test_primary_single_source_has_four_units(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            work = runner.select_work(self.manifest, Path(directory), "program", ["sclc"])
        self.assertEqual(len(work), 4)
        self.assertTrue(all(source == "sclc" for _, source in work))

    def test_completion_markers_are_respected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            run_dir = Path(directory)
            run = next(item for item in self.manifest["runs"] if item["phase"] == "program")
            marker = runner.completion_marker(run_dir, run, "normal")
            marker.parent.mkdir(parents=True)
            marker.write_text("{}\n")
            work = runner.select_work(
                self.manifest, run_dir, "program", ["normal"], items=[run["id"]]
            )
        self.assertEqual(work, [])

    def test_unknown_item_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "Unknown"):
                runner.select_work(
                    self.manifest, Path(directory), "all", ["normal"], items=["not-a-run"]
                )


if __name__ == "__main__":
    unittest.main()

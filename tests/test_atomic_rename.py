from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ATOMIC_RENAME = ROOT / "scripts" / "atomic_rename.py"


class AtomicRenameTests(unittest.TestCase):
    def run_helper(
        self,
        source: Path,
        target: Path,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                "-B",
                str(ATOMIC_RENAME),
                str(source),
                str(target),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_renames_directory_when_target_is_absent(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "staged"
            target = root / "installed"
            source.mkdir()
            (source / "marker.txt").write_text("verified", encoding="utf-8")

            completed = self.run_helper(source, target)

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertFalse(source.exists())
            self.assertEqual(
                (target / "marker.txt").read_text(encoding="utf-8"),
                "verified",
            )
            result = json.loads(completed.stdout)
            self.assertTrue(result["ok"])
            self.assertIn("evidence", result)
            self.assertIn("limits", result)

    def test_refuses_existing_target_without_touching_source(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "staged"
            target = root / "installed"
            source.mkdir()
            target.mkdir()
            (source / "source.txt").write_text("source", encoding="utf-8")
            (target / "target.txt").write_text("target", encoding="utf-8")

            completed = self.run_helper(source, target)

            self.assertEqual(completed.returncode, 1)
            self.assertTrue((source / "source.txt").is_file())
            self.assertTrue((target / "target.txt").is_file())
            result = json.loads(completed.stdout)
            self.assertFalse(result["ok"])
            self.assertFalse(result["limits"]["overwrite_allowed"])


if __name__ == "__main__":
    unittest.main()

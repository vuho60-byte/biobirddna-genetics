from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERIFY = ROOT / "scripts" / "verify_package.py"


def run_verifier(root: Path = ROOT) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-B", str(VERIFY), "--root", str(root)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


class PackageIntegrityTests(unittest.TestCase):
    def test_release_tree_matches_checksum_manifests(self) -> None:
        completed = run_verifier()

        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(completed.stdout)
        self.assertTrue(summary["ok"])
        self.assertEqual(summary["missing_files"], [])
        self.assertEqual(summary["unexpected_files"], [])
        self.assertEqual(summary["checksum_errors"], [])
        self.assertEqual(summary["version_errors"], [])
        self.assertIn("limits", summary)
        self.assertIn("evidence", summary)

    def test_integrity_checker_rejects_unexpected_bytecode(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = Path(temp_dir) / "skill"
            shutil.copytree(
                ROOT,
                copied_root,
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
            )
            rogue = copied_root / "src" / "rogue.pyc"
            rogue.write_bytes(b"not-release-content")

            completed = run_verifier(copied_root)

        self.assertEqual(completed.returncode, 1)
        summary = json.loads(completed.stdout)
        self.assertEqual(summary["unexpected_files"], ["src/rogue.pyc"])

    def test_non_object_manifest_returns_structured_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = Path(temp_dir) / "skill"
            shutil.copytree(
                ROOT,
                copied_root,
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
            )
            (copied_root / "PACKAGE_CHECKSUMS.json").write_text(
                "[]\n",
                encoding="utf-8",
                newline="\n",
            )

            completed = run_verifier(copied_root)

        self.assertEqual(completed.returncode, 1)
        summary = json.loads(completed.stdout)
        self.assertFalse(summary["ok"])
        self.assertTrue(
            any("root must be an object" in error for error in summary["checksum_errors"])
        )

    def test_integrity_checker_rejects_declared_file_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = Path(temp_dir) / "skill"
            shutil.copytree(
                ROOT,
                copied_root,
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
            )
            readme = copied_root / "README.md"
            readme.write_bytes(readme.read_bytes() + b"\n")

            completed = run_verifier(copied_root)

        self.assertEqual(completed.returncode, 1)
        summary = json.loads(completed.stdout)
        self.assertTrue(
            any(
                error == "README.md: SHA-256 mismatch"
                for error in summary["checksum_errors"]
            )
        )

    def test_integrity_checker_rejects_npm_version_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = Path(temp_dir) / "skill"
            shutil.copytree(
                ROOT,
                copied_root,
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
            )
            package_path = copied_root / "npm" / "package.json"
            package = json.loads(package_path.read_text(encoding="utf-8"))
            package["version"] = "9.9.9"
            package_path.write_text(
                json.dumps(package, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
                newline="\n",
            )

            completed = run_verifier(copied_root)

        self.assertEqual(completed.returncode, 1)
        summary = json.loads(completed.stdout)
        self.assertIn(
            "npm/package.json: expected version '0.2.0', got '9.9.9'",
            summary["version_errors"],
        )


if __name__ == "__main__":
    unittest.main()

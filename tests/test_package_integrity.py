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
REBUILD = ROOT / "scripts" / "rebuild_release_inventory.py"
INSTALL_SH = ROOT / "scripts" / "install_skill.sh"


def run_verifier(root: Path = ROOT) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-B", str(VERIFY), "--root", str(root)],
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=False,
    )


def make_consistent_copy(temp_dir: str) -> Path:
    """Copy the tree without dev caches and rebuild its checksum inventories.

    The result is a self-consistent release tree that verifies cleanly, so
    tests can then pollute it (pytest caches, .git, rogue bytecode) and assert
    exactly how the verifier and installer react.
    """
    copied_root = Path(temp_dir) / "skill"
    shutil.copytree(
        ROOT,
        copied_root,
        ignore=shutil.ignore_patterns(
            "__pycache__", "*.pyc", "*.pyo", ".pytest_cache", ".git"
        ),
    )
    subprocess.run(
        [sys.executable, "-B", str(REBUILD), "--root", str(copied_root)],
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=True,
    )
    return copied_root


def pollute_with_pytest_cache_only(root: Path) -> None:
    """Simulate a pytest run with bytecode suppressed (PYTHONDONTWRITEBYTECODE).

    Leaves only .pytest_cache directories behind — the exact pollution that
    the manifest builder must refuse to bake into release inventories.
    """
    for cache_root in (root / ".pytest_cache", root / "tests" / ".pytest_cache"):
        (cache_root / "v" / "cache").mkdir(parents=True)
        (cache_root / "CACHEDIR.TAG").write_text(
            "Signature: 8a477f597d28d172789f06886806bc55\n", encoding="utf-8"
        )
        (cache_root / "v" / "cache" / "nodeids").write_text("[]\n", encoding="utf-8")


def pollute_with_dev_caches(root: Path) -> None:
    """Simulate what plain `pytest` and `git init` leave behind in a tree.

    Includes __pycache__/*.pyc because any pytest/python run without -B writes
    interpreter bytecode caches into scripts/, src/, and tests/.
    """
    for pycache_root in (
        root / "tests" / "__pycache__",
        root / "src" / "biobirddna_genetics" / "__pycache__",
    ):
        pycache_root.mkdir(parents=True)
        (pycache_root / "x.cpython-311.pyc").write_bytes(b"dev-bytecode")
    pollute_with_pytest_cache_only(root)
    git_dir = root / ".git"
    (git_dir / "refs" / "heads").mkdir(parents=True)
    (git_dir / "HEAD").write_text("ref: refs/heads/main\n", encoding="utf-8")
    (git_dir / "config").write_text("[core]\n", encoding="utf-8")


class PackageIntegrityTests(unittest.TestCase):
    def test_private_metadata_drift_fails_even_after_rehash(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = make_consistent_copy(temp_dir)
            path = root / "MANIFEST.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["npm_distribution"]["private"] = False
            path.write_text(json.dumps(manifest), encoding="utf-8")
            subprocess.run([sys.executable, "-B", str(REBUILD), "--root", str(root)],
                           check=True, capture_output=True)
            completed = run_verifier(root)
            self.assertEqual(completed.returncode, 1)
            self.assertIn("npm_distribution.private: metadata mismatch",
                          json.loads(completed.stdout)["metadata_errors"])

    def test_nested_sync_note_is_not_allowlisted(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = make_consistent_copy(temp_dir)
            (root / "tests/SYNC-NOTE.md").write_text("rogue", encoding="utf-8")
            result = json.loads(run_verifier(root).stdout)
            self.assertFalse(result["ok"])
            self.assertIn("tests/SYNC-NOTE.md", result["unexpected_files"])

    def test_local_sync_note_is_reported_separately_and_not_bundled(self) -> None:
        import hashlib
        with tempfile.TemporaryDirectory() as temp_dir:
            root = make_consistent_copy(temp_dir)
            note = b"Local receipt only\n"
            (root / "SYNC-NOTE.md").write_bytes(note)
            completed = run_verifier(root)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertEqual(json.loads(completed.stdout)["local_notes"], [
                {"path": "SYNC-NOTE.md", "bytes": len(note),
                 "sha256": hashlib.sha256(note).hexdigest(),
                 "trust": "local_untrusted_not_release_content"}
            ])
            subprocess.run([sys.executable, "-B", str(REBUILD), "--root", str(root)],
                           check=True, capture_output=True)
            self.assertNotIn("SYNC-NOTE.md", (root / "SHA256SUMS.txt").read_text())
            (root / "OTHER-NOTE.md").write_text("rogue", encoding="utf-8")
            self.assertEqual(run_verifier(root).returncode, 1)

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
            copied_root = make_consistent_copy(temp_dir)
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
                ignore=shutil.ignore_patterns(
                    "__pycache__", "*.pyc", "*.pyo", ".pytest_cache", ".git"
                ),
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
                ignore=shutil.ignore_patterns(
                    "__pycache__", "*.pyc", "*.pyo", ".pytest_cache", ".git"
                ),
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

    def test_verifier_ignores_pytest_cache_and_git_directories(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = make_consistent_copy(temp_dir)
            pollute_with_dev_caches(copied_root)

            completed = run_verifier(copied_root)

        self.assertEqual(completed.returncode, 0, completed.stdout)
        summary = json.loads(completed.stdout)
        self.assertTrue(summary["ok"])
        self.assertEqual(summary["unexpected_files"], [])
        self.assertEqual(summary["missing_files"], [])
        self.assertEqual(summary["checksum_errors"], [])

    def test_verifier_still_rejects_bytecode_when_caches_present(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = make_consistent_copy(temp_dir)
            pollute_with_dev_caches(copied_root)
            (copied_root / "src" / "rogue.pyc").write_bytes(b"not-release-content")

            completed = run_verifier(copied_root)

        self.assertEqual(completed.returncode, 1)
        summary = json.loads(completed.stdout)
        self.assertEqual(summary["unexpected_files"], ["src/rogue.pyc"])

    def test_manifest_rebuild_never_bakes_pytest_cache_into_inventories(self) -> None:
        """Regression: the builder must be cache-immune, not just the verifier.

        On a tree where pytest ran (only .pytest_cache left behind — bytecode
        suppressed), rebuild_release_inventory.py used to exit 0 and bake
        .pytest_cache/CACHEDIR.TAG and .pytest_cache/v/cache/nodeids into
        SHA256SUMS.txt / PACKAGE_CHECKSUMS.json, after which verify_package.py
        failed. Now the builder refuses (release_temporary_file_forbidden), the
        manifests stay cache-free, and the tree still verifies cleanly.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = make_consistent_copy(temp_dir)
            pollute_with_pytest_cache_only(copied_root)

            rebuilt = subprocess.run(
                [sys.executable, "-B", str(REBUILD), "--root", str(copied_root)],
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )

            self.assertNotEqual(rebuilt.returncode, 0, rebuilt.stdout)
            self.assertIn("release_temporary_file_forbidden", rebuilt.stderr)
            self.assertIn(".pytest_cache", rebuilt.stderr)

            sums_text = (copied_root / "SHA256SUMS.txt").read_text(encoding="utf-8")
            self.assertNotIn(".pytest_cache", sums_text)
            manifest = json.loads(
                (copied_root / "PACKAGE_CHECKSUMS.json").read_text(encoding="utf-8")
            )
            self.assertFalse(
                [
                    entry["path"]
                    for entry in manifest["files"]
                    if ".pytest_cache" in entry["path"]
                ]
            )

            completed = run_verifier(copied_root)

        self.assertEqual(completed.returncode, 0, completed.stdout)
        summary = json.loads(completed.stdout)
        self.assertTrue(summary["ok"])

    @unittest.skipUnless(
        shutil.which("sh"),
        "POSIX sh is not available (Windows uses scripts/install_skill.ps1)",
    )
    def test_installer_succeeds_on_pytest_polluted_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = make_consistent_copy(temp_dir)
            pollute_with_dev_caches(copied_root)

            destination = Path(temp_dir) / "dest"
            completed = subprocess.run(
                ["sh", str(copied_root / "scripts" / "install_skill.sh"),
                 str(destination)],
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )

            self.assertEqual(
                completed.returncode, 0, completed.stdout + completed.stderr
            )
            target = destination / "biobirddna-genetics"
            self.assertTrue((target / "VERSION").is_file())
            self.assertFalse((target / ".pytest_cache").exists())
            self.assertFalse((target / "tests" / ".pytest_cache").exists())
            self.assertFalse((target / ".git").exists())
            self.assertEqual(list(target.rglob("__pycache__")), [])
            self.assertEqual(list(target.rglob("*.pyc")), [])

    def test_integrity_checker_rejects_npm_version_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = Path(temp_dir) / "skill"
            shutil.copytree(
                ROOT,
                copied_root,
                ignore=shutil.ignore_patterns(
                    "__pycache__", "*.pyc", "*.pyo", ".pytest_cache", ".git"
                ),
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
            "npm/package.json: expected version '0.3.0', got '9.9.9'",
            summary["version_errors"],
        )


if __name__ == "__main__":
    unittest.main()

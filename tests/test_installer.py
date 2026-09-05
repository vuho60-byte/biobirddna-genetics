from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts" / "install_skill.ps1"
POWERSHELL = shutil.which("pwsh") or shutil.which("powershell")


@unittest.skipUnless(POWERSHELL, "PowerShell is not installed")
class InstallerTests(unittest.TestCase):
    def test_powershell_installer_creates_verified_clean_target(self) -> None:
        if os.environ.get("BIOBIRDDNA_INSTALLER_CHILD") == "1":
            self.skipTest("installer recursion guard")
        with tempfile.TemporaryDirectory() as temp_dir:
            destination_root = Path(temp_dir) / "Codex skills with spaces"
            child_env = os.environ.copy()
            child_env["BIOBIRDDNA_INSTALLER_CHILD"] = "1"
            completed = subprocess.run(
                [
                    POWERSHELL,
                    "-NoProfile",
                    "-File",
                    str(INSTALLER),
                    "-DestinationRoot",
                    str(destination_root),
                ],
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
                env=child_env,
            )
            target = destination_root / "biobirddna-genetics"

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(target.is_dir())
            self.assertEqual(list(target.rglob("*.pyc")), [])
            self.assertEqual(list(target.rglob("__pycache__")), [])

            verified = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(target / "scripts" / "verify_package.py"),
                    "--root",
                    str(target),
                ],
                cwd=target,
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )

        self.assertEqual(verified.returncode, 0, verified.stderr)
        self.assertTrue(json.loads(verified.stdout)["ok"])

    def test_powershell_installer_fails_if_target_appears_before_finalize(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_root = Path(temp_dir)
            destination_root = temp_root / "Codex skills with spaces"
            signal = temp_root / "verification-complete.signal"
            proceed = temp_root / "continue.signal"
            wrapper = temp_root / "verified-python.cmd"
            wrapper.write_text(
                "\r\n".join(
                    [
                        "@echo off",
                        "python %*",
                        "if errorlevel 1 exit /b %errorlevel%",
                        f'type nul > "{signal}"',
                        ":wait",
                        f'if exist "{proceed}" exit /b 0',
                        "ping 127.0.0.1 -n 1 -w 50 > nul",
                        "goto wait",
                        "",
                    ]
                ),
                encoding="utf-8",
                newline="",
            )
            process = subprocess.Popen(
                [
                    POWERSHELL,
                    "-NoProfile",
                    "-File",
                    str(INSTALLER),
                    "-DestinationRoot",
                    str(destination_root),
                    "-PythonExecutable",
                    str(wrapper),
                ],
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            deadline = time.monotonic() + 10.0
            while not signal.is_file() and time.monotonic() < deadline:
                if process.poll() is not None:
                    break
                time.sleep(0.01)
            if not signal.is_file():
                stdout, stderr = process.communicate(timeout=30)
                self.fail(
                    "installer did not reach finalize gate: "
                    f"returncode={process.returncode}, stdout={stdout!r}, "
                    f"stderr={stderr!r}"
                )

            target = destination_root / "biobirddna-genetics"
            target.mkdir()
            proceed.write_text("continue", encoding="utf-8")
            stdout, stderr = process.communicate(timeout=30)

            self.assertNotEqual(process.returncode, 0, stdout)
            self.assertEqual(list(target.iterdir()), [], stderr)


if __name__ == "__main__":
    unittest.main()

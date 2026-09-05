from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    "SKILL.md",
    "README.md",
    "agents/openai.yaml",
    "pyproject.toml",
    "references/architecture.md",
    "references/data-contracts.md",
    "references/scientific-boundaries.md",
    "references/autoresearch-feedback.md",
    "planning/MODEL_CODER_HANDOFF.md",
    "assets/schemas/feedback-input.schema.json",
    "assets/schemas/feedback-artifact.schema.json",
    "src/biobirddna_genetics/cli.py",
    "npm/package.json",
    "npm/bin/biobirddna-genetics.mjs",
    "npm/bin/biobirddna-feedback.mjs",
    "scripts/atomic_rename.py",
    "scripts/verify_examples.py",
    "scripts/verify_package.py",
    "standalone/biobirddna_genetics_standalone.py",
    "standalone/biobirddna_genetics_node.mjs",
    "knowledge/T02_knowledge_base_v1.0.json",
]


def main() -> int:
    for _stream in (sys.stdout, sys.stderr):
        if hasattr(_stream, "reconfigure"):
            _stream.reconfigure(encoding="utf-8")
    child_env = os.environ.copy()
    child_env["PYTHONDONTWRITEBYTECODE"] = "1"
    missing = [path for path in REQUIRED if not (ROOT / path).is_file()]
    result = {
        "root": str(ROOT),
        "missing_files": missing,
        "tests_passed": False,
        "examples_verified": False,
        "npm_tests_passed": False,
        "npm_pack_verified": False,
        "package_verified": False,
    }
    if not missing:
        completed = subprocess.run(
            [
                sys.executable,
                "-B",
                "-m",
                "unittest",
                "discover",
                "-s",
                str(ROOT / "tests"),
                "-v",
            ],
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            capture_output=True,
            check=False,
            env=child_env,
        )
        result["tests_passed"] = completed.returncode == 0
        result["test_stdout"] = completed.stdout
        result["test_stderr"] = completed.stderr
        examples = subprocess.run(
            [sys.executable, "-B", str(ROOT / "scripts" / "verify_examples.py")],
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            capture_output=True,
            check=False,
            env=child_env,
        )
        result["examples_verified"] = examples.returncode == 0
        result["examples_stdout"] = examples.stdout
        result["examples_stderr"] = examples.stderr
        npm_command = shutil.which("npm")
        if npm_command:
            npm_tests = subprocess.run(
                [npm_command, "test"],
                cwd=ROOT / "npm",
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
                env=child_env,
            )
            result["npm_tests_passed"] = npm_tests.returncode == 0
            result["npm_test_stdout"] = npm_tests.stdout
            result["npm_test_stderr"] = npm_tests.stderr
            npm_pack = subprocess.run(
                [npm_command, "pack", "--dry-run", "--json"],
                cwd=ROOT / "npm",
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
                env=child_env,
            )
            result["npm_pack_verified"] = npm_pack.returncode == 0
            result["npm_pack_stdout"] = npm_pack.stdout
            result["npm_pack_stderr"] = npm_pack.stderr
        else:
            result["npm_test_stderr"] = "npm_not_found"
            result["npm_pack_stderr"] = "npm_not_found"
        package = subprocess.run(
            [
                sys.executable,
                "-B",
                str(ROOT / "scripts" / "verify_package.py"),
                "--root",
                str(ROOT),
            ],
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            capture_output=True,
            check=False,
            env=child_env,
        )
        result["package_verified"] = package.returncode == 0
        result["package_stdout"] = package.stdout
        result["package_stderr"] = package.stderr
    result["evidence"] = {
        "tests": "python -B -m unittest discover -s tests -v",
        "examples": "scripts/verify_examples.py",
        "npm_tests": "npm test",
        "npm_pack": "npm pack --dry-run --json",
        "package": "scripts/verify_package.py",
    }
    result["limits"] = {
        "purpose": "development_and_release_validation",
        "scientific_claim": "none",
        "gene_status": "not_assessed",
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return (
        0
        if (
            not missing
            and result["tests_passed"]
            and result["examples_verified"]
            and result["npm_tests_passed"]
            and result["npm_pack_verified"]
            and result["package_verified"]
        )
        else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tomllib
from pathlib import Path
from typing import Any

SELF_EXCLUDED = "PACKAGE_CHECKSUMS.json"
SUMS_FILE = "SHA256SUMS.txt"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative_files(root: Path) -> dict[str, Path]:
    return {
        path.relative_to(root).as_posix(): path
        for path in root.rglob("*")
        if path.is_file()
    }


def _read_sha256sums(path: Path) -> tuple[dict[str, str], list[str]]:
    entries: dict[str, str] = {}
    errors: list[str] = []
    if not path.is_file():
        return entries, [f"{SUMS_FILE}: missing"]
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        match = re.fullmatch(r"([0-9a-fA-F]{64})  \./(.+)", line)
        if match is None:
            errors.append(f"{SUMS_FILE}:{line_number}: invalid entry")
            continue
        relative_path = Path(match.group(2)).as_posix()
        if relative_path in entries:
            errors.append(f"{SUMS_FILE}:{line_number}: duplicate {relative_path}")
            continue
        entries[relative_path] = match.group(1).lower()
    return entries, errors


def _version_errors(
    root: Path,
    checksum_manifest: dict[str, Any],
) -> list[str]:
    errors: list[str] = []
    try:
        expected = (root / "VERSION").read_text(encoding="utf-8").strip()
    except OSError as exc:
        return [f"VERSION: {exc}"]

    candidates: dict[str, str | None] = {
        "PACKAGE_CHECKSUMS.json": checksum_manifest.get("version"),
    }
    try:
        candidates["MANIFEST.json"] = json.loads(
            (root / "MANIFEST.json").read_text(encoding="utf-8")
        ).get("version")
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"MANIFEST.json: {exc}")
    try:
        candidates["pyproject.toml"] = tomllib.loads(
            (root / "pyproject.toml").read_text(encoding="utf-8")
        )["project"]["version"]
    except (OSError, KeyError, tomllib.TOMLDecodeError) as exc:
        errors.append(f"pyproject.toml: {exc}")
    try:
        candidates["npm/package.json"] = json.loads(
            (root / "npm" / "package.json").read_text(encoding="utf-8")
        ).get("version")
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"npm/package.json: {exc}")

    version_sources = {
        "src/biobirddna_genetics/__init__.py": r'__version__\s*=\s*"([^"]+)"',
        "src/biobirddna_genetics/cli.py": r'"engine_version":\s*"([^"]+)"',
    }
    for relative_path, pattern in version_sources.items():
        try:
            source = (root / relative_path).read_text(encoding="utf-8")
        except OSError as exc:
            errors.append(f"{relative_path}: {exc}")
            continue
        match = re.search(pattern, source)
        if match is None:
            errors.append(f"{relative_path}: version marker missing")
        else:
            candidates[relative_path] = match.group(1)

    for source, version in candidates.items():
        if version != expected:
            errors.append(f"{source}: expected version {expected!r}, got {version!r}")
    return errors


def verify(root: Path) -> dict[str, Any]:
    root = root.resolve()
    actual_files = _relative_files(root)
    checksum_errors: list[str] = []
    version_errors: list[str] = []

    manifest_path = root / SELF_EXCLUDED
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        manifest = {"files": []}
        checksum_errors.append(f"{SELF_EXCLUDED}: {exc}")
    if not isinstance(manifest, dict):
        manifest = {"files": []}
        checksum_errors.append(f"{SELF_EXCLUDED}: root must be an object")

    declared: dict[str, dict[str, Any]] = {}
    raw_entries = manifest.get("files", [])
    if not isinstance(raw_entries, list):
        checksum_errors.append(f"{SELF_EXCLUDED}: files must be an array")
        raw_entries = []
    for index, entry in enumerate(raw_entries):
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            checksum_errors.append(f"{SELF_EXCLUDED}: files[{index}] is invalid")
            continue
        relative_path = Path(entry["path"]).as_posix()
        if (
            Path(relative_path).is_absolute()
            or ".." in Path(relative_path).parts
            or relative_path == SELF_EXCLUDED
        ):
            checksum_errors.append(
                f"{SELF_EXCLUDED}: invalid path {entry['path']!r}"
            )
            continue
        if relative_path in declared:
            checksum_errors.append(f"{SELF_EXCLUDED}: duplicate {relative_path}")
            continue
        declared[relative_path] = entry

    expected_actual = set(actual_files) - {SELF_EXCLUDED}
    declared_paths = set(declared)
    missing_files = sorted(declared_paths - set(actual_files))
    unexpected_files = sorted(expected_actual - declared_paths)

    for relative_path in sorted(declared_paths & set(actual_files)):
        path = actual_files[relative_path]
        entry = declared[relative_path]
        actual_bytes = path.stat().st_size
        actual_hash = _sha256(path)
        if entry.get("bytes") != actual_bytes:
            checksum_errors.append(
                f"{relative_path}: expected {entry.get('bytes')!r} bytes, "
                f"got {actual_bytes}"
            )
        if str(entry.get("sha256", "")).lower() != actual_hash:
            checksum_errors.append(f"{relative_path}: SHA-256 mismatch")

    declared_total = sum(
        entry.get("bytes", 0)
        for entry in declared.values()
        if isinstance(entry.get("bytes"), int)
    )
    if manifest.get("file_count") != len(declared):
        checksum_errors.append(
            f"{SELF_EXCLUDED}: file_count does not match files array"
        )
    if manifest.get("total_bytes") != declared_total:
        checksum_errors.append(
            f"{SELF_EXCLUDED}: total_bytes does not match files array"
        )
    if manifest.get("algorithm") != "SHA-256":
        checksum_errors.append(f"{SELF_EXCLUDED}: algorithm must be SHA-256")
    if manifest.get("self_excluded") != SELF_EXCLUDED:
        checksum_errors.append(
            f"{SELF_EXCLUDED}: self_excluded must name {SELF_EXCLUDED}"
        )

    sums, sums_errors = _read_sha256sums(root / SUMS_FILE)
    checksum_errors.extend(sums_errors)
    expected_sums = set(actual_files) - {SELF_EXCLUDED, SUMS_FILE}
    if set(sums) != expected_sums:
        for relative_path in sorted(expected_sums - set(sums)):
            checksum_errors.append(f"{SUMS_FILE}: missing {relative_path}")
        for relative_path in sorted(set(sums) - expected_sums):
            checksum_errors.append(f"{SUMS_FILE}: unexpected {relative_path}")
    for relative_path in sorted(set(sums) & set(actual_files)):
        if sums[relative_path] != _sha256(actual_files[relative_path]):
            checksum_errors.append(f"{SUMS_FILE}: SHA-256 mismatch for {relative_path}")

    if manifest_path.is_file():
        version_errors = _version_errors(root, manifest)

    ok = not (missing_files or unexpected_files or checksum_errors or version_errors)
    return {
        "ok": ok,
        "root": str(root),
        "files_verified": len(declared_paths & set(actual_files)),
        "missing_files": missing_files,
        "unexpected_files": unexpected_files,
        "checksum_errors": checksum_errors,
        "version_errors": version_errors,
        "evidence": {
            "checksum_manifest": SELF_EXCLUDED,
            "checksum_list": SUMS_FILE,
            "algorithm": "SHA-256",
        },
        "limits": {
            "purpose": "package_integrity_only",
            "scientific_claim": "none",
            "gene_status": "not_assessed",
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify the exact BioBirdDNA Genetics release tree."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    result = verify(build_parser().parse_args(argv).root)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

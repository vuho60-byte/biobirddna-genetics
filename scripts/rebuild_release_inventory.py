from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path


CHECKSUM_MANIFEST = "PACKAGE_CHECKSUMS.json"
SUMS_FILE = "SHA256SUMS.txt"
GENERATED = {CHECKSUM_MANIFEST, SUMS_FILE}
FORBIDDEN_PARTS = {
    ".biobirddna-cloud",
    ".cache",
    ".git",
    ".pytest_cache",
    "__pycache__",
    "node_modules",
}
FORBIDDEN_SUFFIXES = {".bak", ".lock", ".pyc", ".pyo", ".temp", ".tmp"}


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _release_files(root: Path) -> list[tuple[str, Path]]:
    files: list[tuple[str, Path]] = []
    for path in root.rglob("*"):
        # VCS state and the single local receipt are not release content.
        relative = path.relative_to(root).as_posix()
        if ".git" in path.relative_to(root).parts or relative == "SYNC-NOTE.md":
            continue
        if path.is_symlink():
            raise ValueError(f"release_symlink_forbidden:{path}")
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        parts = Path(relative).parts
        if (
            any(part in FORBIDDEN_PARTS or part.startswith("stage-") for part in parts)
            or path.suffix.lower() in FORBIDDEN_SUFFIXES
            or path.name.endswith("~")
        ):
            raise ValueError(f"release_temporary_file_forbidden:{relative}")
        files.append((relative, path))
    return sorted(files, key=lambda item: item[0].encode("utf-8"))


def _render(root: Path, package_name: str) -> tuple[bytes, bytes, dict]:
    source_files = [
        (relative, path.read_bytes())
        for relative, path in _release_files(root)
        if relative not in GENERATED
    ]
    sums_bytes = b"".join(
        f"{_sha256(data)}  ./{relative}\n".encode("utf-8")
        for relative, data in source_files
    )
    manifest_files = sorted(
        [*source_files, (SUMS_FILE, sums_bytes)],
        key=lambda item: item[0].encode("utf-8"),
    )
    entries = [
        {
            "path": relative,
            "bytes": len(data),
            "sha256": _sha256(data),
        }
        for relative, data in manifest_files
    ]
    manifest = {
        "package": package_name,
        "version": (root / "VERSION").read_text(encoding="utf-8").strip(),
        "algorithm": "SHA-256",
        "self_excluded": CHECKSUM_MANIFEST,
        "file_count": len(entries),
        "total_bytes": sum(entry["bytes"] for entry in entries),
        "files": entries,
    }
    manifest_bytes = (
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    ).encode("utf-8")
    return sums_bytes, manifest_bytes, manifest


def _write_atomic(path: Path, data: bytes) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def rebuild(root: Path, package_name: str | None = None) -> dict:
    root = root.resolve()
    version = (root / "VERSION").read_text(encoding="utf-8").strip()
    resolved_name = (
        package_name
        if package_name is not None
        else f"BioBirdDNA_Genetics_Skill_v{version}"
    )
    first = _render(root, resolved_name)
    second = _render(root, resolved_name)
    if first[:2] != second[:2]:
        raise RuntimeError("release_inventory_nondeterministic")
    _write_atomic(root / SUMS_FILE, first[0])
    _write_atomic(root / CHECKSUM_MANIFEST, first[1])
    verified = _render(root, resolved_name)
    if first[:2] != verified[:2]:
        raise RuntimeError("release_inventory_postwrite_drift")
    return first[2]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild deterministic BioBirdDNA release inventories."
    )
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--package-name")
    args = parser.parse_args()
    manifest = rebuild(args.root, package_name=args.package_name)
    print(json.dumps(manifest, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

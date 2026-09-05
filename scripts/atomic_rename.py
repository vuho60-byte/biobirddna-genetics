from __future__ import annotations

import argparse
import ctypes
import errno
import json
import os
import sys
from pathlib import Path

AT_FDCWD = -100
RENAME_NOREPLACE = 1
RENAME_EXCL = 0x00000004


def _raise_last_os_error(operation: str) -> None:
    error_number = ctypes.get_errno()
    if error_number == 0:
        error_number = errno.EIO
    raise OSError(error_number, f"{operation}: {os.strerror(error_number)}")


def mechanism() -> str:
    if os.name == "nt":
        return "windows_rename_noreplace"
    if sys.platform.startswith("linux"):
        return "linux_renameat2_noreplace"
    if sys.platform == "darwin":
        return "darwin_renamex_np_excl"
    return "unsupported"


def atomic_rename_noreplace(source: Path, target: Path) -> None:
    source_bytes = os.fsencode(source)
    target_bytes = os.fsencode(target)

    if os.name == "nt":
        os.rename(source, target)
        return

    libc = ctypes.CDLL(None, use_errno=True)
    if sys.platform.startswith("linux"):
        renameat2 = getattr(libc, "renameat2", None)
        if renameat2 is None:
            raise OSError(
                errno.ENOSYS,
                "renameat2(RENAME_NOREPLACE) is unavailable; refusing unsafe fallback",
            )
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        ctypes.set_errno(0)
        if (
            renameat2(
                AT_FDCWD,
                source_bytes,
                AT_FDCWD,
                target_bytes,
                RENAME_NOREPLACE,
            )
            != 0
        ):
            _raise_last_os_error("renameat2(RENAME_NOREPLACE)")
        return

    if sys.platform == "darwin":
        renamex_np = getattr(libc, "renamex_np", None)
        if renamex_np is None:
            raise OSError(
                errno.ENOSYS,
                "renamex_np(RENAME_EXCL) is unavailable; refusing unsafe fallback",
            )
        renamex_np.argtypes = [
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renamex_np.restype = ctypes.c_int
        ctypes.set_errno(0)
        if renamex_np(source_bytes, target_bytes, RENAME_EXCL) != 0:
            _raise_last_os_error("renamex_np(RENAME_EXCL)")
        return

    raise OSError(
        errno.ENOTSUP,
        f"atomic no-replace rename is unsupported on {sys.platform}",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Atomically rename a verified directory without replacement."
    )
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    source = Path(os.path.abspath(args.source))
    target = Path(os.path.abspath(args.target))
    result = {
        "ok": False,
        "source": str(source),
        "target": str(target),
        "evidence": {
            "mechanism": mechanism(),
            "platform": sys.platform,
        },
        "limits": {
            "purpose": "installer_finalization_only",
            "overwrite_allowed": False,
            "same_filesystem_required": True,
            "scientific_claim": "none",
            "gene_status": "not_assessed",
        },
    }
    try:
        if not source.is_dir():
            raise FileNotFoundError(f"staged directory missing: {source}")
        atomic_rename_noreplace(source, target)
    except (OSError, ValueError) as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 1

    result["ok"] = True
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

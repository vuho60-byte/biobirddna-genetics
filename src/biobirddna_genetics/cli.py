from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

from .errors import BioBirdDNAError
from .reporting import write_json
from .service import dispatch


def _reject_json_constant(value: str) -> None:
    raise BioBirdDNAError(f"input_json_non_finite:{value}")


def _reject_non_finite_numbers(value: Any) -> None:
    if isinstance(value, dict):
        for child in value.values():
            _reject_non_finite_numbers(child)
    elif isinstance(value, list):
        for child in value:
            _reject_non_finite_numbers(child)
    elif (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and not math.isfinite(float(value))
    ):
        raise BioBirdDNAError("input_json_non_finite:number")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=_reject_json_constant,
        )
    except FileNotFoundError as exc:
        raise BioBirdDNAError(f"input_file_missing:{path}") from exc
    except json.JSONDecodeError as exc:
        raise BioBirdDNAError(f"input_json_invalid:{path}:{exc.lineno}") from exc
    if not isinstance(value, dict):
        raise BioBirdDNAError("input_root_must_be_object")
    _reject_non_finite_numbers(value)
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="biobirddna-genetics",
        description=(
            "BioBirdDNA cross simulation, inheritance-model inference, "
            "mate ranking, breeding-record normalization, and phenotype-only "
            "pigment screening."
        ),
    )
    parser.add_argument(
        "operation",
        choices=[
            "cross",
            "infer",
            "pigment",
            "rank-mates",
            "normalize-records",
            "simulate-recovery",
        ],
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON to stdout when --output is omitted.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        request = _read_json(args.input)
        result = dispatch(args.operation, request)
    except (BioBirdDNAError, KeyError, TypeError, ValueError) as exc:
        error = {
            "ok": False,
            "operation": args.operation,
            "error": str(exc),
        }
        print(
            json.dumps(
                error,
                ensure_ascii=False,
                indent=2,
                allow_nan=False,
            ),
            file=sys.stderr,
        )
        return 2

    envelope = {
        "ok": True,
        "operation": args.operation,
        "engine_version": "0.2.0",
        "result": result,
    }
    if args.output:
        write_json(args.output, envelope)
        print(str(args.output))
    else:
        indent = 2 if args.pretty else None
        print(
            json.dumps(
                envelope,
                ensure_ascii=False,
                indent=indent,
                allow_nan=False,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

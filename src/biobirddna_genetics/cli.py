from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
from pathlib import Path
from typing import Any

from .errors import BioBirdDNAError
from .reporting import write_json
from .service import dispatch

EXIT_INPUT_ERROR = 2
EXIT_INTERNAL_ERROR = 70

# Human-readable Vietnamese messages for the most common machine error slugs.
# Keys match the slug before the first ":". The machine `error` field stays
# unchanged; `message_vi` is additive. Mirrored verbatim in
# npm/bin/biobirddna-genetics.mjs — keep both maps in sync.
MESSAGES_VI = {
    "input_file_missing": (
        "Không tìm thấy file input — kiểm tra lại đường dẫn."
    ),
    "input_json_invalid": (
        "File input không phải JSON hợp lệ — kiểm tra lại cú pháp."
    ),
    "input_json_non_finite": (
        "Input chứa số không hữu hạn (NaN/Infinity) — chỉ dùng số thường."
    ),
    "input_root_must_be_object": (
        "Nội dung JSON phải là một object ở cấp cao nhất."
    ),
    "operation_unknown": (
        "Thao tác không được hỗ trợ — dùng: cross, infer, pigment, "
        "rank-mates, normalize-records, simulate-recovery."
    ),
    "male_and_female_genotype_required": (
        "Thiếu kiểu gene bố hoặc mẹ — cần cả male_genotype và "
        "female_genotype."
    ),
    "offspring_count_invalid": (
        "Số lượng chim con không hợp lệ — dùng số nguyên không âm."
    ),
    "offspring_sex_invalid": (
        "Giới tính chim con không hợp lệ — dùng male, female hoặc unknown."
    ),
    "internal_error": (
        "Lỗi nội bộ của công cụ, không phải do dữ liệu bạn nhập — báo lại "
        "cho người phát triển kèm traceback."
    ),
}
MESSAGE_VI_FALLBACK = (
    "Dữ liệu đầu vào không hợp lệ — xem mã lỗi ở trường error và đối chiếu "
    "references/data-contracts.md."
)


def message_vi(code: str) -> str:
    return MESSAGES_VI.get(code.split(":", 1)[0], MESSAGE_VI_FALLBACK)


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
    for _stream in (sys.stdout, sys.stderr):
        if hasattr(_stream, "reconfigure"):
            _stream.reconfigure(encoding="utf-8")
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        request = _read_json(args.input)
        result = dispatch(args.operation, request)
    except BioBirdDNAError as exc:
        code = str(exc)
        error = {
            "ok": False,
            "operation": args.operation,
            "error": code,
            "message_vi": message_vi(code),
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
        return EXIT_INPUT_ERROR
    except Exception as exc:
        # Internal bug (TypeError/KeyError/...): NOT an input problem.
        # Preserve the traceback and use a distinct exit code so wrappers
        # never present engine bugs as user input errors.
        traceback.print_exc(file=sys.stderr)
        code = f"internal_error:{type(exc).__name__}"
        error = {
            "ok": False,
            "operation": args.operation,
            "error": code,
            "message_vi": message_vi(code),
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
        return EXIT_INTERNAL_ERROR

    envelope = {
        "ok": True,
        "operation": args.operation,
        "engine_version": "0.4.0",
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

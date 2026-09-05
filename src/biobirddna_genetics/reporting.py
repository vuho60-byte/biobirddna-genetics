from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def write_json(path: str | Path, value: Any) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(
            value,
            ensure_ascii=False,
            indent=2,
            sort_keys=False,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )


def markdown_cross_report(result: dict[str, Any]) -> str:
    lines = [
        f"# Kết quả mô phỏng ghép — {result['model_id']}",
        "",
        f"- Mô hình: `{result['mode']}`",
        f"- Kiểu gene trống: `{result.get('male_genotype', 'uncertain')}`",
        f"- Kiểu gene mái: `{result.get('female_genotype', 'uncertain')}`",
        "",
        "## Tỷ lệ kiểu hình đời con",
        "",
        "| Nhóm | Xác suất |",
        "| --- | ---: |",
    ]
    for key, probability in result["phenotype_probabilities"].items():
        lines.append(f"| `{key}` | {probability:.2%} |")
    lines.extend(
        [
            "",
            "## Giới hạn",
            "",
            "Đây là mô phỏng theo giả định mô hình và kiểu gene đầu vào; không phải xác nhận gene hoặc cam kết tỷ lệ của một lứa nhỏ.",
            "",
        ]
    )
    return "\n".join(lines)

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "scripts" / "biobirddna_cli.py"
CASES = [
    ("cross", "cross_autosomal_recessive.json"),
    ("cross", "cross_z_linked_recessive.json"),
    ("infer", "infer_z_linked_case.json"),
    ("pigment", "pigment_progressive_greying.json"),
    ("rank-mates", "mate_ranking_case.json"),
    ("normalize-records", "normalize_breeding_records.json"),
    ("simulate-recovery", "simulate_recovery.json"),
]


def main() -> int:
    for operation, filename in CASES:
        input_path = ROOT / "examples" / "input" / filename
        output_path = ROOT / "examples" / "output" / filename
        print(f"== {operation}: {filename} ==")
        completed = subprocess.run(
            [
                sys.executable,
                "-B",
                str(CLI),
                operation,
                "--input",
                str(input_path),
                "--output",
                str(output_path),
            ],
            cwd=ROOT,
            check=False,
        )
        if completed.returncode != 0:
            return completed.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

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
TOLERANCE = 1e-9
GOLDEN_REL_TOL = 1e-9
GOLDEN_ABS_TOL = 1e-12


def _trees_equivalent(actual: Any, golden: Any) -> bool:
    if isinstance(actual, bool) or isinstance(golden, bool):
        return actual is golden
    if isinstance(actual, (int, float)) and isinstance(golden, (int, float)):
        return math.isclose(
            actual, golden, rel_tol=GOLDEN_REL_TOL, abs_tol=GOLDEN_ABS_TOL
        )
    if type(actual) is not type(golden):
        return False
    if isinstance(actual, dict):
        if actual.keys() != golden.keys():
            return False
        return all(_trees_equivalent(actual[key], golden[key]) for key in actual)
    if isinstance(actual, list):
        if len(actual) != len(golden):
            return False
        return all(
            _trees_equivalent(actual_item, golden_item)
            for actual_item, golden_item in zip(actual, golden)
        )
    return actual == golden


def golden_matches(actual_bytes: bytes, golden_path: Path) -> bool:
    if not golden_path.is_file():
        return False
    golden_bytes = golden_path.read_bytes()
    if actual_bytes == golden_bytes:
        return True
    try:
        golden_tree = json.loads(golden_bytes)
        actual_tree = json.loads(actual_bytes)
    except (json.JSONDecodeError, UnicodeDecodeError):
        return False
    return _trees_equivalent(actual_tree, golden_tree)


def _distribution_error(path: str, values: list[Any]) -> str | None:
    if not values:
        return f"{path}: probability distribution must not be empty"
    if any(
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not math.isfinite(value)
        or value < 0.0
        for value in values
    ):
        return f"{path}: probabilities must be finite, non-negative numbers"
    total = math.fsum(values)
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=TOLERANCE):
        return f"{path}: probability total {total!r} differs from 1"
    return None


def _probability_checks(value: Any, path: str = "$") -> list[tuple[str, str | None]]:
    checks: list[tuple[str, str | None]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_path = f"{path}.{key}"
            if (
                isinstance(child, dict)
                and (
                    key.endswith("_probabilities")
                    or key.endswith("_posterior")
                    or key == "model_posteriors"
                    or key == "scores"
                )
            ):
                checks.append(
                    (child_path, _distribution_error(child_path, list(child.values())))
                )
            if key == "ranking" and isinstance(child, list) and child:
                posteriors = [
                    item["posterior"]
                    for item in child
                    if isinstance(item, dict) and "posterior" in item
                ]
                if len(posteriors) == len(child):
                    checks.append(
                        (
                            f"{child_path}[*].posterior",
                            _distribution_error(
                                f"{child_path}[*].posterior",
                                posteriors,
                            ),
                        )
                    )
            if key == "configurations" and isinstance(child, list) and child:
                posteriors = [
                    item["posterior_within_model"]
                    for item in child
                    if isinstance(item, dict) and "posterior_within_model" in item
                ]
                if len(posteriors) == len(child):
                    checks.append(
                        (
                            f"{child_path}[*].posterior_within_model",
                            _distribution_error(
                                f"{child_path}[*].posterior_within_model",
                                posteriors,
                            ),
                        )
                    )
            checks.extend(_probability_checks(child, child_path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            checks.extend(_probability_checks(child, f"{path}[{index}]"))
    return checks


def find_probability_errors(value: Any) -> list[str]:
    return [error for _, error in _probability_checks(value) if error is not None]


def find_non_finite_errors(value: Any, path: str = "$") -> list[str]:
    errors: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            errors.extend(find_non_finite_errors(child, f"{path}.{key}"))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            errors.extend(find_non_finite_errors(child, f"{path}[{index}]"))
    elif (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and not math.isfinite(value)
    ):
        errors.append(f"{path}: value must be a finite JSON number")
    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify example goldens, probability invariants, and byte stability."
    )
    parser.add_argument(
        "--golden-dir",
        type=Path,
        default=ROOT / "examples" / "output",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    for _stream in (sys.stdout, sys.stderr):
        if hasattr(_stream, "reconfigure"):
            _stream.reconfigure(encoding="utf-8")
    args = build_parser().parse_args(argv)
    golden_dir = args.golden_dir.resolve()
    golden_mismatches: list[str] = []
    unstable_examples: list[str] = []
    execution_errors: list[str] = []
    probability_errors: list[str] = []
    numeric_errors: list[str] = []
    distributions_checked = 0
    verified_examples = 0

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_root = Path(temp_dir)
        first_dir = temp_root / "first"
        second_dir = temp_root / "second"
        first_dir.mkdir()
        second_dir.mkdir()

        for operation, filename in CASES:
            generated: list[Path] = []
            case_failed = False
            for run_dir in (first_dir, second_dir):
                output_path = run_dir / filename
                completed = subprocess.run(
                    [
                        sys.executable,
                        "-B",
                        str(CLI),
                        operation,
                        "--input",
                        str(ROOT / "examples" / "input" / filename),
                        "--output",
                        str(output_path),
                    ],
                    cwd=ROOT,
                    text=True,
                    encoding="utf-8",
                    capture_output=True,
                    check=False,
                )
                if completed.returncode != 0:
                    execution_errors.append(
                        f"{filename}: CLI returned {completed.returncode}: "
                        f"{completed.stderr.strip()}"
                    )
                    case_failed = True
                    break
                generated.append(output_path)

            if case_failed:
                continue

            first_bytes = generated[0].read_bytes()
            second_bytes = generated[1].read_bytes()
            if first_bytes != second_bytes:
                unstable_examples.append(filename)

            golden_path = golden_dir / filename
            if not golden_matches(first_bytes, golden_path):
                golden_mismatches.append(filename)

            payload = json.loads(first_bytes)
            checks = _probability_checks(payload)
            distributions_checked += len(checks)
            case_probability_errors = [
                f"{filename}:{error}"
                for _, error in checks
                if error is not None
            ]
            probability_errors.extend(case_probability_errors)
            case_numeric_errors = [
                f"{filename}:{error}"
                for error in find_non_finite_errors(payload)
            ]
            numeric_errors.extend(case_numeric_errors)

            if (
                filename not in unstable_examples
                and filename not in golden_mismatches
                and not case_probability_errors
                and not case_numeric_errors
            ):
                verified_examples += 1

    summary = {
        "ok": not (
            execution_errors
            or golden_mismatches
            or unstable_examples
            or probability_errors
            or numeric_errors
        ),
        "examples_checked": len(CASES),
        "verified_examples": verified_examples,
        "distributions_checked": distributions_checked,
        "byte_stable": not unstable_examples,
        "golden_mismatches": golden_mismatches,
        "unstable_examples": unstable_examples,
        "probability_errors": probability_errors,
        "numeric_errors": numeric_errors,
        "execution_errors": execution_errors,
        "evidence": {
            "generator": "scripts/biobirddna_cli.py",
            "input_directory": "examples/input",
            "golden_directory": str(golden_dir),
            "examples": [filename for _, filename in CASES],
            "golden_comparison": {
                "mode": "parsed_json_tree",
                "fast_path": "byte_equality",
                "numeric_rel_tol": GOLDEN_REL_TOL,
                "numeric_abs_tol": GOLDEN_ABS_TOL,
                "bool_comparison": "identity_before_numeric",
                "exact_fields": [
                    "key_sets",
                    "types",
                    "string_values",
                    "null",
                    "list_lengths_and_order",
                ],
                "byte_stability": "same_machine_double_run_byte_exact",
            },
        },
        "limits": {
            "purpose": "release_verification_only",
            "scientific_claim": "none",
            "gene_status": "not_assessed",
        },
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

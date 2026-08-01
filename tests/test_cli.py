from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "scripts" / "biobirddna_cli.py"


class CliTests(unittest.TestCase):
    def test_cross_cli_writes_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "result.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(CLI),
                    "cross",
                    "--input",
                    str(ROOT / "examples" / "input" / "cross_autosomal_recessive.json"),
                    "--output",
                    str(output),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(payload["ok"])
            self.assertEqual(payload["operation"], "cross")

    def test_normalize_records_cli_preserves_observed_zero(self) -> None:
        request = {
            "schema_version": "1.0.0",
            "case_id": "CM-CASE-CLI",
            "pairings": [
                {
                    "pairing_id": "CM-MAT-CLI",
                    "male_parent_id": "CM-IND-M",
                    "female_parent_id": "CM-IND-F",
                    "male_parent_phenotype": "normal",
                    "female_parent_phenotype": "normal",
                    "male_parent_sex_method": "dna",
                    "male_parent_sex_confidence": 0.99,
                    "female_parent_sex_method": "dna",
                    "female_parent_sex_confidence": 0.99,
                }
            ],
            "clutches": [
                {
                    "clutch_id": "CM-CLU-CLI",
                    "pairing_id": "CM-MAT-CLI",
                    "offspring_records_coverage": "complete",
                    "stage_counts": {
                        stage: {
                            "status": "observed",
                            "coverage": "complete",
                            "count": 0,
                        }
                        for stage in (
                            "eggs_laid",
                            "fertile_eggs",
                            "hatched",
                            "survived_to_first_molt",
                        )
                    },
                }
            ],
            "offspring": [],
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "records.json"
            output_path = Path(temp_dir) / "result.json"
            input_path.write_text(
                json.dumps(request, ensure_ascii=False),
                encoding="utf-8",
            )

            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(CLI),
                    "normalize-records",
                    "--input",
                    str(input_path),
                    "--output",
                    str(output_path),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["operation"], "normalize-records")
            self.assertEqual(
                payload["result"]["clutches"][0]["stage_counts"]["eggs_laid"],
                {"status": "observed", "coverage": "complete", "count": 0},
            )

    def test_simulation_recovery_cli_writes_seeded_result(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "simulation.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(CLI),
                    "simulate-recovery",
                    "--input",
                    str(
                        ROOT
                        / "examples"
                        / "input"
                        / "simulate_recovery.json"
                    ),
                    "--output",
                    str(output),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["operation"], "simulate-recovery")
            self.assertEqual(payload["result"]["base_seed"], 20260730)
            self.assertTrue(
                payload["result"]["limits"]["decision_support_only"]
            )

    def test_cli_rejects_nonstandard_nan_json_literal(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "non-finite.json"
            input_path.write_text(
                (
                    '{"model":{"model_id":"M1_AR",'
                    '"mode":"autosomal_recessive","penetrance":NaN},'
                    '"male_genotype":"N/N","female_genotype":"N/N"}'
                ),
                encoding="utf-8",
            )

            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(CLI),
                    "cross",
                    "--input",
                    str(input_path),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 2)
            self.assertEqual(completed.stdout, "")
            error = json.loads(completed.stderr)
            self.assertEqual(error["error"], "input_json_non_finite:NaN")

    def test_cli_rejects_json_number_that_overflows_to_infinity(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "overflow.json"
            input_path.write_text(
                (
                    '{"model":{"model_id":"M1_AR",'
                    '"mode":"autosomal_recessive"},'
                    '"male_genotype":"N/N","female_genotype":"N/N",'
                    '"unused":1e309}'
                ),
                encoding="utf-8",
            )

            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(CLI),
                    "cross",
                    "--input",
                    str(input_path),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 2)
            self.assertEqual(completed.stdout, "")
            error = json.loads(completed.stderr)
            self.assertEqual(error["error"], "input_json_non_finite:number")

    def test_script_entrypoint_does_not_write_bytecode_into_skill_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            copied_root = Path(temp_dir) / "skill with spaces"
            (copied_root / "scripts").mkdir(parents=True)
            (copied_root / "examples" / "input").mkdir(parents=True)
            shutil.copy2(CLI, copied_root / "scripts" / CLI.name)
            shutil.copytree(
                ROOT / "src",
                copied_root / "src",
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
            )
            input_path = (
                copied_root
                / "examples"
                / "input"
                / "cross_autosomal_recessive.json"
            )
            shutil.copy2(
                ROOT / "examples" / "input" / input_path.name,
                input_path,
            )

            completed = subprocess.run(
                [
                    sys.executable,
                    str(copied_root / "scripts" / CLI.name),
                    "cross",
                    "--input",
                    str(input_path),
                    "--output",
                    str(copied_root / "result.json"),
                ],
                cwd=copied_root,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(list(copied_root.rglob("*.pyc")), [])
            self.assertEqual(list(copied_root.rglob("__pycache__")), [])


if __name__ == "__main__":
    unittest.main()

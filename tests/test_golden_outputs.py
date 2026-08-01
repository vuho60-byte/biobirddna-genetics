from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERIFY = ROOT / "scripts" / "verify_examples.py"


def run_verifier(*extra_args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-B", str(VERIFY), *extra_args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


class GoldenOutputTests(unittest.TestCase):
    def test_all_goldens_are_locked_and_probability_invariants_hold(self) -> None:
        completed = run_verifier()

        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(completed.stdout)
        self.assertEqual(summary["examples_checked"], 7)
        self.assertEqual(summary["verified_examples"], 7)
        self.assertGreater(summary["distributions_checked"], 0)
        self.assertTrue(summary["byte_stable"])
        self.assertEqual(summary["golden_mismatches"], [])
        self.assertEqual(summary["probability_errors"], [])
        self.assertIn("limits", summary)
        self.assertIn("evidence", summary)

    def test_verifier_detects_golden_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            golden_dir = Path(temp_dir) / "output"
            shutil.copytree(ROOT / "examples" / "output", golden_dir)
            drifted = golden_dir / "cross_autosomal_recessive.json"
            drifted.write_bytes(drifted.read_bytes() + b" ")

            completed = run_verifier("--golden-dir", str(golden_dir))

        self.assertEqual(completed.returncode, 1)
        summary = json.loads(completed.stdout)
        self.assertEqual(
            summary["golden_mismatches"],
            ["cross_autosomal_recessive.json"],
        )

    def test_probability_guard_rejects_distribution_drift(self) -> None:
        self.assertTrue(VERIFY.is_file(), "verify_examples.py must exist")
        spec = importlib.util.spec_from_file_location("verify_examples", VERIFY)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        errors = module.find_probability_errors(
            {
                "result": {
                    "phenotype_probabilities": {
                        "normal": 0.4,
                        "target_colour": 0.5,
                    }
                }
            }
        )

        self.assertEqual(len(errors), 1)
        self.assertIn("phenotype_probabilities", errors[0])

    def test_probability_guard_rejects_empty_distribution(self) -> None:
        self.assertTrue(VERIFY.is_file(), "verify_examples.py must exist")
        spec = importlib.util.spec_from_file_location("verify_examples", VERIFY)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        errors = module.find_probability_errors(
            {"result": {"phenotype_probabilities": {}}}
        )

        self.assertEqual(len(errors), 1)
        self.assertIn("must not be empty", errors[0])

    def test_probability_guard_checks_simulation_model_posteriors(self) -> None:
        self.assertTrue(VERIFY.is_file(), "verify_examples.py must exist")
        spec = importlib.util.spec_from_file_location("verify_examples", VERIFY)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        errors = module.find_probability_errors(
            {"model_posteriors": {"M1_AR": 0.6, "M0_NONGENETIC": 0.3}}
        )

        self.assertEqual(len(errors), 1)
        self.assertIn("model_posteriors", errors[0])

    def test_json_number_guard_rejects_non_finite_scalar_metrics(self) -> None:
        self.assertTrue(VERIFY.is_file(), "verify_examples.py must exist")
        spec = importlib.util.spec_from_file_location("verify_examples", VERIFY)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        errors = module.find_non_finite_errors(
            {
                "uncertainty": {
                    "entropy_bits": float("nan"),
                    "multiclass_brier": float("inf"),
                }
            }
        )

        self.assertEqual(len(errors), 2)
        self.assertTrue(
            all("finite JSON number" in error for error in errors)
        )


if __name__ == "__main__":
    unittest.main()

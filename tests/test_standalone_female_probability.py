"""female_probability contract for the standalone reference implementations.

references/data-contracts.md section 3 promises fail-closed rejection of a
Z-linked model with female_probability != 0.5 in every runtime, with the
shared machine-readable code ``z_linked_female_probability_must_be_0_5``.
These tests pin that promise onto standalone/biobirddna_genetics_standalone.py
and standalone/biobirddna_genetics_node.mjs, and prove the autosomal path
honors a non-0.5 value identically in both files (no silent bias).
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from standalone.biobirddna_genetics_standalone import (  # noqa: E402
    Model,
    cross_known,
)

NODE = shutil.which("node")
NODE_ENTRYPOINT = ROOT / "standalone" / "biobirddna_genetics_node.mjs"


def _run_node_cross(request: dict) -> subprocess.CompletedProcess[str]:
    with tempfile.NamedTemporaryFile(
        "w",
        suffix=".json",
        delete=False,
        encoding="utf-8",
    ) as handle:
        json.dump(request, handle)
        input_path = handle.name
    try:
        return subprocess.run(
            [NODE, str(NODE_ENTRYPOINT), "cross", input_path],
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            capture_output=True,
            check=False,
        )
    finally:
        Path(input_path).unlink(missing_ok=True)


class PythonStandaloneFemaleProbabilityTests(unittest.TestCase):
    def test_z_linked_non_half_female_probability_rejected(self) -> None:
        for mode in ("z_linked_recessive", "z_linked_dominant"):
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError) as context:
                    Model.from_dict({"mode": mode, "female_probability": 0.9})
                self.assertEqual(
                    str(context.exception),
                    "z_linked_female_probability_must_be_0_5",
                )

    def test_z_linked_accepts_explicit_and_default_half(self) -> None:
        explicit = Model.from_dict(
            {"mode": "z_linked_recessive", "female_probability": 0.5}
        )
        defaulted = Model.from_dict({"mode": "z_linked_recessive"})
        self.assertEqual(explicit.female_probability, 0.5)
        self.assertEqual(defaulted.female_probability, 0.5)

    def test_out_of_range_female_probability_rejected(self) -> None:
        for bad in (-0.1, 1.5, float("nan")):
            with self.subTest(value=bad):
                with self.assertRaises(ValueError) as context:
                    Model.from_dict(
                        {
                            "mode": "autosomal_recessive",
                            "female_probability": bad,
                        }
                    )
                self.assertEqual(
                    str(context.exception),
                    "female_probability_must_be_between_0_and_1",
                )

    def test_autosomal_cross_honors_non_half_female_probability(self) -> None:
        model = Model.from_dict(
            {"mode": "autosomal_recessive", "female_probability": 0.7}
        )
        result = cross_known(model, "N/m", "N/m")
        phenotypes = result["phenotype_probabilities"]
        female_mass = sum(
            value
            for key, value in phenotypes.items()
            if key.startswith("female:")
        )
        self.assertTrue(math.isclose(female_mass, 0.7, rel_tol=0.0, abs_tol=1e-12))
        self.assertTrue(
            math.isclose(
                phenotypes["female:target_colour"],
                0.7 * 0.25,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        )

    def test_non_genetic_cross_honors_non_half_female_probability(self) -> None:
        model = Model.from_dict(
            {"mode": "non_genetic", "female_probability": 0.7}
        )
        result = cross_known(model, "NA", "NA")
        phenotypes = result["phenotype_probabilities"]
        female_mass = sum(
            value
            for key, value in phenotypes.items()
            if key.startswith("female:")
        )
        self.assertTrue(math.isclose(female_mass, 0.7, rel_tol=0.0, abs_tol=1e-12))


@unittest.skipUnless(NODE, "Node.js is not installed")
class NodeStandaloneFemaleProbabilityTests(unittest.TestCase):
    def test_z_linked_non_half_female_probability_rejected(self) -> None:
        completed = _run_node_cross(
            {
                "model": {
                    "mode": "z_linked_recessive",
                    "female_probability": 0.9,
                },
                "male_genotype": "N/m",
                "female_genotype": "m/W",
            }
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(
            "z_linked_female_probability_must_be_0_5",
            completed.stderr,
        )

    def test_z_linked_accepts_explicit_half(self) -> None:
        completed = _run_node_cross(
            {
                "model": {
                    "mode": "z_linked_recessive",
                    "female_probability": 0.5,
                },
                "male_genotype": "N/m",
                "female_genotype": "m/W",
            }
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_out_of_range_female_probability_rejected(self) -> None:
        completed = _run_node_cross(
            {
                "model": {
                    "mode": "autosomal_recessive",
                    "female_probability": 1.5,
                },
                "male_genotype": "N/m",
                "female_genotype": "N/m",
            }
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(
            "female_probability_must_be_between_0_and_1",
            completed.stderr,
        )

    def test_autosomal_non_half_matches_python_standalone(self) -> None:
        model_dict = {
            "mode": "autosomal_recessive",
            "female_probability": 0.7,
        }
        completed = _run_node_cross(
            {
                "model": model_dict,
                "male_genotype": "N/m",
                "female_genotype": "N/m",
            }
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        node_phenotypes = json.loads(completed.stdout)["phenotype_probabilities"]
        python_phenotypes = cross_known(
            Model.from_dict(model_dict), "N/m", "N/m"
        )["phenotype_probabilities"]
        self.assertEqual(set(node_phenotypes), set(python_phenotypes))
        for key, value in python_phenotypes.items():
            self.assertTrue(
                math.isclose(
                    node_phenotypes[key],
                    value,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ),
                key,
            )


if __name__ == "__main__":
    unittest.main()

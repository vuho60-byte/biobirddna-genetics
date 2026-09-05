from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.enums import Sex
from biobirddna_genetics.errors import ValidationError
from biobirddna_genetics.inference import (
    CONFIGURATION_PRIOR_WARNING,
    default_candidate_models,
    infer_candidate_models,
)
from biobirddna_genetics.models import OffspringObservation


class InferenceTests(unittest.TestCase):
    def test_candidate_set_requires_m0_and_multiple_models(self) -> None:
        observations = [
            OffspringObservation("normal", 1, Sex.UNKNOWN),
        ]
        models = default_candidate_models(target_label="target")

        with self.assertRaisesRegex(
            ValidationError,
            "^candidate_models_require_exactly_one_M0_NONGENETIC$",
        ):
            infer_candidate_models(
                [models[0]],
                male_parent_phenotype="normal",
                female_parent_phenotype="normal",
                observations=observations,
            )

        with self.assertRaisesRegex(
            ValidationError,
            "^candidate_models_require_multiple_models$",
        ):
            infer_candidate_models(
                [models[-1]],
                male_parent_phenotype="unknown",
                female_parent_phenotype="unknown",
                observations=observations,
            )

    def test_duplicate_candidate_model_ids_fail_closed(self) -> None:
        models = default_candidate_models(target_label="target")
        duplicate = type(models[0])(
            **{
                **models[1].to_dict(),
                "model_id": models[0].model_id,
                "mode": models[1].mode,
            }
        )

        with self.assertRaisesRegex(
            ValidationError,
            "^candidate_model_id_duplicate:M1_AR$",
        ):
            infer_candidate_models(
                [models[0], duplicate, models[-1]],
                male_parent_phenotype="normal",
                female_parent_phenotype="normal",
                observations=[
                    OffspringObservation("normal", 1, Sex.UNKNOWN),
                ],
            )

    def test_equal_posterior_order_is_deterministic_by_model_id(self) -> None:
        models = default_candidate_models(target_label="target")
        model_type = type(models[0])
        shared = {
            **models[0].to_dict(),
            "mode": models[0].mode,
        }
        z_model = model_type(**{**shared, "model_id": "Z_AR"})
        a_model = model_type(**{**shared, "model_id": "A_AR"})
        observations = [
            OffspringObservation("target", 2, Sex.UNKNOWN),
        ]

        first = infer_candidate_models(
            [z_model, a_model, models[-1]],
            male_parent_phenotype="unknown",
            female_parent_phenotype="unknown",
            observations=observations,
        )
        second = infer_candidate_models(
            [a_model, z_model, models[-1]],
            male_parent_phenotype="unknown",
            female_parent_phenotype="unknown",
            observations=observations,
        )

        self.assertEqual(first, second)
        self.assertEqual(first["best_model"], "A_AR")
        self.assertEqual(
            [item["model_id"] for item in first["ranking"][:2]],
            ["A_AR", "Z_AR"],
        )

    def test_sex_skewed_data_prefers_z_linked_recessive(self) -> None:
        observations = [
            OffspringObservation("normal", 5, Sex.MALE),
            OffspringObservation("target", 0, Sex.MALE),
            OffspringObservation("normal", 0, Sex.FEMALE),
            OffspringObservation("target", 6, Sex.FEMALE),
        ]
        result = infer_candidate_models(
            default_candidate_models(target_label="target"),
            male_parent_phenotype="target",
            female_parent_phenotype="normal",
            observations=observations,
        )
        self.assertEqual(result["best_model"], "M4_ZR")
        self.assertGreater(result["best_model_posterior"], 0.8)

    def test_collapsed_sex_observations_are_supported(self) -> None:
        observations = [
            OffspringObservation("normal", 8, Sex.UNKNOWN),
            OffspringObservation("target", 3, Sex.UNKNOWN),
        ]
        result = infer_candidate_models(
            default_candidate_models(target_label="target"),
            male_parent_phenotype="normal",
            female_parent_phenotype="normal",
            observations=observations,
        )
        self.assertTrue(result["ranking"])

    def test_z_linked_candidate_set_rejects_non_default_female_probability(
        self,
    ) -> None:
        # A non-0.5 female_probability on any candidate would bias model
        # comparison against Z-linked candidates (which always predict a 0.5
        # sex ratio), so the whole set is rejected fail-closed.
        models = default_candidate_models(target_label="target")
        model_type = type(models[0])
        biased_autosomal = model_type(
            **{
                **models[0].to_dict(),
                "mode": models[0].mode,
                "female_probability": 0.7,
            }
        )

        with self.assertRaisesRegex(
            ValidationError,
            "^female_probability_must_be_0_5_when_z_linked_model_present$",
        ):
            infer_candidate_models(
                [biased_autosomal, *models[1:]],
                male_parent_phenotype="normal",
                female_parent_phenotype="normal",
                observations=[
                    OffspringObservation("normal", 1, Sex.UNKNOWN),
                ],
            )

    def test_autosomal_only_candidate_set_accepts_custom_female_probability(
        self,
    ) -> None:
        models = default_candidate_models(target_label="target")
        model_type = type(models[0])
        biased_autosomal = model_type(
            **{
                **models[0].to_dict(),
                "mode": models[0].mode,
                "female_probability": 0.7,
            }
        )

        result = infer_candidate_models(
            [biased_autosomal, models[-1]],
            male_parent_phenotype="normal",
            female_parent_phenotype="normal",
            observations=[
                OffspringObservation("normal", 5, Sex.MALE),
                OffspringObservation("target", 2, Sex.FEMALE),
            ],
        )
        self.assertTrue(result["ranking"])

    def test_default_female_probability_candidate_set_accepted(self) -> None:
        # All default models carry female_probability 0.5, including the
        # Z-linked candidates: the shared default must keep working.
        result = infer_candidate_models(
            default_candidate_models(target_label="target"),
            male_parent_phenotype="normal",
            female_parent_phenotype="normal",
            observations=[
                OffspringObservation("normal", 3, Sex.UNKNOWN),
                OffspringObservation("target", 1, Sex.UNKNOWN),
            ],
        )
        self.assertTrue(result["ranking"])

    def test_model_priors_reject_non_finite_and_negative_values(self) -> None:
        observations = [
            OffspringObservation("normal", 1, Sex.UNKNOWN),
        ]
        cases = [
            (float("nan"), "model_prior_non_finite:M1_AR"),
            (float("inf"), "model_prior_non_finite:M1_AR"),
            (-0.1, "model_prior_negative:M1_AR"),
        ]

        for prior, error in cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(
                    ValidationError,
                    f"^{error}$",
                ):
                    infer_candidate_models(
                        default_candidate_models(target_label="target"),
                        male_parent_phenotype="normal",
                        female_parent_phenotype="normal",
                        observations=observations,
                        model_priors={"M1_AR": prior},
                    )

    def test_configuration_prior_ambiguity_emits_warning(self) -> None:
        # Sex-unknown 6 normal : 2 target from normal x normal parents:
        # M4_ZR and M1_AR share the same best-configuration likelihood, and
        # the posterior gap (0.527 vs 0.265) comes from M4_ZR having 2
        # compatible parent configurations while M1_AR has 4 under the
        # uniform per-configuration prior. The ranking must carry the
        # statistical caution.
        observations = [
            OffspringObservation("normal", 6, Sex.UNKNOWN),
            OffspringObservation("target", 2, Sex.UNKNOWN),
        ]
        result = infer_candidate_models(
            default_candidate_models(target_label="target"),
            male_parent_phenotype="normal",
            female_parent_phenotype="normal",
            observations=observations,
        )
        self.assertEqual(result["best_model"], "M4_ZR")
        self.assertEqual(result["warnings"], [CONFIGURATION_PRIOR_WARNING])

    def test_decisive_evidence_does_not_emit_configuration_prior_warning(
        self,
    ) -> None:
        # Sexed offspring (criss-cross pattern, the infer golden example):
        # log-evidence gap between M4_ZR and every runner-up far exceeds the
        # log configuration-count ratio, so no warning may appear.
        observations = [
            OffspringObservation("normal", 5, Sex.MALE),
            OffspringObservation("target", 0, Sex.MALE),
            OffspringObservation("normal", 0, Sex.FEMALE),
            OffspringObservation("target", 6, Sex.FEMALE),
        ]
        result = infer_candidate_models(
            default_candidate_models(target_label="target"),
            male_parent_phenotype="target",
            female_parent_phenotype="normal",
            observations=observations,
        )
        self.assertEqual(result["best_model"], "M4_ZR")
        self.assertGreater(result["best_model_posterior"], 0.99)
        self.assertEqual(result["warnings"], [])

    def test_configuration_prior_warning_visible_in_cli_stdout_envelope(
        self,
    ) -> None:
        request = {
            "target_label": "target",
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "offspring": [
                {"sex": "unknown", "phenotype": "normal", "count": 6},
                {"sex": "unknown", "phenotype": "target", "count": 2},
            ],
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "request.json"
            input_path.write_text(
                json.dumps(request, ensure_ascii=False),
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(ROOT / "scripts" / "biobirddna_cli.py"),
                    "infer",
                    "--input",
                    str(input_path),
                ],
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn(CONFIGURATION_PRIOR_WARNING, completed.stdout)
        envelope = json.loads(completed.stdout)
        self.assertEqual(
            envelope["result"]["warnings"],
            [CONFIGURATION_PRIOR_WARNING],
        )

    def test_incompatible_parent_models_use_json_null_not_infinity(self) -> None:
        result = infer_candidate_models(
            default_candidate_models(target_label="target"),
            male_parent_phenotype="unsupported",
            female_parent_phenotype="unsupported",
            observations=[
                OffspringObservation("normal", 1, Sex.UNKNOWN),
            ],
        )

        incompatible = [
            score
            for score in result["ranking"]
            if score.get("reason") == "parent_phenotype_incompatible_with_model"
        ]
        self.assertTrue(incompatible)
        self.assertTrue(
            all(score["log_evidence"] is None for score in incompatible)
        )
        json.dumps(result, allow_nan=False)


if __name__ == "__main__":
    unittest.main()

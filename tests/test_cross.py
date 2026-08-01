from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.cross import cross_known, cross_uncertain
from biobirddna_genetics.enums import InheritanceMode
from biobirddna_genetics.errors import ValidationError
from biobirddna_genetics.genotype import normalize_distribution
from biobirddna_genetics.models import InheritanceModel
from biobirddna_genetics.service import dispatch


class CrossTests(unittest.TestCase):
    def test_autosomal_recessive_heterozygote_by_affected(self) -> None:
        model = InheritanceModel(
            model_id="AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label="target",
        )
        result = cross_known(model, "N/m", "m/m")
        self.assertAlmostEqual(result["phenotype_probabilities"]["male:target"], 0.25)
        self.assertAlmostEqual(result["phenotype_probabilities"]["female:target"], 0.25)
        target_total = sum(
            value
            for key, value in result["phenotype_probabilities"].items()
            if key.endswith(":target")
        )
        self.assertAlmostEqual(target_total, 0.5)

    def test_z_linked_recessive_affected_male_by_normal_female(self) -> None:
        model = InheritanceModel(
            model_id="ZR",
            mode=InheritanceMode.Z_LINKED_RECESSIVE,
            target_label="target",
        )
        result = cross_known(model, "m/m", "N/W")
        self.assertAlmostEqual(result["phenotype_probabilities"]["female:target"], 0.5)
        self.assertAlmostEqual(result["phenotype_probabilities"]["male:normal"], 0.5)
        self.assertNotIn("male:target", result["phenotype_probabilities"])

    def test_canonical_duplicate_genotype_probabilities_are_combined(self) -> None:
        model = InheritanceModel(
            model_id="AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label="target",
        )
        result = cross_uncertain(
            model,
            {"N/m": 0.5, "m/N": 0.5},
            {"m/m": 1.0},
        )
        self.assertEqual(result["male_genotype_posterior"], {"N/m": 1.0})

    def test_uncertain_parent_genotype_marginalization(self) -> None:
        model = InheritanceModel(
            model_id="AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label="target",
        )
        result = cross_uncertain(
            model,
            {"N/N": 0.5, "N/m": 0.5},
            {"m/m": 1.0},
        )
        target_total = sum(
            value
            for key, value in result["phenotype_probabilities"].items()
            if key.endswith(":target")
        )
        self.assertAlmostEqual(target_total, 0.25)

    def test_uncertain_distribution_accepts_relative_weights_above_one(
        self,
    ) -> None:
        model = InheritanceModel(
            model_id="AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label="target",
        )
        result = cross_uncertain(
            model,
            {"N/N": 2.0, "N/m": 1.0},
            {"m/m": 3.0},
        )

        self.assertAlmostEqual(
            result["male_genotype_posterior"]["N/N"],
            2.0 / 3.0,
        )
        self.assertEqual(result["female_genotype_posterior"], {"m/m": 1.0})

    def test_probability_distribution_rejects_non_finite_values(self) -> None:
        for non_finite in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(non_finite=non_finite):
                with self.assertRaisesRegex(
                    ValidationError,
                    "probability_distribution_non_finite",
                ):
                    normalize_distribution({"N/N": non_finite, "N/m": 1.0})

    def test_probability_distribution_rejects_non_numeric_values(self) -> None:
        for invalid in (None, True, "1"):
            with self.subTest(invalid=invalid):
                with self.assertRaisesRegex(
                    ValidationError,
                    "^probability_distribution_invalid$",
                ):
                    normalize_distribution({"N/N": invalid, "N/m": 1.0})

    def test_model_rejects_colliding_phenotype_labels(self) -> None:
        with self.assertRaisesRegex(
            ValidationError,
            "^phenotype_labels_must_be_distinct$",
        ):
            InheritanceModel(
                model_id="AR",
                mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
                normal_label="same",
                target_label="same",
            )

        with self.assertRaisesRegex(
            ValidationError,
            "^phenotype_labels_must_be_distinct$",
        ):
            InheritanceModel(
                model_id="INC",
                mode=InheritanceMode.INCOMPLETE_DOMINANT,
                normal_label="normal",
                target_label="target",
                intermediate_label="target",
            )

    def test_non_genetic_cross_requires_exact_na_genotypes(self) -> None:
        model = InheritanceModel(
            model_id="M0_NONGENETIC",
            mode=InheritanceMode.NON_GENETIC,
            target_label="target",
        )

        with self.assertRaisesRegex(
            ValidationError,
            "^non_genetic_model_requires_NA_genotype$",
        ):
            cross_known(model, "garbage", "NA")

        with self.assertRaisesRegex(
            ValidationError,
            "^non_genetic_model_requires_NA_genotype$",
        ):
            cross_uncertain(
                model,
                {"garbage": 1.0},
                {"NA": 1.0},
            )

    def test_cross_request_rejects_mixed_known_and_uncertain_modes(self) -> None:
        with self.assertRaisesRegex(
            ValidationError,
            "^cross_parent_input_modes_mutually_exclusive$",
        ):
            dispatch(
                "cross",
                {
                    "model": {
                        "model_id": "M1_AR",
                        "mode": "autosomal_recessive",
                    },
                    "male_genotype": "N/m",
                    "female_genotype": "m/m",
                    "male_genotypes": {"N/m": 1.0},
                    "female_genotypes": {"m/m": 1.0},
                },
            )

    def test_dispatch_rejects_non_finite_numbers_before_routing(self) -> None:
        with self.assertRaisesRegex(
            ValidationError,
            "^input_json_non_finite:number$",
        ):
            dispatch(
                "cross",
                {
                    "model": {
                        "model_id": "M1_AR",
                        "mode": "autosomal_recessive",
                    },
                    "male_genotype": "N/m",
                    "female_genotype": "m/m",
                    "unused": float("inf"),
                },
            )

    def test_cross_request_rejects_unknown_fields(self) -> None:
        with self.assertRaisesRegex(
            ValidationError,
            "^cross_request_unknown_field:unknown$",
        ):
            dispatch(
                "cross",
                {
                    "model": {
                        "model_id": "M1_AR",
                        "mode": "autosomal_recessive",
                    },
                    "male_genotype": "N/m",
                    "female_genotype": "m/m",
                    "unknown": True,
                },
            )


if __name__ == "__main__":
    unittest.main()

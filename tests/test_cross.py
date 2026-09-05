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

    def test_z_linked_model_rejects_non_default_female_probability(self) -> None:
        for mode in (
            InheritanceMode.Z_LINKED_RECESSIVE,
            InheritanceMode.Z_LINKED_DOMINANT,
        ):
            with self.subTest(mode=mode):
                with self.assertRaisesRegex(
                    ValidationError,
                    "^z_linked_female_probability_must_be_0_5$",
                ):
                    InheritanceModel(
                        model_id="Z",
                        mode=mode,
                        target_label="target",
                        female_probability=0.6,
                    )

        explicit_default = InheritanceModel(
            model_id="ZR",
            mode=InheritanceMode.Z_LINKED_RECESSIVE,
            target_label="target",
            female_probability=0.5,
        )
        result = cross_known(explicit_default, "m/m", "N/W")
        self.assertAlmostEqual(
            result["phenotype_probabilities"]["female:target"],
            0.5,
        )

    def test_dispatch_cross_rejects_z_linked_female_probability(self) -> None:
        with self.assertRaisesRegex(
            ValidationError,
            "^z_linked_female_probability_must_be_0_5$",
        ):
            dispatch(
                "cross",
                {
                    "model": {
                        "model_id": "M4_ZR",
                        "mode": "z_linked_recessive",
                        "target_label": "target",
                        "female_probability": 0.6,
                    },
                    "male_genotype": "N/m",
                    "female_genotype": "N/W",
                },
            )

    def test_female_probability_is_honored_or_rejected_never_ignored(self) -> None:
        # Contract guard: for every mode, a non-default female_probability is
        # either reflected exactly in the female phenotype marginal or
        # rejected fail-closed. Reverting the Z-linked rejection makes this
        # fail because Z-linked crosses would silently emit a 0.5 marginal.
        p_female = 0.7
        genotypes = {
            InheritanceMode.AUTOSOMAL_RECESSIVE: ("N/m", "N/m"),
            InheritanceMode.AUTOSOMAL_DOMINANT: ("N/m", "N/m"),
            InheritanceMode.INCOMPLETE_DOMINANT: ("N/m", "N/m"),
            InheritanceMode.Z_LINKED_RECESSIVE: ("N/m", "N/W"),
            InheritanceMode.Z_LINKED_DOMINANT: ("N/m", "N/W"),
            InheritanceMode.NON_GENETIC: ("NA", "NA"),
        }
        for mode, (male, female) in genotypes.items():
            with self.subTest(mode=mode):
                try:
                    model = InheritanceModel(
                        model_id="M",
                        mode=mode,
                        target_label="target",
                        female_probability=p_female,
                    )
                except ValidationError as exc:
                    self.assertEqual(
                        str(exc),
                        "z_linked_female_probability_must_be_0_5",
                    )
                    continue
                result = cross_known(model, male, female)
                female_marginal = sum(
                    value
                    for key, value in result[
                        "phenotype_probabilities"
                    ].items()
                    if key.startswith("female:")
                )
                self.assertAlmostEqual(female_marginal, p_female)

    def test_female_probability_range_is_validated(self) -> None:
        for out_of_range in (-0.1, 1.5):
            with self.subTest(out_of_range=out_of_range):
                with self.assertRaisesRegex(
                    ValidationError,
                    "^female_probability_invalid$",
                ):
                    InheritanceModel.from_dict(
                        {
                            "model_id": "M1_AR",
                            "mode": "autosomal_recessive",
                            "female_probability": out_of_range,
                        }
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

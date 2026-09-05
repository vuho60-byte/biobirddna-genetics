from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.enums import InheritanceMode, Sex
from biobirddna_genetics.errors import ValidationError
from biobirddna_genetics.models import (
    GeneticProfile,
    InheritanceModel,
    ModelHypothesis,
    PairingWeights,
    PedigreeIndividual,
)
from biobirddna_genetics.optimizer import rank_mating_pairs
from biobirddna_genetics.pedigree import Pedigree


class OptimizerTests(unittest.TestCase):
    def test_same_individual_pair_is_always_blocked(self) -> None:
        pedigree = Pedigree([PedigreeIndividual("SELF", Sex.MALE)])
        model = InheritanceModel(
            model_id="AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label="target",
        )

        result = rank_mating_pairs(
            [
                GeneticProfile(
                    "SELF",
                    Sex.MALE,
                    genotype_posteriors={"AR": {"N/m": 1.0}},
                )
            ],
            [
                GeneticProfile(
                    "SELF",
                    Sex.FEMALE,
                    genotype_posteriors={"AR": {"N/m": 1.0}},
                )
            ],
            [ModelHypothesis(model, 1.0)],
            pedigree,
            max_offspring_inbreeding=1.0,
            block_direct_ancestor_and_full_sibling=False,
        )

        pair = result["pair_ranking"][0]
        self.assertFalse(pair["eligible"])
        self.assertIn("same_individual_pair", pair["block_reasons"])

    def test_full_sibling_pair_is_blocked_and_unrelated_pair_is_eligible(self) -> None:
        pedigree = Pedigree(
            [
                PedigreeIndividual("S", Sex.MALE),
                PedigreeIndividual("D", Sex.FEMALE),
                PedigreeIndividual("M", Sex.MALE, "S", "D"),
                PedigreeIndividual("F1", Sex.FEMALE, "S", "D"),
                PedigreeIndividual("F2", Sex.FEMALE),
            ]
        )
        model = InheritanceModel(
            model_id="AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label="target",
        )
        male = GeneticProfile(
            "M",
            Sex.MALE,
            genotype_posteriors={"AR": {"N/m": 1.0}},
        )
        females = [
            GeneticProfile(
                "F1",
                Sex.FEMALE,
                genotype_posteriors={"AR": {"m/m": 1.0}},
            ),
            GeneticProfile(
                "F2",
                Sex.FEMALE,
                genotype_posteriors={"AR": {"N/N": 1.0}},
            ),
        ]
        result = rank_mating_pairs(
            [male],
            females,
            [ModelHypothesis(model, 1.0)],
            pedigree,
            weights=PairingWeights(
                target=100.0,
                information=0.0,
                inbreeding=0.0,
                health=0.0,
                uncertainty=0.0,
            ),
        )
        by_female = {item["female_id"]: item for item in result["pair_ranking"]}
        self.assertFalse(by_female["F1"]["eligible"])
        self.assertIn("full_sibling_pair", by_female["F1"]["block_reasons"])
        self.assertTrue(by_female["F2"]["eligible"])
        self.assertGreater(by_female["F1"]["score"], by_female["F2"]["score"])
        self.assertEqual(result["pair_ranking"][0]["female_id"], "F2")

    def test_relative_hypothesis_weights_above_one_are_normalized(self) -> None:
        pedigree = Pedigree(
            [
                PedigreeIndividual("M", Sex.MALE),
                PedigreeIndividual("F", Sex.FEMALE),
            ]
        )
        model = InheritanceModel(
            model_id="AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label="target",
        )
        result = rank_mating_pairs(
            [
                GeneticProfile(
                    "M",
                    Sex.MALE,
                    genotype_posteriors={"AR": {"N/m": 1.0}},
                )
            ],
            [
                GeneticProfile(
                    "F",
                    Sex.FEMALE,
                    genotype_posteriors={"AR": {"m/m": 1.0}},
                )
            ],
            [ModelHypothesis(model, 10.0)],
            pedigree,
        )

        self.assertEqual(len(result["pair_ranking"]), 1)
        self.assertTrue(result["pair_ranking"][0]["eligible"])

    def test_pairing_weights_and_scores_must_remain_finite(self) -> None:
        for kwargs, error in (
            ({"target": float("inf")}, "target_weight_invalid"),
            ({"health": -1.0}, "health_weight_invalid"),
        ):
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValidationError, f"^{error}$"):
                    PairingWeights(**kwargs)

        pedigree = Pedigree(
            [
                PedigreeIndividual("M", Sex.MALE),
                PedigreeIndividual("F", Sex.FEMALE),
            ]
        )
        ar = InheritanceModel(
            model_id="AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label="target",
        )
        ad = InheritanceModel(
            model_id="AD",
            mode=InheritanceMode.AUTOSOMAL_DOMINANT,
            target_label="target",
        )
        male = GeneticProfile(
            "M",
            Sex.MALE,
            genotype_posteriors={
                "AR": {"m/m": 1.0},
                "AD": {"N/N": 1.0},
            },
        )
        female = GeneticProfile(
            "F",
            Sex.FEMALE,
            genotype_posteriors={
                "AR": {"m/m": 1.0},
                "AD": {"N/N": 1.0},
            },
        )

        with self.assertRaisesRegex(
            ValidationError,
            "^pair_score_non_finite:M:F$",
        ):
            rank_mating_pairs(
                [male],
                [female],
                [ModelHypothesis(ar, 10.0), ModelHypothesis(ad, 10.0)],
                pedigree,
                weights=PairingWeights(
                    target=1.7e308,
                    information=1.7e308,
                    inbreeding=0.0,
                    health=0.0,
                    uncertainty=0.0,
                ),
            )


if __name__ == "__main__":
    unittest.main()

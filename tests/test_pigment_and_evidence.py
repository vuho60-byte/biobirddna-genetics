from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.enums import EvidenceGrade
from biobirddna_genetics.errors import ValidationError
from biobirddna_genetics.evidence import claim_allowed
from biobirddna_genetics.models import PigmentObservation
from biobirddna_genetics.pigment import rank_pigment_hypotheses
from standalone.biobirddna_genetics_standalone import pigment as standalone_pigment


class PigmentAndEvidenceTests(unittest.TestCase):
    def test_progressive_white_pattern_prefers_pgr(self) -> None:
        result = rank_pigment_hypotheses(
            PigmentObservation(
                white_area_ratio=0.2,
                progressive_across_molts=True,
                present_in_first_plumage=False,
                image_quality="Q3",
            )
        )
        self.assertEqual(result["primary_hypothesis"], "MEL-PGR")
        self.assertEqual(result["limits"]["gene_status"], "not_assessed")

    def test_alternative_hypotheses_exclude_zero_score_labels(self) -> None:
        # Golden-style case: only MEL-PGR and MEL-LEU accrue score; the other
        # slots in ranking[1:4] hold score-0.0 labels and must be dropped.
        result = rank_pigment_hypotheses(
            PigmentObservation(
                white_area_ratio=0.18,
                progressive_across_molts=True,
                present_in_first_plumage=False,
                image_quality="Q3",
            )
        )
        self.assertEqual(result["primary_hypothesis"], "MEL-PGR")
        self.assertEqual(result["alternative_hypotheses"], ["MEL-LEU"])
        for label in result["alternative_hypotheses"]:
            self.assertGreater(result["scores"][label], 0.0)

    def test_alternative_hypotheses_keep_three_when_all_positive(self) -> None:
        # Control: enough rules fire that ranking[1:4] is all positive-score
        # labels, so the three-alternative shape is unchanged.
        result = rank_pigment_hypotheses(
            PigmentObservation(
                brown_shift_score=0.6,
                grey_shift_score=0.5,
                dark_expansion_score=0.4,
                red_patch_change_score=0.3,
                image_quality="Q3",
            )
        )
        self.assertEqual(result["primary_hypothesis"], "MEL-BRN")
        self.assertEqual(
            result["alternative_hypotheses"],
            ["MEL-DIL", "MEL-MEL", "CAR-ALT"],
        )
        for label in result["alternative_hypotheses"]:
            self.assertGreater(result["scores"][label], 0.0)

    def test_q0_and_q1_force_recapture_despite_strong_visible_features(self) -> None:
        for image_quality in ("Q0", "Q1"):
            with self.subTest(image_quality=image_quality):
                result = rank_pigment_hypotheses(
                    PigmentObservation(
                        white_area_ratio=0.99,
                        brown_shift_score=1.0,
                        grey_shift_score=1.0,
                        dark_expansion_score=1.0,
                        red_patch_change_score=1.0,
                        pattern_retained=True,
                        iris_red_confirmed=True,
                        progressive_across_molts=True,
                        image_quality=image_quality,
                    )
                )

                self.assertEqual(result["primary_hypothesis"], "UNK")
                self.assertIs(result["needs_recapture"], True)
                self.assertEqual(
                    result["limits"]["measurement_domain"],
                    "visible_phenotype_only",
                )
                self.assertEqual(result["limits"]["gene_status"], "not_assessed")
                self.assertEqual(
                    result["limits"]["inheritance_status"],
                    "not_assessed",
                )

    def test_standalone_q0_and_q1_also_force_recapture(self) -> None:
        for image_quality in ("Q0", "Q1"):
            with self.subTest(image_quality=image_quality):
                result = standalone_pigment(
                    {
                        "white_area_ratio": 0.99,
                        "brown_shift_score": 1.0,
                        "grey_shift_score": 1.0,
                        "dark_expansion_score": 1.0,
                        "red_patch_change_score": 1.0,
                        "iris_red_confirmed": True,
                        "progressive_across_molts": True,
                        "image_quality": image_quality,
                    }
                )

                self.assertEqual(result["primary_hypothesis"], "UNK")
                self.assertIs(result["needs_recapture"], True)
                self.assertEqual(result["limits"]["gene_status"], "not_assessed")

    def test_gene_claim_requires_g5(self) -> None:
        self.assertFalse(claim_allowed("molecular_locus_supported", EvidenceGrade.G4))
        self.assertTrue(claim_allowed("molecular_locus_supported", EvidenceGrade.G5))

    def test_domain_object_rejects_invalid_boolean_and_quality_types(self) -> None:
        for kwargs, error in (
            ({"iris_red_confirmed": "false"}, "iris_red_confirmed_invalid"),
            ({"pattern_retained": "false"}, "pattern_retained_invalid"),
            ({"image_quality": "BOGUS"}, "image_quality_invalid"),
        ):
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValidationError, f"^{error}$"):
                    PigmentObservation(**kwargs)


if __name__ == "__main__":
    unittest.main()

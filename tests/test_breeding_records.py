from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.errors import ValidationError
from biobirddna_genetics.service import dispatch


def complete_record_request() -> dict:
    return {
        "schema_version": "1.0.0",
        "case_id": "CM-CASE-000001",
        "pairings": [
            {
                "pairing_id": "CM-MAT-000002",
                "male_parent_id": "CM-IND-B",
                "female_parent_id": "CM-IND-A",
                "male_parent_phenotype": "target_colour",
                "female_parent_phenotype": "normal",
                "male_parent_sex_method": "dna",
                "male_parent_sex_confidence": 0.99,
                "female_parent_sex_method": "dna",
                "female_parent_sex_confidence": 0.99,
                "reciprocal_group_id": "CM-RECIP-000001",
            },
            {
                "pairing_id": "CM-MAT-000001",
                "male_parent_id": "CM-IND-A",
                "female_parent_id": "CM-IND-B",
                "male_parent_phenotype": "normal",
                "female_parent_phenotype": "target_colour",
                "male_parent_sex_method": "dna",
                "male_parent_sex_confidence": 0.99,
                "female_parent_sex_method": "dna",
                "female_parent_sex_confidence": 0.99,
                "reciprocal_group_id": "CM-RECIP-000001",
            },
        ],
        "clutches": [
            {
                "clutch_id": "CM-CLU-000001",
                "pairing_id": "CM-MAT-000001",
                "offspring_records_coverage": "complete",
                "stage_counts": {
                    "eggs_laid": {
                        "status": "observed",
                        "coverage": "complete",
                        "count": 4,
                    },
                    "fertile_eggs": {
                        "status": "observed",
                        "coverage": "complete",
                        "count": 3,
                    },
                    "hatched": {
                        "status": "observed",
                        "coverage": "complete",
                        "count": 2,
                    },
                    "survived_to_first_molt": {
                        "status": "missing",
                        "coverage": "unknown",
                    },
                },
            },
            {
                "clutch_id": "CM-CLU-000002",
                "pairing_id": "CM-MAT-000002",
                "offspring_records_coverage": "complete",
                "stage_counts": {
                    "eggs_laid": {
                        "status": "observed",
                        "coverage": "complete",
                        "count": 0,
                    },
                    "fertile_eggs": {
                        "status": "observed",
                        "coverage": "complete",
                        "count": 0,
                    },
                    "hatched": {
                        "status": "observed",
                        "coverage": "complete",
                        "count": 0,
                    },
                    "survived_to_first_molt": {
                        "status": "observed",
                        "coverage": "complete",
                        "count": 0,
                    },
                },
            },
        ],
        "offspring": [
            {
                "offspring_record_id": "CM-OFF-000001",
                "individual_id": "CM-IND-OFF-000001",
                "clutch_id": "CM-CLU-000001",
                "outcome": "alive",
                "fertility_status": "fertile",
                "hatch_status": "hatched",
                "survival_status": "survived",
                "sex": "male",
                "sex_method": "dna",
                "sex_confidence": 0.99,
                "phenotype_timeline": [
                    {
                        "observation_id": "CM-PHENO-000001",
                        "sequence_index": 1,
                        "stage": "first_plumage",
                        "status": "observed",
                        "phenotype": "normal",
                    }
                ],
            },
            {
                "offspring_record_id": "CM-OFF-000002",
                "individual_id": "CM-IND-OFF-000002",
                "clutch_id": "CM-CLU-000001",
                "outcome": "post_hatch_death",
                "fertility_status": "fertile",
                "hatch_status": "hatched",
                "survival_status": "died_before_first_molt",
                "sex": "unknown",
                "sex_method": "unknown",
                "sex_confidence": None,
                "phenotype_timeline": [
                    {
                        "observation_id": "CM-PHENO-000002",
                        "sequence_index": 0,
                        "stage": "hatch",
                        "status": "missing",
                    }
                ],
            },
            {
                "offspring_record_id": "CM-OFF-000003",
                "individual_id": None,
                "clutch_id": "CM-CLU-000001",
                "outcome": "embryo_death",
                "fertility_status": "fertile",
                "hatch_status": "not_hatched",
                "survival_status": "not_applicable",
                "sex": "unknown",
                "sex_method": "unknown",
                "sex_confidence": None,
                "phenotype_timeline": [],
            },
            {
                "offspring_record_id": "CM-OFF-000004",
                "individual_id": None,
                "clutch_id": "CM-CLU-000001",
                "outcome": "unfertilized",
                "fertility_status": "unfertilized",
                "hatch_status": "not_applicable",
                "survival_status": "not_applicable",
                "sex": "unknown",
                "sex_method": "unknown",
                "sex_confidence": None,
                "phenotype_timeline": [],
            },
        ],
    }


class BreedingRecordTests(unittest.TestCase):
    def normalize(self, request: dict) -> dict:
        try:
            return dispatch("normalize-records", request)
        except ValidationError as exc:
            if str(exc) == "operation_unknown:normalize-records":
                self.fail("GEN-010 normalize-records operation is not implemented")
            self.fail(f"valid GEN-010 record was rejected: {exc}")

    def test_normalizes_losses_and_normal_phenotype_without_mutating_input(self) -> None:
        request = complete_record_request()
        original = copy.deepcopy(request)

        result = self.normalize(request)

        self.assertEqual(request, original)
        self.assertEqual(
            [item["outcome"] for item in result["offspring"]],
            ["alive", "post_hatch_death", "embryo_death", "unfertilized"],
        )
        self.assertEqual(
            result["offspring"][0]["phenotype_timeline"][0]["phenotype"],
            "normal",
        )
        self.assertEqual(result["evidence"]["genetics_grade"], "not_assessed")
        self.assertTrue(result["limits"]["decision_support_only"])
        self.assertEqual(result["limits"]["gene_status"], "not_assessed")

    def test_missing_stage_remains_distinct_from_observed_zero(self) -> None:
        result = self.normalize(complete_record_request())

        missing = result["clutches"][0]["stage_counts"]["survived_to_first_molt"]
        observed_zero = result["clutches"][1]["stage_counts"]["eggs_laid"]

        self.assertEqual(
            missing,
            {"status": "missing", "coverage": "unknown", "count": None},
        )
        self.assertEqual(
            observed_zero,
            {"status": "observed", "coverage": "complete", "count": 0},
        )
        self.assertIn(
            {
                "field": (
                    "clutches/CM-CLU-000001/stage_counts/"
                    "survived_to_first_molt"
                ),
                "impact": "survival_stage_not_observed",
            },
            result["missing_data"],
        )

    def test_partial_stage_coverage_is_reported_as_incomplete(self) -> None:
        request = complete_record_request()
        request["clutches"][0]["stage_counts"]["hatched"]["coverage"] = "partial"

        result = self.normalize(request)

        self.assertEqual(
            result["clutches"][0]["stage_counts"]["hatched"],
            {"status": "observed", "coverage": "partial", "count": 2},
        )
        self.assertIn(
            {
                "field": "clutches/CM-CLU-000001/stage_counts/hatched",
                "impact": "hatch_stage_incomplete_coverage",
            },
            result["missing_data"],
        )

    def test_preserves_unfertilized_and_hatch_failure_records(self) -> None:
        request = complete_record_request()
        request["clutches"][0]["stage_counts"]["eggs_laid"]["count"] = 5
        request["clutches"][0]["stage_counts"]["fertile_eggs"]["count"] = 4
        request["offspring"].append(
            {
                "offspring_record_id": "CM-OFF-000005",
                "individual_id": None,
                "clutch_id": "CM-CLU-000001",
                "outcome": "hatch_failure",
                "fertility_status": "fertile",
                "hatch_status": "not_hatched",
                "survival_status": "not_applicable",
                "sex": "unknown",
                "sex_method": "unknown",
                "sex_confidence": None,
                "phenotype_timeline": [],
            }
        )

        result = self.normalize(request)

        self.assertEqual(
            [item["outcome"] for item in result["offspring"][-2:]],
            ["unfertilized", "hatch_failure"],
        )

    def test_requires_parent_sex_provenance(self) -> None:
        request = complete_record_request()
        request["pairings"][0]["male_parent_sex_method"] = "unknown"
        request["pairings"][0]["male_parent_sex_confidence"] = None

        with self.assertRaisesRegex(
            ValidationError,
            "^parent_sex_provenance_required:CM-MAT-000002:male$",
        ):
            dispatch("normalize-records", request)

    def test_complete_stage_counts_must_match_offspring_records(self) -> None:
        request = complete_record_request()
        request["clutches"][0]["stage_counts"]["eggs_laid"]["count"] = 3

        with self.assertRaisesRegex(
            ValidationError,
            "^stage_offspring_count_conflict:"
            "CM-CLU-000001:eggs_laid:3:4$",
        ):
            dispatch("normalize-records", request)

    def test_requires_every_stage_key_explicitly(self) -> None:
        request = complete_record_request()
        del request["clutches"][0]["stage_counts"]["hatched"]

        with self.assertRaisesRegex(
            ValidationError,
            "^stage_count_required:CM-CLU-000001:hatched$",
        ):
            dispatch("normalize-records", request)

    def test_sorts_phenotype_timeline_by_sequence_index(self) -> None:
        request = complete_record_request()
        request["offspring"][0]["phenotype_timeline"].insert(
            0,
            {
                "observation_id": "CM-PHENO-000003",
                "sequence_index": 2,
                "stage": "adult",
                "status": "observed",
                "phenotype": "normal",
            },
        )
        request["offspring"][0]["phenotype_timeline"].append(
            {
                "observation_id": "CM-PHENO-000004",
                "sequence_index": 0,
                "stage": "hatch",
                "status": "missing",
            },
        )

        result = self.normalize(request)

        self.assertEqual(
            [
                item["sequence_index"]
                for item in result["offspring"][0]["phenotype_timeline"]
            ],
            [0, 1, 2],
        )

    def test_reciprocal_crosses_keep_parent_direction(self) -> None:
        result = self.normalize(complete_record_request())

        self.assertEqual(
            [item["pairing_id"] for item in result["pairings"]],
            ["CM-MAT-000001", "CM-MAT-000002"],
        )
        self.assertEqual(
            [item["cross_direction_key"] for item in result["pairings"]],
            ["CM-IND-A>CM-IND-B", "CM-IND-B>CM-IND-A"],
        )
        self.assertEqual(
            {item["reciprocal_group_id"] for item in result["pairings"]},
            {"CM-RECIP-000001"},
        )

    def test_rejects_duplicate_record_ids(self) -> None:
        request = complete_record_request()
        request["offspring"][1]["offspring_record_id"] = "CM-OFF-000001"

        with self.assertRaisesRegex(
            ValidationError,
            "^offspring_record_id_duplicate:CM-OFF-000001$",
        ):
            dispatch("normalize-records", request)

    def test_rejects_increasing_observed_stage_counts(self) -> None:
        request = complete_record_request()
        request["clutches"][0]["stage_counts"]["hatched"]["count"] = 4

        with self.assertRaisesRegex(
            ValidationError,
            "^clutch_stage_count_increases:CM-CLU-000001:"
            "fertile_eggs:3:hatched:4$",
        ):
            dispatch("normalize-records", request)

    def test_requires_method_and_confidence_for_known_sex(self) -> None:
        request = complete_record_request()
        request["offspring"][0]["sex_method"] = "unknown"
        request["offspring"][0]["sex_confidence"] = None

        with self.assertRaisesRegex(
            ValidationError,
            "^known_sex_requires_method_and_confidence:CM-OFF-000001$",
        ):
            dispatch("normalize-records", request)

    def test_rejects_empty_pairings_like_the_versioned_schema(self) -> None:
        request = complete_record_request()
        request["pairings"] = []

        with self.assertRaisesRegex(ValidationError, "^pairings_must_not_be_empty$"):
            dispatch("normalize-records", request)

    def test_rejects_unknown_fields_like_the_versioned_schema(self) -> None:
        request = complete_record_request()
        request["offspring"][0]["inferred_gene"] = "forbidden"

        with self.assertRaisesRegex(
            ValidationError,
            "^offspring_record_unknown_field:inferred_gene$",
        ):
            dispatch("normalize-records", request)

    def test_rejects_outcome_that_contradicts_stage_path(self) -> None:
        request = complete_record_request()
        request["offspring"][0].update(
            {
                "fertility_status": "unfertilized",
                "hatch_status": "not_applicable",
                "survival_status": "not_applicable",
            }
        )
        request["offspring"][3].update(
            {
                "outcome": "alive",
                "fertility_status": "fertile",
                "hatch_status": "hatched",
                "survival_status": "survived",
            }
        )

        with self.assertRaisesRegex(
            ValidationError,
            "^offspring_outcome_stage_conflict:CM-OFF-000001:alive$",
        ):
            dispatch("normalize-records", request)

    def test_rejects_absent_schema_required_offspring_fields(self) -> None:
        for field in ("sex_confidence", "phenotype_timeline"):
            with self.subTest(field=field):
                request = complete_record_request()
                del request["offspring"][1][field]

                with self.assertRaisesRegex(
                    ValidationError,
                    f"^offspring_field_required:CM-OFF-000002:{field}$",
                ):
                    dispatch("normalize-records", request)


if __name__ == "__main__":
    unittest.main()

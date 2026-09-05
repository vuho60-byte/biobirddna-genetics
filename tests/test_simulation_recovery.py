from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.errors import ValidationError
from biobirddna_genetics.service import dispatch
from biobirddna_genetics import simulation


ALL_CANDIDATES = [
    "M1_AR",
    "M2_AD",
    "M3_INC",
    "M4_ZR",
    "M5_ZD",
    "M0_NONGENETIC",
]


def ar_experiment(
    *,
    brood_sizes: list[int],
    replicates: int,
) -> dict:
    return {
        "experiment_id": "SIM-AR",
        "generating_model": {
            "model_id": "M1_AR",
            "mode": "autosomal_recessive",
            "target_label": "target",
        },
        "male_genotype": "N/m",
        "female_genotype": "N/m",
        "male_parent_phenotype": "normal",
        "female_parent_phenotype": "normal",
        "brood_sizes": brood_sizes,
        "replicates": replicates,
        "candidate_model_ids": list(ALL_CANDIDATES),
    }


def zr_experiment(
    *,
    brood_sizes: list[int],
    replicates: int,
) -> dict:
    return {
        "experiment_id": "SIM-ZR",
        "generating_model": {
            "model_id": "M4_ZR",
            "mode": "z_linked_recessive",
            "target_label": "target",
        },
        "male_genotype": "m/m",
        "female_genotype": "N/W",
        "male_parent_phenotype": "target",
        "female_parent_phenotype": "normal",
        "brood_sizes": brood_sizes,
        "replicates": replicates,
        "candidate_model_ids": list(ALL_CANDIDATES),
    }


def recovery_request(experiments: list[dict]) -> dict:
    return {
        "schema_version": "1.0.0",
        "base_seed": 20260730,
        "goodness_of_fit_replicates": 40,
        "misspecification_alpha": 0.05,
        "experiments": experiments,
    }


class SimulationRecoveryTests(unittest.TestCase):
    def simulate(self, request: dict) -> dict:
        try:
            return dispatch("simulate-recovery", request)
        except ValidationError as exc:
            if str(exc) == "operation_unknown:simulate-recovery":
                self.fail("GEN-011 simulate-recovery operation is not implemented")
            self.fail(f"valid GEN-011 request was rejected: {exc}")

    def test_same_seed_is_object_deterministic_and_does_not_mutate_input(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[8, 24], replicates=5)]
        )
        original = copy.deepcopy(request)

        first = self.simulate(request)
        second = self.simulate(request)

        self.assertEqual(request, original)
        self.assertEqual(first, second)
        self.assertEqual(first["base_seed"], 20260730)
        self.assertTrue(
            all(
                replicate["seed"] <= 9007199254740991
                for summary in first["experiments"][0]["by_brood_size"]
                for replicate in summary["replicate_results"]
            )
        )
        self.assertEqual(first["limits"]["gene_status"], "not_assessed")
        self.assertTrue(first["limits"]["synthetic_only"])
        self.assertTrue(first["limits"]["decision_support_only"])
        self.assertEqual(
            first["limits"]["posterior_scope"],
            "relative_to_candidate_model_set",
        )

    def test_known_ar_and_zr_models_recover_at_large_brood_size(self) -> None:
        request = recovery_request(
            [
                ar_experiment(brood_sizes=[96], replicates=40),
                zr_experiment(brood_sizes=[96], replicates=40),
            ]
        )
        request["goodness_of_fit_replicates"] = 100

        result = self.simulate(request)
        summaries = {
            item["experiment_id"]: item["by_brood_size"][0]
            for item in result["experiments"]
        }

        self.assertGreaterEqual(
            summaries["SIM-AR"]["recovery"]["top1_rate"],
            0.75,
        )
        self.assertGreaterEqual(
            summaries["SIM-ZR"]["recovery"]["top1_rate"],
            0.75,
        )
        self.assertNotIn(
            "candidate_model_set_misspecified_suspected",
            summaries["SIM-AR"]["warnings"],
        )
        self.assertNotIn(
            "candidate_model_set_misspecified_suspected",
            summaries["SIM-ZR"]["warnings"],
        )

    def test_small_sample_uncertainty_is_explicit(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[4], replicates=24)]
        )

        summary = self.simulate(request)["experiments"][0]["by_brood_size"][0]
        interval = summary["recovery"]["wilson_95"]

        self.assertTrue(summary["small_sample_warning"])
        self.assertIn("small_sample_high_uncertainty", summary["warnings"])
        self.assertLessEqual(interval["low"], interval["high"])
        self.assertGreaterEqual(interval["low"], 0.0)
        self.assertLessEqual(interval["high"], 1.0)
        self.assertIn(
            "true_model_posterior_quantiles",
            summary["uncertainty"],
        )
        self.assertTrue(summary["posterior_discrimination_bins"])
        self.assertNotIn(
            "observed_true_model_rate",
            summary["posterior_discrimination_bins"][0],
        )

    def test_misspecified_candidate_set_uses_absolute_fit_warning(self) -> None:
        experiment = {
            "experiment_id": "SIM-MISSPEC",
            "generating_model": {
                "model_id": "M3_INC",
                "mode": "incomplete_dominant",
                "target_label": "target",
                "intermediate_label": "intermediate",
            },
            "male_genotype": "N/N",
            "female_genotype": "m/m",
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "target",
            "brood_sizes": [64],
            "replicates": 24,
            "candidate_model_ids": [
                "M1_AR",
                "M2_AD",
                "M4_ZR",
                "M5_ZD",
                "M0_NONGENETIC",
            ],
        }
        request = recovery_request([experiment])

        summary = self.simulate(request)["experiments"][0]["by_brood_size"][0]

        self.assertIn(
            "candidate_model_set_misspecified_suspected",
            summary["warnings"],
        )
        self.assertTrue(
            summary["goodness_of_fit"]["all_candidates_below_alpha"]
        )
        self.assertTrue(
            any(
                candidate["status"] == "evaluated"
                and candidate["model_id"] != "M0_NONGENETIC"
                for candidate in summary["goodness_of_fit"]["candidates"]
            )
        )
        self.assertIsNone(summary["recovery"]["top1_rate"])

    def test_bootstrap_reselects_parent_configuration_each_time(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[32], replicates=3)]
        )
        request["goodness_of_fit_replicates"] = 20

        with patch.object(
            simulation,
            "score_model",
            wraps=simulation.score_model,
        ) as tracked_score:
            self.simulate(request)

        self.assertGreater(tracked_score.call_count, len(ALL_CANDIDATES))

    def test_calibration_uses_realized_true_model_across_design(self) -> None:
        request = recovery_request(
            [
                ar_experiment(brood_sizes=[32], replicates=16),
                zr_experiment(brood_sizes=[32], replicates=16),
            ]
        )

        calibration = self.simulate(request)["calibration_by_design"][0]

        self.assertEqual(calibration["status"], "estimated")
        self.assertEqual(calibration["brood_size"], 32)
        self.assertEqual(
            calibration["method"],
            "empirical_reliability_under_declared_simulation_design",
        )
        bins = next(
            item["bins"]
            for item in calibration["models"]
            if item["model_id"] == "M1_AR"
        )
        self.assertTrue(bins)
        self.assertIn("observed_true_model_rate", bins[0])
        self.assertNotIn("observed_top1_rate", bins[0])
        self.assertEqual(sum(item["count"] for item in bins), 32)
        self.assertAlmostEqual(
            sum(
                item["observed_true_model_rate"] * item["count"]
                for item in bins
            ),
            16.0,
        )

    def test_brier_score_matches_multiclass_definition(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[16], replicates=2)]
        )
        replicate = self.simulate(request)["experiments"][0][
            "by_brood_size"
        ][0]["replicate_results"][0]

        expected = sum(
            (
                probability
                - (1.0 if model_id == "M1_AR" else 0.0)
            )
            ** 2
            for model_id, probability in replicate[
                "model_posteriors"
            ].items()
        )

        self.assertAlmostEqual(
            replicate["multiclass_brier"],
            expected,
            places=15,
        )

    def test_wilson_interval_matches_fixed_numerical_oracle(self) -> None:
        interval = simulation._wilson_interval(5, 10)

        self.assertIsNotNone(interval)
        self.assertAlmostEqual(interval["low"], 0.236593090512564, places=15)
        self.assertAlmostEqual(interval["high"], 0.7634069094874361, places=15)

    def test_m0_non_genetic_candidate_is_mandatory(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[32], replicates=5)]
        )
        request["experiments"][0]["candidate_model_ids"].remove(
            "M0_NONGENETIC"
        )

        with self.assertRaisesRegex(
            ValidationError,
            "^simulation_candidate_set_requires_exactly_one_M0_NONGENETIC:"
            "SIM-AR$",
        ):
            dispatch("simulate-recovery", request)

    def test_candidate_set_requires_multiple_models(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[32], replicates=3)]
        )
        request["experiments"][0]["generating_model"] = {
            "model_id": "M0_NONGENETIC",
            "mode": "non_genetic",
            "target_label": "target",
        }
        request["experiments"][0]["male_genotype"] = "NA"
        request["experiments"][0]["female_genotype"] = "NA"
        request["experiments"][0]["candidate_model_ids"] = [
            "M0_NONGENETIC"
        ]

        with self.assertRaisesRegex(
            ValidationError,
            "^simulation_candidate_set_requires_multiple_models:SIM-AR$",
        ):
            dispatch("simulate-recovery", request)

    def test_bootstrap_resolution_must_support_selected_alpha(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[32], replicates=3)]
        )
        request["goodness_of_fit_replicates"] = 10

        with self.assertRaisesRegex(
            ValidationError,
            "^goodness_of_fit_replicates_insufficient_for_alpha:10:0.05$",
        ):
            dispatch("simulate-recovery", request)

    def test_generating_model_requires_schema_model_id(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[32], replicates=3)]
        )
        del request["experiments"][0]["generating_model"]["model_id"]

        with self.assertRaisesRegex(
            ValidationError,
            "^simulation_generating_model_id_required$",
        ):
            dispatch("simulate-recovery", request)

    def test_matching_model_id_requires_matching_candidate_snapshot(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[32], replicates=3)]
        )
        request["experiments"][0]["generating_model"]["penetrance"] = 0.0

        with self.assertRaisesRegex(
            ValidationError,
            "^simulation_generating_model_snapshot_mismatch:SIM-AR:M1_AR$",
        ):
            dispatch("simulate-recovery", request)

    def test_known_generating_id_keeps_canonical_semantics_when_omitted(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[32], replicates=3)]
        )
        experiment = request["experiments"][0]
        experiment["generating_model"]["mode"] = "autosomal_dominant"
        experiment["candidate_model_ids"].remove("M1_AR")

        with self.assertRaisesRegex(
            ValidationError,
            "^simulation_generating_model_snapshot_mismatch:SIM-AR:M1_AR$",
        ):
            dispatch("simulate-recovery", request)

    def test_semantic_text_fields_reject_surrounding_whitespace(self) -> None:
        cases = [
            (
                lambda experiment: experiment.__setitem__(
                    "experiment_id",
                    " SIM-AR ",
                ),
                "simulation_experiment_id_whitespace_forbidden",
            ),
            (
                lambda experiment: experiment["generating_model"].__setitem__(
                    "model_id",
                    " M1_AR ",
                ),
                "simulation_generating_model_id_whitespace_forbidden",
            ),
            (
                lambda experiment: experiment["generating_model"].__setitem__(
                    "mode",
                    " autosomal_recessive ",
                ),
                "simulation_generating_model_mode_whitespace_forbidden",
            ),
            (
                lambda experiment: experiment["generating_model"].__setitem__(
                    "normal_allele",
                    " N ",
                ),
                "simulation_generating_model_normal_allele_whitespace_forbidden",
            ),
            (
                lambda experiment: experiment.__setitem__(
                    "male_genotype",
                    " N/m ",
                ),
                "simulation_male_genotype_whitespace_forbidden",
            ),
            (
                lambda experiment: experiment.__setitem__(
                    "male_parent_phenotype",
                    " normal ",
                ),
                "simulation_male_parent_phenotype_whitespace_forbidden",
            ),
            (
                lambda experiment: experiment["candidate_model_ids"].__setitem__(
                    0,
                    " M1_AR ",
                ),
                "simulation_candidate_model_id_whitespace_forbidden",
            ),
            (
                lambda experiment: experiment["generating_model"].__setitem__(
                    "description",
                    " padded description ",
                ),
                "simulation_generating_model_description_whitespace_forbidden",
            ),
        ]

        for mutate, error in cases:
            with self.subTest(error=error):
                request = recovery_request(
                    [ar_experiment(brood_sizes=[32], replicates=3)]
                )
                mutate(request["experiments"][0])
                with self.assertRaisesRegex(
                    ValidationError,
                    f"^{error}$",
                ):
                    dispatch("simulate-recovery", request)

    def test_non_genetic_generating_model_requires_na_genotypes(self) -> None:
        experiment = {
            "experiment_id": "SIM-M0",
            "generating_model": {
                "model_id": "M0_NONGENETIC",
                "mode": "non_genetic",
                "target_label": "target",
            },
            "male_genotype": "NA",
            "female_genotype": "NA",
            "male_parent_phenotype": "unknown",
            "female_parent_phenotype": "unknown",
            "brood_sizes": [4],
            "replicates": 1,
            "candidate_model_ids": ["M1_AR", "M0_NONGENETIC"],
        }

        for role, invalid in (("male", "garbage"), ("female", "also-garbage")):
            with self.subTest(role=role):
                request = recovery_request([copy.deepcopy(experiment)])
                request["goodness_of_fit_replicates"] = 10
                request["misspecification_alpha"] = 0.1
                request["experiments"][0][f"{role}_genotype"] = invalid
                with self.assertRaisesRegex(
                    ValidationError,
                    "^simulation_non_genetic_genotype_must_be_NA:"
                    f"SIM-M0:{role}$",
                ):
                    dispatch("simulate-recovery", request)

    def test_parent_phenotype_must_match_generating_genotype(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[32], replicates=3)]
        )
        request["experiments"][0]["male_parent_phenotype"] = "target"

        with self.assertRaisesRegex(
            ValidationError,
            "^simulation_parent_phenotype_genotype_incompatible:"
            "SIM-AR:male$",
        ):
            dispatch("simulate-recovery", request)

    def test_bootstrap_draw_budget_fails_before_sampling(self) -> None:
        request = recovery_request(
            [ar_experiment(brood_sizes=[100], replicates=1)]
        )
        request["goodness_of_fit_replicates"] = 10_000

        with patch.object(
            simulation,
            "_goodness_of_fit",
            side_effect=AssertionError("budget guard was bypassed"),
        ):
            with self.assertRaisesRegex(
                ValidationError,
                "^simulation_bootstrap_work_budget_exceeded$",
            ):
                dispatch("simulate-recovery", request)

    def test_all_parent_genotypes_validate_before_any_sampling(self) -> None:
        valid = ar_experiment(brood_sizes=[16], replicates=1)
        invalid = ar_experiment(brood_sizes=[16], replicates=1)
        invalid["experiment_id"] = "SIM-LATE-INVALID"
        invalid["male_genotype"] = "invalid"
        invalid["male_parent_phenotype"] = "unknown"
        invalid["female_parent_phenotype"] = "unknown"
        request = recovery_request([valid, invalid])

        with patch.object(
            simulation,
            "_multinomial_counts",
            side_effect=AssertionError("sampling started before preflight ended"),
        ):
            with self.assertRaisesRegex(
                ValidationError,
                "^autosomal_genotype_invalid$",
            ):
                dispatch("simulate-recovery", request)

    def test_recovery_is_order_independent_when_models_tie(self) -> None:
        experiment = {
            "experiment_id": "SIM-TIE",
            "generating_model": {
                "model_id": "M2_AD",
                "mode": "autosomal_dominant",
                "target_label": "target",
            },
            "male_genotype": "N/N",
            "female_genotype": "N/N",
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "brood_sizes": [32],
            "replicates": 4,
            "candidate_model_ids": list(ALL_CANDIDATES),
        }
        reordered = copy.deepcopy(experiment)
        reordered["candidate_model_ids"] = [
            "M3_INC",
            "M5_ZD",
            "M2_AD",
            "M1_AR",
            "M4_ZR",
            "M0_NONGENETIC",
        ]

        first = self.simulate(recovery_request([experiment]))[
            "experiments"
        ][0]["by_brood_size"][0]
        second = self.simulate(recovery_request([reordered]))[
            "experiments"
        ][0]["by_brood_size"][0]

        self.assertEqual(first["recovery"], second["recovery"])
        self.assertEqual(
            [
                replicate["best_models"]
                for replicate in first["replicate_results"]
            ],
            [
                replicate["best_models"]
                for replicate in second["replicate_results"]
            ],
        )
        self.assertTrue(
            any(
                len(replicate["best_models"]) > 1
                and replicate["best_model"] is None
                for replicate in first["replicate_results"]
            )
        )


if __name__ == "__main__":
    unittest.main()

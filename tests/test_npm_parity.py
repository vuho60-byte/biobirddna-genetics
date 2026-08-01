from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.errors import ValidationError
from biobirddna_genetics.service import dispatch


CASES = [
    ("cross", "cross_autosomal_recessive.json"),
    ("cross", "cross_z_linked_recessive.json"),
    ("infer", "infer_z_linked_case.json"),
    ("pigment", "pigment_progressive_greying.json"),
    ("rank-mates", "mate_ranking_case.json"),
    ("normalize-records", "normalize_breeding_records.json"),
]


def assert_json_parity(
    case: unittest.TestCase,
    expected: Any,
    actual: Any,
    path: str = "$",
) -> None:
    if (
        isinstance(expected, (int, float))
        and not isinstance(expected, bool)
        and isinstance(actual, (int, float))
        and not isinstance(actual, bool)
    ):
        case.assertTrue(
            math.isclose(
                float(expected),
                float(actual),
                rel_tol=1e-12,
                abs_tol=1e-12,
            ),
            f"{path}: expected {expected!r}, got {actual!r}",
        )
        return
    case.assertIs(
        type(actual),
        type(expected),
        f"{path}: expected type {type(expected).__name__}, "
        f"got {type(actual).__name__}",
    )
    if isinstance(expected, dict):
        case.assertEqual(
            set(actual),
            set(expected),
            f"{path}: object fields differ",
        )
        for key in expected:
            assert_json_parity(
                case,
                expected[key],
                actual[key],
                f"{path}.{key}",
            )
    elif isinstance(expected, list):
        case.assertEqual(
            len(actual),
            len(expected),
            f"{path}: array length differs",
        )
        for index, (expected_item, actual_item) in enumerate(
            zip(expected, actual, strict=True)
        ):
            assert_json_parity(
                case,
                expected_item,
                actual_item,
                f"{path}[{index}]",
            )
    else:
        case.assertEqual(actual, expected, f"{path}: values differ")


class NpmParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.node = shutil.which("node")
        if cls.node is None:
            raise AssertionError("node runtime is required for npm parity tests")

    def node_dispatch(self, operation: str, request: dict[str, Any]) -> Any:
        node_cli = ROOT / "npm" / "bin" / "biobirddna-genetics.mjs"
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "request.json"
            input_path.write_text(
                json.dumps(request, ensure_ascii=False),
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    self.node,
                    str(node_cli),
                    operation,
                    "--input",
                    str(input_path),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        envelope = json.loads(completed.stdout)
        self.assertEqual(envelope["engine_version"], "0.2.0")
        self.assertEqual(envelope["operation"], operation)
        self.assertTrue(envelope["ok"])
        return envelope["result"]

    def node_error(self, operation: str, request: dict[str, Any]) -> str:
        node_cli = ROOT / "npm" / "bin" / "biobirddna-genetics.mjs"
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "request.json"
            input_path.write_text(
                json.dumps(request, ensure_ascii=False),
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    self.node,
                    str(node_cli),
                    operation,
                    "--input",
                    str(input_path),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
        self.assertEqual(completed.returncode, 2, completed.stdout)
        return str(json.loads(completed.stderr)["error"])

    def test_deterministic_public_operations_match_python_contract(self) -> None:
        for operation, filename in CASES:
            with self.subTest(operation=operation, filename=filename):
                input_path = ROOT / "examples" / "input" / filename
                request = json.loads(input_path.read_text(encoding="utf-8"))
                expected = dispatch(operation, request)
                actual = self.node_dispatch(operation, request)
                assert_json_parity(self, expected, actual)

    def test_simulation_recovery_seed_stream_and_candidate_order_match(self) -> None:
        request = {
            "schema_version": "1.0.0",
            "base_seed": 7,
            "goodness_of_fit_replicates": 10,
            "misspecification_alpha": 0.1,
            "experiments": [
                {
                    "experiment_id": "SIM-AUDIT",
                    "generating_model": {
                        "model_id": "M1_AR",
                        "mode": "autosomal_recessive",
                        "target_label": "target",
                    },
                    "male_genotype": "N/m",
                    "female_genotype": "N/m",
                    "male_parent_phenotype": "normal",
                    "female_parent_phenotype": "normal",
                    "brood_sizes": [2],
                    "replicates": 1,
                    "candidate_model_ids": [
                        "M1_AR",
                        "M0_NONGENETIC",
                    ],
                }
            ],
        }

        expected = dispatch("simulate-recovery", request)
        actual = self.node_dispatch("simulate-recovery", request)
        self.assertEqual(actual, expected)

    def test_unknown_sex_inference_matches_python_collapse(self) -> None:
        request = {
            "target_label": "target",
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "offspring": [
                {"sex": "unknown", "phenotype": "target", "count": 10},
                {"sex": "unknown", "phenotype": "normal", "count": 2},
            ],
            "model_priors": {
                "M1_AR": 1,
                "M2_AD": 1,
                "M3_INC": 1,
                "M4_ZR": 1,
                "M5_ZD": 1,
                "M0_NONGENETIC": 1,
            },
        }

        assert_json_parity(
            self,
            dispatch("infer", request),
            self.node_dispatch("infer", request),
        )

    def test_inference_accepts_nonnegative_relative_prior_weights(self) -> None:
        request = json.loads(
            (
                ROOT
                / "examples"
                / "input"
                / "infer_z_linked_case.json"
            ).read_text(encoding="utf-8")
        )
        request["model_priors"]["M4_ZR"] = 2.5

        assert_json_parity(
            self,
            dispatch("infer", request),
            self.node_dispatch("infer", request),
        )

    def test_non_genetic_cross_and_relative_weights_match_python(self) -> None:
        requests = [
            {
                "model": {
                    "model_id": "M0_NONGENETIC",
                    "mode": "non_genetic",
                    "target_label": "target",
                },
                "male_genotype": "NA",
                "female_genotype": "NA",
            },
            {
                "model": {
                    "model_id": "M1_AR",
                    "mode": "autosomal_recessive",
                    "target_label": "target",
                },
                "male_genotypes": {"N/N": 2, "N/m": 1},
                "female_genotypes": {"m/m": 3},
            },
        ]

        for request in requests:
            with self.subTest(request=request):
                assert_json_parity(
                    self,
                    dispatch("cross", request),
                    self.node_dispatch("cross", request),
                )

    def test_equal_inference_ties_match_deterministic_model_order(self) -> None:
        shared = {
            "mode": "autosomal_recessive",
            "target_label": "target",
        }
        request = {
            "models": [
                {"model_id": "Z_AR", **shared},
                {"model_id": "A_AR", **shared},
                {
                    "model_id": "M0_NONGENETIC",
                    "mode": "non_genetic",
                    "target_label": "target",
                },
            ],
            "male_parent_phenotype": "unknown",
            "female_parent_phenotype": "unknown",
            "offspring": [
                {"sex": "unknown", "phenotype": "target", "count": 2}
            ],
        }

        assert_json_parity(
            self,
            dispatch("infer", request),
            self.node_dispatch("infer", request),
        )

    def test_inference_request_validation_matches_schema_and_node(self) -> None:
        base = {
            "target_label": "target",
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "offspring": [{"phenotype": "normal", "count": 1}],
        }
        assert_json_parity(
            self,
            dispatch("infer", base),
            self.node_dispatch("infer", base),
        )

        invalid_cases = [
            (
                {**base, "unknown": True},
                "inference_request_unknown_field:unknown",
            ),
            ({**base, "target_label": 7}, "target_label_required"),
            (
                {**base, "male_parent_phenotype": 7},
                "male_parent_phenotype_required",
            ),
            (
                {
                    **base,
                    "offspring": [
                        {"phenotype": "normal", "count": "1"}
                    ],
                },
                "offspring_count_invalid",
            ),
            (
                {
                    **base,
                    "offspring": [
                        {"phenotype": "normal", "count": True}
                    ],
                },
                "offspring_count_invalid",
            ),
            (
                {
                    **base,
                    "offspring": [
                        {
                            "phenotype": "normal",
                            "count": 1,
                            "unknown": True,
                        }
                    ],
                },
                "offspring_observation_unknown_field:unknown",
            ),
        ]

        for request, error in invalid_cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValidationError, f"^{error}$"):
                    dispatch("infer", request)
                self.assertEqual(self.node_error("infer", request), error)

    def test_incompatible_model_log_evidence_is_json_null_in_both_runtimes(
        self,
    ) -> None:
        request = {
            "target_label": "target",
            "male_parent_phenotype": "unsupported",
            "female_parent_phenotype": "unsupported",
            "offspring": [
                {"sex": "unknown", "phenotype": "normal", "count": 1}
            ],
        }

        assert_json_parity(
            self,
            dispatch("infer", request),
            self.node_dispatch("infer", request),
        )

    def test_zero_likelihood_configuration_does_not_discard_model(self) -> None:
        request = {
            "models": [
                {
                    "model_id": "M1_AR",
                    "mode": "autosomal_recessive",
                    "target_label": "target",
                    "observation_error": 0,
                },
                {
                    "model_id": "M0_NONGENETIC",
                    "mode": "non_genetic",
                    "target_label": "target",
                },
            ],
            "male_parent_phenotype": "unknown",
            "female_parent_phenotype": "unknown",
            "offspring": [
                {"sex": "unknown", "phenotype": "target", "count": 1}
            ],
        }

        expected = dispatch("infer", request)
        ar_score = next(
            item
            for item in expected["ranking"]
            if item["model_id"] == "M1_AR"
        )
        self.assertIsNotNone(ar_score["log_evidence"])
        self.assertEqual(len(ar_score["configurations"]), 9)
        assert_json_parity(
            self,
            expected,
            self.node_dispatch("infer", request),
        )

    def test_brown_shift_pigment_case_matches_full_python_rules(self) -> None:
        request = {
            "white_area_ratio": 0.0,
            "brown_shift_score": 1.0,
            "grey_shift_score": 0.0,
            "dark_expansion_score": 0.0,
            "red_patch_change_score": 0.0,
            "pattern_retained": True,
            "iris_red_confirmed": False,
            "flash_red_eye_risk": False,
            "progressive_across_molts": False,
            "present_in_first_plumage": True,
            "localized_after_injury": False,
            "bilateral_symmetry": 1.0,
            "image_quality": "Q4",
        }

        assert_json_parity(
            self,
            dispatch("pigment", request),
            self.node_dispatch("pigment", request),
        )

    def test_q0_and_q1_pigment_abstention_matches_node(self) -> None:
        for image_quality in ("Q0", "Q1"):
            with self.subTest(image_quality=image_quality):
                request = {
                    "white_area_ratio": 0.99,
                    "brown_shift_score": 1.0,
                    "grey_shift_score": 1.0,
                    "dark_expansion_score": 1.0,
                    "red_patch_change_score": 1.0,
                    "pattern_retained": True,
                    "iris_red_confirmed": True,
                    "progressive_across_molts": True,
                    "image_quality": image_quality,
                }
                python_result = dispatch("pigment", request)
                node_result = self.node_dispatch("pigment", request)

                self.assertEqual(python_result["primary_hypothesis"], "UNK")
                self.assertIs(python_result["needs_recapture"], True)
                assert_json_parity(self, python_result, node_result)

    def test_pigment_request_rejects_schema_invalid_coercions(self) -> None:
        invalid_cases = [
            ({"image_quality": "BOGUS"}, "image_quality_invalid"),
            (
                {"white_area_ratio": "0.5"},
                "white_area_ratio_invalid",
            ),
            (
                {"iris_red_confirmed": "false"},
                "iris_red_confirmed_invalid",
            ),
            (
                {"pattern_retained": "false"},
                "pattern_retained_invalid",
            ),
            (
                {"unknown": True},
                "pigment_request_unknown_field:unknown",
            ),
        ]

        for request, error in invalid_cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValidationError, f"^{error}$"):
                    dispatch("pigment", request)
                self.assertEqual(self.node_error("pigment", request), error)

    def test_rank_mates_preserves_pedigree_and_welfare_hard_blocks(self) -> None:
        request = {
            "pedigree": [
                {"individual_id": "S", "sex": "male"},
                {
                    "individual_id": "D",
                    "sex": "female",
                    "sire_id": "S",
                },
            ],
            "males": [
                {
                    "individual_id": "S",
                    "sex": "male",
                    "welfare_blocked": True,
                    "genotype_posteriors": {},
                }
            ],
            "females": [
                {
                    "individual_id": "D",
                    "sex": "female",
                    "genotype_posteriors": {},
                }
            ],
            "hypotheses": [
                {
                    "weight": 1.0,
                    "model": {
                        "model_id": "M1_AR",
                        "mode": "autosomal_recessive",
                        "normal_allele": "N",
                        "variant_allele": "m",
                        "normal_label": "normal",
                        "target_label": "target",
                    },
                }
            ],
            "max_offspring_inbreeding": 0.0625,
            "block_direct_ancestor_and_full_sibling": True,
        }

        assert_json_parity(
            self,
            dispatch("rank-mates", request),
            self.node_dispatch("rank-mates", request),
        )

    def test_rank_mates_accepts_relative_hypothesis_weights_above_one(self) -> None:
        request = json.loads(
            (
                ROOT / "examples" / "input" / "mate_ranking_case.json"
            ).read_text(encoding="utf-8")
        )
        for hypothesis in request["hypotheses"]:
            hypothesis["weight"] *= 10

        assert_json_parity(
            self,
            dispatch("rank-mates", request),
            self.node_dispatch("rank-mates", request),
        )

    def test_cross_validation_rejects_label_collisions_and_null_weights(
        self,
    ) -> None:
        cases = [
            (
                {
                    "model": {
                        "model_id": "M1_AR",
                        "mode": "autosomal_recessive",
                        "normal_label": "same",
                        "target_label": "same",
                    },
                    "male_genotype": "N/m",
                    "female_genotype": "m/m",
                },
                "phenotype_labels_must_be_distinct",
            ),
            (
                {
                    "model": {
                        "model_id": "M1_AR",
                        "mode": "autosomal_recessive",
                        "target_label": "target",
                    },
                    "male_genotypes": {"N/N": None, "N/m": 1.0},
                    "female_genotypes": {"m/m": 1.0},
                },
                "probability_distribution_invalid",
            ),
        ]
        for request, error in cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValidationError, f"^{error}$"):
                    dispatch("cross", request)
                self.assertEqual(self.node_error("cross", request), error)

    def test_inference_count_has_operational_work_budget(self) -> None:
        request = {
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "offspring": [
                {
                    "sex": "unknown",
                    "phenotype": "normal",
                    "count": 1_000_001,
                }
            ],
        }
        with self.assertRaisesRegex(
            ValidationError,
            "^offspring_count_invalid$",
        ):
            dispatch("infer", request)
        self.assertEqual(
            self.node_error("infer", request),
            "offspring_count_invalid",
        )

    def test_inference_log_evidence_matches_at_work_budget_boundary(self) -> None:
        request = {
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "offspring": [
                {
                    "sex": "unknown",
                    "phenotype": "normal",
                    "count": 1_000_000,
                }
            ],
        }
        python_ranking = {
            item["model_id"]: item
            for item in dispatch("infer", request)["ranking"]
        }
        node_ranking = {
            item["model_id"]: item
            for item in self.node_dispatch("infer", request)["ranking"]
        }
        self.assertEqual(set(node_ranking), set(python_ranking))
        for model_id, expected in python_ranking.items():
            with self.subTest(model_id=model_id):
                actual = node_ranking[model_id]
                if expected["log_evidence"] is None:
                    self.assertIsNone(actual["log_evidence"])
                else:
                    self.assertLessEqual(
                        abs(
                            expected["log_evidence"]
                            - actual["log_evidence"]
                        ),
                        1e-12,
                    )

    def test_cross_model_rejects_schema_invalid_coercions(self) -> None:
        base = {
            "model": {
                "model_id": "M1_AR",
                "mode": "autosomal_recessive",
                "target_label": "target",
            },
            "male_genotype": "N/m",
            "female_genotype": "m/m",
        }
        cases = [
            ({"target_label": None}, "target_label_required"),
            ({"model_id": None}, "model_id_required"),
            ({"penetrance": "0.5"}, "penetrance_invalid"),
        ]
        for change, error in cases:
            request = {
                **base,
                "model": {**base["model"], **change},
            }
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValidationError, f"^{error}$"):
                    dispatch("cross", request)
                self.assertEqual(self.node_error("cross", request), error)

        missing_id = {
            **base,
            "model": {
                key: value
                for key, value in base["model"].items()
                if key != "model_id"
            },
        }
        with self.assertRaisesRegex(ValidationError, "^model_id_required$"):
            dispatch("cross", missing_id)
        self.assertEqual(
            self.node_error("cross", missing_id),
            "model_id_required",
        )

        unknown_field = {
            **base,
            "model": {**base["model"], "unknown": True},
        }
        with self.assertRaisesRegex(
            ValidationError,
            "^model_unknown_field:unknown$",
        ):
            dispatch("cross", unknown_field)
        self.assertEqual(
            self.node_error("cross", unknown_field),
            "model_unknown_field:unknown",
        )

    def test_rank_mates_rejects_null_and_coercible_controls(self) -> None:
        base = json.loads(
            (
                ROOT / "examples" / "input" / "mate_ranking_case.json"
            ).read_text(encoding="utf-8")
        )
        cases = [
            (
                {"max_offspring_inbreeding": None},
                "max_offspring_inbreeding_invalid",
            ),
            (
                {"block_direct_ancestor_and_full_sibling": None},
                "block_direct_ancestor_and_full_sibling_invalid",
            ),
            (
                {
                    "hypotheses": [
                        {**base["hypotheses"][0], "weight": "1"}
                    ]
                },
                "model_weight_invalid",
            ),
        ]
        for change, error in cases:
            request = {**base, **change}
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValidationError, f"^{error}$"):
                    dispatch("rank-mates", request)
                self.assertEqual(self.node_error("rank-mates", request), error)

    def test_rank_mates_rejects_schema_invalid_profiles_and_unknowns(
        self,
    ) -> None:
        base = json.loads(
            (
                ROOT / "examples" / "input" / "mate_ranking_case.json"
            ).read_text(encoding="utf-8")
        )
        cases: list[tuple[dict[str, Any], str]] = []

        request = json.loads(json.dumps(base))
        request["males"][0]["individual_id"] = None
        cases.append((request, "individual_id_required"))

        request = json.loads(json.dumps(base))
        request["males"][0]["health_risk"] = "0.1"
        cases.append((request, "health_risk_invalid"))

        request = json.loads(json.dumps(base))
        request["males"][0]["welfare_blocked"] = None
        cases.append((request, "welfare_blocked_invalid"))

        request = json.loads(json.dumps(base))
        request["unknown"] = True
        cases.append(
            (request, "rank_mates_request_unknown_field:unknown")
        )

        request = json.loads(json.dumps(base))
        request["males"][0]["unknown"] = True
        cases.append((request, "genetic_profile_unknown_field:unknown"))

        request = json.loads(json.dumps(base))
        request["weights"]["unknown"] = 1
        cases.append((request, "weights_unknown_field:unknown"))

        request = json.loads(json.dumps(base))
        request["hypotheses"][0]["unknown"] = True
        cases.append((request, "model_hypothesis_unknown_field:unknown"))

        for request, error in cases:
            with self.subTest(error=error):
                with self.assertRaisesRegex(ValidationError, f"^{error}$"):
                    dispatch("rank-mates", request)
                self.assertEqual(self.node_error("rank-mates", request), error)


if __name__ == "__main__":
    unittest.main()

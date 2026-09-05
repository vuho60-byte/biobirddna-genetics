from __future__ import annotations

import json
import math
import random
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.errors import ValidationError
from biobirddna_genetics.inference import (
    _log_gamma,
    _logsumexp,
    _neumaier_sum,
    _portable_exp,
    _portable_log,
    _portable_log2,
)
from biobirddna_genetics.service import dispatch
from biobirddna_genetics.simulation import (
    _entropy_bits,
    _order_independent_mean,
    _quantile,
)


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
                encoding="utf-8",
                capture_output=True,
                check=False,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        envelope = json.loads(completed.stdout)
        self.assertEqual(envelope["engine_version"], "0.3.0")
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
                encoding="utf-8",
                capture_output=True,
                check=False,
            )
        self.assertEqual(completed.returncode, 2, completed.stdout)
        envelope = json.loads(completed.stderr)
        self.assertIs(envelope["ok"], False)
        self.assertEqual(envelope["operation"], operation)
        self.assertIsInstance(envelope["message_vi"], str)
        self.assertTrue(envelope["message_vi"])
        return str(envelope["error"])

    def test_error_envelope_matches_python_cli_including_message_vi(
        self,
    ) -> None:
        request = {
            "model": {"model_id": "M1_AR", "mode": "autosomal_recessive"},
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "request.json"
            input_path.write_text(
                json.dumps(request, ensure_ascii=False),
                encoding="utf-8",
            )
            python_completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(ROOT / "scripts" / "biobirddna_cli.py"),
                    "cross",
                    "--input",
                    str(input_path),
                ],
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )
            node_completed = subprocess.run(
                [
                    self.node,
                    str(ROOT / "npm" / "bin" / "biobirddna-genetics.mjs"),
                    "cross",
                    "--input",
                    str(input_path),
                ],
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )
        self.assertEqual(python_completed.returncode, 2)
        self.assertEqual(node_completed.returncode, 2)
        python_envelope = json.loads(python_completed.stderr)
        node_envelope = json.loads(node_completed.stderr)
        self.assertEqual(
            set(python_envelope),
            {"ok", "operation", "error", "message_vi"},
        )
        self.assertEqual(set(node_envelope), set(python_envelope))
        for key in ("ok", "operation", "error", "message_vi"):
            self.assertEqual(
                node_envelope[key],
                python_envelope[key],
                f"error envelope field {key} differs across runtimes",
            )
        self.assertEqual(
            python_envelope["error"],
            "male_and_female_genotype_required",
        )

    def test_internal_engine_fault_exits_70_in_both_runtimes(self) -> None:
        # An injected engine bug (TypeError from dispatch) must NOT be
        # reported as an input error: both CLIs must exit 70 with the
        # internal_error:<Type> slug and the traceback/stack preserved.
        # The Node side runs the real shipped bin file against a shim
        # src/index.mjs whose dispatch throws, so no production test hook
        # is needed in either runtime.
        example = ROOT / "examples" / "input" / "cross_autosomal_recessive.json"
        python_driver = (
            "import sys\n"
            f"sys.path.insert(0, {str(ROOT / 'src')!r})\n"
            "import biobirddna_genetics.cli as cli\n"
            "def boom(operation, request):\n"
            "    raise TypeError('injected internal bug')\n"
            "cli.dispatch = boom\n"
            "sys.exit(cli.main(['cross', '--input', sys.argv[1]]))\n"
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            (temp / "bin").mkdir()
            (temp / "src").mkdir()
            shutil.copy2(
                ROOT / "npm" / "bin" / "biobirddna-genetics.mjs",
                temp / "bin" / "biobirddna-genetics.mjs",
            )
            real_index = (ROOT / "npm" / "src" / "index.mjs").resolve()
            (temp / "src" / "index.mjs").write_text(
                "export { ENGINE_VERSION, BioBirdDNAError } from "
                + json.dumps(real_index.as_uri())
                + ";\n"
                "export function dispatch() {"
                " throw new TypeError('injected internal bug'); }\n",
                encoding="utf-8",
            )
            python_completed = subprocess.run(
                [sys.executable, "-B", "-c", python_driver, str(example)],
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )
            node_completed = subprocess.run(
                [
                    self.node,
                    str(temp / "bin" / "biobirddna-genetics.mjs"),
                    "cross",
                    "--input",
                    str(example),
                ],
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )
        self.assertEqual(python_completed.returncode, 70, python_completed.stderr)
        self.assertEqual(node_completed.returncode, 70, node_completed.stderr)
        self.assertIn(
            "Traceback (most recent call last)", python_completed.stderr
        )
        self.assertIn("TypeError: injected internal bug", python_completed.stderr)
        self.assertIn("TypeError: injected internal bug", node_completed.stderr)
        python_envelope = json.loads(
            python_completed.stderr[
                python_completed.stderr.index('{\n  "ok"'):
            ]
        )
        node_envelope = json.loads(
            node_completed.stderr.strip().splitlines()[-1]
        )
        self.assertEqual(set(python_envelope), set(node_envelope))
        for envelope in (python_envelope, node_envelope):
            self.assertIs(envelope["ok"], False)
            self.assertEqual(envelope["operation"], "cross")
            self.assertEqual(envelope["error"], "internal_error:TypeError")
        self.assertEqual(
            python_envelope["message_vi"], node_envelope["message_vi"]
        )

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

    def test_configuration_prior_warning_matches_python(self) -> None:
        # Sex-unknown 6:2 from normal x normal parents: ranking gap between
        # M4_ZR and M1_AR is attributable to configuration counts alone, so
        # both runtimes must emit the identical statistical caution.
        request = {
            "target_label": "target",
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "offspring": [
                {"sex": "unknown", "phenotype": "normal", "count": 6},
                {"sex": "unknown", "phenotype": "target", "count": 2},
            ],
        }
        expected = dispatch("infer", request)
        self.assertEqual(
            expected["warnings"],
            ["model_ranking_may_reflect_uniform_configuration_prior_not_data"],
        )
        actual = self.node_dispatch("infer", request)
        self.assertEqual(actual["warnings"], expected["warnings"])
        assert_json_parity(self, expected, actual)

        decisive = json.loads(
            (
                ROOT / "examples" / "input" / "infer_z_linked_case.json"
            ).read_text(encoding="utf-8")
        )
        expected_decisive = dispatch("infer", decisive)
        self.assertEqual(expected_decisive["warnings"], [])
        actual_decisive = self.node_dispatch("infer", decisive)
        self.assertEqual(actual_decisive["warnings"], [])
        assert_json_parity(self, expected_decisive, actual_decisive)

    def test_dunder_target_labels_match_python_in_both_runtimes(self) -> None:
        # JS prototype-chain keys must behave like any other label
        # (regression for silent Node corruption with "__proto__").
        for label in ("__proto__", "constructor", "prototype"):
            with self.subTest(label=label):
                request = {
                    "target_label": label,
                    "male_parent_phenotype": "normal",
                    "female_parent_phenotype": "normal",
                    "offspring": [
                        {"sex": "unknown", "phenotype": "normal", "count": 6},
                        {"sex": "unknown", "phenotype": label, "count": 2},
                    ],
                }
                expected = dispatch("infer", request)
                self.assertEqual(expected["best_model"], "M4_ZR")
                assert_json_parity(
                    self,
                    expected,
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

    def test_z_linked_female_probability_rejection_matches_node(self) -> None:
        cross_request = {
            "model": {
                "model_id": "M4_ZR",
                "mode": "z_linked_recessive",
                "target_label": "target",
                "female_probability": 0.6,
            },
            "male_genotype": "N/m",
            "female_genotype": "N/W",
        }
        with self.assertRaisesRegex(
            ValidationError,
            "^z_linked_female_probability_must_be_0_5$",
        ):
            dispatch("cross", cross_request)
        self.assertEqual(
            self.node_error("cross", cross_request),
            "z_linked_female_probability_must_be_0_5",
        )

        infer_request = {
            "target_label": "target",
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "offspring": [
                {"sex": "unknown", "phenotype": "target", "count": 2}
            ],
            "models": [
                {
                    "model_id": "M1_AR",
                    "mode": "autosomal_recessive",
                    "target_label": "target",
                    "female_probability": 0.7,
                },
                {
                    "model_id": "M4_ZR",
                    "mode": "z_linked_recessive",
                    "target_label": "target",
                },
                {
                    "model_id": "M0_NONGENETIC",
                    "mode": "non_genetic",
                    "target_label": "target",
                },
            ],
        }
        with self.assertRaisesRegex(
            ValidationError,
            "^female_probability_must_be_0_5_when_z_linked_model_present$",
        ):
            dispatch("infer", infer_request)
        self.assertEqual(
            self.node_error("infer", infer_request),
            "female_probability_must_be_0_5_when_z_linked_model_present",
        )

        autosomal_only = {
            **infer_request,
            "models": [
                model
                for model in infer_request["models"]
                if model["model_id"] != "M4_ZR"
            ],
        }
        assert_json_parity(
            self,
            dispatch("infer", autosomal_only),
            self.node_dispatch("infer", autosomal_only),
        )

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


NODE_PRIMITIVE_RUNNER = """
import { __numericPrimitivesForTest as P } from "%INDEX_URL%";
import { readFileSync } from "node:fs";
const inp = JSON.parse(readFileSync(0, "utf8"));
const enc = (v) => Number.isFinite(v) ? v : String(v);
console.log(JSON.stringify({
  portable_log: inp.portable_log.map((v) => enc(P.portableLog(v))),
  portable_exp: inp.portable_exp.map((v) => enc(P.portableExp(v))),
  portable_log2: inp.portable_log.map((v) => enc(P.portableLog2(v))),
  log_gamma: inp.log_gamma.map((v) => enc(P.logGamma(v))),
  logsumexp: inp.logsumexp.map((vs) => enc(P.logSumExp(
    vs.map((x) => x === "-Infinity" ? -Infinity : x)
  ))),
  entropy: inp.entropy.map((vs) => enc(P.entropy(
    Object.fromEntries(vs.map((x, i) => [String(i), x]))
  ))),
  quantile: inp.quantile.map(([vs, p]) => enc(P.quantile(vs, p))),
  mean: inp.mean.map((vs) => enc(P.orderIndependentMean(vs))),
  neumaier: inp.neumaier.map((vs) => enc(P.neumaierSum(vs))),
}));
"""


class NumericPrimitiveParityTests(unittest.TestCase):
    """Cross-runtime bit-identity sweep for the shared numeric primitives.

    Both runtimes intentionally run the same fdlibm-derived portable
    log/exp (npm/src/index.mjs portableLog/portableExp mirrored in
    src/biobirddna_genetics/inference.py), because glibc and V8 libm
    disagree in the last bit for roughly 3-10% of inputs. Threshold here
    is therefore exact bit equality (max abs diff 0.0) — do not loosen.
    """

    maxDiff = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.node = shutil.which("node")
        if cls.node is None:
            raise AssertionError("node runtime is required for npm parity tests")

    @staticmethod
    def _sweep_inputs() -> dict[str, Any]:
        rng = random.Random(20260807)
        log_inputs = [
            5e-324,
            2.2250738585072014e-308,
            1e-300,
            1e-12,
            0.001,
            0.05,
            1.0 / 3.0,
            0.5,
            0.9999999999999999,
            1.0,
            1.0000000000000002,
            2.0,
            math.pi,
            1e6,
            1e300,
            1.7976931348623157e308,
        ]
        log_inputs += [rng.random() for _ in range(400)]
        log_inputs += [rng.random() * 1e6 for _ in range(400)]
        exp_inputs = [
            -746.0,
            -745.2,
            -709.0,
            -700.0,
            -100.0,
            -1.0,
            -1e-15,
            0.0,
            1e-15,
            0.5,
            1.0,
            100.0,
            709.0,
        ]
        exp_inputs += [-rng.random() * 50.0 for _ in range(400)]
        exp_inputs += [(rng.random() - 0.5) * 2.0 for _ in range(400)]
        log_gamma_inputs = [
            0.5,
            0.75,
            1.0,
            1.5,
            2.0,
            3.0,
            7.5,
            41.0,
            101.0,
            1001.0,
            100001.0,
            1000001.0,
        ]
        log_gamma_inputs += [1.0 + index / 7.0 for index in range(1, 100)]
        logsumexp_inputs = [
            [-1.2, -3.4, -0.5],
            [-700.0, -701.5, -702.3, -690.1],
            [0.0, 0.0, 0.0],
            ["-Infinity", -2.0, -1.0],
            ["-Infinity"],
            [-1e-15, -2e-15, -3.3e-13, -0.9999999],
            [math.log(1.0 / 6.0)] * 6,
            [-(index * 0.037) % 5 - index * 1e-13 for index in range(50)],
        ]
        entropy_inputs = [
            [0.5, 0.5],
            [0.25, 0.25, 0.25, 0.25],
            [0.9999999999, 1e-10],
            [1e-300, 1.0 - 1e-300],
            [0.527, 0.265, 0.1, 0.05, 0.033, 0.025],
            [1.0 / 3.0] * 3,
            [1.0 / 7.0, 2.0 / 7.0, 4.0 / 7.0],
        ]
        quantile_inputs = [
            ([0.1], 0.5),
            ([0.3, 0.9], 0.05),
            ([0.3, 0.9], 0.95),
            ([0.1, 0.2, 0.7, 0.9], 0.5),
            ([1.0 / 3.0, 1.0 / 7.0, 2.0 / 3.0, 0.123456789, 0.9], 0.05),
            ([1.0 / 3.0, 1.0 / 7.0, 2.0 / 3.0, 0.123456789, 0.9], 0.95),
            ([5.5e-13, 1.0, 2.0], 0.33),
            ([0.0, 1.0], 0.0),
            ([0.0, 1.0], 1.0),
        ]
        mean_inputs = [
            [0.1, 0.2, 0.30000000000000004],
            [1e16, 1.0, -1e16, 2.0],
            [rng.random() for _ in range(31)],
        ]
        neumaier_inputs = [
            [1e100, 1.0, -1e100, 1.0],
            [0.1] * 10,
            [rng.random() - 0.5 for _ in range(64)],
        ]
        return {
            "portable_log": log_inputs,
            "portable_exp": exp_inputs,
            "log_gamma": log_gamma_inputs,
            "logsumexp": logsumexp_inputs,
            "entropy": entropy_inputs,
            "quantile": quantile_inputs,
            "mean": mean_inputs,
            "neumaier": neumaier_inputs,
        }

    def _node_primitives(self, payload: dict[str, Any]) -> dict[str, Any]:
        index_url = (ROOT / "npm" / "src" / "index.mjs").as_uri()
        with tempfile.TemporaryDirectory() as temp_dir:
            runner = Path(temp_dir) / "runner.mjs"
            runner.write_text(
                NODE_PRIMITIVE_RUNNER.replace("%INDEX_URL%", index_url),
                encoding="utf-8",
            )
            completed = subprocess.run(
                [self.node, str(runner)],
                input=json.dumps(payload),
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                capture_output=True,
                check=False,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(completed.stdout)

    @staticmethod
    def _bits(value: float) -> bytes:
        return struct.pack("<d", value)

    def test_numeric_primitives_are_bit_identical_across_runtimes(self) -> None:
        payload = self._sweep_inputs()
        node = self._node_primitives(payload)
        decode = lambda value: (  # noqa: E731
            float(value) if isinstance(value, str) else value
        )
        python_results: dict[str, list[float]] = {
            "portable_log": [_portable_log(v) for v in payload["portable_log"]],
            "portable_exp": [_portable_exp(v) for v in payload["portable_exp"]],
            "portable_log2": [
                _portable_log2(v) for v in payload["portable_log"]
            ],
            "log_gamma": [_log_gamma(v) for v in payload["log_gamma"]],
            "logsumexp": [
                _logsumexp(
                    [
                        float("-inf") if item == "-Infinity" else item
                        for item in values
                    ]
                )
                for values in payload["logsumexp"]
            ],
            "entropy": [_entropy_bits(v) for v in payload["entropy"]],
            "quantile": [
                _quantile(values, p) for values, p in payload["quantile"]
            ],
            "mean": [_order_independent_mean(v) for v in payload["mean"]],
            "neumaier": [_neumaier_sum(v) for v in payload["neumaier"]],
        }
        threshold = 0.0  # exact bit identity — do not loosen
        report_lines = [
            "primitive parity sweep (threshold: exact bit identity, "
            "max_abs_diff must be 0.0):"
        ]
        failures: list[str] = []
        for name, python_values in python_results.items():
            node_values = [decode(value) for value in node[name]]
            self.assertEqual(len(node_values), len(python_values), name)
            mismatches = 0
            max_abs_diff = 0.0
            for index, (python_value, node_value) in enumerate(
                zip(python_values, node_values)
            ):
                if math.isnan(python_value) and math.isnan(node_value):
                    continue
                if self._bits(python_value) != self._bits(node_value):
                    mismatches += 1
                    if math.isfinite(python_value) and math.isfinite(
                        node_value
                    ):
                        max_abs_diff = max(
                            max_abs_diff,
                            abs(python_value - node_value),
                        )
                    else:
                        max_abs_diff = float("inf")
                    failures.append(
                        f"{name}[{index}]: python={python_value!r} "
                        f"node={node_value!r}"
                    )
            report_lines.append(
                f"  {name}: cases={len(python_values)} "
                f"bit_mismatches={mismatches} max_abs_diff={max_abs_diff}"
            )
        print("\n".join(report_lines))
        self.assertEqual(
            failures,
            [],
            "cross-runtime primitive divergence beyond threshold "
            f"{threshold}:\n" + "\n".join(failures[:20]),
        )


class NumericOrderAndScaleParityTests(unittest.TestCase):
    """Larger fixtures + order independence at exact (bit-level) equality."""

    maxDiff = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.node = shutil.which("node")
        if cls.node is None:
            raise AssertionError("node runtime is required for npm parity tests")

    def node_dispatch(self, operation: str, request: dict[str, Any]) -> Any:
        return NpmParityTests.node_dispatch(self, operation, request)  # type: ignore[arg-type]

    def test_large_brood_simulation_recovery_matches_node_exactly(self) -> None:
        # Big broods + several replicates exercise logGamma at large counts,
        # quantile interpolation, entropy, and the compensated logsumexp —
        # exactly the primitives that diverged on the pre-portable-math
        # numerics. Equality is exact (==), matching the existing
        # seed-stream parity test; do not loosen to a tolerance.
        request = {
            "schema_version": "1.0.0",
            "base_seed": 20260807,
            "goodness_of_fit_replicates": 25,
            "misspecification_alpha": 0.05,
            "experiments": [
                {
                    "experiment_id": "SIM-LARGE-ZR",
                    "generating_model": {
                        "model_id": "M4_ZR",
                        "mode": "z_linked_recessive",
                        "target_label": "target",
                    },
                    "male_genotype": "N/m",
                    "female_genotype": "N/W",
                    "male_parent_phenotype": "normal",
                    "female_parent_phenotype": "normal",
                    "brood_sizes": [24, 48],
                    "replicates": 6,
                    "candidate_model_ids": [
                        "M1_AR",
                        "M2_AD",
                        "M3_INC",
                        "M4_ZR",
                        "M5_ZD",
                        "M0_NONGENETIC",
                    ],
                }
            ],
        }

        expected = dispatch("simulate-recovery", request)
        actual = self.node_dispatch("simulate-recovery", request)
        self.assertEqual(actual, expected)

    def test_infer_model_order_is_irrelevant_in_both_runtimes(self) -> None:
        # Reversing the candidate model order must yield IDENTICAL results
        # (bit level) in each runtime, and both runtimes must agree exactly.
        shared = {"target_label": "target"}
        models = [
            {"model_id": "M1_AR", "mode": "autosomal_recessive", **shared},
            {"model_id": "M2_AD", "mode": "autosomal_dominant", **shared},
            {"model_id": "M3_INC", "mode": "incomplete_dominant", **shared},
            {"model_id": "M4_ZR", "mode": "z_linked_recessive", **shared},
            {"model_id": "M5_ZD", "mode": "z_linked_dominant", **shared},
            {"model_id": "M0_NONGENETIC", "mode": "non_genetic", **shared},
        ]
        base = {
            "target_label": "target",
            "male_parent_phenotype": "normal",
            "female_parent_phenotype": "normal",
            "offspring": [
                {"sex": "unknown", "phenotype": "normal", "count": 6},
                {"sex": "unknown", "phenotype": "target", "count": 2},
            ],
        }
        forward = {**base, "models": models}
        reversed_request = {**base, "models": list(reversed(models))}

        python_forward = dispatch("infer", forward)
        python_reversed = dispatch("infer", reversed_request)
        self.assertEqual(python_forward, python_reversed)

        node_forward = self.node_dispatch("infer", forward)
        node_reversed = self.node_dispatch("infer", reversed_request)
        self.assertEqual(node_forward, node_reversed)
        self.assertEqual(node_forward, python_forward)


if __name__ == "__main__":
    unittest.main()

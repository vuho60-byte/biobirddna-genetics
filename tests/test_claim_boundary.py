from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics import service
from biobirddna_genetics.enums import EvidenceGrade
from biobirddna_genetics.errors import EvidenceBoundaryError
from biobirddna_genetics.evidence import claim_allowed, safe_claim_language
from biobirddna_genetics.service import dispatch

CLI = ROOT / "scripts" / "biobirddna_cli.py"

EXAMPLE_BY_OPERATION = {
    "cross": "cross_autosomal_recessive.json",
    "infer": "infer_z_linked_case.json",
    "pigment": "pigment_progressive_greying.json",
    "rank-mates": "mate_ranking_case.json",
    "normalize-records": "normalize_breeding_records.json",
    "simulate-recovery": "simulate_recovery.json",
}


def load_example(filename: str) -> dict:
    path = ROOT / "examples" / "input" / filename
    return json.loads(path.read_text(encoding="utf-8"))


class ClaimBoundaryWiringTests(unittest.TestCase):
    """The evidence boundary must sit on the real output path (dispatch),
    not only in the test suite (review 3.1)."""

    def test_every_public_operation_emits_safe_claim_language(self) -> None:
        self.assertEqual(
            set(service._OPERATION_CLAIM_CEILINGS),
            set(EXAMPLE_BY_OPERATION),
        )
        for operation, filename in EXAMPLE_BY_OPERATION.items():
            with self.subTest(operation=operation):
                claim_type, grade = service._OPERATION_CLAIM_CEILINGS[
                    operation
                ]
                # The static ceiling itself must be a permitted claim: the
                # boundary never emits language above the grade.
                self.assertTrue(claim_allowed(claim_type, grade))
                result = dispatch(operation, load_example(filename))
                self.assertEqual(
                    result["limits"]["claim_language"],
                    safe_claim_language(claim_type, grade),
                )
                self.assertEqual(
                    result["limits"]["claim_language"],
                    claim_type,
                )

    def test_disallowed_claim_cannot_leave_dispatch(self) -> None:
        # If anyone remaps an operation to a claim above its evidence grade,
        # dispatch must raise before emitting a result. This test also fails
        # if the boundary wiring in dispatch() is reverted (the illegal
        # ceiling would then go unnoticed and no error would be raised).
        illegal_ceilings = {
            "infer": ("supported_inheritance_pattern", EvidenceGrade.G2),
            "pigment": ("candidate_inheritance_model", EvidenceGrade.G0),
            "cross": ("molecular_locus_supported", EvidenceGrade.G0),
        }
        for operation, ceiling in illegal_ceilings.items():
            with self.subTest(operation=operation, ceiling=ceiling):
                original = service._OPERATION_CLAIM_CEILINGS[operation]
                service._OPERATION_CLAIM_CEILINGS[operation] = ceiling
                try:
                    with self.assertRaises(EvidenceBoundaryError):
                        dispatch(
                            operation,
                            load_example(EXAMPLE_BY_OPERATION[operation]),
                        )
                finally:
                    service._OPERATION_CLAIM_CEILINGS[operation] = original

    def test_unmapped_operation_fails_closed(self) -> None:
        with self.assertRaisesRegex(
            EvidenceBoundaryError,
            "^operation_claim_ceiling_missing:bogus-op$",
        ):
            service._apply_claim_boundary("bogus-op", {"limits": {}})

    def test_valid_claims_leave_existing_limits_untouched(self) -> None:
        cross_result = dispatch("cross", load_example(EXAMPLE_BY_OPERATION["cross"]))
        self.assertEqual(
            cross_result["limits"]["gene_status"],
            "model_assumption_only",
        )
        self.assertEqual(
            cross_result["limits"]["inheritance_status"],
            "simulation_not_confirmation",
        )
        self.assertEqual(
            cross_result["limits"]["claim_language"],
            "simulate_cross",
        )
        infer_result = dispatch("infer", load_example(EXAMPLE_BY_OPERATION["infer"]))
        self.assertEqual(
            infer_result["limits"]["claim"],
            "candidate_inheritance_model_only",
        )
        self.assertEqual(infer_result["limits"]["gene_status"], "not_assessed")
        self.assertEqual(
            infer_result["limits"]["claim_language"],
            "candidate_inheritance_model",
        )

    def test_cli_output_carries_claim_language(self) -> None:
        # The public CLI envelope must carry the boundary-computed language,
        # so an agent wrapper is no longer the only protection.
        for operation, expected in (
            ("cross", "simulate_cross"),
            ("pigment", "phenotype_hypothesis"),
        ):
            with self.subTest(operation=operation):
                completed = subprocess.run(
                    [
                        sys.executable,
                        "-B",
                        str(CLI),
                        operation,
                        "--input",
                        str(
                            ROOT
                            / "examples"
                            / "input"
                            / EXAMPLE_BY_OPERATION[operation]
                        ),
                    ],
                    cwd=ROOT,
                    text=True,
                    encoding="utf-8",
                    capture_output=True,
                    check=False,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                envelope = json.loads(completed.stdout)
                self.assertTrue(envelope["ok"])
                self.assertEqual(
                    envelope["result"]["limits"]["claim_language"],
                    expected,
                )


if __name__ == "__main__":
    unittest.main()

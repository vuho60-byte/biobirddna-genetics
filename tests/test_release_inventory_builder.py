from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.rebuild_release_inventory import rebuild


ROOT = Path(__file__).resolve().parents[1]


class ReleaseInventoryBuilderTests(unittest.TestCase):
    def test_rebuild_is_deterministic_and_respects_self_exclusions(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "VERSION").write_text("0.2.0\n", encoding="utf-8")
            (root / "sample.txt").write_text("stable\n", encoding="utf-8")

            rebuild(root, package_name="BioBirdDNA_Genetics_Skill_v0.2.0")
            first_checksums = (root / "PACKAGE_CHECKSUMS.json").read_bytes()
            first_sums = (root / "SHA256SUMS.txt").read_bytes()

            rebuild(root, package_name="BioBirdDNA_Genetics_Skill_v0.2.0")

            self.assertEqual(
                (root / "PACKAGE_CHECKSUMS.json").read_bytes(),
                first_checksums,
            )
            self.assertEqual((root / "SHA256SUMS.txt").read_bytes(), first_sums)
            manifest = json.loads(first_checksums)
            self.assertNotIn(
                "PACKAGE_CHECKSUMS.json",
                {entry["path"] for entry in manifest["files"]},
            )
            self.assertIn(
                "SHA256SUMS.txt",
                {entry["path"] for entry in manifest["files"]},
            )

    def test_gen_015_depends_only_on_completed_transport_prerequisites(self) -> None:
        backlog = json.loads(
            (ROOT / "planning" / "backlog.json").read_text(encoding="utf-8")
        )
        issues = {issue["id"]: issue for issue in backlog["issues"]}

        issue = issues["GEN-015"]
        self.assertEqual(issue["depends_on"], ["GEN-013", "GEN-014"])
        self.assertNotIn("GEN-050", issue["depends_on"])
        self.assertNotIn("GEN-051", issue["depends_on"])

        acceptance = " ".join(issue["acceptance_criteria"]).lower()
        for boundary in (
            "oauth",
            "edge",
            "changed deployment",
            "unchanged no-op",
            "no table storage queue model-run persistence",
            "pending human review",
        ):
            self.assertIn(boundary, acceptance)


if __name__ == "__main__":
    unittest.main()

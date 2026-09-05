from __future__ import annotations

import json
import shutil
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NODE = shutil.which("node")
ENTRYPOINT = ROOT / "standalone" / "biobirddna_genetics_node.mjs"
INPUT = ROOT / "examples" / "input" / "cross_autosomal_recessive.json"


@unittest.skipUnless(NODE, "Node.js is not installed")
class NodeStandaloneTests(unittest.TestCase):
    def test_documented_entrypoint_runs_from_windows_path_with_spaces(self) -> None:
        completed = subprocess.run(
            [NODE, str(ENTRYPOINT), "cross", str(INPUT)],
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            capture_output=True,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(completed.stdout.strip(), "Node entrypoint produced no output")
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["model_id"], "M1_AR")
        self.assertEqual(payload["limits"]["gene_status"], "model_assumption_only")


if __name__ == "__main__":
    unittest.main()

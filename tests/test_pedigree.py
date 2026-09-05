from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biobirddna_genetics.enums import Sex
from biobirddna_genetics.models import PedigreeIndividual
from biobirddna_genetics.pedigree import Pedigree


class PedigreeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.pedigree = Pedigree(
            [
                PedigreeIndividual("S", Sex.MALE),
                PedigreeIndividual("D", Sex.FEMALE),
                PedigreeIndividual("B1", Sex.MALE, "S", "D"),
                PedigreeIndividual("B2", Sex.FEMALE, "S", "D"),
                PedigreeIndividual("U", Sex.FEMALE),
            ]
        )

    def test_full_sibling_relationship(self) -> None:
        self.assertAlmostEqual(self.pedigree.relationship("B1", "B2"), 0.5)
        self.assertAlmostEqual(self.pedigree.offspring_inbreeding("B1", "B2"), 0.25)
        self.assertTrue(self.pedigree.are_full_siblings("B1", "B2"))

    def test_unrelated_founder(self) -> None:
        self.assertAlmostEqual(self.pedigree.relationship("B1", "U"), 0.0)
        self.assertAlmostEqual(self.pedigree.offspring_inbreeding("B1", "U"), 0.0)


if __name__ == "__main__":
    unittest.main()

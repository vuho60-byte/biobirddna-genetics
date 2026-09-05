from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class ResearchCoverageTests(unittest.TestCase):
    def test_exact_tracker_ids_and_snapshot_digest(self):
        coverage = json.loads((ROOT / 'knowledge/research-coverage-v3.json').read_text(encoding='utf-8'))
        snapshot = ROOT / coverage['source']['path']
        self.assertEqual(hashlib.sha256(snapshot.read_bytes()).hexdigest(), coverage['source']['sha256'])
        rows = json.loads(snapshot.read_text(encoding='utf-8'))['values'][1:]
        self.assertEqual({r[0] for r in rows if r}, {t['topic_id'] for t in coverage['topics']})
        self.assertEqual(len(coverage['topics']), 27)


def _topic_test(topic_id):
    def check(self):
        coverage = json.loads((ROOT / 'knowledge/research-coverage-v3.json').read_text(encoding='utf-8'))
        topic = next(t for t in coverage['topics'] if t['topic_id'] == topic_id)
        rows = json.loads((ROOT / coverage['source']['path']).read_text(encoding='utf-8'))['values']
        source_row = rows[topic['source']['row'] - 1]
        self.assertEqual(source_row[0], topic_id)
        self.assertEqual(source_row[4], topic['title'])
        self.assertEqual(source_row[5], topic['tracker_status'])
        self.assertEqual(source_row[17], topic['tracker_evidence_status'])
        self.assertEqual(topic['evidence_grade'], 'not_assessed')
        self.assertEqual(topic['scientific_confidence'], 'not_assessed')
        self.assertIn(topic['evidence'], ('missing', 'local_ontology_only'))
        if topic['evidence'] == 'missing':
            self.assertEqual(topic['research_findings'], [])
        self.assertTrue(topic['limits'])
        for path in topic['implemented_in']:
            self.assertTrue((ROOT / path).is_file(), path)
        for entry in topic['behavior_tests']:
            path, name = entry.split('::')
            self.assertIn('def ' + name + '(', (ROOT / path).read_text(encoding='utf-8'))
        for candidate in topic['candidate_models']:
            self.assertIn(candidate, source_row[25] + source_row[28])
        self.assertIn(f'tests/test_research_coverage.py::test_{topic_id}_source_and_boundary', topic['test_ids'])
    return check


for _id in range(1, 28):
    setattr(ResearchCoverageTests, f'test_T{_id:02d}_source_and_boundary', _topic_test(f'T{_id:02d}'))

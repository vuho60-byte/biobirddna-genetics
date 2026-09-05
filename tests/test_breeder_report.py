from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))


class BreederReportTests(unittest.TestCase):
    def test_missing_zero_unknown_sex_and_loss_are_preserved(self):
        from biobirddna_genetics.breeder_report import build_breeder_report
        request = self.request()
        child = request['records']['offspring'][0]
        child.update(outcome='lost_to_followup', survival_status='unknown', sex='unknown',
                     sex_method='unknown', sex_confidence=None)
        result = build_breeder_report(request)
        self.assertIn('Mất theo dõi', result['markdown'])
        self.assertIn('chưa rõ', result['markdown'])
        self.assertIsNone(result['records']['clutches'][0]['stage_counts']['survived_to_first_molt']['count'])
        child.update(outcome='post_hatch_death', survival_status='died_before_first_molt')
        request['records']['clutches'][0]['stage_counts']['survived_to_first_molt'] = {
            'status': 'observed', 'coverage': 'complete', 'count': 0}
        result = build_breeder_report(request)
        self.assertIn('Chết sau nở', result['markdown'])
        self.assertEqual(result['records']['clutches'][0]['stage_counts']['survived_to_first_molt']['count'], 0)
        self.assertIn('Sống tới thay lông đầu: 0;', result['markdown'])

    def test_bad_provenance_stage_and_unknown_pairing_fail_closed(self):
        from biobirddna_genetics.breeder_report import build_breeder_report
        from biobirddna_genetics.errors import ValidationError
        for field, value in [('source', ''), ('synthetic', 'true')]:
            request = self.request()
            request['provenance'][field] = value
            with self.subTest(field=field), self.assertRaises(ValidationError):
                build_breeder_report(request)
        request = self.request()
        request['records']['clutches'][0]['stage_counts']['eggs_laid']['count'] = 4
        with self.assertRaisesRegex(ValidationError, 'stage_offspring_count_conflict'):
            build_breeder_report(request)
        request = self.request()
        request['simulations'] = [{'pairing_id': 'wrong', 'cross_request': {}}]
        with self.assertRaisesRegex(ValidationError, 'unknown_or_duplicate'):
            build_breeder_report(request)

    def test_bad_pedigree_and_markdown_labels(self):
        from biobirddna_genetics.breeder_report import build_breeder_report
        from biobirddna_genetics.errors import ValidationError
        request = self.request()
        request['pedigree'] = [{'individual_id': 'CM-IND-A', 'sex': 'female'}]
        with self.assertRaisesRegex(ValidationError, 'pedigree_sex_conflict'):
            build_breeder_report(request)
        request = self.request()
        request['records']['pairings'][0]['male_parent_phenotype'] = '<script>x</script>|[click](evil)'
        result = build_breeder_report(request)
        table = result['markdown'].split('## Mô phỏng')[0]
        self.assertNotIn('<script>', table)
        self.assertIn('&lt;script&gt;', table)
        self.assertIn('\\|\\[click\\]', table)

    def test_cli_rejects_input_overwrite_and_invalid_input_without_touching_output(self):
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory() as folder:
            source, output = Path(folder) / 'input.json', Path(folder) / 'report.md'
            source.write_text('{}', encoding='utf-8')
            output.write_text('preserve', encoding='utf-8')
            cli = ROOT / 'scripts/breeder_report.py'
            bad = subprocess.run([sys.executable, '-B', str(cli), '--input', str(source), '--output', str(output)], capture_output=True)
            self.assertEqual(bad.returncode, 2)
            self.assertEqual(output.read_text(), 'preserve')
            same = subprocess.run([sys.executable, '-B', str(cli), '--input', str(source), '--output', str(source)], capture_output=True)
            self.assertEqual(same.returncode, 2)
            self.assertEqual(source.read_text(), '{}')

    def test_new_golden_and_cli_are_byte_stable(self):
        import subprocess
        import tempfile
        from biobirddna_genetics.breeder_report import build_breeder_report
        source = ROOT / 'examples/input/breeder_report.json'
        request = json.loads(source.read_text(encoding='utf-8'))
        expected = (ROOT / 'examples/output/breeder_report.md').read_bytes()
        self.assertEqual(build_breeder_report(request)['markdown'].encode('utf-8'), expected)
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'report.md'
            run = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/breeder_report.py'),
                                  '--input', str(source), '--output', str(output)], capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertEqual(output.read_bytes(), expected)

    def test_explicit_cross_keeps_core_probabilities_and_all_limits(self):
        from biobirddna_genetics.breeder_report import build_breeder_report
        from biobirddna_genetics.service import dispatch
        request = self.request()
        cross = json.loads((ROOT / 'examples/input/cross_autosomal_recessive.json').read_text())
        request['simulations'] = [{'pairing_id': request['records']['pairings'][0]['pairing_id'], 'cross_request': cross}]
        result = build_breeder_report(request)
        self.assertEqual(result['simulations'][0]['result'], dispatch('cross', cross))
        self.assertEqual(sum(result['simulations'][0]['result']['phenotype_probabilities'].values()), 1)
        self.assertIn('25.00%', result['markdown'])
        self.assertIn('simulation_not_confirmation', result['markdown'])

    def test_full_siblings_and_welfare_are_blocked_even_with_zero_limit_override_attempt(self):
        from biobirddna_genetics.breeder_report import build_breeder_report
        from biobirddna_genetics.errors import ValidationError
        request = self.request()
        request['pedigree'] = [
            {'individual_id': 'CM-IND-A', 'sex': 'male', 'sire_id': 'S', 'dam_id': 'D'},
            {'individual_id': 'CM-IND-B', 'sex': 'female', 'sire_id': 'S', 'dam_id': 'D'},
            {'individual_id': 'S', 'sex': 'male'}, {'individual_id': 'D', 'sex': 'female'}]
        request['welfare_blocked_pairing_ids'] = [request['records']['pairings'][0]['pairing_id']]
        result = build_breeder_report(request)
        review = result['pairing_reviews'][0]
        self.assertEqual(review['offspring_inbreeding'], .25)
        self.assertEqual(review['status'], 'blocked')
        self.assertIn('full_sibling_pair', review['blocked_reasons'])
        self.assertIn('welfare_blocked', review['blocked_reasons'])
        self.assertIn('DỪNG', result['markdown'])
        request['block_direct_ancestor_and_full_sibling'] = False
        with self.assertRaises(ValidationError):
            build_breeder_report(request)

    def test_direct_parent_is_blocked_when_only_child_has_a_pedigree_row(self):
        from biobirddna_genetics.breeder_report import build_breeder_report
        request = self.request()
        request['pedigree'] = [{'individual_id': 'CM-IND-B', 'sex': 'female', 'sire_id': 'CM-IND-A'}]
        result = build_breeder_report(request)
        self.assertIn('direct_ancestor_pair', result['pairing_reviews'][0]['blocked_reasons'])

    def request(self):
        records = json.loads((ROOT / 'examples/input/normalize_breeding_records.json').read_text(encoding='utf-8'))
        records['pairings'] = records['pairings'][:1]
        return {'schema_version': '1.0.0', 'records': records,
                'provenance': {'source': 'synthetic fixture', 'synthetic': True}}

    def test_records_remain_exactly_the_core_output_and_missing_is_visible(self):
        from biobirddna_genetics.breeder_report import build_breeder_report
        from biobirddna_genetics.service import dispatch
        request = self.request()
        before = copy.deepcopy(request)
        result = build_breeder_report(request)
        self.assertEqual(result['records'], dispatch('normalize-records', request['records']))
        self.assertEqual(request, before)
        self.assertIn('chưa ghi nhận', result['markdown'])
        self.assertIn('Không thụ tinh', result['markdown'])
        self.assertIn('Nở thất bại', result['markdown'])
        self.assertIn('not_assessed', result['markdown'])
        self.assertEqual(result['pairing_reviews'][0]['offspring_inbreeding'], None)
        self.assertEqual(result['simulations'], [])

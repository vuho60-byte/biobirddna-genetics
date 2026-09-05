"""Map a preserved tracker snapshot to local software coverage, never gene claims."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'knowledge/sources/tracker-v3-20260905.json'
# Deliberate, reviewed links to existing primitives, not completion of topic exit gates.
COVERAGE = {
    'T02': ('pigment', 'pigment_and_evidence'),
    'T04': ('pedigree', 'pedigree'),
    'T08': ('breeding_records', 'breeding_records'),
    'T11': ('breeding_records', 'breeding_records'),
    'T12': ('breeding_records', 'breeding_records'),
    'T13': ('breeding_records', 'breeding_records'),
    'T14': ('pigment', 'pigment_and_evidence'),
    'T15': ('pigment', 'pigment_and_evidence'),
    'T16': ('pigment', 'pigment_and_evidence'),
    'T20': ('simulation', 'simulation_recovery'),
    'T21': ('inference', 'inference'),
    'T23': ('pedigree', 'pedigree'),
}
CANDIDATES = {
    'T14': ['albinism', 'leucism', 'progressive greying'],
    'T15': ['autosomal', 'liên kết Z', 'đa gene', 'không di truyền'],
    'T21': ['mô hình cạnh tranh'],
}


def build() -> dict:
    path = ROOT / SOURCE
    snapshot = json.loads(path.read_text(encoding='utf-8'))
    topics = []
    for row_number, values in enumerate(snapshot['values'][1:], 2):
        if not values:
            continue
        row = dict(zip(snapshot['values'][0], values))
        topic_id = row['Topic ID']
        if topic_id not in {f'T{i:02d}' for i in range(1, 28)}:
            raise ValueError('unexpected_tracker_topic_id')
        implemented, behavior_tests = [], []
        if topic_id in COVERAGE:
            module, test_module = COVERAGE[topic_id]
            implemented.append(f'src/biobirddna_genetics/{module}.py')
            test_path = f'tests/test_{test_module}.py'
            tree = ast.parse((ROOT / test_path).read_text(encoding='utf-8'))
            behavior_tests = [f'{test_path}::{node.name}' for node in ast.walk(tree)
                              if isinstance(node, ast.FunctionDef) and node.name.startswith('test_')]
        if topic_id == 'T02':
            implemented.append('knowledge/T02_knowledge_base_v1.0.json')
        topics.append({
            'topic_id': topic_id, 'title': row['Topic'],
            'source': {'path': SOURCE, 'row': row_number, 'version': 'v3',
                       'referenced_document': row.get('Source', '')},
            'tracker_status': row.get('Status', ''),
            'tracker_evidence_status': row.get('Evidence Status', ''),
            'tracker_evidence_grade': row.get('Evidence Grade', ''),
            'tracker_allowed_conclusion': row.get('Allowed Conclusion', ''),
            'tracker_evidence_gap': row.get('Evidence Gap', ''),
            'tracker_source_review_file': row.get('Source Review File', ''),
            'tracker_file_availability': row.get('File Availability', ''),
            'candidate_models': CANDIDATES.get(topic_id, []),
            'candidate_status': 'tracker_questions_only_not_supported_models',
            'evidence': 'local_ontology_only' if topic_id == 'T02' else 'missing',
            'evidence_grade': 'not_assessed', 'scientific_confidence': 'not_assessed',
            'research_findings': [],
            'status': 'partial_software_support' if implemented else 'not_implemented',
            'implemented_in': implemented, 'behavior_tests': behavior_tests,
            'test_ids': [f'tests/test_research_coverage.py::test_{topic_id}_source_and_boundary'],
            'exit_gate': row.get('Đầu ra / Exit Gate', ''),
            'limits': [
                'Tracker là lộ trình và câu hỏi nghiên cứu, không chứng minh cơ chế di truyền.',
                'Test mapping kiểm nguồn và ranh giới; test core kiểm phần mềm, không xác minh kết quả nghiên cứu.',
                'Không tự nâng trạng thái tracker hoặc xác nhận gene; thiếu tài liệu nghiên cứu thì evidence: missing.',
                'Mô hình đa gene, survival và cấu trúc quần thể chưa được triển khai; không suy rộng từ primitive một locus.',
            ],
        })
    if len(topics) != 27 or len({t['topic_id'] for t in topics}) != 27:
        raise ValueError('tracker_must_contain_27_unique_topics')
    return {
        'schema_version': '1.0.0', 'tracker_version': 'v3',
        'retrieved_on': '2026-09-05',
        'source': {'path': SOURCE, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                   'spreadsheet_id': '1Av2TyP3XvPp_BPsxr8HdJkGzMh-kpYIYFnZlOmDKv90',
                   'sheet': 'Untitled', 'range': 'A1:AL30',
                   'title': 'Lộ trình nghiên cứu màu lông chào mào – Tracker v3 (27 topic) — BẢN ĐANG DÙNG'},
        'topics': topics,
    }


if __name__ == '__main__':
    (ROOT / 'knowledge/research-coverage-v3.json').write_text(
        json.dumps(build(), ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

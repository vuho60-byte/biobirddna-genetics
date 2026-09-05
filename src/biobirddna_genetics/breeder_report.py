"""Local Vietnamese breeder reports over the existing research primitives."""
from __future__ import annotations

import hashlib
import html
import json
from typing import Any

from . import __version__
from .enums import Sex
from .errors import ValidationError
from .models import PedigreeIndividual
from .pedigree import Pedigree
from .service import dispatch

MAX_OFFSPRING_INBREEDING = 0.0625  # Existing rank-mates default; no owner override here.
SEX = {'male': 'trống', 'female': 'mái', 'unknown': 'chưa rõ'}
OUTCOME = {'alive': 'Còn sống', 'unfertilized': 'Không thụ tinh',
           'embryo_death': 'Chết phôi', 'unhatched': 'Không nở',
           'hatch_failure': 'Nở thất bại', 'post_hatch_death': 'Chết sau nở',
           'lost_to_followup': 'Mất theo dõi', 'unknown': 'Chưa rõ'}
STAGE = {'eggs_laid': 'Trứng đẻ', 'fertile_eggs': 'Trứng có phôi',
         'hatched': 'Đã nở', 'survived_to_first_molt': 'Sống tới thay lông đầu'}
REASON = {'direct_ancestor_pair': 'Quan hệ trực hệ (cha/mẹ–con hoặc tổ tiên–hậu duệ)',
          'full_sibling_pair': 'Anh chị em ruột', 'welfare_blocked': 'Có cờ dừng phúc lợi',
          'offspring_inbreeding_exceeds_limit': 'Hệ số F vượt ngưỡng chặn mặc định'}


def _object(value, allowed, path):
    if not isinstance(value, dict):
        raise ValidationError(f'{path}:must_be_object')
    unknown = set(value) - set(allowed)
    if unknown:
        raise ValidationError(f'{path}:unknown_field:{sorted(unknown)[0]}')
    return value


def _text(value, path):
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f'{path}:required_text')
    return value


def _array(value, path):
    if not isinstance(value, list):
        raise ValidationError(f'{path}:must_be_array')
    return value


def _safe(value):
    text = html.escape(str(value), quote=False).replace('\n', ' ').replace('\r', ' ')
    for char in ('\\', '`', '*', '_', '[', ']', '|', '#'):
        text = text.replace(char, '\\' + char)
    return text


def _pedigree(request, pairings):
    raw = _array(request.get('pedigree', []), 'pedigree')
    entries, by_id = [], {}
    for i, item in enumerate(raw):
        where = f'pedigree/{i}'
        _object(item, {'individual_id', 'sex', 'sire_id', 'dam_id'}, where)
        identifier = _text(item.get('individual_id'), where + '/individual_id')
        if identifier in by_id:
            raise ValidationError('pedigree_duplicate_individual_id')
        try:
            sex = Sex(item.get('sex', 'unknown'))
        except (ValueError, TypeError) as exc:
            raise ValidationError(where + '/sex:invalid') from exc
        for key in ('sire_id', 'dam_id'):
            if item.get(key) is not None:
                _text(item[key], where + '/' + key)
        entries.append(PedigreeIndividual(identifier, sex, item.get('sire_id'), item.get('dam_id')))
        by_id[identifier] = item
    for pair in pairings:
        for key, expected in (('male_parent_id', 'male'), ('female_parent_id', 'female')):
            item = by_id.get(pair[key])
            if item and item.get('sex', 'unknown') not in (expected, 'unknown'):
                raise ValidationError(f'pairings/{pair["pairing_id"]}/{key}:pedigree_sex_conflict')
    # Materialize referenced unknown founders only for relationship calculations.
    # A child naming its sire is already enough to block a direct pairing,
    # even if the sire has no standalone pedigree row.
    if entries:
        referenced = {p[key] for p in pairings for key in ('male_parent_id', 'female_parent_id')}
        referenced.update(e.sire_id for e in entries if e.sire_id)
        referenced.update(e.dam_id for e in entries if e.dam_id)
        entries.extend(PedigreeIndividual(i) for i in sorted(referenced - set(by_id)))
    return Pedigree(entries) if entries else None, by_id


def build_breeder_report(request: dict[str, Any]) -> dict[str, Any]:
    """Normalize records without mutation and render only explicit assumptions.

    No genotype is derived from parental labels or offspring records. The
    separate simulations array must identify a pairing and supply a cross request.
    """
    _object(request, {'schema_version', 'records', 'provenance', 'pedigree',
                      'simulations', 'welfare_blocked_pairing_ids'}, 'report')
    if request.get('schema_version') != '1.0.0':
        raise ValidationError('report/schema_version:unsupported')
    provenance = _object(request.get('provenance'), {'source', 'synthetic'}, 'provenance')
    _text(provenance.get('source'), 'provenance/source')
    if not isinstance(provenance.get('synthetic'), bool):
        raise ValidationError('provenance/synthetic:must_be_boolean')
    records = dispatch('normalize-records', request.get('records'))
    pairs = {pair['pairing_id']: pair for pair in records['pairings']}
    pedigree, pedigree_rows = _pedigree(request, records['pairings'])
    blocked_ids = _array(request.get('welfare_blocked_pairing_ids', []), 'welfare_blocked_pairing_ids')
    if any(not isinstance(i, str) or i not in pairs for i in blocked_ids):
        raise ValidationError('welfare_blocked_pairing_ids:unknown_pairing')
    reviews = []
    for identifier, pair in pairs.items():
        reasons, coefficient = [], None
        male, female = pair['male_parent_id'], pair['female_parent_id']
        if pedigree and (male in pedigree_rows or female in pedigree_rows):
            coefficient = pedigree.offspring_inbreeding(male, female)
            relation = pedigree.hard_relation_block(male, female)
            if relation:
                reasons.append(relation)
            if coefficient > MAX_OFFSPRING_INBREEDING:
                reasons.append('offspring_inbreeding_exceeds_limit')
        if identifier in blocked_ids:
            reasons.append('welfare_blocked')
        reviews.append({'pairing_id': identifier, 'offspring_inbreeding': coefficient,
                        'blocked_reasons': reasons,
                        'status': 'blocked' if reasons else 'human_review_required',
                        'ancestry_status': 'declared_not_dna_verified',
                        'welfare_status': 'blocked' if identifier in blocked_ids else 'not_assessed'})
    simulations, seen = [], set()
    for index, item in enumerate(_array(request.get('simulations', []), 'simulations')):
        where = f'simulations/{index}'
        _object(item, {'pairing_id', 'cross_request'}, where)
        identifier = _text(item.get('pairing_id'), where + '/pairing_id')
        if identifier not in pairs or identifier in seen:
            raise ValidationError(where + '/pairing_id:unknown_or_duplicate')
        seen.add(identifier)
        cross_request = item.get('cross_request')
        result = dispatch('cross', cross_request)
        simulations.append({'pairing_id': identifier, 'assumptions': cross_request, 'result': result})
    encoded = json.dumps(request, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)
    report = {'schema_version': '1.0.0', 'engine_version': __version__,
              'input_sha256': hashlib.sha256(encoded.encode('utf-8')).hexdigest(),
              'provenance': dict(provenance), 'records': records, 'simulations': simulations,
              'pairing_reviews': reviews, 'evidence_grade': 'not_assessed',
              'limits': {'decision_support_only': True, 'gene_status': 'not_assessed',
                         'human_review_status': 'pending', 'pedigree_declared_not_verified': True,
                         'unknown_founders_assumed_unrelated_for_calculation_only': True,
                         'max_offspring_inbreeding': MAX_OFFSPRING_INBREEDING,
                         'health_safety_not_assessed': True}}
    report['markdown'] = render_breeder_report(report)
    return report


def render_breeder_report(report: dict[str, Any]) -> str:
    records = report['records']
    lines = ['# Báo cáo hồ sơ cặp ghép', '',
             f"Hồ sơ: {_safe(records['case_id'])}. BioBirdDNA {_safe(report['engine_version'])}.", '',
             f"Nguồn khai báo: {_safe(report['provenance']['source'])}.",
             ('Dữ liệu minh họa tổng hợp (synthetic), không phải đàn chim thật.'
              if report['provenance']['synthetic'] else 'Dữ liệu do người dùng cung cấp, chưa được xác minh độc lập.'), '',
             '**Chỉ hỗ trợ xem xét; cần người có trách nhiệm duyệt trước khi ghép thật.**', '',
             '## Bố mẹ và đời con', '',
             '| Cặp | Bố: mã / kiểu hình | Mẹ: mã / kiểu hình | Nguồn giới tính bố / mẹ |',
             '|---|---|---|---|']
    for p in records['pairings']:
        lines.append(f"| {_safe(p['pairing_id'])} | {_safe(p['male_parent_id'])} / {_safe(p['male_parent_phenotype'])} | {_safe(p['female_parent_id'])} / {_safe(p['female_parent_phenotype'])} | {_safe(p['male_parent_sex_method'])} / {_safe(p['female_parent_sex_method'])} |")
    for clutch in records['clutches']:
        lines += ['', f"### Lứa {_safe(clutch['clutch_id'])} — cặp {_safe(clutch['pairing_id'])}", '',
                  f"Độ đầy đủ hồ sơ con: {_safe(clutch['offspring_records_coverage'])}.", '']
        for stage, value in clutch['stage_counts'].items():
            count = 'chưa ghi nhận' if value['status'] == 'missing' else str(value['count'])
            lines.append(f"- {STAGE[stage]}: {count}; độ đầy đủ: {_safe(value['coverage'])}.")
    lines += ['', '| Mã trứng/con | Lứa | Kết quả | Giới tính | Diễn biến kiểu hình |', '|---|---|---|---|---|']
    for child in records['offspring']:
        timeline = '; '.join(f"{o['sequence_index']} / {o['stage']}: {o['phenotype'] if o['status'] == 'observed' else 'chưa ghi nhận'}" for o in child['phenotype_timeline']) or 'chưa ghi nhận'
        lines.append(f"| {_safe(child['offspring_record_id'])} | {_safe(child['clutch_id'])} | {OUTCOME[child['outcome']]} | {SEX[child['sex']]} | {_safe(timeline)} |")
    lines += ['', '## Mô phỏng và xác suất kiểu hình', '',
              'Tên màu và số con quan sát không tự cho biết kiểu gene. Chỉ tính khi có mô hình và kiểu gene giả định được nhập riêng.',
              'Xác suất dưới đây tính trên toàn bộ con, gộp cả giới tính (không phải tỷ lệ có điều kiện trong riêng trống hoặc mái).',
              'Ví dụ 25% mái màu mục tiêu nghĩa là 25% tổng số con; một lứa nhỏ không được bảo đảm đúng tỷ lệ này.', '']
    if not report['simulations']:
        lines += ['Chưa đủ giả định mô hình/kiểu gene để tính xác suất. Không suy gene từ kiểu hình bố mẹ.', '']
    for simulation in report['simulations']:
        result = simulation['result']
        lines += [f"### Cặp {_safe(simulation['pairing_id'])} — mô hình giả định {_safe(result['model_id'])}", '',
                  'Giả định đầu vào (allele ký hiệu, không phải tên gene đã xác nhận):', '', '```json',
                  json.dumps(simulation['assumptions'], ensure_ascii=False, indent=2), '```', '',
                  '| Giới tính và kiểu hình | Xác suất trên tổng số con |', '|---|---:|']
        for label, probability in sorted(result['phenotype_probabilities'].items()):
            sex, phenotype = label.split(':', 1)
            lines.append(f'| {SEX.get(sex, _safe(sex))}: {_safe(phenotype)} | {100 * probability:.2f}% |')
    lines += ['', '## Cận huyết và khuyến nghị xem xét cặp', '',
              'F tính từ phả hệ khai báo. Thiếu tổ tiên không đồng nghĩa không có quan hệ; F = 0 không chứng minh không cận huyết.',
              'Ngưỡng chặn mặc định của engine là F > 6,25%; đây không phải ngưỡng an toàn sức khỏe được chứng nhận.',
              'Nếu chưa có cờ phúc lợi, sức khỏe vẫn là chưa đánh giá. Ưu tiên hoàn thiện hồ sơ và để người phụ trách/thú y xem xét.', '',
              '| Cặp | F kỳ vọng của con | Khuyến nghị |', '|---|---:|---|']
    for review in report['pairing_reviews']:
        coefficient = review['offspring_inbreeding']
        shown = 'chưa đủ phả hệ' if coefficient is None else f'{100 * coefficient:.2f}%'
        advice = ('DỪNG xem xét ghép: ' + '; '.join(REASON[r] for r in review['blocked_reasons'])
                  if review['blocked_reasons'] else 'Chưa phê duyệt ghép; bổ sung/phúc tra phả hệ, sức khỏe và chờ người phụ trách duyệt')
        lines.append(f"| {_safe(review['pairing_id'])} | {shown} | {advice} |")
    lines += ['', '## Dữ liệu thiếu và giới hạn bằng chứng', '',
              'Evidence grade: `not_assessed`. Không tự cấp G-level. Posterior/mô phỏng không xác nhận gene.',
              'Số 0 đã quan sát khác dữ liệu chưa ghi nhận. Giữ cả trứng không nở, chết phôi, chết sau nở và mất theo dõi.',
              'Mô phỏng chưa mô hình hóa sống sót theo genotype, cấu trúc quần thể hoặc nhiều locus.', '']
    for missing in records['missing_data']:
        lines.append(f"- {_safe(missing['field'])}: {_safe(missing['impact'])}.")
    if not records['missing_data']:
        lines.append('- Bộ chuẩn hóa không ghi nhận trường thiếu; điều này không chứng minh hồ sơ đúng hoặc đầy đủ về khoa học.')
    lines += ['', '### Toàn bộ Limits từ kết quả core và báo cáo', '', '```json',
              json.dumps({'records': records['limits'], 'simulations': [s['result']['limits'] for s in report['simulations']],
                          'report': report['limits']}, ensure_ascii=False, indent=2), '```', '',
              '### Hồ sơ chuẩn hóa để đối chiếu (giữ nguyên kết quả core)', '', '```json',
              json.dumps(records, ensure_ascii=False, indent=2), '```', '',
              f"Mã SHA-256 input JSON chuẩn hóa: `{report['input_sha256']}`.", '']
    return '\n'.join(lines)
